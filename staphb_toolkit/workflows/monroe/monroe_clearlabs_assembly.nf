#!/usr/bin/env nextflow

//Description: Adaptation of ARTIC Network nCoV-2019 Bioinformatics SOP
//Available: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Authors of this Nextflow: Kevin Libuit; based on ont_assembly workflow written by Kelsey Florek and Abigail Shockey
//Email: kevin@libuitsci.com

//starting parameters
params.reads = ""
params.outdir = ""
params.run_prefix = "artic_ncov19"
params.primers = "V3"
params.pipe = ""

// Input channels
Channel
    .value( "${params.primers}")
    .ifEmpty { exit 1, "Primers used must be included." }
    .set { polish_primers }

Channel
    .fromPath( "${params.reads}/*.fastq*")
    .ifEmpty { exit 1, "Cannot find any fastq files in: ${params.fastq_dir} Path must not end with /" }
    .set { fastq_reads }

// Run artic pipeline using medaka
process artic_medaka_pipeline {
  publishDir "${params.outdir}/pipeline_medaka", mode: 'copy'
  publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: '*.fasta'
  errorStrategy 'ignore'

  input:
    val primers from polish_primers
    file(fastq) from fastq_reads

  output:
    file "*.primertrimmed.rg.sorted.bam" into alignment_file
    file "*{.primertrimmed.rg,.primers.vcf,.vcf.gz,.trimmed.rg,.fail.vcf}*"
    file "*.consensus.fasta" into consensus_fasta

  script:
    """
    # get samplename by dropping file extension
    filename=${fastq}
    samplename=\${filename%.*}
    artic minion --medaka --normalise ${params.normalise} --threads ${task.cpus} --scheme-directory /artic-ncov2019/primer_schemes --read-file ${fastq} nCoV-2019/${primers} \$samplename
    """
}

//QC of read data
process samtools {
  tag "$name"

  input:
  file(alignment) from alignment_file

  output:
  file "*_samtoolscoverage.tsv" into alignment_qc

  shell:
  """
  # get samplename by dropping extension
  filename=${alignment}
  samplename=\$(echo \${filename} | cut -d "." -f 1)
  samtools coverage ${alignment} -o \${samplename}_samtoolscoverage.tsv
  """
}
//Typing of SC2 assemblies
process pangolin_typing {
  tag "$name"

  publishDir "${params.outdir}/pangolin_reports", mode: 'copy', pattern: "*_lineage_report.csv"

  input:
  file(assembly) from consensus_fasta

  output:
  file "*_lineage_report.csv" into pangolin_lineages

  shell:
  """
  # get samplename by dropping extension
  filename=${assembly}
  samplename=\$(echo \${filename} | cut -d "." -f 1)

  pangolin ${assembly} --outfile \${samplename}_lineage_report.csv
  """
}

//Collect and format report
process assembly_results{
  publishDir "${params.outdir}/", mode: 'copy', pattern: "monroe_summary*.csv"

  echo true

  input:
  file(cg_pipeline_results) from alignment_qc.collect()
  file(pangolin_lineage) from pangolin_lineages.collect()

  output:
  file "monroe_summary*.csv"


  script:
"""
#!/usr/bin/env python3
import os, sys
import glob, csv
import xml.etree.ElementTree as ET
from datetime import datetime

today = datetime.today()
today = today.strftime("%m%d%y")

class result_values:
    def __init__(self,id):
        self.id = id
        self.aligned_bases = "NA"
        self.percent_cvg = "NA"
        self.mean_depth = "NA"
        self.mean_base_q = "NA"
        self.mean_map_q = "NA"
        self.monroe_qc = "NA"
        self.pangolin_lineage = "NA"
        self.pangolin_probability = "NA"
        self.pangolin_notes = "NA"


#get list of result files
samtools_results = glob.glob("*_samtoolscoverage.tsv")
pangolin_results = glob.glob("*_lineage_report.csv")
results = {}

# collect samtools results
for file in samtools_results:
    id = file.split("_samtoolscoverage.tsv")[0]
    result = result_values(id)
    monroe_qc = []
    with open(file,'r') as tsv_file:
        tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
        for line in tsv_reader:
            result.aligned_bases = line["covbases"]
            result.percent_cvg = line["coverage"]
            if float(line["coverage"]) < 98:
                monroe_qc.append("coverage <98%")
            result.mean_depth = line["meandepth"]
            result.mean_base_q = line["meanbaseq"]
            if float(line["meanbaseq"]) < 30:
                monroe_qc.append("meanbaseq < 30")
            result.mean_map_q = line["meanmapq"]
            if float(line["meanmapq"]) < 30:
                monroe_qc.append("meanmapq < 30")
        if len(monroe_qc) == 0:
            result.monroe_qc = "PASS"
        else:
            result.monroe_qc ="WARNING: " + '; '.join(monroe_qc)

    file = (id + "_lineage_report.csv")
    with open(file,'r') as csv_file:
        csv_reader = list(csv.DictReader(csv_file, delimiter=","))
        for line in csv_reader:
            if line["status"] == "fail":
                result.pangolin_lineage = "failed pangolin qc"
            else:
                result.pangolin_lineage = line["lineage"]
                result.pangolin_probability = line["probability"]
                result.pangolin_notes = line["note"]

    results[id] = result


#create output file
with open(f"monroe_summary_{today}.csv",'w') as csvout:
    writer = csv.writer(csvout,delimiter=',')
    writer.writerow(["sample","aligned_bases","percent_cvg", "mean_depth", "mean_base_q", "mean_map_q", "monroe_qc", "pangolin_lineage", "pangolin_probability", "pangolin_notes"])
    for id in results:
        result = results[id]
        writer.writerow([result.id,result.aligned_bases,result.percent_cvg,result.mean_depth,result.mean_base_q,result.mean_map_q,result.monroe_qc,result.pangolin_lineage,result.pangolin_probability,result.pangolin_notes])
"""

}
