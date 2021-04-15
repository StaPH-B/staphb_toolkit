#!/usr/bin/env nextflow

//Description: Adaptation of ARTIC Network nCoV-2019 Bioinformatics SOP
//Available: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Authors of this Nextflow: Kelsey Florek and Abigail Shockey
//Email: kelsey.florek@slh.wisc.edu

//starting parameters
params.fast5_dir = ""
params.fastq_dir = ""
params.sequencing_summary = ""
params.outdir = ""
params.primers = ""
params.run_prefix = "artic_ncov19"
params.pipe = ""

// Input channels
Channel
    .value( "${params.primers}")
    .ifEmpty { exit 1, "Primers used must be included." }
    .set { polish_primers }

// If we have fast5 files then start with basecalling
if(params.basecalling){
  Channel
      .fromPath( "${params.fast5_dir}")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .into { raw_fast5; polish_fast5 }

  process guppy_basecalling {
    input:
      file(fast5s) from raw_fast5.collect()

    output:
      file "fastq/*.fastq" into fastq_reads
      file "fastq/sequencing_summary.txt" into sequencing_summary

    script:
      if(params.basecalling_mode == "fast"){
        """
        guppy_basecaller --chunk_size ${params.chunk_size} --chunks_per_runner ${params.chunks_per_runner} --gpu_runners_per_device ${params.gpu_runners_per_device} ${params.basecalling_params} -i . -s fastq -x auto -r
        """
      }else{
        """
        guppy_basecaller --chunk_size ${params.chunk_size} --chunks_per_runner ${params.chunks_per_runner} --gpu_runners_per_device ${params.gpu_runners_per_device} ${params.basecalling_params} -i . -s fastq -x auto -r
        """
      }
  }
}

// If we have already basecalled get fastqs and fast5s for polishing
else {
  Channel
      .fromPath( "${params.fastq_dir}/*.fastq*")
      .ifEmpty { exit 1, "Cannot find any fastq files in: ${params.fastq_dir} Path must not end with /" }
      .set { fastq_reads }

  Channel
      .fromPath( "${params.fast5_dir}")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .set { polish_fast5 }

  if(params.polishing == "nanopolish"){
    Channel
        .fromPath( "${params.sequencing_summary}")
        .ifEmpty { exit 1, "Cannot find sequencing summary in: ${params.sequencing_summary}" }
        .set { sequencing_summary }
    }
}

// Demultiplex fastqs
process guppy_demultiplexing {
  publishDir "${params.outdir}/demultiplexing", mode: 'copy'

  input:
    file(fastqs) from fastq_reads.collect()

  output:
    path("barcode*",type:'dir') into demultiplexed_reads

  script:
    """
      guppy_barcoder -t ${task.cpus} --require_barcodes_both_ends -i . -s . ${params.demultiplexing_params} -q 0 -r
    """
}

// Run artic gupplyplex
process artic_guppyplex {
  publishDir "${params.outdir}/guppyplex", mode: 'copy'
  errorStrategy 'ignore'

  input:
    file(reads) from demultiplexed_reads.flatten()

  output:
    file "${params.run_prefix}_barcode*.fastq" into polish_files

  script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --directory barcode* \
    --prefix ${params.run_prefix}
    """
}

// Run artic pipeline using nanopolish
if(params.polishing=="nanopolish"){
  process artic_nanopolish_pipeline {
    publishDir "${params.outdir}/pipeline_nanopolish", mode: 'copy'
    errorStrategy 'ignore'

    input:
      val primers from polish_primers
      tuple file(fastq), path(fast5path), file(sequencing_summary) from polish_files .combine(polish_fast5) .combine(sequencing_summary)

    output:
      file *.primrtrimmed.bam into alignment_file
      file "*{.primertrimmed.bam,.vcf,.variants.tab}"
      file "*.consensus.fasta" into consensus_fasta

    script:
      """
      # get samplename by dropping file extension
      filename=${fastq}
      samplename=\${filename%.*}

      artic minion --normalise ${params.normalise} --threads ${task.cpus} --scheme-directory /artic-ncov2019/primer_schemes --fast5-directory ${fast5path}  --sequencing-summary ${sequencing_summary} --read-file ${fastq} nCoV-2019/${primers} \$samplename
      """
  }
}

// Run artic pipeline using medaka
else {
  process artic_medaka_pipeline {
    publishDir "${params.outdir}/pipeline_medaka", mode: 'copy'
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: '*.fasta'
    errorStrategy 'ignore'

    input:
      val primers from polish_primers
      file(fastq) from polish_files

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
