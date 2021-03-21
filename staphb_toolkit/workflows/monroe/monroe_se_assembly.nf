#!/usr/bin/env nextflow

//Description: Workflow for quality control of raw illumina reads
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""
params.primerSet = ""
params.primerPath = workflow.projectDir + params.primerSet
params.report = ""
params.pipe = ""

//setup channel to read in and pair the fastq files
Channel
    .fromPath( "${params.reads}/*.{fastq,fq}.*")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

Channel
  .fromPath(params.primerPath, type:'file')
  .ifEmpty{
    println("A bedfile for primers is required. Set with 'params.primerPath'.")
    exit 1
  }
  .view { "Primer BedFile : $it"}
  .set { primer_bed }

//Step0: Preprocess reads - change name to end at first underscore if specified
process preProcess {
  input:
  file(read) from raw_reads

  output:
  file(read) into raw_reads_trim
  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    """
  }else{
  """
  """
  }
}

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"

  input:
  file(read) from raw_reads_trim

  output:
  file("*_trimmed.fastq.gz") into trimmed_reads

  script:
  """
  # get samplename by dropping file extension
  filename=${read}
  samplename=\$(echo \${filename} | cut -d "." -f 1)

  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads ${task.cpus} ${read} \${samplename}_trimmed.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} > \${samplename}.trim.stats.txt
  """
}
//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"

  input:
  file(read) from trimmed_reads

  output:
  file("*_1.clean.fastq.gz") into cleaned_reads
  script:
  """
  # get samplename by dropping string after first underscore
  filename=${read}
  samplename=\$(echo \${filename} | cut -d "_" -f 1)

  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${read} out1=\${samplename}.rmadpt_1.fastq.gz ref=/bbmap/resources/adapters.fa stats=\${samplename}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=\${samplename}.rmadpt_1.fastq.gz out1=\${samplename}_1.clean.fastq.gz outm=\${samplename}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=\${samplename}.phix.stats.txt
  """
}

cleaned_reads
  .combine(primer_bed)
  .set { cleaned_reads_with_primer_bed }

//Assemble cleaned reads with iVar
process ivar {
  tag "$name"

  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.sorted.bam"
  publishDir "${params.outdir}/SC2_reads", mode: 'copy',pattern:"*_SC2*.fastq.gz"
  publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "*_consensus.fasta"


  input:
  set file(read), file(primer_bed) from cleaned_reads_with_primer_bed

  output:
  file "*_consensus.fasta" into assembled_genomes
  file "*.sorted.bam" into alignment_file
  file "*_SC2*.fastq.gz" into sc2_reads

  shell:
"""
# get samplename by dropping string after first underscore
filename=${read}
samplename=\$(echo \${filename} | cut -d "_" -f 1)

ln -s /reference/nCoV-2019.reference.fasta ./nCoV-2019.reference.fasta
minimap2 -K 20M -x sr -a ./nCoV-2019.reference.fasta ${read} | samtools view -u -h -F 4 - | samtools sort > SC2.bam
samtools index SC2.bam
samtools flagstat SC2.bam
samtools sort -n SC2.bam > SC2_sorted.bam
samtools fastq -f2 -F4 -1 \${samplename}_SC2_R1.fastq.gz SC2_sorted.bam -s singletons.fastq.gz

ivar trim -i SC2.bam -b !{primer_bed} -p ivar -e

samtools sort  ivar.bam > \${samplename}.sorted.bam
samtools index \${samplename}.sorted.bam
samtools flagstat \${samplename}.sorted.bam

samtools mpileup -f ./nCoV-2019.reference.fasta -d 1000000 -A -B -Q 0 \${samplename}.sorted.bam | ivar consensus -p ivar -m ${params.ivar_mindepth} -t ${params.ivar_minfreq} -n N
echo \">\${samplename}\" > \${samplename}_consensus.fasta

seqtk seq -U -l 50 ivar.fa | tail -n +2 >> \${samplename}_consensus.fasta
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
  file(assembly) from assembled_genomes

  output:
  file "*_lineage_report.csv" into pangolin_lineages

  shell:
  """
  # get samplename by dropping string after first underscore
  filename=${assembly}
  samplename=\$(echo \${filename} | cut -d "_" -f 1)

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
