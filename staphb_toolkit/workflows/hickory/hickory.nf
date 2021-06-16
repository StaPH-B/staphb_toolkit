#!/usr/bin/env nextflow

//Description: Workflow for determining an ideal reference genome given a set of fastq files
//Author: Rachael St. Jacques
//email: rachael.stjacques@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
    .set { raw_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(reads) into read_files_fastqc, read_files_trimming

  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
    """
    """
  }
}

//Step1a: FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/logs/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from read_files_fastqc

  output:
  file("*_fastqc.{zip,html}") into fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"
  if(params.savetrimmedreads){
    publishDir "${params.outdir}/trimmed", mode: 'copy'
  }
  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_trimmed{_1,_2}.fastq.gz") into trimmed_reads
  file("${name}.trim.stats.txt") into trimmomatic_stats

  script:
  """
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${task.cpus} ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} 2> ${name}.trim.stats.txt
  mv ${name}*1P.fastq.gz ${name}_trimmed_1.fastq.gz
  mv ${name}*2P.fastq.gz ${name}_trimmed_2.fastq.gz
  """
}
//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"
  publishDir "${params.outdir}/logs/cleanedreads/stats", mode: 'copy',pattern:"*.stats.txt"
  publishDir "${params.outdir}/logs/cleanedreads/reads", mode: 'copy',pattern:"*.fastq.gz"

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads, reads_files_minimap
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats

  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}

//Step3: Assemble trimmed reads with Shovill
process shovill {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/logs/assemblies", mode: 'copy'

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}.fasta") into assembled_genomes_quality, assembled_genomes

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory}  --outdir . --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv contigs.fa ${name}.fasta
  """
}

//Step3a: Assembly Quality Report
process quast {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/logs/quast",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_report

  script:
  """
  quast.py ${assembly} -o .
  mv report.txt ${name}.quast.tsv
  """
}

//Step4: Determine reference genome from assembled raw_reads
process centroid {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  file(assembly) from assembled_genomes.collect()

  output:
  file("*_centroid_ref.fasta") into centroid_out, ref_minimap

  script:
  """
  mkdir assemblies
  mv *.fasta ./assemblies
  centroid.py ./assemblies
  ref=\$(cat centroid_out.txt | awk -F. '{print \$1}')
  ln ./assemblies/\${ref}.fasta ./\${ref}_centroid_ref.fasta
  """
}

//Step5: Map reads to reference genome
process sam_files {
  publishDir "${params.outdir}/logs/sam_files", mode: 'copy'
  tag "$name"

  input:
  set val(name), file(reads) from reads_files_minimap
  file(reference) from ref_minimap


  output:
  tuple name, file("${name}.sam") into sam_percent

  script:
  """
  minimap2 -ax sr *_centroid_ref.fasta ${reads[0]} ${reads[1]} > ${name}.sam
  """
}

//step6: get coverage
process reference_mapping{
  publishDir "${params.outdir}/logs/bam_files", mode: 'copy' ,pattern:'*.sorted.bam'
  publishDir "${params.outdir}/logs/indexed_bam_files", mode: 'copy', pattern:'*.bai'
  publishDir "${params.outdir}/logs/cvg_tsvs", mode: 'copy' ,pattern:'*.tsv'
  tag "$name"

  input:
  set val(name), file(sam_files) from sam_percent

  output:
  file("${name}.samtools.cvg.tsv") into samtools_cvg_tsvs
  file("${name}.sorted.bam")
  file("${name}.sorted.bam.bai")

  script:
  """
  samtools view -S -b ${name}.sam -o ${name}.bam
  samtools sort ${name}.bam > ${name}.sorted.bam
  samtools index ${name}.sorted.bam
  samtools flagstat ${name}.sorted.bam
  samtools coverage ${name}.sorted.bam -o ${name}.samtools.cvg.tsv
  """
}

process generate_report{
  publishDir "${params.outdir}/", mode: 'copy', pattern:"hickory*"
//  tag "$name"

  input:
  file(samtools_coverage) from samtools_cvg_tsvs.collect()

  output:
  file("*hickory_read_map_summary*")

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

#get list of result files
samtools_results = glob.glob("*.samtools.cvg.tsv")
results = {}

# collect samtools results
for file in samtools_results:
    id = file.split(".samtools.cvg.tsv")[0]
    result = result_values(id)
    with open(file,'r') as tsv_file:
        tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
        for line in tsv_reader:
            result.aligned_bases = line["covbases"]
            result.percent_cvg = line["coverage"]
            result.mean_depth = line["meandepth"]
            result.mean_base_q = line["meanbaseq"]
            result.mean_map_q = line["meanmapq"]

    results[id] = result

#create output file
with open(f"hickory_read_map_summary_{today}.csv",'w') as csvout:
    writer = csv.writer(csvout,delimiter=',')
    writer.writerow(["sample","aligned_bases","percent_cvg", "mean_depth", "mean_base_q", "mean_map_q"])
    for id in results:
        result = results[id]
        writer.writerow([result.id,result.aligned_bases,result.percent_cvg,result.mean_depth,result.mean_base_q,result.mean_map_q])
  """
}
