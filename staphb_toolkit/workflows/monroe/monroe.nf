#!/usr/bin/env nextflow

//Description: Workflow for quality control of raw illumina reads
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""
params.primers =""
params.report = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs(  "${params.reads}/*{R1,R2,_1,_2}*.fastq.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

Channel
    .fromPath(params.report)
    .set { report }
    
//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file("*{R1,R2,_1,_2}.fastq.gz") into raw_reads_clean

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

//Trim reads and remove adapters
process seqyclean {
  tag "$name"

  input:
  set val(name), file(reads) from raw_reads_clean

  output:
  tuple val("${name}"), "${name}_clean_PE{1,2}.fastq" into cleaned_reads

  script:
"""
  seqyclean -1 ${name}_R1.fastq.gz -2 ${name}_R2.fastq.gz -minlen ${params.min_read_length} -o ${name}_clean -c ${params.contaminants} ${params.quality_trimming}
"""
}


//Assemble cleaned reads with iVar
process ivar {
  publishDir "${params.outdir}/consensus_assemblies", mode: 'copy',pattern:"*_consensus.fasta"
  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.sorted.bam"

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}_consensus.fasta") into assembled_genomes, assembled_genomes_msa
  tuple name, file("${name}.sorted.bam") into alignment_file

  shell:
"""
minimap2 -K 20M -x sr -a /reference/nCoV-2019.reference.fasta !{reads[0]} !{reads[1]} | samtools view -u -h -F 4 - | samtools sort > SC2.bam
samtools index SC2.bam
samtools flagstat SC2.bam
ivar trim -i SC2.bam -b /reference/ARTIC-${params.primers}.bed -p ivar -e

samtools sort  ivar.bam > ${name}.sorted.bam
samtools index ${name}.sorted.bam
samtools flagstat ${name}.sorted.bam

samtools mpileup -f /reference/nCoV-2019.reference.fasta -d 1000000 -A -B -Q 0 ${name}.sorted.bam | ivar consensus -p ivar -m 1 -t 0 -n N
echo '>${name}' > ${name}_consensus.fasta

seqtk seq -U -l 50 ivar.fa | tail -n +2 >> ${name}_consensus.fasta
"""
}

//QC of read data
process samtools {
  publishDir "${params.outdir}/alignments",mode:'copy',overwrite: false

  input:
  set val(name), file(alignment) from alignment_file

  output:
  file "${name}_samtoolscoverage.tsv" into alignment_qc

  shell:
  """
  samtools coverage ${alignment} -o ${name}_samtoolscoverage.tsv
  """
}

//Collect and format report
process alignment_results{
  publishDir "${params.outdir}", mode: 'copy'
  echo true

  input:
  file(cg_pipeline_results) from alignment_qc.collect()

  output:
  file "consensus_statstics.csv"

  script:
"""
#!/usr/bin/env python3
import os, sys
import glob, csv
import xml.etree.ElementTree as ET
class result_values:
    def __init__(self,id):
        self.id = id
        self.aligned_bases = "NA"
        self.percent_cvg = "NA"
        self.mean_depth = "NA"
        self.mean_base_q = "NA"
        self.mean_map_q = "NA"


#get list of result files
samtools_results = glob.glob("*_samtoolscoverage.tsv")

results = {}

# collect samtools results
for file in samtools_results:
    id = file.split("_samtoolscoverage.tsv")[0]
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
with open("consensus_statstics.csv",'w') as csvout:
    writer = csv.writer(csvout,delimiter=',')
    writer.writerow(["sample","aligned_bases","percent_cvg", "mean_depth", "mean_base_q", "mean_map_q"])
    for id in results:
        result = results[id]
        writer.writerow([result.id,result.aligned_bases,result.percent_cvg,result.mean_depth,result.mean_base_q,result.mean_map_q])
"""

}

// Generate msa from consensus sequences using MAFFT
process msa{
  publishDir "${params.outdir}/msa", mode: 'copy'
  echo true

  input:
  file(assemblies) from assembled_genomes_msa.collect()

  output:
  file "msa.fasta" into msa_snp, msa_vcf, msa_tree

  shell:
  """
  cat *.fasta > assemblies.fasta
  mafft --globalpair --maxiterate 1000 assemblies.fasta > msa.fasta
  """
}

// Generate SNP matrix from MAFFT alignment
process snp_matrix{
  publishDir "${params.outdir}", mode: 'copy'
  echo true

  input:
  file(alignment) from msa_snp

  output:
  file "pairwise_snp_distance_matrix.tsv"

  shell:
  """
  snp-dists ${alignment} > pairwise_snp_distance_matrix.tsv
  """
}

// Generate multi-sample vcf from MAFFT alignment
process vcf{
  publishDir "${params.outdir}", mode: 'copy'
  echo true

  input:
  file(msa) from msa_vcf

  output:
  file "msa.vcf" into vcf

  shell:
  """
  snp-sites -v ${msa} -o msa.vcf
  """
}

//Infer ML tree from MAFFT alignment
process iqtree {
  publishDir "${params.outdir}",mode:'copy'

  input:
  file("msa.fasta") from msa_tree

  output:
  file("msa.tree") optional true

  script:
    """
    numGenomes=`grep -o '>' msa.fasta | wc -l`
    if [ \$numGenomes -gt 3 ]
    then
      iqtree -nt AUTO -s msa.fasta -m 'GTR+G4' -bb 1000
      mv msa.fasta.contree msa.tree
    fi
    """
}

process snp_frequency{
  publishDir "${params.outdir}", mode: 'copy'
  echo true

  input:
  file(vcf) from vcf

  output:
  file "snp_frequencies.tsv"

  script:
"""
#!/usr/bin/env python3
import csv

with open('msa.vcf','r') as vcf:
    with open('snp_frequencies.tsv','w') as tsv:
        writer = csv.writer(tsv, delimiter='\t')
        writer.writerow(['Genomic Position','Reference Allele','Alternative Allele','Frequency'])
        for line in vcf:
          if '#' in line:
              next
          else:
              sline = line.strip().split("\t")
              pos = sline[1]
              ref = sline[3]
              alt = sline[4]
              vals = sline[9:len(sline)]
              vals = list(map(str, vals))
              vals = [sub.replace("2", "1") for sub in vals]
              vals = [sub.replace("3", "1") for sub in vals]
              vals = list(map(int, vals))
              freq = sum(vals)
              writer.writerow([pos,ref,alt,freq])
"""
}

process render{
  publishDir "${params.outdir}/report", mode: 'copy'
  echo true

  input:
  file("pairwise_snp_distance_matrix.tsv") from matrix
  file("msa.tree") from outChannel
  file(rmd) from report

  output:
  file "monroe_report.pdf"                                                                                                                                                                                         
  shell:
"""
Rscript /reports/render.R pairwise_snp_distance_matrix.tsv msa.tree ${rmd}
}     
