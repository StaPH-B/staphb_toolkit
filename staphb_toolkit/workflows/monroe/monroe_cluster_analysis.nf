#!/usr/bin/env nextflow

//Description: Workflow for quality control of raw illumina reads
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.assemblies = ""
params.outdir = ""
params.report = ""

//setup channel to read in and pair the fastq files
Channel
    .fromPath("${params.assemblies}/*.fasta")
    .ifEmpty { exit 1, "Cannot find any assemblies matching: ${params.reads}" }
    .set { assemblies}

Channel
    .fromPath(params.report)
    .set { report }


// Generate msa from consensus sequences using MAFFT
process msa{
  publishDir "${params.outdir}/cluster_analysis/",mode:'copy',overwrite: false

  input:
  file(assemblies) from assemblies.collect()

  output:
  file "*msa.fasta" into msa_snp, msa_tree

  shell:
  """
  cat *.fasta > assemblies.fasta
  mafft --thread -1 assemblies.fasta > msa.fasta
  """
}

// Generate SNP matrix from MAFFT alignment
process snp_matrix{
  echo true

  input:
  file(alignment) from msa_snp

  output:
  file "*pairwise_snp_distance_matrix.tsv" into matrix

  shell:
  """
  snp-dists ${alignment} > pairwise_snp_distance_matrix.tsv
  """
}

//Infer ML tree from MAFFT alignment
process iqtree {
  publishDir "${params.outdir}/cluster_analysis/",mode:'copy', overwrite: false

  input:
  file(msa) from msa_tree

  output:
  file("*msa.tree") into ML_tree

  script:
    """
    numGenomes=`grep -o '>' ${msa} | wc -l`
    if [ \$numGenomes -gt 3 ]
    then
      iqtree -nt AUTO -s ${msa} -m 'GTR+G4' -bb 1000
      mv msa.fasta.contree msa.tree
    fi
    """
}

process render{
  publishDir "${params.outdir}/cluster_analysis/report", mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/cluster_analysis/images", mode: 'copy', pattern: "*.png"
  publishDir "${params.outdir}/cluster_analysis/msa", mode: 'copy', pattern: "*snp_distance_matrix.tsv", overwrite: false
  echo true

  input:
  file(pairwise_snp_distance_matrix) from matrix
  file("msa.tree") from ML_tree
  file(rmd) from report

  output:
  file "monroe_cluster_report.pdf"
  file "ML_tree.png"
  file "SNP_heatmap.png"
  file "snp_distance_matrix.tsv"
  shell:
"""
date=\$(date '+%m%d%y')
cp ${rmd} ./report_template.Rmd
Rscript /reports/render.R ${pairwise_snp_distance_matrix} msa.tree ./report_template.Rmd
mv report.pdf monroe_cluster_report.pdf
"""
}
