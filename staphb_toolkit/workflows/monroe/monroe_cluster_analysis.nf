#!/usr/bin/env nextflow

//Description: Perform multiple sequence alinmment of SC2 assemblies to generate a SNP-distance matrix & ML phylogenetic tree
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.assemblies = ""
params.outdir = ""
params.report = ""
params.pipe = ""

//setup channel to read in and pair the fastq files
Channel
    .fromPath("${params.assemblies}/*.fasta")
    .ifEmpty { exit 1, "Cannot find any fasta assemblies in ${params.assemblies}" }
    .set { assemblies}

Channel
    .fromPath(params.report)
    .set { report }


// Generate msa from consensus sequences using MAFFT
process msa{
  publishDir "${params.outdir}/msa",mode:'copy',overwrite: false

  input:
  file(assemblies) from assemblies.collect()

  output:
  file "*msa.fasta" into msa_snp, msa_tree

  shell:
  """
  date=\$(date '+%m%d%y')
  cat *.fasta > assemblies.fasta
  mafft --thread -1 assemblies.fasta > \${date}_msa.fasta
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
  publishDir "${params.outdir}/msa",mode:'copy', overwrite: false

  input:
  file(msa) from msa_tree

  output:
  file("*msa.tree") into ML_tree

  script:
    """
    date=\$(date '+%m%d%y')
    numGenomes=`grep -o '>' ${msa} | wc -l`
    if [ \$numGenomes -gt 3 ]
    then
      iqtree -nt AUTO -s ${msa} -m ${params.iqtree_model} -bb ${params.iqtree_boostraps}
      mv \${date}_msa.fasta.contree \${date}_msa.tree
    fi
    """
}

process render{
  publishDir "${params.outdir}/", mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/images", mode: 'copy', pattern: "*.png"
  publishDir "${params.outdir}/msa", mode: 'copy', pattern: "*snp_distance_matrix.tsv", overwrite: false
  echo true

  input:
  file(pairwise_snp_distance_matrix) from matrix
  file(ml_tree) from ML_tree
  file(rmd) from report

  output:
  file "*monroe_cluster_report.pdf"
  file "*ML_tree.png"
  file "*SNP_heatmap.png"
  file "*snp_distance_matrix.tsv"
  shell:
"""
date=\$(date '+%m%d%y')
cp ${rmd} ./report_template.Rmd
Rscript /reports/render.R ${pairwise_snp_distance_matrix} ${ml_tree} ./report_template.Rmd
mv report.pdf \${date}_monroe_cluster_report.pdf
mv ML_tree.png \${date}_ML_tree.png
mv SNP_heatmap.png \${date}_SNP_heatmap.png
mv snp_distance_matrix.tsv \${date}_snp_distance_matrix.tsv
"""
}
