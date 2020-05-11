#!/usr/bin/env nextflow

//Description: Workflow for regenerating report
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

params.snp_matrix = ""
params.cg_tree = ""
params.ar_tsv = ""
params.rmd = ""
params.logo = ""


Channel
  .fromPath(params.snp_matrix)
  .set{ snp_mat }

Channel
  .fromPath(params.cg_tree)
  .set{ cgtree }

Channel
  .fromPath(params.rmd)
  .set{ report }

Channel
  .fromPath(params.logo)
  .set{ logo }

if(params.ar_tsv){
  Channel
    .fromPath(params.ar_tsv)
    .set{ ar_tsv }
}else{
  Channel
    .empty()
    .set{ ar_tsv }
}


process render{
  publishDir "${params.outdir}", mode: 'copy'
  stageInMode 'link'

  input:
  file snp from snp_mat
  file tree from cgtree
  file ar from ar_tsv
  file rmd from report
  file dryad_logo from logo

  output:
  file("cluster_report.pdf")
  file(rmd)

  shell:
  """
  Rscript /reports/render.R ${snp} ${tree} ${rmd} ${ar}
  mv report.pdf cluster_report.pdf
  """
}
