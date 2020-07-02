---
title: "StaPH-B ToolKit Workflows"
layout: page
---

Workflows in the ToolKit are run using [Nextflow](https://www.nextflow.io/). The ToolKit incorporates the Nextflow executable and runs the provided workflow internally. Nextflow only requires Java version 8 or later and is extremely flexible across different compute environments. By default all workflows are designed to work with either Docker or Singularity locally or in the cloud using AWS. If you would like to add an additional environment contact the tool's creator.

## Contents
  * [Cutshaw](#cutshaw)
  * [Dryad](#dryad)
  * [Foushee](#foushee)
  * [Monroe](#monroe)
  * [Tredegar](#tredegar)


## Cutshaw
Read the usage guide [here](/staphb_toolkit/workflow_docs/cutshaw).
## Dryad
Dryad is a pipeline to detect antimicrobial resistance genes, construct reference free core-genome, and/or SNP phylogenetic trees for examining prokaryote relatedness in outbreaks. Dryad will perform both a reference free core-genome analysis based off of the approach outlined by [Oakeson et. al](https://www.ncbi.nlm.nih.gov/pubmed/30158193) and/or a SNP analysis using the [CFSAN-SNP](https://snp-pipeline.readthedocs.io/en/latest/readme.html) pipeline.
Read the usage guide [here](/staphb_toolkit/workflow_docs/dryad).

## Foushee
Read the usage guide [here](/staphb_toolkit/workflow_docs/foushee).
## Monroe
Read the usage guide [here](/staphb_toolkit/workflow_docs/monroe).
## Tredegar
Read the usage guide [here](/staphb_toolkit/workflow_docs/tredegar).
