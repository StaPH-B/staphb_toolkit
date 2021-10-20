---
title: "StaPH-B ToolKit Workflows"
layout: page
---

Workflows in the ToolKit are run using [Nextflow](https://www.nextflow.io/). The ToolKit incorporates the Nextflow executable and runs the provided workflow internally. Nextflow only requires Java version 8 or later and is extremely flexible across different compute environments. By default all workflows are designed to work with either Docker or Singularity locally or in the cloud using AWS. If you would like to add an additional environment contact the tool's creator.

## Workflows
  * [Cutshaw](/staphb_toolkit/workflow_docs/cutshaw) - Bioinformatics Pipeline for Instrument Validation and Assessment Technical Proficiency in WGS Protocols.
  * [Dryad](/staphb_toolkit/workflow_docs/dryad) - A pipeline to construct reference free core-genome or SNP phylogenetic trees for examining prokaryote relatedness in outbreaks.
  * [Foushee](/staphb_toolkit/workflow_docs/foushee) - Bioinformatics pipeline for reference-free SNP analysis of Group-A *Streptococcus* (GAS) isolates.
  * [Tredegar](/staphb_toolkit/workflow_docs/tredegar) - Bioinformatics pipeline for infectious disease WGS data QC.
  * [Cecret](/staphb_toolkit/workflow_docs/cecret) - Workflow for generating SARS-CoV-2 consensus sequences from single or paired-end fastq.gz or fastq reads from amplicon prepared Illumina libraries.
  * [Monroe](/staphb_toolkit/workflow_docs/monroe) - Bioinformatics pipeline for SARS-CoV-2 genome assembly and sample cluster detection. (Generally Cecret is the preferred workflow)
