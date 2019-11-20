---
category: Workflows
path: '/Workflows/:id'
title: 'Foushee'

layout: nil
---

# Foushee v1.0
Bioinformatics pipeline for reference-free SNP analysis of Group-A *Streptococcus* (GAS) isolates

## Data workflow:
![Foushee pipeline](../assets/Foushee_v1.0.png)

## Reference-free SNP analysis for GAS isolates of same emm-type
Read quality assessment, genome assembly, and taxonomic predictions (including GAS emm-type) are performed for all input data using [Tredegar](https://staph-b.github.io/staphb_toolkit/#/tredegar-README), a QC workflow within the StaPH-B tookit. Tredegar results are used to group isolates by emm-type before performing reference-free Single-Nucleotide Polymorphism (SNP) analysis on each emm-type group independently using [kSNP3](https://www.ncbi.nlm.nih.gov/pubmed/25913206). The resulting core_SNPs_matrix.fasta file is used to generate a pairwise SNP distance matrix using [snp-dists](https://github.com/tseemann/snp-dists). 