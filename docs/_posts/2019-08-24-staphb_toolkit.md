---
path: '/:id'
title: 'Running staphb_toolkit'

layout: nil
---
**staphb_toolkit** allows the user access to all of the configured containerized programs hosted on the [StaPH-B Docker Repository](https://github.com/StaPH-B/docker-builds). The usage of the executable is shown below:
```     
usage: staphb_toolkit [optional arguments] <application> [application arguments]

optional arguments:
  -h, --help      show this help message and exit
  --overide_path  Overide the automatic path mounting that is performed for
                  docker containers (Not yet operational)

application:

    spades        SPAdes - St. Petersburg genome assembler
    mash          MASH - Fast genome and metagenome distance estimation using
                  MinHash
    seqyclean     SeqyClean - Pre-process and clean NGS data in order to
                  prepare for downstream analysis
    shovill       Shovill - Faster SPAdes assembler
    prokka        Prokka - Rapid prokaryotic genome annotation
    abricate      Abricate - Mass screening of contigs for antimicrobial and
                  virulence genes
    iqtree        IQ-TREE - A fast and effective stochastic algorithm to infer
                  phylogenetic trees by maximum likelihood.
    lyveset       LYVE-SET - a method of using hqSNPs to create a phylogeny.
    quast         Quast - Genome assembly evaluation tool.
    roary         Roary - Rapid large-scale prokaryote pan genome analysis.
    seqsero       SeqSero - Salmonella serotyping from genome sequencing data.
    serotypefinder
                  SerotypeFinder - identifies the serotype in total or partial
                  sequenced isolates of E. coli.
    unicycler     Unicycler - an assembly pipeline for bacterial genomes.
```

**staphb_toolkit_workflows** allows the user access to the workflows and pipelines that have been developed by members of StaPH-B through a single executable. The usage of the executable is shown below:

```     
usage: staphb_toolkit_workflows [optional arguments] <application> [application arguments]
```
