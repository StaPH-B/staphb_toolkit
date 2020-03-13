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

As an example, using the assembler SPAdes would look like this:

```
staphb_toolkit spades -1 sample001_R1.fastq.gz -2 sample001_R2.fastq.gz -o sample001_spades -t 8
```

Workflows that have multiple scripts can be used by calling each script individually. Here is an example calling a few different scripts from Lyve-SET:

```
staphb_toolkit lyveset shuffleSplitReads.pl --numcpus 8 -o interleaved *.fastq.gz
staphb_toolkit lyveset set_manage.pl --create setTest
staphb_toolkit lyveset launch_set.pl --numcpus 8 -ref yourProject/ref/reference.fasta yourProject
```
