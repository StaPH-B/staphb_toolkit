---
title: "Using the ToolKit"
layout: page
---

Using the ToolKit is simply done by running the commands `staphb-tk` for running individual tools or `staphb-wf` for running the different workflows incorporated into the toolkit. If you would like more information about the different workflows available in the ToolKit visit the [workflows](/workflows) page.

## Contents
  * [Running Applications](#running-applications-using-the-toolkit)
    - [Pipes and Paths](#special-note-about-autopathing-and-pipes)
  * [Included Applications](#included-applications)
  * [Running Workflows](#using-the-toolkit-to-run-workflows)

<br>
## Running applications using the ToolKit
Running the toolkit using the the `staphb-tk` command provides a menu of options available for running tools in the toolkit:
```
usage: staphb-tk [optional arguments] <application> [application arguments]

optional arguments:
  -h, --help            show this help message and exit
  --docker_config <path>, -c <path>
                        Configuration file for container images and tags; if
                        none provided, configuration will be set to
                        staphb_toolkit/core/docker_config.json
  --list, -l            List all of the software available in the toolkit.

custom program execution:

    mash_species        MASH Species uses a custom database to identify the
                        isolate species.
```

The typical usage of the ToolKit involves a command structure that calls the toolkit i.e. `staphb-tk` followed by the application i.e. `spades` then the parameters associated with that tool. For example using the command `staphb-tk spades` the output shows the options available to the SPAdes assembly tool:

```
SPAdes genome assembler v3.13.0

Usage: /SPAdes-3.13.0-Linux/bin/spades.py [options] -o <output_dir>

Basic options:
-o	<output_dir>	directory to store all the resulting files (required)
--sc			this flag is required for MDA (single-cell) data
--meta			this flag is required for metagenomic sample data
...
```
If we wanted to run the SPAdes assembler on a pair of fastq files the command would be:
```
staphb-tk spades -1 <path to fwd reads> -2 <path to rev reads> -o <output dir>
```

### Special note about autopathing and pipes
The ToolKit will automatically mount paths in your command from your host file system. This allows the toolkit to interact with docker or singularity containers without needing your input on how to mount things. However, if you wish to use a path to a file contained inside the container the autopathing will still try to find that file on your host system therefore you must use an `$` to indicate the path is located inside the container as shown in the command below:

```
staphb-tk mash dist $/db/RefSeqSketchesDefaults.msh mash_output/sample_sketch.msh
```

In addition, pipes are by default read by the bash interpreter. If you wish to use a pipe in your command and you want that pipe to run inside the container, you must use the bash escape character `\` to signify that you want the pipe run in the container. For example:

```
staphb-tk mash dist $/db/RefSeqSketchesDefaults.msh mash_output/sample_sketch.msh \> mash_output/sample_distances.tab
```
<br>
## Included Applications

To list the available software included with the ToolKit use the command `staphb-tk --list`.

**Note**: Some programs have been customized with additional functionality and are available under the **custom program execution** menu of the ToolKit help.
```
Available programs:
Command                  Description
-------                  -----------
abricate                 Abricate - Mass screening of contigs for antimicrobial and virulence genes
bbtools                  BBTools - Suite of fast, multithreaded bioinformatics tools for DNA and RNA sequence data
bwa                      BWA - mapping low-divergent sequences against a large reference genome
canu-racon               Canu-Racon - Ultrafast consensus module for raw de novo assembly of long, uncorrected reads.
cfsan-snp                CFSAN-SNP - SNP calling pipeline from the FDA CFSAN laboratory
circlator                Circlator - A tool to circularize genome assemblies
clustalo                 ClustalO - A fast multiple sequence alignment program
emm-typing-tool          Emm-typing-tool - Group A streptococci emm typing tool for NGS data
fastani                  FastANI - Fast whole-genome sequence average nucleotide identity (ANI) estimation
fastqc                   FastQC - A quality control tool for high throughput sequence data.
fasttree                 FastTree - Infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.
filtlong                 Filtlong - Quality filtering tool for long reads
flye                     Flye - De novo assembler for single molecule sequencing reads using repeat graphs
iqtree                   IQ-TREE - A fast and effective stochastic algorithm to infer phylogenetic trees by maximum likelihood.
kma                      KMA - Mapping method designed to map raw reads directly against redundant databases, in an ultra-fast manner using seed and extend.
kraken                   Kraken - Taxonomic sequence classification system
kraken-build             Kraken-build - Build a kraken database
kraken2                  Kraken2 - The second version of the Kraken taxonomic sequence classification system
kraken2-build            Karken2-build - Build a kraken2 database
ksnp3                    kSNP2 - Identifies the pan-genome SNPs in a set of genome sequences, and estimates phylogenetic trees based upon those SNPs.
legsta                   Legsta - In silico Legionella pneumophila Sequence Based Typing
lyveset                  LYVE-SET - a method of using hqSNPs to create a phylogeny.
mash                     MASH - Fast genome and metagenome distance estimation using MinHash
mashtree                 MashTree - Create a tree using Mash distances
medaka                   Medaka - Sequence correction provided by ONT Research
minimap2                 Minimap2 - a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database.
mlst                     MLST - Scan contig files against PubMLST typing schemes
mugsy                    Mugsy - A multiple whole genome aligner.
multiqc                  MultiQC - Aggregate results from bioinformatics analyses across many samples into a single report.
nanoplot                 NanoPlot - Plotting scripts for long read sequencing data
ncbi-amrfinder-plus      NCBI AMRFinderPlus - Designed to find acquired antimicrobial resistance genes and some point mutations in protein or assembled nucleotide sequences.
orthofinder              OrthoFinder - Phylogenetic orthology inference for comparative genomics
pilon                    Pilon - Automated genome assembly improvement and variant detection tool
plasmidseeker            PlasmidSeeker - A k-mer based program for the identification of known plasmids from whole-genome sequencing reads
prokka                   Prokka - Rapid prokaryotic genome annotation
quast                    Quast - Genome assembly evaluation tool.
rasusa                   RASUA - Randomly subsample sequencing reads to a specified coverage
raxml                    RAxML -Maximum likelihood tree builder.
roary                    Roary - Rapid large-scale prokaryote pan genome analysis.
salmid                   SalmID - Rapid confirmation of Salmonella spp. and subspp. from sequence data
samtools                 Samtools - A suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories.
seqsero                  SeqSero - Salmonella serotyping from genome sequencing data.
seqsero2                 SeqSero2 - Salmonella serotype prediction from genome sequencing data.
seqyclean                SeqyClean - Pre-process and clean NGS data in order to prepare for downstream analysis
seroba                   Seroba - k-mer based Pipeline to identify the Serotype from Illumina NGS reads
serotypefinder           SerotypeFinder - identifies the serotype in total or partial sequenced isolates of E. coli.
shovill                  Shovill - Faster SPAdes assembler
sistr                    SISTR - Salmonella in silico typing resource command-line tool
skesa                    SKESA - NCBI's de novo genome assemlber
snippy                   Snippy - Rapid haploid variant calling and core genome alignment
snp-dists                SNP-dists - Pairwise SNP distance matrix from a FASTA sequence alignment
spades                   SPAdes - St. Petersburg genome assembler
sra-toolkit              SRA ToolKit - Collection of tools and libraries for using data in the INSDC Sequence Read Archives.
staramr                  StarAMR - Scans genome contigs against the ResFinder, PlasmidFinder, and PointFinder databases.
tiptoft                  TipToft - Predict plasmids from uncorrected long read data
trimmomatic              Trimmoamtic - Flexible read trimming tool for Illumina NGS data
unicycler                Unicycler - an assembly pipeline for bacterial genomes.
wtdbg2                   WTDBG2 - Fuzzy Bruijn graph approach to long noisy reads assembly
```
<br>
## Using the ToolKit to run workflows
The ToolKit also provides the ability to run workflows using the `staphb-wf` command. Information and usage of specific workflows is available on the [workflow page](/workflows).
```
usage: staphb-wf [optional arguments] <workflow> [workflow arguments]

optional arguments:
  -h, --help  show this help message and exit

workflows:

    tredegar  Quality control of WGS read data.
    monroe    Consensus assembly for SARS-CoV-2 from ARTIC + Illumina
              protocols.
    dryad     A comprehensive tree building program.
```

To run the workflow simply use the command and workflow allowing with the workflow specific commands.
