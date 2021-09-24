# CECRET

A workflow for generating consensus sequences from single or paired-end fastq.gz or fastq reads from amplicon prepared Illumina libraries.

# USAGE

```
usage: staphb-wf [optional arguments] <workflow> [workflow arguments] cecret [--reads_type {paired,single}] [--output <output_path>] [--profile {docker,singularity}] [--config CONFIG] [--get_config] [--resume] [reads_path]
```
## Default Usage
```
staphb-wf cecret Sequencing_reads
```
## Getting a config file for User-supplied goals
```
staphb-wf cecret --get_config
```
## Using single-end instead of paired-end reads
```
staphb-wf cecret --reads_type single Sequencing_reads
```
**Cecret** can actually handle both read types at the same time, as long as both types are reads are in separate directories. Specify the single end reads in the config file `params.single_reads = <directry with single reads>` and set the `reads_path` to the directory with paired end reads.

## Annotating a collection of fastas
```
staphb-wf cecret --annotation fastas
```
Note: set `params.relatedness = true` in order to get a multiple sequence alignment, SNP matrix, and newick file for the collection of fastas.

## Required parameters :
- [primer_bed](./configs/artic_V3_nCoV-2019.bed) bedfile for primer sequences 
  - Default is [artic](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3)'s SARS-CoV-2 V3 primer set
  - Change to user-supplied bedfile with `params.primer_bed`
- [reference_genome](./configs/MN908947.3.fasta) fasta file of genome to align to
  - Default is MN908947.3/SARS-CoV-2
  - Change to user-supplied fasta file with `params.reference_genome`

## Optional (recommended) parameters :
- [gff_file](./configs/MN908947.3.gff) for ivar variants (recommended)
  - Default is MN908947.3/SARS-CoV-2
  - Change to user-supplied gff file with `params.gff_file`
  - If not using a gff file, set `params.ivar_variants = false`
- kraken2 (recommended)
  - Default is `false`
  - Set `params.kraken2 = true` and `kraken2_db = <path to kraken2 database>`
- [amplicon file](./configs/nCoV-2019.insert.bed)
  - Default is the amplicons from artic's V3 primers.
  - Change to user-supplied bedfile with `params.amplicon_bed`
  - If not using, set `params.bedtools_multicov = false`  
- Creating a multiple sequencing alignment, SNP matrix, and treefile with mafft, snp-dist, and iqtree
  - Default is `false`
  - If this is desired, set `params.relatedness = true`

This workflow is also available as a standalone repository, [https://github.com/UPHL-BioNGS/Cecret](https://github.com/UPHL-BioNGS/Cecret), with extended documentation.

# Questions Worth Asking

## How is `cecret` different than `monroe`?

It's not all that different. [monroe](../monroe) uses minimap2 for mapping/aligning and cleans reads with bbduk and trimmomatic. Running the aligned reads through ivar for primer trimming and consensus creation is the core for both workflows. 

## What if I am using an amplicon based library that is not SARS-CoV-2?

Change the following relevant paramters:
* `params.reference_genome`
* `params.primer_bed`
* `params.amplicon_bed` or `params.bedtools_multicov = false`
* `params.gff_file` or `params.ivar_variants = false`
* `params.pangolin = false`
* `params.nextclade = false`
* `params.vadr = false` or create a new vadr container with the appropriate build and adjust the parameters of the vadr process in a [config file](./configs/cecret_config_template.config)
* `params.kraken2_organism = "<organism name>"` or keep `params.kraken2 = false`

## How can I tell if certain amplicons are failing?

There are two ways to do this. 

### With bedtools multicov : 
`cecret/bedtools_multicov` has a file for each sample.
This is standard bedtools multicov output, so it doesn't have a header.

- Column 1 : The reference
- Column 2 : Start of amplicon
- Column 3 : End of amplicon
- Column 4 : Amplicon number
- Column 5-6 : version number and strand from bedfile
- Column 7 : (Column G) is the depth observed for that amplicon for that sample.

### With samtools ampliconstats :
`cecret/samtools_ampliconstats` has a file for each sample
Row number 126 (FDEPTH) has a column for each amplicon (also without a header). To get this row for all of your samples, you can grep the keyword "FDEPTH" from each sample.

```
grep "^FDEPTH" cecret/samtools_ampliconstats/* > samtools_ampliconstats_all.tsv
``` 

## This workflow has too many bells and whistles. I really only care about generating a consensus fasta. How do I get rid of all the extras?

Change the parameters in your config file and set most of them to false. 

```
params.fastqc = false
params.ivar_variants = false
params.samtools_stats = false
params.samtools_coverage = false
params.samtools_flagstat = false
params.samtools_depth = false
params.bedtools_multicov = false
params.samtools_ampliconstats = false
params.samtools_plot_ampliconstats = false
params.pangolin = false
params.nextclade = false
params.vadr = false
```

And, yes, this means I added some bells and whistles so you could turn off the bells and whistles. /irony

## Where do I find the rest of the documentation?

Cecret's standalone workflow repository : [https://github.com/UPHL-BioNGS/Cecret](https://github.com/UPHL-BioNGS/Cecret)
