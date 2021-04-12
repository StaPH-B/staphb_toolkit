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
cecret can actually handle both read types at the same time, as long as both types are reads are in separate directories. Specify the single end reads in the config file `params.single_reads = <directry with single reads>` and set the `reads_path` to the directory with paired end reads.

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
  - Set `params.kraken2 = true` and `kraken2_db = <path to kraken2 database>`. Instructions are below for the kraken2 and human database.
- Creating a multiple sequencing alignment, SNP matrix, and treefile with mafft, snp-dist, and iqtree
  - Default is `false`
  - If this is desired, set `params.relatedness = true`

#### Downloading the h+v kraken2 database :
Set `params.kraken2 = true` and `params.kraken2_db = 'kraken2_db'`
```
mkdir -p kraken2_db
cd kraken2_db
wget https://storage.googleapis.com/sars-cov-2/kraken2_h%2Bv_20200319.tar.gz
tar -zxf kraken2_h+v_20200319.tar.gz
rm -rf kraken2_h+v_20200319.tar.gz
cd ../
```

This workflow is for the staphB toolkit. A full list of parameters and documentation can be found `cecret`'s standalone repository : [https://github.com/UPHL-BioNGS/Cecret](https://github.com/UPHL-BioNGS/Cecret)

## Questions Worth Asking

### How is `cecret` different than `monroe`?

It's not all that different. [monroe](../monroe) uses minimap2 for mapping/aligning and cleans reads with bbduk and trimmomatic. Running the aligned reads through ivar for primer trimming and consensus creation is the core for both workflows. 

### Why is bcftools set to 'false' by default?

There's nothing wrong with the bcftools process, and the vcf created by bcftools is rather handy for additional analyses. The `'staphb/bcftools:latest'` container is really popular, and has issues downloading during high traffic times. I don't want to have to handle issues of users not understanding why the container did not download. /Sorry

If you want to get the output from bcftools, set `params.bcftools = true` 

### Can I get images of my SNPs and indels?

Yes. Set `params.bamsnap = true`. This is false by default because of how long it takes. It will work with variants called by `ivar` and `bcftools`, although it is **MUCH** faster with the vcf created by bcftools. 

Warning : will not work on all variants. This is due to how bamsnap runs.

### What if I am using an amplicon based library that is not SARS-CoV-2?

In your config file, set your `params.reference_genome`, `params.primer_bed`, and `params.gff_file` appropriately.

You'll also want to set `params.pangolin = false` and `params.nextclade = false`

### This workflow has too many bells and whistles. I really only care about generating a consensus fasta. How do I get rid of all the extras?

Change the parameters in your config file and set most of them to false. 

```
params.fastqc = false
params.ivar_variants = false
params.samtools_stats = false
params.samtools_coverage = false
params.samtools_flagstat = false
params.bedtools = false
params.samtools_ampliconstats = false
params.pangolin = false
params.nextclade = false
```

And, yes, this means I added some bells and whistles so you could turn off the bells and whistles. /irony

### How do I get `cecret` to get my files ready for [GenBank](https://submit.ncbi.nlm.nih.gov/subs/genbank/) or [GISAID](https://www.gisaid.org/) submission?

Create a comma-delimited file with at least the following headers: `Sample_ID`, `Submission_id`, and `Collection_Date` and set `params.sample_file` to that file. Then set `params.rename = true`. Additional headers can be accepted and used for preparing the GenBank submission fastas such as `BioProject`,`SRA`, and `Isolate`. 

Defaults for fasta header values are 
```
country="USA" 
host="Human"
GenBank isolate = SARS-CoV-2/$host/$country/$submission_id/$year
GISAID isolate = hCoV-19/$country/$submission_id/$year
```

Example `sample_file.csv` to go with fastq files `134568-UT-M04492-191220_S44_L001_R{1,2}_001.fastq.gz`, `456578_{1,2}.fastq.gz`, `248074-UT-M05067-200511_S9_L001_001.fastq.gz`
```
Sample_ID,Submission_ID,Collection_Date,Bioproject
134568,UT-UPHL-1902281085,01/12/2021,PRJNA614995
456578,UT-UPHL-1902853304,2020-10-01,PRJNA614995
248074,UT-UPHL-1902032935,missing,PRJNA614995
```

Sometimes sequencing fails, so there are parameters for how many non-ambiguous bases a fasta needs in order to get incorporated into the final file. This can be set with `params.gisaid_threshold` (Default is '`params.gisaid_threshold = 25000`') and `params.genbank_threshold` (Default is '`params.genbank_threshold = 15000`').

### Where do I find the rest of the documentation?

Cecret's standalone workflow repository : [https://github.com/UPHL-BioNGS/Cecret](https://github.com/UPHL-BioNGS/Cecret)
