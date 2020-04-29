---
title: 'Monroe v1.0.0'
layout: page
---
SARS-CoV-2 genome assembly and sample cluster detection for read data generated using the [ARTIC PCR tiling protocols](https://artic.network/ncov-2019) and either  Illumina-platform (e.g. MiSeq) or Oxford Nanopore Technologies sequencing.

## Data workflow:
Monroe consists of three separate NextFlow pipelines for Illumina paired-end read assembly (`pe_assebly`) Oxford Nanopore Technlogies read assembly (`ont_assembly`) cluster analysis from assembled SC2 genomes (`cluster_analysis`):
![Monroe pipeline](/assets/workflows/monroe/Monroe_v1.0.png)

## Paired-End Read Assembly:
Monroe `pe_assembly` Uses [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [BBDuk](http://seqanswers.com/forums/showthread.php?t=42776) to perform read trimming and adapter removal  and Illumina read data prior to mapping read data to a reference SARS-CoV-2 genome (Wuhan-1; NCBI RefSeq [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)). Paired-fastq files are pulled from the alignment file with [SAMtools](http://www.htslib.org/doc/samtools.html)--these filtered read data (i.e. paired reads that map to the reference genome) can be uploaded to public repositories.

ARTIC primers are trimmed from the alignment file and a consensus assembly is generated using [iVar v1.2.1](https://github.com/andersen-lab/ivar).


### Quick Start

````
$ staphb-wf monroe pe_assembly <input_dir> -o <output_dir> --primers <ARTIC_primer_version>
````
- `<input_dir>` is a positional argument that designates the path to an input directory containing paired-end read data (fastq). <br />
- `-o` specifies the directory to which Monroe will write all results; if an `<output_dir>`is not provided, results will be written to a `monroe_run_<date>` directory. <br />
- `--primers` The ARTIC primer set version used to generate the read data must be designated as `V1`, `V2`, or `V3`.

### Other Options
- `--profile`: Nextflow profile, either Docker or Singularity. Default will try docker first, then singularity if the docker executable cannot be found.
- `--config`, `-c`: Path to a custom Nextflow configureation
- `--resume`: Resume a previous run
- `--get_config`: Create a template config file for pipeline customization

### Output:
Monroe `pe_assembly` will organize all output into four subdirectories under the specified `<output_dir>`:
- `<output_dir>/alignments`: Sorted BAM files after minimap2 mapping and ivar primer trimming
- `<output_dir>/assemblies`: Consensus genome assemblies in fasta format as well as a `<date>_quality_metrics.tsv` file comprising of quality metrics for all genome assemblies
- `<output_dir>/SC2_reads`: Paired read data that have mapped to the Wu-Han-1 reference genome
- `<output_dir>/logs`: NextFlow execution report (`Monroe_execution_report.html`), trace file (`Monroe_trace.txt`), and task work directories.

Sample quality metrics file:

![Monroe pipeline](/assets/workflows/monroe/quality_metrics.png)
- sample: isolate ID pulled from the fastq file
- aligned_bases: number of bases mapped to the SARS-CoV-2 reference genome
- percent_cvg: percent of reference genome with mapped read data
- mean_depth: mean depth-of-coverage
- mean_base_q: average quality of basecalls for read data mapped to the reference genome
- mean_map_q: mean mapping quality
- status: "PASS" if percent_cvg >80%, mean_base_q >30, and mean_map_q >30; "WARNING" if any of these quality thresholds are not met


### Docker Images
The base NextFlow configuration profiles (Docker, Singularity) for Monroe `pe_assebly` incorporate the following StaPH-B Docker Images:

| Monroe `pe_assembly` Process   | Function  | Docker Image  | Comment|
|---|---|---|---|---|
| preProcess  | Renames input read files for downstream processing | staphb/fastqc_container  | Light-weight container for quick text processing  |
| trim  | Quality trimming of input read data  | staphb/trimmomatic:0.39  | |
| cleanreads  | Adapter removal of input read data  | staphb/bbtools:38.76  | |
| ivar  | Read mapping and consensus genome assembly  | staphb/ivar:1.2.1-SC2  |  |
| samtools  | Gathering alignment quality metrics  | staphb/samtools:1.10  | |
| assembly_results  | Curating assembly quality metrics  | staphb/tiptoft:1.0.0 | Light-weight container with python3 |

## Oxford Nanopore Technlogies  Read Assembly:

guppy
guppy artic minion artic medaka

### Quick Start
````
$ staphb-wf monroe ont_assembly <input_dir> <sequencing_summary> -o <output_dir> --primers <ARTIC_primer_version>
````
- `<input_dir>` is a positional argument that designates the path to an input directory containing Oxford Nanopore Technologies read read data; input data can be either in either Fast5 or FastQ format--if Fast5, the `--ont_basecalling` must also be invoked
- `<sequencing_summary>` is a positional argument that designates the path to the location of the sequencing summary
- `-o` specifies the directory to which Monroe will write all results; if an `<output_dir>`is not provided, results will be written to a `monroe_run_<date>` directory. <br />
- `--primers` The ARTIC primer set version used to generate the read data must be designated as `V1`, `V2`, or `V3`.

### Other Options
- `--run_prefix`: Desired run prefix. Default = `artic_ncov19`
- `--ont_basecalling`: perform high accuracy basecalling using GPU (only use if you have setup a GPU compatible device; must be invoked if input data is in Fast5 format.
- `--profile`: Nextflow profile, either Docker or Singularity. Default will try docker first, then singularity if the docker executable cannot be found.
- `--config`, `-c`: Path to a custom Nextflow configureation
- `--resume`: Resume a previous run
- `--get_config`: Create a template config file for pipeline customization

### Output:


## Cluster Analysis:
Monroe `cluster_analysis` uses [Mafft](https://mafft.cbrc.jp/alignment/software/) to perform multiple-sequence alignment of all SARS-CoV-2 genomes provided. The resulting alignment fasta file is used to generate a pairwise-snp distance matrix with [snp-dists](https://github.com/tseemann/snp-dists) and a maximum-likelihood phylogeneitc tree with [IQ-Tree](http://www.iqtree.org/).

Output from snp-dists and IQ-Tree are curated into a single pdf report using the [StaPH-B cluster-report-env](https://hub.docker.com/r/staphb/cluster-report-env)

### Quick Start
````
$ staphb-wf monroe cluster_analysis <input_dir> -o <output_dir>
````
- `<input_dir>` is a positional argument that designates the path to an input directory containing the SARS-CoV-2 assembly files (fasta). <br />
- `-o` specifies the directory to which Monroe will write all results; if an `<output_dir>`is not provided, results will be written to a `monroe_run_<date>` directory. <br />

### Other Options
- `--profile`: Nextflow profile, either Docker or Singularity. Default will try docker first, then singularity if the docker executable cannot be found.
- `--config`, `-c`: Path to a custom Nextflow configureation
- `--resume`: Resume a previous run
- `--get_config`: Create a template config file for pipeline customization

### Output:
Monroe `cluster_analysis` will write the final pdf report to the specified `<output_dir>`. All other output will be organized into three subdirectories:
- `<output_dir>/images`: PNG files of the maximum-likelihood tree and color-coded SNP-distance matrix
- `<output_dir>/msa`: Mafft alignment file (fasta), IQ-Tree newick file, and SNP-dists pairwise-distance matrix
- `<output_dir>/logs`: NextFlow execution report (`Monroe_execution_report.html`), trace file (`Monroe_trace.txt`), and task work directories.


### Docker Images
The base NextFlow configuration profiles (Docker, Singularity) for Monroe `cluster_analysis` incorporate the following StaPH-B Docker Images:

| Monroe `cluster_analysis` Process   | Function  | Docker Image  | Comment|
|---|---|---|---|---|
| msa  | Performing multi-sequence alignment with Mafft | staphb/mafft:7.450  |  |
| snp_matrix  | Generating pairwise snp-distance matrix from Mafft msa  | staphb/snp-dists:0.6.2  | |
| iqtree  | Generating maximum-likelihood phylogenetic tree from Mafft msa  | staphb/iqtree:1.6.7 | |
| render  | Curating all output into a single pdf report   | staphb/cluster-report-env:1.0  |  |

## Version History

<b>Current version: v1.0.0 April 29, 2020</b>

Version 1.0.0 is the first stable version of Monroe

## Author
[Kevin G. Libuit](https://github.com/kevinlibuit), DCLS Bioinformatics Lead Scientist <br />
[Kelsey R Florek](https://github.com/k-florek), WSLH Bioinformatics Scientist  <br />
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatics Fellow
