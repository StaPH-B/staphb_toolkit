---
title: 'Monroe v1.0.0'
layout: page
---
Bioinformatics pipeline for SARS-CoV-2 genome assembly and sample cluster detection.
- SARS-CoV-2 genome assembly can be performed from read data generated using the [ARTIC PCR tiling protocols](https://artic.network/ncov-2019) (V1, V2, or V3) with either an Illumina sequencing platform (e.g. MiSeq) or an Oxford Nanopore Technologies MinIon device
- Cluster detection can be performed from input assembly files (fasta) generated from any sequencing protocol

## Data workflow:
Monroe consists of three separate NextFlow pipelines for Illumina paired-end read assembly (`pe_assebly`) Oxford Nanopore Technlogies read assembly (`ont_assembly`) cluster analysis from assembled SC2 genomes (`cluster_analysis`):
![Monroe pipeline](/assets/workflows/monroe/Monroe_v1.0.png)

Monroe's NextFlow pipelines can be executed using the following command format:
```
$ staphbe-wf monroe <monroe_pipeline> [options]
```
- `<monroe_pipeline>` must be `pe_assembly`, `ont_assembly`, or `cluster_analysis`

## Paired-End Read Assembly:
Monroe `pe_assembly` uses [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [BBDuk](http://seqanswers.com/forums/showthread.php?t=42776) to perform read trimming and adapter/PhiX removal prior to mapping read data to a reference SARS-CoV-2 genome (Wuhan-1; NCBI RefSeq [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)) with [minimap2](https://www.ncbi.nlm.nih.gov/pubmed/29750242). Paired-fastq files are pulled from the alignment file with [SAMtools](http://www.htslib.org/doc/samtools.html)--these filtered read data (i.e. paired reads that map to the reference genome) are stored as separate files that can be uploaded to public repositories such as NCBI SRA.

The minimap2 alignment file is also used to generate a consensus assembly after trimming ARTIC primers using [iVar v1.2.1](https://github.com/andersen-lab/ivar).

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

| Process   | Function  | Docker Image  | Comment|
|---|---|---|---|---|
| preProcess  | Renames input read files for downstream processing | staphb/fastqc_container  | Light-weight container for quick text processing  |
| trim  | Quality trimming of input read data with bbduk  | staphb/trimmomatic:0.39  | `trimmomatic` parameters set to: minlength=75, windowsize=4, & qualitytrimscore=30|
| cleanreads  | Adapter and PhiX removal from input read data  | staphb/bbtools:38.76  |`bbduk` default parameters used |
| ivar  | Read mapping and consensus genome assembly  | staphb/ivar:1.2.1-SC2  | `ivar consensus` parameters set to: minimum frequency (-t)=0, minimum depth (-m)=1 |
| samtools  | Gathering alignment quality metrics  | staphb/samtools:1.10  | `samtools coverage` default parameters used|
| assembly_results  | Curating assembly quality metrics  | staphb/tiptoft:1.0.0 | Light-weight container with python3 |

Default docker images and parameters listed above can be adjusted by:
1. Copying the template `pe_assebly` config file (`$ staphb-wf monroe pe_assembly --get_config`)
2. Using a text editor to change the `<date>_pe_assembly.config` file
3. Specifying your custom config file (i.e. the edited `<date>_pe_assembly.config>` file) when running the pipeline:<br />

```
$ staphb-wf monroe pe_assembly <input_dir> -o <output_dir> --primers <ARTIC_primer_version> -c <custom_config_file> [options]
```

## Oxford Nanopore Technlogies (ONT) Read Assembly:

Monroe `ont_assembly` can accept ONT Fast5 or FastQ read data. If Fast5 files are provided, high accuracy basecalling will be performed using a GPU-optimized environment. ONT FastQ files, in accordance to the [ARTIC bioinformatics protocols](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html), undergo demultiplexing and read filtering prior to genome assembly with either [NanoPolish](https://nanopolish.readthedocs.io/en/latest/) or [Medaka](https://github.com/nanoporetech/medaka).


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
- `--ont_basecalling`: perform high accuracy basecalling using GPU (only use if you have setup a GPU compatible device); must be invoked if input data is in Fast5 format.
- `--profile`: Nextflow profile, either Docker or Singularity. Default will try docker first, then singularity if the docker executable cannot be found.
- `--config`, `-c`: Path to a custom Nextflow configureation
- `--resume`: Resume a previous run
- `--get_config`: Create a template config file for pipeline customization

### Output:

### Docker Images
The base NextFlow configuration profiles (Docker, Singularity) for Monroe `cluster_analysis` incorporate the following StaPH-B Docker Images:

| Process   | Function  | Docker Image  | Comments|
|---|---|---|---|---|
| guppy_basecalling  | Performing high accuracy basecalling using GPU | genomicpariscentre/guppy-gpu| Used if `--ont_basecalling` is invoked; must have setup a GPU compatible device|
| guppy_demultiplexing  | Demultiplexing samples by barcodes present  | genomicpariscentre/guppy  | `guppy_barcoder` default parameters used |
| artic_guppyplex  | Filtering read data by read-length to remove chimeric reads  |genomicpariscentre/guppy | `artic guppyplex` parameters set to: --min-length=400, --max-length=700  |
| artic_nanopolish_pipeline  | Performing genome assembly with NanoPolish   | staphb/artic-ncov2019-nanopolish  | `artic minion` parameters set to: --normalize=200 |
| artic_medaka_pipeline  | Performing genome assembly with Medaka   | staphb/artic-ncov2019-nanopolish  | `artic medaka` parameters set to: --normalize=200 |

Default docker images and parameters listed above can be adjusted by:
1. Copying the template `ont_assembly` config file (`$ staphb-wf monroe ont_assembly --get_config`)
2. Using a text editor to change the `<date>_ont_assembly.config` file
3. Specifying your custom config file (i.e. the edited `<date>_ont_assembly.config>` file) when running the pipeline: < br/>

```
$ staphb-wf monroe ont_assembly <input_dir> -o <output_dir> --primers <ARTIC_primer_version> -c <custom_config_file> [options]
```

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

| Process   | Function  | Docker Image  | Comments |
|---|---|---|---|---|
| msa  | Performing multi-sequence alignment with Mafft | staphb/mafft:7.450  | `mafft` default parameters used|
| snp_matrix  | Generating pairwise snp-distance matrix from Mafft msa  | staphb/snp-dists:0.6.2  | `snp-dists` default parameters used|
| iqtree  | Generating maximum-likelihood phylogenetic tree from Mafft msa  | staphb/iqtree:1.6.7 | `iqtree` set to: substitution model (-m)=GTR+4, bootstrap replicates (-bb)=1000|
| render  | Curating all output into a single pdf report   | staphb/cluster-report-env:1.0  |   |

Default docker images and parameters listed above can be adjusted by:
1. Copying the template `cluster_analysis` config file (`$ staphb-wf monroe cluster_analysis --get_config`)
2. Using a text editor to change the `<date>_cluster_analysis.config` file
3. Specifying your custom config file (i.e. the edited `<date>_cluster_analysis.config>` file) when running the pipeline:<br />

```
staphb-wf monroe cluster_analysis <input_dir> -o <output_dir> -c <custom_config_file> [options]
```

## Version History

<b>Current version: v1.0.0 April 29, 2020</b>

Version 1.0.0 is the first stable version of Monroe

## Authors
[Kevin G. Libuit](https://github.com/kevinlibuit), DCLS Bioinformatics Lead Scientist <br />
[Kelsey R Florek](https://github.com/k-florek), WSLH Bioinformatics Scientist  <br />
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatics Fellow
