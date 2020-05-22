---
title: 'Tredegar'
layout: page
---

# Tredegar v2.1.0
Bioinformatics pipeline for infectious disease WGS data QC

## Data workflow:
![Tredegar pipeline](/staphb_toolkit/assets/workflows/tredegar/Tredegar_v2.1.png)

## Taxonomic Prediction:
Tredegar uses [Mash](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x), a fast genome distance estimation algorithm, to predict the taxonomic identity of all isolates by querying input read data aginst the entire RefSeq databse. This method allows for accurate taxonomic identification up to the species level.
For *Salmonella enterica*, *Escherichia coli*, and Group A *Streptoccus* isolates, subspecies characterization is determined using [SeqSero](http://jcm.asm.org/content/early/2015/03/05/JCM.00323-15), [SerotypeFinder](http://jcm.asm.org/content/53/8/2410.full.pdf+html), and [the emm typer tool from Public Health England (PHE)](https://github.com/phe-bioinformatics/emm-typing-tool) respectively.

## Read Quality
De novo assemblies are performed using [Shovill](https://github.com/tseemann/shovill) and assembly metrics are collected using [Quast](https://github.com/ablab/quast). Quast calculations for assembly length are used to gauge genome coverage using the CDC's CG Pipeline--Q-score averages for the forward and reverse reads are also gauged using the [CG Pipeline](https://github.com/lskatz/CG-Pipeline).

---

## Quick Start:

````
$ staphb-wf tredegar.py <input_dir> -o <output_dir>
````

`<input_dir>` can be the path to an input directory containing paired-end fastq read data.
If an `<output_dir>` is not provided, results will be written to a `tredegar_run_<date>` directory.


## Other Options
- `--profile`: Nextflow profile, either Docker or Singularity. Default will try docker first, then singularity if the docker executable cannot be found.
- `--config`,`-c`, Path to a custom NextFlow configuration file
- `--get_config`: Get a Nextflow configuration template for Tredegar
- `--resume`: Resume a previous run

## Output:
The final `Tredegar_results.csv` file will be writtein directly to the `<output_dir>`. Mash, Quast, SeqSero, SerotypeFinder, Shovill and Emm-typer results will be writtein to `<tool>` subdirectories within an `<output_dir>/logs` directory, where `<tool>` is the name of the analytical tool that generated the results. The Nexttflow execution report (`Tredegar_execution_report.html`), trace file (`Tredegar_trace.txt`), and task work directories, will also be written to the  `<output_dir>/logs`. E.g. The output of a Tredegar run after analyzing four samples (`ecoli`, `GAS`, `listeria`, and `salmonella`) will be organized as:

`````
tredegar_results/
├── logs
│   ├── 20_05_21_10_58_27_Tredegar_execution_report.html
│   ├── 20_05_21_10_58_27_Tredegar_trace.txt
│   ├── cg_pipeline
│   │   ├── ecoli_readMeterics.tsv
│   │   ├── GAS_readMeterics.tsv
│   │   ├── listeria_readMeterics.tsv
│   │   └── salmonella_readMeterics.tsv
│   ├── emmtyper
│   │   └── GAS_R2.results.xml
│   ├── mash
│   │   └── mash_species.tsv
│   ├── quast
│   │   ├── ecoli_report.tsv
│   │   ├── GAS_report.tsv
│   │   ├── listeria_report.tsv
│   │   └── salmonella_report.tsv
│   ├── seqsero
│   │   └── salmonella_seqsero.txt
│   ├── serotypefinder
│   │   └── ecoli_serotypefinder.txt
│   ├── shovill
│   │   ├── ecoli.contigs.fa
│   │   ├── GAS.contigs.fa
│   │   ├── listeria.contigs.fa
│   │   └── salmonella.contigs.fa
│   ├── work
│   │   |
|   | [...]
|   |
| [...]
└── Tredegar_results.tsv
`````

**logs/Tredegar_execution_report.html** - dryad_workflow_2.<br />
**logs/Tredegar_trace.txt** -. <br />
**logs/cg_pipeline/\<sample\>_readMetrics.csv** - <br />
**logs/emmtyper/\<sample\>_R2.results.xml** - <br />
**logs/mash/mash_species.tsv** - <br />
**logs/quast/\<sample\>_report.tsv** - <br />
**logs/seqsero/\<sample\>_seqsero.txt** - <br />
**logs/serotypefinder/\<sample\>._serotypefinder.txt** - <br />
**logs/shovill/\<sample\>.contigs.fa** - <br />
**logs/work/** - <br />
**Tredegar_results.tsv** -



## Sample Tredegar Report
![Sample_tredegar report](/staphb_toolkit/assets/workflows/tredegar/tred2_sample_out.png)
- sample: isolate ID pulled from the fastq file
- r1_q and r2_q: Average Q-score for the forward and reverse reads, respectively, calculated by CG Pipeline
- est_genome_length: Shovill assembly length calculated by Quast
- est_cvg: Calculated genome coverage calculated by CG Pipeline
- number_contigs: Total number of contigs in assembly calculated by Quast
- predicted_species: Genus and species prediction by Mash
- predicted_serotype: Serotype predictions by SeqSero, SerotypeFinder, and PHE's emm typer for *Salmonella enterica* and *Escherichia coli* and Group A *Streptoccus* isolates, respectively. If a subspecies prediction cannot be made this field will be populated with "NA".  


## Docker Images
The base NextFlow configuration profiles (Docker, Singularity, and AWS) incorporate the following StaPH-B Docker Images:

| Tredegar Process   | Function  | Docker Image  | Comment|
|---|---|---|---|---|
| preProcess  | Renames input read files for downstream processing | staphb/fastqc_container  | Light-weight container for quick text processing  |
| mash_sketch  | Creates mash sketch files for each sample  | staphb/mash:2.1  | `mash sketch` default parameters used |
| mash_dist  | Calculates mash distances for each sketch against a pre-skecteched RefSeq database  | staphb/mash:2.1  | `mash dist` default parameters used to query against the pre-sketch RefSeq database available in the staphb/mash:2.1 Docker Image (/db) |
| mash_species  | Makes taxonomic species predictions based on closest mash distances | staphb/tiptoft:1.0.0  | Light-weight container with python3  |
| trim  | Quality trimming of input read data with bbduk  | staphb/trimmomatic:0.39  | `trimmomatic` parameters set to: minlength=100, windowsize=10, & qualitytrimscore=20|
| cleanreads  | Adapter and PhiX removal from input read data  | staphb/bbtools:38.76  |`bbduk` default parameters used |
| shovill | De novo genome assembly  | staphb/shovill:1.0.4  |`shovill.py` default parameters used |
| quast  | Quality assessment of de novo genome assemblies  | staphb/quast:5.0.2  |`quast` default parameters used |
| cg_pipeline  | Read quality assessment and genome coverage calculations  | staphb/lyveset:1.1.4f  | `run_assembly_readMetrics.pl` in fast mode (--fast)  |
| emmtype_finder  | Emm-type predictions for GAS isolates  | staphb/emm-typing-tool:0.0.1  | `emm_typing.py` default parameters used to query against the GAS database available in the staphb/emm-typing-tool:0.0.1 Docker Image (/db)|
| seqsero  | Serotype predictions for Salmonella isolates  | staphb/seqsero:1.0.1  | `SeqSero.py` default parameters used|
| serotypefinder  | Serotype predictions for E.coli isolates  | staphb/serotypefinder:1.1  |`SeqSero.py` parameters set to: nucleotide agreement (-k)=95.00, coverage (-l)=10; querying against the serotypefinder database available in the staphb/serotypefinder:1.1 Docker Image (/serotypefinder/database)|
| results | Curate all output into single report   |  staphb/tiptoft:1.0.0  | Light-weight container with python3  |

Default docker images and parameters listed above can be adjusted by:
1. Copying the template tredegar config file (`$ staphb-wf tredegar --get_config`)
2. Using a text editor to change the `<date>_tredegar.config` file
3. Specifying your custom config file (i.e. the edited `<date>_tredegar.config>` file) when running the pipeline:<br />

```
$ staphb-wf tredegar <input_dir> -o <output_dir> -c <custom_config_file> [optoins]
```

## Version History

<b>Current version: v2.1 April 21, 2020</b>

Updates in v2.1
1. Written into the StaPH-B ToolKit
2. Report building using standard python libraries
3. Replace SeqyClean with Trimmomatic and BBDuk

Updates in [v2.0](https://github.com/kevinlibuit/Tredegar)
1. Use StaPH-B docker images rather than local downloads
2. Incorporate CG-Pipeline for Q-score and coverage calculations
3. Incorporate Shovill and Quast to get accurate genome length estimations
4. Adjust E.coli serotyping to account for similar alleles
5. Use only O-type predictions from SerotypeFinder

Version 1.0 available [here](https://github.com/kevinlibuit/Tredegar/releases)

## Author
[Kevin G. Libuit](https://github.com/kevinlibuit), DCLS Bioinformatics Lead Scientist
