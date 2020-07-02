---
title: 'Cutshaw'
layout: page
---

# Cutshaw v1.0
Bioinformatics Pipeline for Instrument Validation and Assessment Technical Proficiency in WGS Protocols

## Data workflow:
![Cutshaw pipeline](/staphb_toolkit/assets/workflows/cutshaw/Cutshaw_v1.0.png)

Cutshaw, a DCLS-developed workflow based on the [U.S. Food and Drug Administration’s (FDA) Center for Food Safety and Applied Nutrition (CFSAN) GenomeTrakr Proficiency Assessment workflow](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000185), was developed to assess a scientist’s proficiency on the Illumina MiSeq sequencing platform. CutShaw was developed in collaboration with FDA CFSAN and derives minimum quality thresholds from the [2018 Genomics for Food Safety (Gen-FS) Proficiency Testing exercise](https://github.com/CFSAN-Biostatistics/wgs_competency).


---

## Quick Start:

````
$ staphb-wf cutshaw <input_dir> --isolate_key <isolate_key_file> --report_title "<report_title>"
````

`<input_dir>` can be the path to an input directory containing paired-end fastq read data.
`<report_title>` can be any string designating the title of the cutshaw report (e.g. name of instrument or scientist being assessed); must be enclosed in quotes
`isolate_key` must be a csv file indicating the 2018 FDA Gen-FS PT isolate (isolate_id) associated with each sample (sample_id) in the project directory. e.g.

```
$ ls .
input_dir
isolate_key.csv

$ ls input_dir
sample_001_R1.fastq.gz
sample_001_R2.fastq.gz
sample_002_R1.fastq.gz
sample_002_R2.fastq.gz

$ cat isolate_key.csv
sample_id,isolate_id
sample_001,SAP18-8729
sample_002,SAP18-H9654
```
PT isolates currently compatable  with the Cutshaw pipeline:
- SAP18-0432    Salmonella enterica subsp. enterica serovar Enteritidis
- SAP18-H9654   Salmonella enterica subsp. enterica serovar Enteritidis
- SAP18-6199    Salmonella enterica subsp. enterica serovar Typhimurium
- SAP18-8729    Salmonella enterica subsp. enterica serovar Newport
- LMP18-H2446   Listeria monocytogenes
- LMP18-H8393   Listeria monocytogenes

If an `<output_dir>` is not provided, results will be written to a `tredegar_run_<date>` directory.


## Other Options
- `-o`: Output directory
- `--profile`: Nextflow profile, either Docker or Singularity. Default will try docker first, then singularity if the docker executable cannot be found.
- `--config`,`-c`, Path to a custom NextFlow configuration file
- `--get_config`: Get a Nextflow configuration template for Tredegar
- `--resume`: Resume a previous run

## Output:
The final report file will be written directly to the `<output_dir>`.

![Sample report](/staphb_toolkit/assets/workflows/cutshaw/sample_report.html)


## Version History

<b>Current version: v1.0.0 July 03, 2020</b>

Version 1.0.0 is the first stable version of Cutshaw

## Authors
[Kevin G. Libuit](https://github.com/kevinlibuit), DCLS Bioinformatics Lead Scientist <br />
[Joseph Baugher](https://github.com/jdbaugher), FDA Bioinformatics Scientist
