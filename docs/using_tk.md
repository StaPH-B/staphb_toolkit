---
title: "Using the ToolKit"
layout: page
---

Once the ToolKit has been installed using the ToolKit is simply done by running the command `staphb-tk`. If you would like more information about the different workflows available in the ToolKit visit the [workflows](/staphb_toolkit/workflows) page.

# Contents
  * [Running Applications](#running-applications-using-the-toolkit)
    - [Pipes and Paths](#special-note-about-autopathing-and-pipes)
    - [Included Applications](#included-applications)
  * [Running Workflows](#using-the-toolkit-to-run-workflows)
    - [Controlling Versions](#controlling-versions)
    - [Configuring Workflows](#configuring-workflows)

---

# Running applications using the ToolKit
Running the toolkit using the the `staphb-tk` command provides a menu of options available for running tools in the toolkit:
````
$ staphb-tk
     _______.___________.    ___      .______    __    __         .______   
    /       |           |   /   \     |   _  \  |  |  |  |        |   _  \  
   |   (----`---|  |----`  /  ^  \    |  |_)  | |  |__|  |  ______|  |_)  |
    \   \       |  |      /  /_\  \   |   ___/  |   __   | |______|   _  <  
.----)   |      |  |     /  _____  \  |  |      |  |  |  |        |  |_)  |
|_______/       |__|    /__/     \__\ | _|      |__|  |__|        |______/  

.___________.  ______     ______    __       __  ___  __  .___________.
|           | /  __  \   /  __  \  |  |     |  |/  / |  | |           |
`---|  |----`|  |  |  | |  |  |  | |  |     |  '  /  |  | `---|  |----`
    |  |     |  |  |  | |  |  |  | |  |     |    <   |  |     |  |     
    |  |     |  `--'  | |  `--'  | |  `----.|  .  \  |  |     |  |     
    |__|      \______/   \______/  |_______||__|\__\ |__|     |__|     


StaPH-B ToolKit
Version: 2.0.0
usage: staphb-tk [optional arguments] <application/workflow> [application/workflow arguments]

optional arguments:
-h, --help            show this help message and exit
-l, --list_tools      List all tools in the toolkit.
-w, --list_workflows  List all workflows in the toolkit.
-wv <version>, --workflow_version <version>
                   Version of tool or workflow to run. Default: latest
-c <config_file>, --configuration <config_file>
                   Specify a custom workflow configuration file.
-gc, --get_configuration
                   Get the configuration file for the specified workflow.
                   Note: You may need to specify a version for the
                   workflow using -wv to get the correct configuration
                   file.
-nv [<version>], --nextflow_version [<version>]
                   Get or set the version of nextflow.
--update              Check for and install a ToolKit update.
--auto_update         Toggle automatic ToolKit updates. Default is off.

application or workflow name:
<application/workflow>

```

The typical usage of the ToolKit involves a command structure that calls the toolkit i.e. `staphb-tk` followed by the application i.e. `bwa` then the parameters associated with that tool. For example using the command `staphb-tk bwa` the output shows the options available to the bwa alignment tool:

```
$ staphb-tk bwa   
Status: Downloaded newer image for staphb/bwa:latest
Pulling staphb/bwa:latest ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00

Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
Contact: Heng Li <lh3@sanger.ac.uk>

Usage:   bwa <command> 

Command: index         index sequences in the FASTA format
mem           BWA-MEM algorithm
fastmap       identify super-maximal exact matches
pemerge       merge overlapping paired ends (EXPERIMENTAL)
aln           gapped/ungapped alignment
samse         generate alignment (single ended)
sampe         generate alignment (paired ended)
bwasw         BWA-SW for long queries

shm           manage indices in shared memory
fa2pac        convert FASTA to PAC format
pac2bwt       generate BWT from PAC
pac2bwtgen    alternative algorithm for generating BWT
bwtupdate     update .bwt to the new format
bwt2sa        generate SA from BWT and Occ

Note: To use BWA, you need to first index the genome with `bwa index'.
There are three alignment algorithms in BWA: `mem', `bwasw', and
`aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
first. Please `man ./bwa.1' for the manual.

```
<br>
#### Special note about autopathing and pipes
The ToolKit will automatically mount paths in your command from your host file system. This allows the toolkit to interact with docker or singularity containers without needing your input on how to mount file systems. However, if you wish to use a path to a file contained inside the container the autopathing will still try to find that file on your host system therefore you must use an `$` to indicate the path is located inside the container as shown in the command below:

```
staphb-tk mash dist $/db/RefSeqSketchesDefaults.msh mash_output/sample_sketch.msh
```

In addition, pipes (`<`,`>`,`|`) are by default read by the bash interpreter. If you wish to use a pipe in your command and you want that pipe to run inside the container, you must use the bash escape character `\` to signify that you want the pipe run in the container. For example:

```
staphb-tk mash dist $/db/RefSeqSketchesDefaults.msh mash_output/sample_sketch.msh \> mash_output/sample_distances.tab
```
<br>
### Included Applications

To list the available software included with the ToolKit use the command `staphb-tk -l` or `staphb-tk --list_tools`.

---

# Using the ToolKit to run workflows
The ToolKit also provides the ability to run workflows. To list the available workflows using the `staphb-tk -w` or `staphb-tk --list_workflows` command.

```
$ staphb-tk -w
Available workflows:                                                                                                              
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃          ┃                                                                                                                     ┃
┃Command   ┃Description                                                                                                          ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│          │                                                                                                                     │
│bactopia  │Bactopia is a flexible pipeline for complete analysis of bacterial genomes. The goal of Bactopia is to process your  │
│          │data with a broad set of tools, so that you can get to the fun part of analyses quicker!                             │
│          │                                                                                                                     │
│cecret    │Cecret is a workflow developed by @erinyoung at the Utah Public Health Laborotory for SARS-COV-2 sequencing with the │
│          │artic/Illumina hybrid library prep workflow for MiSeq data.                                                          │
│          │                                                                                                                     │
│dryad     │Dryad is a pipeline to construct reference free core-genome or SNP phylogenetic trees for examining prokaryote       │
│          │relatedness in outbreaks. Dryad will performs both a reference free core-genome and/or a SNP analysis using the      │
│          │CFSAN-SNP pipeline.                                                                                                  │
│          │                                                                                                                     │
│mycosnp   │MycoSNP is a portable workflow for performing whole genome sequencing analysis of fungal organisms, including Candida│
│          │auris.                                                                                                               │
│          │                                                                                                                     │
│spriggan  │Spriggan is a pipeline used for assembly of bacterial whole genome sequence data and identification of antibiotic    │
│          │resistance genes.                                                                                                    │
│          │                                                                                                                     │
│viralrecon│Viralrecon is a pipeline used to perform assembly and intra-host/low-frequency variant calling for viral samples.    │
└──────────┴─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

```

Additional help is available for each workflow by using the workflow name and the `-h` flag i.e. `staphb-tk bactopia -h`. This will list all of the parameters available for the workflow. For further information about each of the workflows visit their documentation.

### Controlling Versions
The toolkit allows changing the version of nextflow used to run the workflows as well as the version of the workflow. In order to specify a specific version of nextflow use the `-nv` or `--nextflow_version` flag. Likewise to specify a workflow version use the `-wv` or `--workflow_version` flags.

```
staphb-tk -nv 22.04.5 -wv 3.3.20220810 cecret
```

### Configuring Workflows 
Workflow configurations can be obtained using the `-gc` or `--get_configuration` flags. This will either download the most recent configuration file or if a version is specified using the `-wv` flag then the appropriate version of the configuration will be downloaded. For further information on how Nextflow utilizes the configuration file visit the nextflow documentation [here](https://www.nextflow.io/docs/latest/config.html) or visit the workflows documentation page.