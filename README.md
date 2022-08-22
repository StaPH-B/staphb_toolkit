# StaPH-B ToolKit
![Latest Release](https://img.shields.io/github/v/release/StaPH-B/staphb_toolkit)  
[![Build Status](https://travis-ci.org/StaPH-B/staphb_toolkit.svg?branch=master)](https://travis-ci.org/StaPH-B/staphb_toolkit)  

A python library designed to make bioinformatics piplines and applications more accessible to public health scientists.

### [staphb.org/staphb_toolkit/](https://staphb.org/staphb_toolkit/)

## Summary
The StaPH-B ToolKit is a Python library of commonly used bioinformatics tools that help to inform public health action. The StaPH-B ToolKit utilizes the [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) to enable easy access of open-source software without the need of local installation and/or dependency maintenance.

## Motivation
Public health bioinformatics is dependent on open-source software that require carefully curated computational environments and various software dependencies. Setting up and maintaining such environments requires a skill set and expertise absent in most public health laboratories. The [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) have helped generate reproducible computational environments through the use of containerization. However, access to these images is dependent on a working understanding of containerization, which is not available in most laboratories. The ToolKit addresses this issue through the handling of the StaPH-B docker images allowing users to interact with bioinformatis programs without needing to interact directly with mounted file systems and running containers. The goal of the Toolkit is it increase usability while mirroring the functionality of a locally-installed tool.

## Installing and Usage
The ToolKit requires **either** singularity or docker, Python 3.7 or greater, and Java version 8 or later.
The documentation for installing the dependencies can be found here: [https://staph-b.github.io/staphb_toolkit/install](https://staph-b.github.io/staphb_toolkit/install).  
The ToolKit itself can be installed using pip or by cloning the repository from git:

To install using pip:
```
$ pip install staphb_toolkit
```

To install using git:
```
$ git clone https://github.com/StaPH-B/staphb_toolkit.git
$ cd staphb_toolkit/packaging/
$ ./setup.py install
```

Test the pipeline with the following command and ensure you see the same usage output:  
```
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
