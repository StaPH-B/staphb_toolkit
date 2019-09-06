# StaPH-B ToolKit
A python library designed to make programs held within the StaPH-B Docker repository more accessible to public health scientists.

## Summary
The StaPH-B ToolKit is a Python library (Python 3.7) of commonly used bioinformatics tools that help to inform public health action. The StaPH-B ToolKit utilizes the [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) to enable easy access of open-source software without the need of local installation and/or dependency maintenance.

## Motivation
Public health bioinformatics is dependent on open-source software that require carefully curated computational environments and various software dependencies. Setting up and maintaining such environments requires a skill set and expertise absent in most public health laboratories. The [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) have helped generate reproducible computational environments through the use of containerization. However, access to these images is dependent on a working understanding of containerization, which is not available in most laboratories. The ToolKit addresses this issue through the handling of the StaPH-B docker images allowing users to interact with bioinformatis programs without needing to interact directly with mounted file systems and running containers. The goal of the Toolkit is it increase usability while mirroring the functionality of a locally-installed tool.

## Installing and Usage
The full documentation can be found here: [https://staph-b.github.io/staphb_toolkit](https://staph-b.github.io/staphb_toolkit).

Installing the pipeline is done with two simple commands. Note this does not include the installation of Docker or Singularity. One of these container engines must be present for the Toolkit to function. The documentation above has instructions for installing the container engine.

First download the Toolkit using git:  
`git clone https://github.com/StaPH-B/staphb_toolkit.git`

Then install the python dependencies using pip:  
`pip install -r staphb_toolkit/requirements.txt`  
There will only ever be two python dependencies, which are the libraries to interact with either [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/).

Run the pipeline with either of the following commands:  
`staphb_toolkit`  
`staphb_toolkit_workflows`
