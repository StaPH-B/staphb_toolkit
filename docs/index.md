---
layout: page
---

<a href="https://staph-b.github.io/staphb_toolkit/"><img src="assets/staphb-tk_logo.png" style="display:block;margin-left:auto;margin-right:auto;width:600px"></a>

#### Introduction
Installing bioinformatics software can be challenging as applications can often have specific and conflicting dependency requirements. StaPH-B's Docker project [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) has helped by containerizing and distributing commonly used open-source bioinformatics software. The Docker project supports the implementation of reproducible computational environments through docker images. Accessing and using Docker images can however sometimes be challenging. The ToolKit addresses this limitation by allowing users to interact with bioinformatics programs without needing to interact directly with mounted file systems and running containers. The goal of the Toolkit is to increase usability of StaPH-B dockerized applications by mirroring the functionality of a locally-installed tool. 

Additionally, the StaPH-B ToolKit provides access to a variety of custom analytical workflows that have been developed by and for public health bioinformatics scientists. For a full list of workflows visit the workflow page [here]({{ "/workflows" | prepend: site.baseurl }}).

#### Installing
The StaPH-B ToolKit is a Python application (Python 3.7) that utilizes modern workflow approach [Nextflow](https://www.nextflow.io/) alongside containerization via [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io) of commonly used bioinformatics tools that help to inform public health action. The StaPH-B ToolKit utilizes the StaPH-B Docker Images to enable easy access of open-source software without the need of local installation and/or dependency maintenance.

This tool was developed and tested for use on an Ubuntu/Debian OS. However, the tool could also be used on other Linux/Unix and MacOS systems as long as the proper dependencies are met.

#### Contributing
If you are interested in contributing a workflow to the StaPH-B Toolkit please reach out on our slack channel or create an issue on [GitHub](https://github.com/StaPH-B/staphb_toolkit). If you would like to contribute to a containerized application in the toolkit, you can do so through our [StaPH-B Docker Project](https://github.com/StaPH-B/docker-builds).
