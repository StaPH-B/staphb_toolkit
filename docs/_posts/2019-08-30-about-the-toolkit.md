---
path: '/:id'
title: 'About the Toolkit'

layout: nil
---

## Summary
The StaPH-B ToolKit is a Python library (Python 3.7) of commonly used bioinformatics tools that help to inform public health action. The StaPH-B ToolKit utilizes the [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) to enable easy access of open-source software without the need of local installation and/or dependency maintenance.

## Motivation
Public health bioinformatics is dependent on open-source software that require carefully curated computational environments and various software dependencies. Setting up and maintaining such environments requires a skill set and expertise absent in most public health laboratories. The [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) have helped generate reproducible computational environments through the use of containerization. However, access to these images is dependent on a working understanding of containerization, which is not available in most laboratories. The ToolKit addresses this issue through the handling of the StaPH-B docker images allowing users to interact with bioinformatis programs without needing to interact directly with mounted file systems and running containers. The goal of the Toolkit is it increase usability while mirroring the functionality of a locally-installed tool.

Additionally, the StaPH-B ToolKit will provide a single repository for public health bioinformatics scientists to package and distribute custom analytical pipelines for specific public health use-cases.

## Dependency
The Toolkit has been designed to minimize the amount of needed dependencies. The ToolKit has been built using Python 3.7 which can easily be installed on any Unix/Linux operating system. The design philosophy is centered around usability including the installation of dependencies. Using containerization the Toolkit is able to access numerous bioinformatics applications with out the necessity of installing various often conflicting dependencies. However, because of the use of containerization there are several dependencies that are unavoidable. The system must have either [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/) installed. The instructions for installing them can be found [here](#container-installation-instructions). Additionally, the toolkit requires two python libraries to interact with the container engines *spython* and *docker_py*. These can be installed by using *pip* in the following way:
<pre><code class="bash">pip install -r staphb_toolkit/requirements.txt
</code></pre>

Workflows added to the ToolKit must also follow these standards to prevent decreasing the usability of the overall ToolKit.
