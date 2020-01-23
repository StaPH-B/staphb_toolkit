# StaPH-B ToolKit
[![Build Status](https://travis-ci.org/StaPH-B/staphb_toolkit.svg?branch=master)](https://travis-ci.org/StaPH-B/staphb_toolkit)  
A python library designed to make programs held within the StaPH-B Docker repository more accessible to public health scientists.

## Summary
The StaPH-B ToolKit is a Python library (Python 3.7) of commonly used bioinformatics tools that help to inform public health action. The StaPH-B ToolKit utilizes the [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) to enable easy access of open-source software without the need of local installation and/or dependency maintenance.

## Motivation
Public health bioinformatics is dependent on open-source software that require carefully curated computational environments and various software dependencies. Setting up and maintaining such environments requires a skill set and expertise absent in most public health laboratories. The [StaPH-B Docker Images](https://github.com/StaPH-B/docker-builds) have helped generate reproducible computational environments through the use of containerization. However, access to these images is dependent on a working understanding of containerization, which is not available in most laboratories. The ToolKit addresses this issue through the handling of the StaPH-B docker images allowing users to interact with bioinformatis programs without needing to interact directly with mounted file systems and running containers. The goal of the Toolkit is it increase usability while mirroring the functionality of a locally-installed tool.

## Installing and Usage
The full documentation can be found here: [https://staph-b.github.io/staphb_toolkit](https://staph-b.github.io/staphb_toolkit).

### Installing dependencies
The toolkit requires **either** singularity or docker. To install these on a debian based system use the following commands:

#### Docker
```
sudo apt-get update
sudo apt install apt-transport-https ca-certificates curl gnupg-agent software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
sudo usermod -aG docker ${USER}
sudo su - ${USER}
```
#### Singularity
```
sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
sudo apt-get update
sudo apt-get install singularity-container
```
#### Python 3.7 or greater
Python can be installed a number of ways but we recommend using Anaconda:
```
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
source ~/.bashrc
```
When prompted with “Do you wish the installer to initialize Anaconda3 by running conda init?” We recommend “yes”.

### Installing the ToolKit
First download the Toolkit using git:  
`git clone https://github.com/StaPH-B/staphb_toolkit.git`

Then install the python dependencies using pip:  
`pip install -r staphb_toolkit/requirements.txt`  

Run the pipeline with either of the following commands:  
```
staphb_toolkit
staphb_toolkit_workflows
```
