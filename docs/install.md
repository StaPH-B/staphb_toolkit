---
title: "Installing the ToolKit"
layout: page
---

# Installing Dependencies
The Toolkit has been designed to minimize the amount of needed dependencies. The ToolKit has been built using Python 3.6 which can easily be installed on any Unix/Linux operating system following the instructions below. The design philosophy is centered around usability including the installation of dependencies. Using containerization the Toolkit is able to access numerous bioinformatics applications with out the necessity of installing various often conflicting dependencies. However, because of the use of containerization there are several dependencies that are unavoidable.

### Python Version >= 3.6
The Toolkit uses python to interact with the containers and workflows. Most Linux systems come with an adequate version of python.  
To check your version of python use the following command:
```
python -V
```

If you do not have Python 3.6 or greater the best way to install python is using the package installer [Anaconda](https://www.anaconda.com/). This can be done on an Debian/Ubuntu Linux system with the following commands:
```
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
source ~/.bashrc
```
**Note: When prompted with “Do you wish the installer to initialize Anaconda3 by running conda init?” We recommend “yes”.**

For a different system download it through the [Anaconda website](https://www.anaconda.com/distribution/).

### Docker or Singularity
The system must have either Docker or Singularity installed in order to use any of the tools in the Toolkit.  
To check if you have docker installed and working use the following command:  
```
docker --version
docker run hello-world
```
To check if you have singularity installed and working use the following command:  
```
singularity --version
singularity run library://sylabsed/examples/lolcow
```

##### Docker
Installation instructions for installing Docker on an Ubuntu Linux system:  
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

##### Singularity
Installation instructions for installing Singularity on an Debian/Ubuntu Linux system:
```
sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
sudo apt-get update
sudo apt-get install singularity-container
```
