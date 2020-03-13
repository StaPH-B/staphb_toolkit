---
layout: page
---

<img src="/docker-builds/assets/user_guide.png" style="display:block;margin-left:auto;margin-right:auto;width:400px">

# Installing Dependencies
The Toolkit has been designed to minimize the amount of needed dependencies. The ToolKit has been built using Python 3.6 which can easily be installed on any Unix/Linux operating system following the instructions below. The design philosophy is centered around usability including the installation of dependencies. Using containerization the Toolkit is able to access numerous bioinformatics applications with out the necessity of installing various often conflicting dependencies. However, because of the use of containerization there are several dependencies that are unavoidable.

### Python Version >= 3.6
The Toolkit uses python to interact with the containers and workflows. Most Linux systems come with an adequate version of python.  
To check your version of python use the following command:
```
python -V
```

If you do not have Python 3.6 or greater the best way to install python is using the package installer [Anaconda](https://www.anaconda.com/). This can be done with the following commands:
```
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
source ~/.bashrc
```
**Note: When prompted with “Do you wish the installer to initialize Anaconda3 by running conda init?” We recommend “yes”.**




The system must have either Docker or Singularity installed.
