# StaPH-B ToolKit
A python library designed to make programs held within the StaPH-B Docker repository more accessible to public health laboratorians.

### How to use the ToolKit
The ToolKit has been designed to be used by Bioinformaticians with varying levels of experience. The ToolKit has core functions to work with sequencing data files and Docker containers. Commonly used programs have also been incorporated allowing users to run programs while abstracting the Docker container interaction. Finally, some pipelines have been included to provide entire data analysis workflows within this ToolKit. See below for more information on each area.

#### core scripts
Within the ToolKit is a core folder containing several scripts that can be used for various lower level tasks needed in working with the programs and sequence data.
**basemount.py** - provides access to BaseSpace sequencing data storage
**calldocker.py** - accesses the local Docker environment, pulls,starts, and interacts with Docker containers
**docker_config.json** - default docker containers used for the programs and pipelines
**fileparserpy** - provides a data structure to work with sequencing data throughout the pipelines
**sb_libs.py** - provides a common structure for running Bioinformatic programs contained within the **lib** directory
