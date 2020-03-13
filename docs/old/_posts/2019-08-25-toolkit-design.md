---
path: '/:id'
title: 'Toolkit Design'

layout: nil
---
The StaPH-B ToolKit is designed in a hierarchical directory structure. At the top level the executable **staphb_toolkit** will allow the user access to all of the configured containerized programs and various workflows developed by StaPH-B members. Several sub-directories exist that contain the core functionality of the toolkit.

## core
Foundational python files that represent the core functionality the StaPH-B ToolKit.  
**calldocker**: methods to interact with StaPH-B docker images  
**callsing**: methods to interact with StaPH-B docker images using Singularity  
**autopath**: method to handle altering command file paths to work with in a containerized envrionment  
**docker_config**: a configuration file that describes the versions (tags) of StaPH-B docker images utilized by the toolkit  
**fileparser**: methods to parse and identify files storing sequencing reads and store relevant meta-data (e.g. sample names, PE v. SE, etc.)  
**basemount**: methods to access sequence data stored in a basespace basemount file system  
**sb_programs**: methods surrounding running commands in containerized envrionments. Handles process in both a serial and multiprocessing approach.  

## lib
The lib files are designed to hold tools with additional functionality or wrappers around the tool. These are primarily user developed for use in workflows. These are often more narrowly focused applications that don't extend the core functionality of the StaPH-B Toolkit.

## workflows
The scrips here define custom pipelines that have been developed by StaPH-B members for use across muliple state public health laboratories. These workflows use the core functionality of the Toolkit and adhere to the portability standards of the Toolkit.
