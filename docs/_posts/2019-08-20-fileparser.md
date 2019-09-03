---
category: Core
path: '/Core/:id'
title: 'core.fileparser'

layout: nil
---
The **core.fileparser** contains methods for parsing directories and storing file information. These objects help with the handling of paired files that are regularly generated during sequencing.  
<br>
## *class* core.fileparser.**ProcessFastqs**(*path,basespace_output_dir=""*)  
#### The ProcessFastqs object performs an inventory of sequencing read files in the fastq format and identifies paired read files.  
<br>
#### The ProcessFastqs object should always be called with a *path* pointing to the directory containing the fastq files. The fastq files can either all be located in the directory or split into sub-directories. In the case of sub-directories pointing to the top level directory will allow the capture of all fastq files in the sub-directories.
<br>
#### The object also has an optional *basespace_output_dir* parameter. If this parameter is given a path the object will check to see if the files are mounted via basespace basemount. If the conditional is met the files will be automatically copied to the path provided in the *basespace_output_dir*. If the files are mounted via basespace basemount and the *basespace_output_dir* is an empty string. The files will be copied to the current working directory.
<pre><code class="python">import core.fileparser as fp

read_path = '/path/to/read/files'
fastqs = fp.ProcessFastqs(read_path)
</code></pre>

### **id_list()**  
#### Returns a nested list where the structure of the list is [['ID','fwd read path','rev read path'],...]. The ID is parsed from the fastq file name incorporating everything up until the first underscore.
<pre><code class="python">import core.fileparser as fp

read_path = '/path/to/read/files'
fastqs = fp.ProcessFastqs(read_path)
list_of_fastqs = fastqs.id_list()
</code></pre>
### **fastq_paths()**  
#### Returns an ordered list of paths to all fastq files.
<br>
### **id_dict()**  
#### Returns a python dictionary where the key is the ID parsed from the fastq file name incorporating everything up until the first underscore and the value is an object describing the fastq files outlined below:
<br>
## *class* core.fileparser.ProcessFastqs.**Fastqs**(*id='',fwd='',rev='',path='',interleaved=False*)  
#### This object stores the information for the read files. This object is created for every set of reads during the initialization of the ProcessFastqs object.  
#### *id* is the ID parsed from the fastq file name incorporating everything up until the first underscore.
#### *fwd* is the path to the forward read.
#### *rev* is the path to the reverse read.
#### *path* is the path to the read if the read file is unpaired or interleaved.
#### *interleaved* is a Boolean indicator that the read file is interleaved.
<pre><code class="python">import core.fileparser as fp

read_path = '/path/to/read/files'
fastqs = fp.ProcessFastqs(read_path)
fastq_dict = fastqs.id_dict()
fwd_read_path = fastq_dict['id'].fwd
rev_read_path = fastq_dict['id'].rev
</code></pre>
