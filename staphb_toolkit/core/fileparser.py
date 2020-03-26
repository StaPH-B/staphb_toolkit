#!/usr/bin/env python3

import os
import re
from staphb_toolkit.core.basemount import Basemount

#class for containing all sequencing run information
class ProcessFastqs:
    #list contianing all isolate ids
    ids = []

    #class containing all of the read information for an isolate
    class Fastqs:
        #isolate id
        id = ''
        #path to forward reads fastq
        fwd = ''
        #path to reverse reads fastq
        rev = ''
        #path to reads fastq if unpaired or interleaved
        path = ''
        #switch for interleaved
        inteleaved = False
        #switch for paired or unpaired
        paired = True
        def __init__(self,id='',fwd='',rev='',path='',interleaved=False):
            #set all variables to definition arguments
            self.id = id
            self.fwd = fwd
            self.rev = rev
            self.path = path
            self.interleaved = interleaved
            if path:
                self.paired = False

    #data dictonary containing all of the read objects
    reads = {}

    def __init__(self,path, output_dir=''):
        #Ensure input path exisits
        if not os.path.isdir(path):
            raise ValueError(path + " " + "not found.")

        if not output_dir:
            output_dir = os.getcwd()

        if os.path.isdir(path+"/AppResults"):
            basemount_project = Basemount(path, output_dir)
            basemount_project.copy_reads()
            path = os.path.join(output_dir, "input_reads")


        for root,dirs,files in os.walk(path):
            #scan path and look for fastq files then gather ids and store in temp lists
            for file in files:
                if '.fastq' in file or '.fastq.gz' in file:
                    #get id and check if we have seen this id before by adding to id list and creating a new read object
                    id = file.split('_')[0]
                    if id not in self.ids:
                        self.ids.append(id)
                        self.reads[id] = self.Fastqs(id)

                    #if fastq file is foward reads add path to .fwd
                    if '_R1' in file or '_1' in file:
                        if not self.reads[id].fwd:
                            self.reads[id].fwd = root + '/' + file
                    #if fastq file is reverese reads add path to .rev
                    elif '_R2' in file or '_2' in file:
                        if not self.reads[id].rev:
                            self.reads[id].rev = root + '/' + file
                    #if fastq file is unpaired or interleaved add path to .path
                    #TODO consider the impact of this and determine best method
                    else:
                        if not self.reads[id].path:
                            self.reads[id].paired = False
                            self.reads[id].path = root + '/' + file

        #notify user if no fastq files were found
        if len(self.ids) == 0:
            raise ValueError("No fastq files found in " + path)

    #return a list of id with each fastq path
    def id_list(self):
        output_list = []
        for id in self.ids:
            if self.reads[id].paired:
                output_list.append([id,self.reads[id].fwd,self.reads[id].rev])
            else:
                output_list.append([id,self.reads[id].path])
        return output_list

    #return a dictonary of id with each fastq path
    def id_dict(self):
        return self.reads

    #return a list of all reads
    def fastq_paths(self):
        output_list = []
        for id in self.ids:
            if self.reads[id].paired:
                output_list.append(self.reads[id].fwd)
                output_list.append(self.reads[id].rev)
            else:
                output_list.append(self.reads[id].path)
        return output_list

    #create symbolic links to raw reads
    def link_reads(self, output_dir):
        #get output directory
        out_path = os.path.abspath(output_dir)
        input_reads_dir = os.path.join(out_path + "/input_reads")

        #for each fastq file create a symbolic link
        for fastq in self.fastq_paths():

            dest = os.path.join(input_reads_dir, os.path.basename(fastq).split('_')[0], re.sub('S\d+_L\d+_R', "", os.path.basename(fastq)))
            dest = dest.replace("_001","")
            dest = dest.replace("_R1", "_1")
            dest = dest.replace("_R2", "_2")

            # If dest dir doesn't exists, create it
            if not os.path.isdir(os.path.dirname(dest)):
                os.makedirs(os.path.dirname(dest))

            if not os.path.exists(dest):
                os.link(fastq, dest)
                print("Hard link for", fastq, "made at", dest)

        # reset fwd/rev paths
        for root, dirs, files in os.walk(output_dir):
            # scan path and look for fastq files then gather ids and store in temp lists
            for file in files:
                if '.fastq' in file or '.fastq.gz' in file:
                    # get id and check if we have seen this id before by adding to id list and creating a new read object
                    id = file.split('_')[0]
                    # if fastq file is foward reads add path to .fwd
                    if '_R1' in file or '_1' in file:
                        del self.reads[id].fwd
                        self.reads[id].fwd = root + '/' + file
                    # if fastq file is reverese reads add path to .rev
                    elif '_R2' in file or '_2' in file:
                        del self.reads[id].rev
                        self.reads[id].rev = root + '/' + file

        return None

    def inputSubdomain(self,config_object={}):
        #see if we got a config object
        if not config_object:
            config_object = {}

        #create file_io if it doesn't exists
        if 'file_io' not in config_object:
            config_object['file_io'] = {}
        if 'input_files' not in config_object['file_io'] or config_object['file_io']['input_files'] == None :
            config_object['file_io']['input_files'] = {}

        for fastqs in self.reads:
            config_object['file_io']['input_files'][self.reads[fastqs].id] = [self.reads[fastqs].fwd,self.reads[fastqs].rev]

        return config_object
