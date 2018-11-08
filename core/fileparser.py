#!/usr/bin/env python3

import os
import argparse
import glob
import re
import sys
from staphB_ToolKit.core import basemount


#class for containing all sequencing run information
class RunFiles:
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

    #data dictonary for containing runtime information
    runtime = {}

    def __init__(self,path, output_dir=''):
        #Ensure input path exisits
        if not os.path.isdir(path):
            raise ValueError(path + " " + "not found.")

        if not output_dir:
            output_dir = os.getcwd()

        if os.path.isdir(path+"/AppResults"):
            basemount_project = basemount.Basemount(path, output_dir)
            basemount_project.copy_reads()
            path = output_dir


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
    def return_id_list(self):
        output_list = []
        for id in self.ids:
            if self.reads[id].paired:
                output_list.append([id,self.reads[id].fwd,self.reads[id].rev])
            else:
                output_list.append([id,self.reads[id].path])
        return output_list

    #return a list of all reads
    def return_fastq_list(self):
        output_list = []
        for id in self.ids:
            if self.reads[id].paired:
                output_list.append(self.reads[id].fwd)
                output_list.append(self.reads[id].rev)
            else:
                output_list.append(self.reads[id].path)
        return output_list

    #function to add runtime data to class **may get depreciated**
    def add_runtime(self,data_type,id,*path):
        if data_type not in self.runtime:
            self.runtime[data_type] = {}
        self.runtime[data_type][id] = path
        return None

    #create symbolic links to raw reads
    def link_reads(self, output_dir):
        #get output directory
        out_path = os.path.abspath(output_dir)
        raw_reads_dir = os.path.join(out_path + "/raw_reads")

        #check if output directory exists, if not create it
        if not os.path.isdir(raw_reads_dir):
            os.makedirs(raw_reads_dir)
            print("Directory for raw reads made: ", raw_reads_dir)

        #for each fastq file create a symbolic link
        for fastq in self.return_fastq_list():
            dest = os.path.join(raw_reads_dir, re.sub('S\d+_L\d+_R', "", os.path.basename(fastq)))
            dest = dest.replace("_001","")

            if not os.path.isfile(dest):
                os.symlink(fastq, dest)
                print("Symbolic link for", fastq, "made at", dest)
            else:
                print("Symbolic link for", fastq, "already exists at", dest)
        return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="fileparser.py <input> [options]")
    parser.add_argument("input", type=str, help="Path to fastq files.")
    parser.add_argument("-o", default="", type=str, help="Path to output directory.")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.o

    run = RunFiles(input_dir)

    if output_dir:
        pass
    else:
        output_dir = input_dir

    print("Isolates found in " + input_dir + ":" + run.ids )
