#!/usr/bin/env python3

import os
import argparse
import glob
import re
import sys
#class for containing read information
class SBReads:
    #list contianing all isolate ids
    idList = []

    #data dictonary containing all of the run information
    #formatted id: ["read1 path","read2 path"]
    rawdata = {}

    #data dictonary for containing runtime information
    data = {}

    #mapping progress tracker
    data['mapProgress'] = {}
    data['mapProgress']['normalized'] = []
    data['mapProgress']['trimmed'] = []
    data['mapProgress']['consensus'] = []

    #trimming progress tracker
    data['trimmed'] = {}

    #remove duplicates tracker
    data['rmDuplicates'] = {}

    #normalize reads tracker
    data['normalized'] = {}

    #consensus tracker
    data['consensus'] = {}

    def __init__(self,path):
        #Ensure input path exisits
        global id
        if not os.path.isdir(path):
            raise ValueError(path + " " + "not found.")

        #list of reads
        self.readList = []
        for root,dirs,files in os.walk(path):
            #scan path and look for fastq files and record ids
            for file in files:
                if '.fastq' in file:
                    if '_R1' in file or '_1' in file:
                        id = file.split('_')[0]
                        if id not in self.idList:
                            self.idList.append(id)
                    if '_R2' in file or '_2' in file:
                        id = file.split('_')[0]
                        if id not in self.idList:
                            self.idList.append(id)
                    self.readList.append(root+'/'+file)

        if len(self.readList) == 0:
            raise ValueError("No read fastq read files found in" + " " + path)

        self.readList.sort()
        for id in self.idList:
            self.rawdata[id] = []
            for read in self.readList:
                if id in read:
                    self.rawdata[id].append(read)


    #return a list of id with each paired path as a 3 item sublist
    def retList(self,):
        l = []
        for id in self.idList:
            l.append([id,self.rawdata[id][0],self.rawdata[id][1]])
        return l

    def addData(self,dataType,id,*path):
        if dataType not in self.data:
            self.data[dataType] = {}
        self.data[dataType][id] = path

    def link_reads(self, output_dir):
        raw_reads_dir = os.getcwd() + "/" + output_dir + "/raw_reads/"

        if not os.path.isdir(raw_reads_dir):
            os.makedirs(raw_reads_dir)
            print("Directory for raw reads made: ", raw_reads_dir)

        for read in self.readList:
            dest = raw_reads_dir + re.sub('S\d+_L\d+_R', "", os.path.basename(read))
            dest = dest.replace("_001","")
            print(dest)
            if not os.path.isfile(dest):
                os.symlink(read, dest)
                print("Sym link for", read, "made at", dest)
            else:
                print("Sym link for", read, "already exists at", dest)
        return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_read_data.py <input> [options]")
    parser.add_argument("input", type=str, help="BS Project or path to directory containing read files")
    parser.add_argument("-o", default="", type=str, help="Name of an output_dir")
    parser.add_argument("-link_reads", action='store_true', help="Link reads to project_dir/raw_reads_dir")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.o
    link_reads = args.link_reads

    reads = SBReads(input_dir)

    if len(output_dir) > 0:
        pass
    else:
        output_dir = input_dir

    if link_reads:
        print("linking raw reads. . .")
        reads.link_reads(output_dir)








