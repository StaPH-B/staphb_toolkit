#!/usr/bin/env python3

import glob
import re
import os,sys
import argparse
import subprocess
import shutil
#class for containing all sequencing run information
class Basemount:
    # isolates
    ids = []
    #path to project
    path = None
    #output directory
    output_dir = None
    # read files
    reads = []

    def __init__(self,path=None, output_dir = None):
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        self.path = path
        self.out = output_dir

        isolates = os.listdir(path+"/Samples/")

        for isolate in isolates:
            if not isolate.startswith('.'):
                self.ids.append(isolate)

        self.reads = glob.glob(self.path + "/Samples/*/Files/*.fastq.gz")

    def copy_reads(self):
        input_reads_dir = self.out + "/input_reads/"
        if not self.reads:
            raise ValueError("No reads found in %s"%self.path)
        else:
            if not os.path.isdir(input_reads_dir):
                os.makedirs(input_reads_dir)
                print("Directory made for raw read files: " + input_reads_dir)

            for read in self.reads:
                dest = os.path.join(input_reads_dir, re.sub('_(.*)', '', os.path.basename(read)), re.sub('S\d+_L\d+_R', "", os.path.basename(read)))
                dest = dest.replace("_001","")

                #If dest dir doesn't exists, create it
                if not os.path.isdir(os.path.dirname(dest)):
                    os.makedirs(os.path.dirname(dest))

                if not os.path.isfile(dest):
                    print("Copying " + read + " to: " + dest)
                    shutil.copyfile(read, dest)
