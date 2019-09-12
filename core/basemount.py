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
        raw_reads_dir = self.out + "/raw_reads/"
        if not self.reads:
            raise ValueError("No reads found in %s"%self.path)
        else:
            if not os.path.isdir(raw_reads_dir):
                os.makedirs(raw_reads_dir)
                print("Directory made for raw read files: " + raw_reads_dir)

            for read in self.reads:
                dest = os.path.join(raw_reads_dir, re.sub('_(.*)', '', os.path.basename(read)), re.sub('S\d+_L\d+_R', "", os.path.basename(read)))
                dest = dest.replace("_001","")

                if os.path.isfile((dest)):
                    pass
                else:
                    print("Copying " + read + " to: " + dest)
                    shutil.copyfile(read, dest)
