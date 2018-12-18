#!/usr/bin/env python3

import glob
import re
import os,sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

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
                dest = raw_reads_dir + re.sub('S\d+_L\d+_R', "", os.path.basename(read))
                dest = dest.replace("_001","")

                if os.path.isfile((dest)):
                    pass
                else:
                    print("Copying " + read + " to: " + dest)
                    shutil.copyfile(read, dest)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="fileparser.py <input> [options]")
    parser.add_argument("input", type=str, help="Path to fastq files.")
    parser.add_argument("-o", default="", type=str, help="Path to output directory.")
    parser.add_argument("-copy_reads", action='store_true', help="copy read files to output dir.")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.o
    copy_reads = args.copy_reads

    if output_dir:
        pass
    else:
        output_dir = os.getcwd()

    project = Basemount(input_dir, output_dir=output_dir)

    print("Isolates in BaseSpace Project " + os.path.basename(input_dir) + ": " + str(project.ids))
    print("Read files: "+ str(project.reads))

    if copy_reads:
        print("copying fasatq files to " + output_dir)
        project.copy_reads()