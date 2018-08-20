#!/usr/bin/env python3

import os
import glob
import re
import sys
import argparse
import datetime
from time import gmtime, strftime

class LocalReads(object):

    def __init__(self, input_dir):
        self.name = input_dir  # path to directory containing read files
        self.path = os.path.expanduser("~") + "/" + self.name  # path to BS project
        self.reads = glob.glob(self.path + "/*.fastq*")  # Read files for all samples within input dir

        if not os.path.isdir(self.name):
            raise ValueError(self.name + " " + "not found.")

        if len(self.reads) == 0:
            raise ValueError("No read files in" + " " + self.input)

        return

    def link_raw_reads(self, output_dir=None, samples=None):
        raw_reads_dir = os.getcwd() + "/" + output_dir + "/raw_reads/"

        if not os.path.isdir(raw_reads_dir):
            os.makedirs(raw_reads_dir)
            print("Directory for raw reads made: ", raw_reads_dir)

        if samples is None:
            reads = self.reads
        else:
            samples = samples.split(", ")
            reads = list()
            for sample in samples:
                reads.extend(glob.glob(self.path + sample + "*.fastq*"))

        for read in reads:
            dest = raw_reads_dir + re.sub('S\d+_L\d+_R', "", os.path.basename(read))
            dest = dest.replace("_001","")
            print(read)
            print(dest)
            if not os.path.isfile(dest):
                os.symlink(read, dest)
                print("Sym link for", read, "made at", dest)
            else:
                print("Sym link for", read, "already exists at", dest)

        return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_local_reads.py <input> [options]")
    parser.add_argument("input", type=str, help="path to directory containing read files")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir. Default will store output within the "
                                                         "<input> directory")
    parser.add_argument("-link_raw_reads", action='store_true', help="Will link reads from <input> BaseSpace Project "
                                                                     "to <output_dir>/raw_reads/")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.o
    link_raw_reads = args.link_raw_reads

    if len(output_dir) > 0:
        pass
    else:
        output_dir = input_dir

    project = LocalReads(input_dir)

    print(project.name, project.reads, output_dir)

    if link_raw_reads:
        print("linking raw reads. . .")
        project.link_raw_reads(output_dir)




