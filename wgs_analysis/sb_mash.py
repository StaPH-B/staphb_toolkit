#!/usr/bin/env python3

import os
import glob
import re
import sys
import argparse
import subprocess
from time import gmtime, strftime

class SBMash(object):

    def __init__(self, project):
        self.name = input_dir  # name of wgs project
        self.path = os.path.expanduser("~") + "/" + self.name  # path to project
        self.reads = glob.glob(self.path + "/raw_reads/*.fastq*")  # Read files

        if not os.path.isdir(self.name):
            raise ValueError(self.name + " " + "not found.")

        if len(self.reads) == 0:
            raise ValueError("No read files in" + " " + self.input)

        return

    def mash_raw_reads(self, samples=None):
        mash_output_dir = self.path + "mash_output/"

        if not os.path.isdir(mash_output_dir):
            os.makedirs(mash_output_dir)
            print("Directory for raw reads made: ", mash_output_dir)

        if samples is None:
            reads = self.reads
        else:
            samples = samples.split(", ")
            reads = list()
            for sample in samples:
                reads.extend(glob.glob(self.path + sample + "*.fastq*"))

        r = re.compile(".*_1.fastq")
        fwd_reads = list(filter(r.match, reads))

        rev_reads = list()
        for read in fwd_reads:
            rev_read = read.replace("_1.fastq", "_2.fastq")
            rev_reads.append(rev_read)

        isolates = []

        for read in fwd_reads:
            File = read
            File_split = File.split('/')
            read = File_split[-1]
            read_split = read.split("_")
            isolate = read_split[0]
            isolates.append(isolate)

        count = 0
        while count < len(isolates):
            isolate = isolates[count]
            r1 = fwd_reads[count]
            r2= rev_reads[count]
            cat_reads = mash_output_dir + isolate + ".fastq"
            mash_db = "/opt/genomes/RefSeqSketchesDefaults.msh"
            mash_sketch = mash_output_dir + isolate + ".fastq.msh"
            mash_result = mash_output_dir + isolate + "_distance.tab"

            subprocess.call("cat" + " " + r1 + " " + r2 + " > " + cat_reads +
                            " && mash sketch -m 2 " + cat_reads, shell = True)
            subprocess.call("mash dist" + " " + mash_db + " " + mash_sketch + " > " + mash_result, shell = True)
            subprocess.call("sort -gk3" + " " + mash_result + " " + "-o" + " " + mash_result, shell = True)

            count= count + 1





if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_local_reads.py <input> [options]")
    parser.add_argument("input", type=str, help="path to directory containing read files")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    input_dir = args.input



    project = SBMash(input_dir)

    #print(project.name, project.reads)
    project.mash_raw_reads()


