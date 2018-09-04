#!/usr/bin/env python3

import os
import glob
import re
import sys
import argparse
import subprocess
from time import gmtime, strftime

class SBTaxon:

    def __init__(self, project_dir, read_dir = None):
        self.name = project_dir  # name of wgs project
        self.path = os.path.expanduser("~") + "/" + self.name  # path to project
        if read_dir is None:
            read_dir = "/raw_reads"
        self.reads = glob.glob(self.path + read_dir + "/*.fastq*")  # Read files
        test = glob.glob("")

        if not os.path.isdir(self.name):
            raise ValueError(self.name + " " + "not found.")

        if len(self.reads) == 0:
            print(self.path + read_dir)
            print(self.reads)
            raise ValueError("No read files in" + " " + self.name)

        return

    def sb_mash(self, samples=None):
        mash_output_dir = self.path + "/mash_output/"

        if not os.path.isdir(mash_output_dir):
            os.makedirs(mash_output_dir)
            print("Directory for mash output made: ", mash_output_dir)

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
    parser = argparse.ArgumentParser(usage="sb_taxon.py <input> [options]")
    parser.add_argument("project", type=str, help="path to Project containing read files")
    parser.add_argument("-reads_dir", default="", type=str, help="path to read files; default = project_dir/raw_read")
    parser.add_argument("-mash", action='store_true', help="Perform taxon prediction using MASH against RefSeq db")


    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    project = args.project
    mash = args.mash

    reads_dir = args.reads_dir
    if len(reads_dir) > 0:
        pass
    else:
        reads_dir =  "/raw_reads/"

    project = SBTaxon(project, reads_dir)
    print("Project selected: " + project.name)


    #print(project.name, project.reads)
    if mash:
        project.sb_mash()

