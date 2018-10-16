#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from staphB_ToolKit.core import fileparser

class Taxon:
    #class object to contain fastq file information
    runfiles = None
    #path to fastq files
    path = None
    #output directory
    output_dir = None
    #mash db location
    mash_db = "/opt/genomes/RefSeqSketchesDefaults.msh"

    def __init__(self,runfiles=None, path, output_dir = None,mash_db =None):
        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.RunFiles(self.path)

        if output_dir:
            self.output_dir = output_dir
        else:
            self.output_dir = os.getcwd()

        if mash_db:
            self.mash_db = mash_db

    def mash(self):
        #create output directory
        mash_output_dir = os.path.join(self.output_dir,"mash_output")
        if not os.path.isdir(mash_output_dir):
            os.makedirs(mash_output_dir)
            print("Directory for mash output made: ", mash_output_dir)

        #data dictionary for storing mash species prediction
        #formatted id: predicted_species
        mash_species = {}
        # run MASH and store top hit in mash_species
        for id in self.runfile.ids:
            mash_species[id] = ""

            fwd = self.runfile.reads[id].fwd
            rev = self.runfile.reads[id].rev
            cat_reads = mash_output_dir + id + ".fastq"
            mash_sketch = mash_output_dir + id + ".fastq.msh"
            mash_result = mash_output_dir + id + "_distance.tab"

            if not os.path.isfile(cat_reads):
                subprocess.call("cat" + " " + fwd + " " + rev + " > " + cat_reads +
                                " && mash sketch -m 2 " + cat_reads, shell = True)

            if not os.path.isfile(mash_result):
                subprocess.call("mash dist" + " " + mash_db + " " + mash_sketch + " > " + mash_result, shell = True)

            subprocess.call("sort -gk3" + " " + mash_result + " " + "-o" + " " + mash_result, shell = True)

            taxon = subprocess.check_output("head -1 " + mash_result +
                                    "| sed 's/.*-\.-//' | awk -F '\t' '{print $1}' | "
                                    "awk -F '_' '{print $1 \"_\" $2}'", shell=True)

            mash_species[id] = taxon.decode('ASCII').rstrip()

        return(mash_species)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_taxon.py <input> [options]")
    parser.add_argument("input", type=str, help="path to dir containing read files")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir")
    parser.add_argument("-mash", action='store_true', help="Perform taxon prediction using MASH against RefSeq db")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = args.input
    mash = args.mash
    output_dir = args.o

    if not output_dir:
        output_dir = path

    taxons = Taxon(path, output_dir)
    print("Project selected: " + taxons.path)

    if mash:
        print(taxons.runfiles.return_fastq_list())
        print(taxons.mash())
