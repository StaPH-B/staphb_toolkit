#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from staphB_ToolKit.fileHandling import sb_read_data

class SBTaxon:

    def __init__(self, path, output_dir = ""):
        self.name = path
        self.reads = sb_read_data.SBReads(self.name)
        if not output_dir:
            output_dir = os.getcwd()
        self.output = output_dir

        if len(self.reads.idList) == 0:
            raise ValueError("No read files in" + " " + self.name)


    def sb_mash(self):
        mash_output_dir = self.output + "/mash_output/"
        mash_db = "/opt/genomes/RefSeqSketchesDefaults.msh"

        if not os.path.isdir(mash_output_dir):
            os.makedirs(mash_output_dir)
            print("Directory for mash output made: ", mash_output_dir)

        #data dictionary for storing mash species prediction
        #formatted id: predicted_species
        mash_species = {}

        # run MASH and store top hit in mash_species
        for id in self.reads.idList:
            mash_species[id] = ""

            r1 = self.reads.rawdata[id][0]
            r2 = self.reads.rawdata[id][1]

            cat_reads = mash_output_dir + id + ".fastq"
            mash_sketch = mash_output_dir + id + ".fastq.msh"
            mash_result = mash_output_dir + id + "_distance.tab"

            if not cat_reads:
                subprocess.call("cat" + " " + r1 + " " + r2 + " > " + cat_reads +
                                " && mash sketch -m 2 " + cat_reads, shell = True)

            if not mash_result:
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

    reads = SBTaxon(path, output_dir)
    print("Project selected: " + reads.name)

    if mash:
        print(reads.sb_mash())