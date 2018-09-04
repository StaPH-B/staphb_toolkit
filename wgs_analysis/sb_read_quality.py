#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import csv
from staphB_ToolKit.fileHandling import sb_read_data
from staphB_ToolKit.wgs_analysis import sb_taxon

class SBReadQuality:

    def __init__(self, path, output_dir = ""):
        self.name = path
        self.reads = sb_read_data.SBReads(self.name)
        if not output_dir:
            output_dir = self.name
        self.output = output_dir

        if len(self.reads.idList) == 0:
            raise ValueError("No read files in" + " " + self.name)


    def sb_cg_pipeline(self, mash=True):
        cgp_output_dir = self.output + "/cg_pipeline_output/"

        if not os.path.isdir(cgp_output_dir):
            os.makedirs(cgp_output_dir)
            print("Directory for cg_pipeline output made: ", cgp_output_dir)


        taxon = sb_taxon.SBTaxon(self.name, self.output)
        mash_species = taxon.sb_mash()

        #data dictrionary containing quality metrics for each isolate
        #formatted id: {r1_q: X, r2_q: X, est_cvg: X}
        isolate_qual = {}

        for id in self.reads.idList:
            isolate_qual[id] = {"r1_q": None, "r2_q": None, "est_cvg": None}

            if 'Salmonella' in mash_species[id] or 'Escherichia' in mash_species[id]:
                genome_length = 5000000
            elif 'Campylobacter' in mash_species[id]:
                genome_length = 1600000
            elif 'Listeria' in mash_species[id]:
                genome_length = 3000000
            elif 'Vibrio' in mash_species[id]:
                genome_length = 4000000
            else:
                genome_length = input("What is the expected genome size of " + id + "?")
                try:
                    float(genome_length)
                except ValueError:
                   print("A number was not entered")

            cgp_out = cgp_output_dir + id + "_readMetrics.tsv"

            fwd = self.reads.rawdata[id][0]
            rev = self.reads.rawdata[id][1]
            if "R1" in fwd:
                reads = fwd.replace("R1", "*")
            else:
                reads = fwd.replace("_1", "*")

            if not os.path.isfile(cgp_out):
                subprocess.call("run_assembly_readMetrics.pl --fast" + " " + reads + " " + "-e" + " "
                                + str(genome_length) + " > " + cgp_out + " ", shell=True)


            with open(cgp_out, "r") as tsv_file:
                tsv_reader = list(csv.DictReader(tsv_file, delimiter='\t'))

                for line in tsv_reader:
                    if fwd in line["File"]:
                        isolate_qual[id]["r1_q"] = line["avgQuality"]
                        isolate_qual[id]["est_cvg"] = float(line["coverage"])
                    if rev in line["File"]:
                        isolate_qual[id]["r2_q"] = line["avgQuality"]
                        isolate_qual[id]["est_cvg"] += float(line["coverage"])

        print(isolate_qual)
        return(isolate_qual)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_read_quality.py <input> [options]")
    parser.add_argument("input", type=str, help="path to dir containing read files")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir")
    parser.add_argument("-cg_pipeline", action='store_true', help="Perform quality assessment using CG_Pipeline")


    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = args.input
    cg_pipeline = args.cg_pipeline
    output_dir = args.o

    if not output_dir:
        output_dir = path

    quality = SBReadQuality(path, output_dir)
    print("Project selected: " + quality.name)

    if cg_pipeline:
        print(quality.sb_cg_pipeline())
