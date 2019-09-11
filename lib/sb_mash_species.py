#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import argparse
import re
import pathlib
import datetime
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

from core import fileparser
from core import sb_programs


class MashSpecies():
    def __init__(self, runfiles=None, path=None, output_dir = None, db=None):
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.ProcessFastqs(self.path, output_dir=output_dir)

        if db:
            self.db=db
        else:
            self.db = "/db/RefSeqSketchesDefaults.msh"

        self.mash_out_dir = os.path.join(self.output_dir,"mash_output")

    def trim_reads(self):
        reads_dict = self.runfiles.id_dict()
        pathlib.Path(os.path.join(self.output_dir, "trimmomatic_output")).mkdir(parents=True, exist_ok=True)

        for id in reads_dict:
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir
                fwd_path = os.path.join("raw_reads", os.path.basename(reads_dict[id].fwd))
                rev_path = os.path.join("raw_reads", os.path.basename(reads_dict[id].rev))
            else:
                self.runfiles.link_reads(output_dir=self.output_dir)
                fwd_path = os.path.join(os.path.basename(reads_dict[id].fwd))
                rev_path = os.path.join(os.path.basename(reads_dict[id].rev))

            #create trimmomatic output directory
            pathlib.Path(os.path.join(self.output_dir,"trimmomatic_output")).mkdir(parents=True, exist_ok=True)

            #docker mounting dictonary
            trimmomatic_mounting = {self.path: '/datain', os.path.join(self.output_dir,"trimmomatic_output"):'/dataout'}

            #command for creating the mash sketch
            trimmomatic_command = f"bash -c 'mkdir -p /dataout/{id}; cd /datain/ && trimmomatic PE /datain/{fwd_path} /datain/{rev_path} -baseout {id}.fq.gz SLIDINGWINDOW:4:30; mv /datain/{id}*.fq.gz /dataout/{id}'"
            #create and run mash sketch object if it doesn't already exist
            if not os.path.isfile(os.path.join(*[self.output_dir,"trimmomatic_output",id, id+"_2P.fq.gz"])):
                #create trimmomatic object
                trimmomatic_obj = sb_programs.Run(command=trimmomatic_command, path=trimmomatic_mounting, docker_image="trimmomatic")
                #run trimmomatic
                trimmomatic_obj.run()

    def run(self, trim=False):
        #Run MASH to get distances
        mash_species = {}
        #dictonary of each set of reads found
        reads_dict = self.runfiles.id_dict()

        if trim:
            print("Trimming reads before running MASH. . .")
            self.trim_reads()

        print("Performing taxonomic predictions with MASH. . .")
        for id in reads_dict:
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir
                reads_dir = "raw_reads"
            else:
                self.runfiles.link_reads(output_dir=self.output_dir)
                reads_dir=self.path

            #Capture read file and path names
            fwd_read = os.path.basename(reads_dict[id].fwd)
            rev_read = os.path.basename(reads_dict[id].fwd)
            read_mount=self.path

            if trim:
                fwd_read = f"{id}_1P.fq.gz"
                rev_read = f"{id}_2P.fq.gz"
                reads_dir = os.path.join("trimmomatic_output", id)
                read_mount=self.output_dir

            fwd_path = os.path.join(reads_dir, fwd_read)
            rev_path = os.path.join(reads_dir, rev_read)

            #create mash output directory
            pathlib.Path(os.path.join(self.output_dir,"mash_output")).mkdir(parents=True, exist_ok=True)

            #docker mounting dictonary
            mash_mounting = {read_mount: '/datain', os.path.join(self.output_dir,"mash_output"):'/dataout'}

            #mash result and sketch file name
            mash_result="{id}_distance.tab".format(id=id)
            sketch_name = id+"_sketch"

            #command for creating the mash sketch
            mash_sketch_command = f"bash -c 'mkdir -p /dataout/{id} && mash sketch -r -m 2 -o /dataout/{id}/{sketch_name} /datain/{fwd_path} /datain/{rev_path}'"

            #command for running mash distance
            sketch_name = id+"_sketch.msh"
            db = "RefSeqSketchesDefaults.msh"
            mash_dist_command = f"bash -c 'mash dist /db/{db} /dataout/{id}/{sketch_name} > /dataout/{id}/{mash_result}'"

            print(id)
            #create and run mash sketch object if it doesn't already exist
            if not os.path.isfile(os.path.join(*[self.output_dir,"mash_output",id, sketch_name])):
                #create mash sketch object
                mash_sketch = sb_programs.Run(command=mash_sketch_command, path=mash_mounting, docker_image="mash")
                #run mash sketch
                mash_sketch.run()

            #create and run mash dist object if it doesn't already exist
            if not os.path.isfile(os.path.join(*[self.output_dir,"mash_output",id,mash_result])):
                #create mash distance object
                mash_dist = sb_programs.Run(command=mash_dist_command, path=mash_mounting, docker_image="mash")
                #run mash distance
                mash_dist.run()

            #path to sorted distance file
            mash_result_sorted = os.path.join(*[self.mash_out_dir,id,id+"_sorted_distance.tab"])

            #path to mash hits file
            mash_hits = open(os.path.join(*[self.mash_out_dir,id,mash_result]), 'r').readlines()

            #write out the hits sorted to a file
            with open(mash_result_sorted, 'w') as output:
                for line in sorted(mash_hits, key=lambda line: line.split()[2]):
                    output.write(line)

            #open the sorted hits and get the top hit
            with open(mash_result_sorted) as file:
                top_hit = file.readline()
                top_hit = re.sub(r'.*-\.-', '', top_hit)
                top_hit=top_hit.split()
                top_hit=top_hit[0]
                top_hit=re.split(r'^([^_]*_[^_]*)(_|\.).*$', top_hit)[1]


            #specify the top hit as the species for this id
            mash_species[id] = top_hit

        #write out the results to mash_species.csv
        with open(os.path.join(self.mash_out_dir,"mash_species.csv"), 'w') as f:
            f.write("Isolate,Predicted Species\n")
            for key in mash_species.keys():
                f.write("%s,%s\n"%(key,mash_species[key]))

        return mash_species
