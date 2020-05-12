#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import os.path
import sys
import json
import re
import pathlib
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

from staphb_toolkit.core import fileparser
from staphb_toolkit.core import sb_programs

class MashSpecies():
    def __init__(self, runfiles=None, path=None, output_dir = None, configuration=None):

        # set configuration file variables
        self.configuration = configuration
        if self.configuration:
            config_file_path = os.path.abspath(configuration)
        else:
            config_file_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4] + "/core/docker_config.json"

        with open(config_file_path, 'r') as config_file:
            self.config = json.load(config_file)

        # create output dir if it doesn't already exist
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        # set runfiles variable using fileparser if runfiles not provided
        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.ProcessFastqs(self.path, output_dir=output_dir)

        # path to untrimmed reads
        self.path = os.path.join(self.path)
        self.mash_out_dir = os.path.join(self.output_dir, "mash_species")

    def create_mash_sketch(self, id, mash_mounting, fwd_read, rev_read, sketch_name):
        if not os.path.isfile(os.path.join(*[self.output_dir, "mash_species", id, sketch_name])):

            # command for creating the mash sketch
            mash_configuration = self.config["parameters"]["mash"]
            mash_sketch_command = f"bash -c 'mkdir -p /dataout/{id} && mash sketch -r -m 2 -o /dataout/{id}/{sketch_name} /datain/{fwd_read} /datain/{rev_read}'"

            # create mash sketch object
            mash_sketch = sb_programs.Run(command=mash_sketch_command, path=mash_mounting, image=mash_configuration["image"], tag=mash_configuration["tag"])
            # run mash sketch
            mash_sketch.run()

    def calc_mash_dist(self, id, mash_mounting, sketch_name, mash_result):
        if not os.path.isfile(os.path.join(*[self.output_dir, "mash_species", id, mash_result])):

            # command for calculating mash distance
            mash_configuration = self.config["parameters"]["mash"]
            mash_dist_command = f"bash -c 'mash dist /db/RefSeqSketchesDefaults.msh /dataout/{id}/{sketch_name} > /dataout/{id}/{mash_result}'"

            # create mash distance object
            mash_dist = sb_programs.Run(command=mash_dist_command, path=mash_mounting, image=mash_configuration["image"], tag=mash_configuration["tag"])
            # run mash distance
            mash_dist.run()

    def run(self):
        # run MASH to get distances
        mash_species = {}
        # dictionary of each set of reads found
        reads_dict = self.runfiles.id_dict()

        print("Running MASH. . .")
        for id in reads_dict:

            # capture read file and path names
            fwd_read = os.path.basename(reads_dict[id].fwd)
            rev_read = os.path.basename(reads_dict[id].fwd)

            # change read dir if hardlinked/copied to input_reads sub dir
            if "input_reads" in self.path:
                reads_dir = os.path.join(self.path, id)
            else:
                reads_dir = self.path

            # create mash output directory
            pathlib.Path(os.path.join(self.output_dir,"mash_species")).mkdir(parents=True, exist_ok=True)

            # docker mounting dictionary
            mash_mounting = {reads_dir: '/datain', os.path.join(self.output_dir,"mash_species"):'/dataout'}

            # command for running mash distance
            sketch_name = id+"_sketch.msh"
            mash_result="{id}_distance.tab".format(id=id)

            # print id, create mash sketch file, and calculate mash distance
            print(id)
            self.create_mash_sketch(id, mash_mounting, fwd_read, rev_read, sketch_name)
            self.calc_mash_dist(id, mash_mounting, sketch_name, mash_result)

            # path to sorted distance file
            mash_result_sorted = os.path.join(*[self.mash_out_dir,id,id+"_sorted_distance.tab"])

            # path to mash hits file
            mash_hits = open(os.path.join(*[self.mash_out_dir,id,mash_result]), 'r').readlines()

            # write out the hits sorted to a file
            with open(mash_result_sorted, 'w') as output:
                for line in sorted(mash_hits, key=lambda line: line.split()[2]):
                    output.write(line)

            # open the sorted hits and get the top hit
            with open(mash_result_sorted, 'r') as file:
                for line in file:
                    top_hit = line

                    # capture only the genus and species of top hit
                    top_hit = re.sub(r'.*-\.-', '', top_hit)
                    top_hit=top_hit.split()
                    top_hit=top_hit[0]
                    top_hit=re.match(r'^[^_]*_[^_]*', top_hit).group(0)
                    top_hit=re.sub(r'.fna', '', top_hit)

                    # Ensure top hit has a definitive species assignment
                    if "_sp." not in top_hit:
                        break
            # specify the top hit as the species for this id
            mash_species[id] = top_hit

        # write out the results to mash_species.csv
        with open(os.path.join(self.mash_out_dir,"mash_species.csv"), 'w') as f:
            f.write("Isolate,Predicted Species\n")
            for key in mash_species.keys():
                f.write("%s,%s\n"%(key,mash_species[key]))

        print(f"sb_mash_species complete! Output saved to {self.mash_out_dir}/mash_species.csv")

        return mash_species
