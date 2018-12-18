#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import shutil
import argparse
import re
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker


class Mash:
    # class object to contain fastq file information
    runfiles = None
    # path to fastq files
    path = None
    # output directory
    output_dir = None

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
            self.runfiles = fileparser.RunFiles(self.path, output_dir=output_dir)

        if db:
            self.db=db
        else:
            self.db = "/db/RefSeqSketchesDefaults.msh"

        self.mash_out_dir = self.output_dir + "/mash_output/"

    def mash(self):
        # create output directory
        mash_out_dir = self.mash_out_dir
        db_name = os.path.basename(self.db)
        if not os.path.isdir(mash_out_dir):
            os.makedirs(mash_out_dir)
            print("Directory for mash output made: ", mash_out_dir)

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            mash_result = "/" + id + "_distance.tab"

            if os.path.isfile(mash_out_dir + mash_result):
                pass
            else:

                # change self.path to local dir if path is a basemounted dir
                if os.path.isdir(self.path + "/AppResults"):
                    self.path = self.output_dir

                # get paths to fastq files
                if self.runfiles.reads[read].paired:
                    fwd = os.path.abspath(self.runfiles.reads[read].fwd).replace(self.path, "")
                    rev = os.path.abspath(self.runfiles.reads[read].rev).replace(self.path,"")
                else:
                    fastq = os.path.basename(self.runfiles.reads[read].path)

                # create paths for data
                mounting = {self.path:'/datain',mash_out_dir:'/dataout'}
                if db is not "RefSeqSketchesDefaults.msh":
                    mounting.update({os.path.dirname(db):'/db'})
                out_dir = '/dataout'
                in_dir = '/datain'

                # build command for creating sketches and generating mash distance table
                # TODO write elif to catch single read data

                if self.runfiles.reads[read].paired:
                    sketch = "bash -c 'cat {in_dir}/{fwd} {in_dir}/{rev} | mash sketch -m 2 - " \
                             "-o {out_dir}/{sketch}'".format(in_dir=in_dir,out_dir=out_dir,sketch=id + "_sketch",
                                                             fwd=fwd,rev=rev)

                    dist = "bash -c 'mash dist /db/{db} {out_dir}/{sketch} > " \
                           "{out_dir}/{mash_result}'".format(in_dir=in_dir,out_dir=out_dir, sketch=id+ "_sketch.msh",
                                                             mash_result=mash_result, db=db_name)

                # call the docker process
                if not os.path.isfile("%s/%s_sketch.sh"%(mash_out_dir, id)):
                    print("Generating MASH sketch for sample " + id)
                    calldocker.call("staphb/mash",sketch,'/dataout',mounting)
                print("Running MASH for sample " + id)
                calldocker.call("staphb/mash", dist, '/dataout',mounting)

    def mash_species(self):
        # capture predicted genus and species
        mash_out_dir = self.mash_out_dir
        mash_species = {}

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            mash_result = mash_out_dir + id + "_distance.tab"
            mash_result_sorted = mash_out_dir + id + "_sorted_distance.tab"

            if not os.path.isfile(mash_result):
                self.mash()

            if not os.path.isfile(mash_result_sorted):
                mash_hits = open(mash_result, 'r').readlines()
                output = open(mash_result_sorted, 'w')
                for line in sorted(mash_hits, key=lambda line: line.split()[2]):
                    output.write(line)
                output.close()

            with open(mash_result_sorted) as file:
                top_hit = file.readline()
                top_hit = re.sub(r'.*-\.-', '', top_hit)
                top_hit=top_hit.split()
                top_hit=top_hit[0]
                top_hit=re.split(r'^([^_]*_[^_]*)(_|\.).*$', top_hit)[1]

            mash_species[id] = top_hit

            print("Predicted species for isolate %s: %s"%(id,top_hit))

        with open(mash_out_dir + "/mash_species.csv", 'w') as f:
            f.write("Isolate,Predicted Species\n")
            for key in mash_species.keys():
                f.write("%s,%s\n"%(key,mash_species[key]))

        return mash_species


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = argparse.ArgumentParser(usage="mash.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")
    parser.add_argument("-species", nargs='?', type=str2bool, default=True, help="return dictionary of predicted "
                                                                                 "species (i.e. genus and species of "
                                                                                 "top Mash hits). default: "
                                                                                 "-species=True")
    parser.add_argument("-db", default="", nargs="?", type=str, help="path to custom mash db")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    species = args.species
    db = args.db

    if not output_dir:
        output_dir = os.getcwd()

    mash_obj = Mash(path=path,output_dir=output_dir,db=db)
    mash_obj.mash()

    if species:
        mash_obj.mash_species()
