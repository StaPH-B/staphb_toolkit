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

from lib import sb_mash
from core import fileparser


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
            self.runfiles = fileparser.RunFiles(self.path, output_dir=output_dir)

        if db:
            self.db=db
        else:
            self.db = "/db/RefSeqSketchesDefaults.msh"

        self.mash_out_dir = self.output_dir + "/mash_output/"


    def main(self):
        #Run MASH to get distances
        mash_species = {}
        for read in self.runfiles.reads:
            id = self.runfiles.reads[read].id

            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir

            fwd = os.path.abspath(self.runfiles.reads[read].fwd).replace(self.path, "")
            rev = os.path.abspath(self.runfiles.reads[read].rev).replace(self.path,"")

            mash_result="{id}_distance.tab".format(id=id)
            pathlib.Path(self.output_dir + "/mash_output/").mkdir(parents=True, exist_ok=True)
            mash_mounting = {self.path: '/datain', self.output_dir + "/mash_output/":'/dataout'}
            out_dir = '/dataout'
            in_dir = '/datain'

            mash_sketch="'mkdir -p {out_dir}/{id}/ && cat {in_dir}/{fwd} {in_dir}/{rev} | mash sketch -r -m 2 " \
                        "-o {out_dir}/{id}/{sketch} -'".format(id=id, in_dir=in_dir, out_dir=out_dir,
                                                               sketch=id+"_sketch", fwd=fwd, rev=rev)
            mash_dist="'mash dist /db/{db} {out_dir}/{id}/{sketch} > {out_dir}/{id}/{mash_result}'".format(
                id=id,in_dir=in_dir,out_dir=out_dir, sketch=id+ "_sketch.msh", mash_result=mash_result,
                db="RefSeqSketchesDefaults.msh")

            if not os.path.isfile(self.output_dir + "/mash_output/" + id + "/" + id+"_sketch.msh"):
                mash_sketch = sb_mash.Mash(executable="bash -c", parameters=mash_sketch, path=mash_mounting)
                mash_sketch.run_lib()
            if not os.path.isfile(self.output_dir + "/mash_output/" + id + "/" + mash_result):
                mash_dist = sb_mash.Mash(executable="bash -c", parameters=mash_dist, path=mash_mounting)
                mash_dist.run_lib()
            mash_result_sorted = self.mash_out_dir + "/" + id + "/" + id + "_sorted_distance.tab"

            mash_hits = open(self.mash_out_dir+"/"+id+"/"+mash_result, 'r').readlines()
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

            with open(self.mash_out_dir + "/mash_species.csv", 'w') as f:
                f.write("Isolate,Predicted Species\n")
                for key in mash_species.keys():
                    f.write("%s,%s\n"%(key,mash_species[key]))

        return mash_species

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_mash_species.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o

    if not output_dir:
        output_dir = os.getcwd()

    mash_species_obj = MashSpecies(path=path,output_dir=output_dir)
    mash_species_obj.main()
