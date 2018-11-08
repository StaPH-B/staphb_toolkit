#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os,sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

import argparse
import subprocess
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker

class Mash:
    #class object to contain fastq file information
    runfiles = None
    #path to fastq files
    path = None
    #output directory
    output_dir = None

    def __init__(self,threads=None,runfiles=None, path=None, output_dir = ""):
        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.RunFiles(self.path)

        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if threads:
            self.threads = threads
        else:
            self.threads = 1

    def mash(self):
        #create output directory
        mash_out_dir = os.path.join(self.output_dir,"mash_output/")
        if not os.path.isdir(mash_out_dir):
            os.makedirs(mash_out_dir)
            print("Directory for mash output made: ", mash_out_dir)

        for read in self.runfiles.reads:
            #get id
            id = self.runfiles.reads[read].id

            # change self.path to local dir if path is a basemounted dir
            if os.path.isdir(self.path + "/AppResults"):
                self.path = self.output_dir

            #get paths to fastq files
            if self.runfiles.reads[read].paired:
                fwd = os.path.abspath(self.runfiles.reads[read].fwd).replace(self.path, "")
                rev = os.path.basename(self.runfiles.reads[read].rev).replace(self.path,"")
            else:
                fastq = os.path.basename(self.runfiles.reads[read].path)

            #create paths for data
            if self.path == self.output_dir:
                mounting = {self.path:'/data'}
                out_dir = '/data'
                in_dir = '/data'
            else:
                mounting = {self.path:'/datain',self.output_dir:'/dataout'}
                out_dir = '/dataout'
                in_dir = '/datain'

            #build command for creating sketches and generating mash distance table
            #TODO write elif to catch single read data
            mash_result = id + "_distance.tab"
            if self.runfiles.reads[read].paired:
                sketch = "bash -c 'cat {in_dir}/{fwd} {in_dir}/{rev} | mash sketch -m 2 - " \
                         "-o {out_dir}/{sketch}'".format(in_dir=in_dir,out_dir=out_dir,mash_out_dir=mash_out_dir,
                                                         sketch=id + "_sketch",fwd=fwd,rev=rev)
                dist = "bash -c 'mash dist /db/RefSeqSketchesDefaults.msh {in_dir}/{sketch} > " \
                       "{out_dir}/{mash_result}'".format(in_dir=in_dir,out_dir=out_dir,
                                                                        mash_out_dir=mash_out_dir,sketch=id+ "_sketch.msh",
                                                                        mash_result=mash_result)

            #call the docker process
            calldocker.call("staphb/mash",sketch,'/dataout',mounting)
            calldocker.call("staphb/mash", dist, '/dataout',mounting)
            subprocess.Popen(["sort", "-gk3", mash_result, "-o", mash_result])
            subprocess.Popen(["mv", id + "_sketch.msh", mash_result, mash_out_dir])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="mash.py <input> [options]")
    parser.add_argument("input", type=str, help="path to dir containing read files")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir")
    parser.add_argument("-t",default=1,type=int,help="number of threads")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    threads = args.t

    if not output_dir:
        output_dir = os.getcwd()

    mash_obj = Mash(threads=threads,path=path,output_dir=output_dir)
    mash_obj.mash()