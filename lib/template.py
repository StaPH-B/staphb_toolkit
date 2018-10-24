#!/usr/bin/env python3

#author:
#email:

import os,sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

import argparse
import subprocess
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker

class BP:
    #class object to contain fastq file information
    runfiles = None
    #path to fastq files
    path = None
    #output directory
    output_dir = None

    def __init__(self,threads=None,runfiles=None, path=None, extra_params = "", output_dir = ""):
        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = os.path.abspath(path)
            self.runfiles = fileparser.RunFiles(self.path)

        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if threads:
            self.threads = threads
        else:
            self.threads = 1

        #add quality metrics string if missing set default
        if extra_params:
            self.extra_params = extra_params
        else:
            self.extra_params = ""

    def bp(self):
        for read in self.runfiles.reads:
            #get id
            id = self.runfiles.reads[read].id
            #get paths to fastq files
            if self.runfiles.reads[read].paired:
                fwd = os.path.basename(self.runfiles.reads[read].fwd)
                rev = os.path.basename(self.runfiles.reads[read].rev)
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

            #build command for running trimmomatic
            if self.runfiles.reads[read].paired:
                command = ""
            else:
                command = ""

            #call the docker process
            calldocker.call("staphb/container",command,'/dataout',mounting)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="template.py <input> [options]")
    parser.add_argument("input", type=str, help="path to dir containing read files")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir")
    parser.add_argument("-t",default=1,type=int,help="number of threads")
    parser.add_argument("--params",default="", type=str, help="additional parameters")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    parameters = args.params
    output_dir = args.o
    threads = args.t

    if not output_dir:
        output_dir = path

    template_obj = BP(threads=threads,path=path,extra_params=parameters,output_dir=output_dir)
    template_obj.bp()
