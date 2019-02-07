#!/usr/bin/env python3

# author: Kevin Libuit
# email: kevin.libuit@dgs.virginia.gov

import os
import sys
import argparse
import psutil
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker


class Shovill:
    # class object to contain fastq file information
    runfiles = None
    # path to fastq files
    path = None
    # output directory
    output_dir = None

    def __init__(self, threads=None, runfiles=None, path=None, output_dir = None, extra_params=None):
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if threads:
            self.threads = threads
        else:
            self.threads = 1

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.RunFiles(self.path, output_dir=output_dir)

        if extra_params:
            self.extra_params = " ".join(extra_params)
        else:
            self.extra_params = ""

        self.shovill_out_dir = self.output_dir + "/shovill_output/"

    def shovill(self):
        # create output directory
        shovill_out_dir = self.shovill_out_dir
        if not os.path.isdir(shovill_out_dir):
            os.makedirs(shovill_out_dir)
            print("Directory for shovill output made: ", shovill_out_dir)
        mem = psutil.virtual_memory()
        mem = str(mem.total >> 30)
        cpu = str(psutil.cpu_count())
        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            shovill_results = "%s/%s/"%(shovill_out_dir, id)
            if not os.path.isdir(shovill_results):
                os.makedirs(shovill_results)

            if not os.path.isfile(shovill_results + "/contigs.fa"):
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
                mounting = {self.path:'/datain',shovill_results:'/dataout'}
                out_dir = '/dataout'
                in_dir = '/datain'

                # build command for creating sketches and generating mash distance table
                # TODO write elif to catch single read data

                if self.runfiles.reads[read].paired:
                    command = "bash -c 'shovill --outdir {out_dir}/ -R1 {in_dir}/{fwd} -R2 {in_dir}/{rev} --ram " \
                              "{mem} --cpus {cpu} {extra_params} --force'".format(in_dir=in_dir, out_dir=out_dir,
                                                                          threads=self.threads,
                                                                          extra_params=self.extra_params, mem=mem,
                                                                          cpu=cpu, fwd=fwd,rev=rev)

                # call the docker process
                print("Generating shovill assembly for sample " + id)
                calldocker.call("staphb/shovill:1.0.4",command,'/dataout',mounting)

            print("Shovill assembly for isolate %s located at %s"%(id,shovill_results))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="mash.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")
    parser.add_argument("-t",default=16,type=int,help="number of threads")
    parser.add_argument("--trim", action='store_true', help="trim adapter sequences")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    threads = args.t
    extra_params = []

    if args.trim:
        extra_params.append("--trim")

    if not output_dir:
        output_dir = os.getcwd()

    shovill_obj = Shovill(path=path,threads=threads,output_dir=output_dir,extra_params=extra_params)
    shovill_obj.shovill()