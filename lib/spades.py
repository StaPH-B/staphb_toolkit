#!/usr/bin/env python3

# author: Kevin Libuit
# email: kevin.libuit@dgs.virginia.gov

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker


class Spades:
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

        self.spades_out_dir = self.output_dir + "/spades_output/"

    def spades(self):
        # create output directory
        spades_out_dir = self.spades_out_dir
        if not os.path.isdir(spades_out_dir):
            os.makedirs(spades_out_dir)
            print("Directory for spades output made: ", spades_out_dir)

        for read in self.runfiles.reads:
            # get id
            id = self.runfiles.reads[read].id
            spades_results = "%s/%s/"%(spades_out_dir, id)
            if not os.path.isdir(spades_results):
                os.makedirs(spades_results)

            if not os.path.isfile(spades_results + "/contigs.fasta"):
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
                mounting = {self.path:'/datain',spades_results:'/dataout'}
                out_dir = '/dataout'
                in_dir = '/datain'

                # build command for creating sketches and generating mash distance table
                # TODO write elif to catch single read data

                if self.runfiles.reads[read].paired:
                    command = "bash -c 'spades.py -1 {in_dir}/{fwd} -2 {in_dir}/{rev} -o " \
                             "{out_dir}/ -t {threads} {extra_params}'".format(in_dir=in_dir,
                                                                                 spades_results=spades_results,
                                                                                 out_dir=out_dir,
                                                                                 threads=self.threads,
                                                                                 extra_params=self.extra_params,
                                                                                 fwd=fwd,rev=rev)

                # call the docker process
                print("Generating SPAdes assembly for sample " + id)
                calldocker.call("staphb/spades",command,'/dataout',mounting)

            print("SPAdes assembly for isolate %s saved to: %s"%(id,spades_results))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="mash.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")
    parser.add_argument("-t",default=16,type=int,help="number of threads")
    parser.add_argument("--plasmid", action='store_true', help="assemble only plasmids")
    parser.add_argument("--only-assembler", action='store_true', help="Run assembly module only")
    parser.add_argument("--careful", action='store_true', help="Tries to reduce number of mismatches and short "
                                                               "indels.")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    threads = args.t
    extra_params = []

    if args.plasmid:
        extra_params.append(args.plasmid)
    if args.only_assembler:
        extra_params.append(args.only_assembler)
    if args.careful:
        extra_params.append(args.careful)

    if not output_dir:
        output_dir = os.getcwd()

    spades_obj = Spades(path=path,threads=threads,output_dir=output_dir,extra_params=extra_params)
    spades_obj.spades()
