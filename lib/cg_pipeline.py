#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os,sys
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

import argparse
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker
from staphB_ToolKit.lib import mash

class CGPipeline:
    #class object to contain fastq file information
    runfiles = None
    #path to fastq files
    path = None
    #output directory
    output_dir = None

    def __init__(self, threads=None, runfiles=None, path=None, output_dir = ""):
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

        if threads:
            self.threads = threads
        else:
            self.threads = 1

        self.cg_out_dir = self.output_dir + "/cg_pipeline_output/"

    def read_metrics(self, from_mash=True):
        cg_out_dir = self.cg_out_dir
        if not os.path.isdir(cg_out_dir):
            os.makedirs(cg_out_dir)
            print("Directory for CG Pipeline output made: ", cg_out_dir)

        mash_species = {}

        if from_mash:
            mash_samples = mash.Mash(path=self.path, output_dir=self.output_dir, threads=threads)
            mash_species = mash_samples.mash_species()


        for read in self.runfiles.reads:
            #get id
            id = self.runfiles.reads[read].id
            cgp_result = id + "_readMetrics.tsv"


            if os.path.isfile(cg_out_dir + cgp_result):
                pass
            else:
                # change self.path to local dir if path is a basemounted dir
                if os.path.isdir(self.path + "/AppResults"):
                    self.path = self.output_dir

                #get paths to fastq files
                if self.runfiles.reads[read].paired:
                    fwd = os.path.abspath(self.runfiles.reads[read].fwd).replace(self.path, "")
                else:
                    fastq = os.path.abspath(self.runfiles.reads[read].path).replace(self.path, "")

                if "R1" in fwd:
                    reads = fwd.replace("R1", "*")
                else:
                    reads = fwd.replace("_1", "*")

                #create paths for data
                if self.path == cg_out_dir:
                    mounting = {self.path:'/data'}
                    out_dir = '/data'
                    in_dir = '/data'
                else:
                    mounting = {self.path:'/datain', cg_out_dir:'/dataout'}
                    out_dir = '/dataout'
                    in_dir = '/datain'

                if from_mash:
                    # set expected genome lengths according to mash hits
                        if 'Salmonella' in mash_species[id] or 'Escherichia' in mash_species[id]:
                            genome_length = 5.0
                        elif 'Campylobacter' in mash_species[id]:
                            genome_length = 1.6
                        elif 'Listeria' in mash_species[id]:
                            genome_length = 3.0
                        elif 'Vibrio' in mash_species[id]:
                            genome_length = 4.0
                        else:
                            genome_length = input("In Mbp, what is the expected genome size of %s?"%(id))
                            try:
                                float(genome_length)
                            except ValueError:
                                print("A number was not entered")
                else:
                    genome_length = input("In Mbp, what is the expected genome size of %s?"%(id))
                    try:
                        float(genome_length)
                    except ValueError:
                        print("A number was not entered")

                genome_length = genome_length*1000000

                # build command for running run_assembly_readMetrics.pl
                command = "bash -c 'run_assembly_readMetrics.pl --fast {in_dir}/{reads} -e {genome_length} > " \
                          "{out_dir}/{cgp_result}'".format(in_dir=in_dir,out_dir=out_dir,reads=reads,
                                                           genome_length=genome_length,cgp_result=cgp_result)

                # call the docker process
                print("Getting read metrics for isolate %s"%(id))
                calldocker.call("staphb/lyveset",command,'/dataout',mounting)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="cg_pipeline.py <input> [options]")
    parser.add_argument("input", type=str, help="path to dir containing read files")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir")
    parser.add_argument("-t",default=1,type=int,help="number of threads")
    parser.add_argument("-from_mash", type=bool, default=True, help="Set expected genome length according "
                                                                    "to MASH species prediction. "
                                                                    "default: -from_mash=True")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    threads = args.t
    from_mash = args.from_mash

    if not output_dir:
        output_dir = os.getcwd()

    CGPipeline_obj = CGPipeline(threads=threads, path=path, output_dir=output_dir)
    CGPipeline_obj.read_metrics(from_mash=from_mash)
