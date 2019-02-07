#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os,sys, glob
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

import argparse
import shutil
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker
from staphB_ToolKit.lib import mash


class SeqSero:
    #class object to contain fastq file information
    runfiles = None
    #path to fastq files
    path = None
    #output directory
    output_dir = None

    def __init__(self, runfiles=None, path=None, output_dir = ""):
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

        self.seqSero_out_dir = self.output_dir + "/SeqSero_output/"

    def seqSero(self, from_mash=True):
        seqSero_out_dir = self.seqSero_out_dir
        if not os.path.isdir(seqSero_out_dir):
            os.makedirs(seqSero_out_dir)
            print("Directory for seqSero output made: ", seqSero_out_dir)
        mash_species = {}

        if from_mash:
            mash_samples = mash.Mash(path=self.path, output_dir=self.output_dir)
            mash_species = mash_samples.mash_species()

        for read in self.runfiles.reads:
            #get id
            id = self.runfiles.reads[read].id

            seqSero_results = "%s/%s/"%(seqSero_out_dir, id)
            if not os.path.isdir(seqSero_results):
                os.makedirs(seqSero_results)

            if not os.path.isfile(seqSero_results+ "/Seqsero_result.txt"):

                # change self.path to local dir if path is a basemounted dir
                if os.path.isdir(self.path + "/AppResults"):
                    self.path = self.output_dir

                # get paths to fastq files
                if self.runfiles.reads[read].paired:
                    fwd = os.path.abspath(self.runfiles.reads[read].fwd).replace(self.path, "")
                    rev = os.path.abspath(self.runfiles.reads[read].rev).replace(self.path,"")
                else:
                    fastq = os.path.basename(self.runfiles.reads[read].path)

                if "R1" in fwd:
                    reads = fwd.replace("R1", "*")
                else:
                    reads = fwd.replace("_1", "*")

                # create paths for data
                # create paths for data
                mounting = {self.path:'/datain', seqSero_results:'/dataout'}
                out_dir = '/dataout'
                in_dir = '/datain'

                if from_mash:
                    # set expected genome lengths according to mash hits
                    if 'Salmonella' not in mash_species[id]:
                        os.rmdir(seqSero_results)
                        pass
                    else:
                        command = "bash -c 'SeqSero.py -m2 -i {in_dir}/{reads} " \
                                  "-d {out_dir}'".format(in_dir=in_dir,out_dir=out_dir, reads=reads)

                        # call the docker process
                        print("Predicting Salmonella serotype for isolate " + id)
                        calldocker.call("staphb/seqsero:1.0.1", command, '/dataout', mounting)

                else:
                    print("HERE!")
                    # build command for running run_assembly_readMetrics.pl
                    command = "bash -c 'SeqSero.py -m2 -i {in_dir}/{reads} " \
                              "-d {out_dir}'".format(in_dir=in_dir, out_dir=out_dir, reads=reads)

                    # call the docker process
                    print("Predicting Salmonella serotype for isolate " + id)
                    calldocker.call("staphb/seqsero:1.0.1", command, '/dataout', mounting)

            print("SeqSero results for isolate %s saved to: %s%s"%(id,seqSero_out_dir,seqSero_results))


if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argprase.ArgumentTypeError('Boolean value expected.')

    parser = argparse.ArgumentParser(usage="seqSero.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")
    parser.add_argument("-from_mash", nargs='?', type=str2bool, default=True, help="Set expected genome length "
                                                                                   "according to MASH species "
                                                                                   "prediction. default: "
                                                                                   "-from_mash=True")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path = os.path.abspath(args.input)
    output_dir = args.o
    from_mash = args.from_mash

    if not output_dir:
        output_dir = os.getcwd()

    SeqSero_obj = SeqSero(path=path, output_dir=output_dir)
    SeqSero_obj.seqSero(from_mash=from_mash)
