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


class SB_lib:
    # path to fastq files
    path = None
    # output directory
    output_dir = None

    def __init__(self, output_dir = None, params=None):
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if params:
            self.params = params

    def run_docker(self):
        # create paths for data
        mounting = {self.path: '/data'}
        out_dir = '/data'
        in_dir = '/data'

        command = ""

        calldocker.call("staphb/{image}".format(image=image), command, '/dataout',mounting)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="mash.py <input> [options]")
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

    lib_obj = SB_lib(path=path,output_dir=output_dir)
    lib_obj.run_docker()
