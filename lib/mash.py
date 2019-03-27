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


class SB_Libs:
    def __init__(self, parameters=None, path=None, docker_image=None):

        self.path=path
        self.docker_image = ""
        self.parameteres = parameters

    def run_lib(self):
        # create paths for data
        mounting = {self.path: '/data'}
        out_dir = '/data'
        in_dir = '/data'

        command = "mash {parameters}".format(parameters=self.parameteres)
        docker_image = self.docker_image

        print(calldocker.call("staphb/{docker_image}".format(docker_image=docker_image), command, '/dataout',mounting))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="mash.py <input> [options]", add_help=False)
    parser.add_argument("mash_parameters", type=str, nargs='?', help="parameters for mash function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    if args.h:
        mash_obj=Mash(parameters="-h")
        mash_obj.mash()
        sys.exit()

    parameters = args.mash_parameters

    mash_obj = Mash(parameters=parameters, path=os.path.getcwd())
    mash_obj.mash()

