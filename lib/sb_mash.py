#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core.sb_libs import SB_lib


class Mash(SB_lib):
    def __init__(self, parameters=None, path=None, docker_image=None):
        SB_lib.__init__(self,parameters=None, path=None, executable = 'mash', docker_image='mash')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_mash.py <input> [options]", add_help=False)
    parser.add_argument("mash_parameters", type=str, nargs='?', help="parameters for mash function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    if args.h:
        mash_obj=Mash(parameters="-h")
        mash_obj.run_lib()
        sys.exit()

    parameters = args.mash_parameters

    mash_obj = Mash(parameters=parameters, path=os.getcwd())
    mash_obj.run_lib()

