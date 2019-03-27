#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core.sb_libs import SB_lib


class Spades(SB_lib):
    def __init__(self, parameters=None, path=None, docker_image=None):
        SB_lib.__init__(self,parameters=None, path=None, executable = 'spades.py', docker_image='spades')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_spades.py <input> [options]", add_help=False)
    parser.add_argument("spades_parameters", type=str, nargs='?', help="parameters for spades function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    if args.h:
        spades_obj=Spades(parameters="-h")
        spades_obj.run_lib()
        sys.exit()

    parameters = args.spades_parameters

    spades_obj = Spades(parameters=parameters, path=os.getcwd())
    spades_obj.run_lib()
