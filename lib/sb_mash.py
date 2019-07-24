#!/usr/bin/env python3

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core.sb_libs import SB_lib


class Mash(SB_lib):
    def __init__(self, parameters=None, path=None, executable = 'mash', docker_image='mash'):
        SB_lib.__init__(self,parameters=parameters, path=path, executable=executable, docker_image=docker_image)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_mash.py <input> [options]", add_help=False)
    parser.add_argument("mash_parameters", type=str, nargs='?', help="parameters for mash function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args, unknown = parser.parse_known_args()

    if args.h:
        mash_obj=Mash(parameters="-h")
        mash_obj.run_lib()
        sys.exit()

    parameters = args.mash_parameters

    if unknown:
        for arg in unknown:
            parameters += f" {arg}"

    mash_obj = Mash(parameters=parameters, path={os.getcwd(): '/data'})
    mash_obj.run_lib()
