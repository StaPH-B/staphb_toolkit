#!/usr/bin/env python3

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core.sb_libs import SB_lib


class IQTree(SB_lib):
    def __init__(self, parameters=None, path=None, executable='iqtree', docker_image='iqtree'):
        SB_lib.__init__(self, parameters=parameters, path=path, executable=executable, docker_image=docker_image)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_iqtree <input> [options]", add_help=False)
    parser.add_argument("iqtree_parameters", type=str, nargs='?', help="parameters for iqtree function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args, unknown = parser.parse_known_args()

    if args.h:
        iqtree_obj=IQTree(parameters="-h")
        iqtree_obj.run_lib()
        sys.exit()

    parameters = args.iqtree_parameters

    if unknown:
        for arg in unknown:
            parameters += f" {arg}"

    iqtree_obj = IQTree(parameters=parameters, path={os.getcwd(): '/data'})
    iqtree_obj.run_lib()
