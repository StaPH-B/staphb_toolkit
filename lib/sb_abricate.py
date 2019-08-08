#!/usr/bin/env python3

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core.sb_libs import SB_lib


class ABRicate(SB_lib):
    def __init__(self, parameters=None, path=None, executable='abricate', docker_image='abricate'):
        SB_lib.__init__(self, parameters=parameters, path=path, executable=executable, docker_image=docker_image)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_abricate <input> [options]", add_help=False)
    parser.add_argument("abricate_parameters", type=str, nargs='?', help="parameters for ABRicate function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args, unknown = parser.parse_known_args()

    if args.h:
        abricate_obj=ABRicate(parameters="-h")
        abricate_obj.run_lib()
        sys.exit()

    parameters = args.unicycler_parameters

    if unknown:
        for arg in unknown:
            parameters += f" {arg}"

    abricate_obj = ABRicate(parameters=parameters, path={os.getcwd(): '/data'})
    abricate_obj.run_lib()
