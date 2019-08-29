#!/usr/bin/env python3

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphb_toolkit.core.sb_libs import SB_lib


class Serotypefinder(SB_lib):
    def __init__(self, parameters=None, path=None, executable = 'serotypefinder.pl', docker_image='serotypefinder'):
        SB_lib.__init__(self,parameters=parameters, path=path, executable=executable, docker_image=docker_image)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="sb_serotypefinder.py <input> [options]", add_help=False)
    parser.add_argument("serotypefinder_parameters", type=str, nargs='?', help="parameters for serotypefinder function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args, unknown = parser.parse_known_args()

    if args.h:
        stf_obj=Serotypefinder(parameters="-h")
        stf_obj.run_lib()
        sys.exit()

    parameters = args.serotypefinder_parameters

    if unknown:
        for arg in unknown:
            parameters += f" {arg}"

    stf_obj = Serotypefinder(parameters=parameters, path={os.getcwd(): '/data'})
    stf_obj.run_lib()