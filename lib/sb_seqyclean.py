#!/usr/bin/env python3

import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphb_toolkit.core.sb_libs import SB_lib


class Seqyclean(SB_lib):
    def __init__(self, parameters=None, path=None, executable='seqyclean', docker_image='seqyclean'):
        SB_lib.__init__(self, parameters=parameters, path=path, executable=executable, docker_image=docker_image)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="seqyclean <input> [options]", add_help=False)
    parser.add_argument("seqyclean_parameters", type=str, nargs='?', help="parameters for seqyclean function")
    parser.add_argument("-h", action='store_true')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args, unknown = parser.parse_known_args()

    if args.h:
        seqyclean_obj=Seqyclean(parameters="-h")
        seqyclean_obj.run_lib()
        sys.exit()

    parameters = args.seqyclean_parameters

    if unknown:
        for arg in unknown:
            parameters += f" {arg}"

    seqyclean_obj = Seqyclean(parameters=parameters, path={os.getcwd(): '/data'})
    seqyclean_obj.run_lib()
