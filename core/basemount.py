#!/usr/bin/env python3

import os
import argparse
import glob
import re
import sys

#class for containing all sequencing run information
class Basemount:
    #class object to contain fastq file information
    runfiles = None
    #path to fastq files
    path = None
    #output directory
    output_dir = None

    def __init__(self,runfiles=None, path=None, output_dir = ""):
        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = os.path.abspath(path)
            self.runfiles = fileparser.RunFiles(self.path)

        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if threads:
            self.threads = threads
        else:
            self.threads = 1