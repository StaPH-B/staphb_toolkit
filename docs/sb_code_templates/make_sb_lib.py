#!/usr/bin/env python3

# This script will take a list of program executables and output a lib file for
# each according to staphb toolkit specifications. A .csv containing a double
# column list of the docker image followed by the executable is required. For example:
# spades,spades.py
# wtdbg,wtdbg.pl
# mash,mash

import sys
import os
import argparse
from string import Template


def get_args():
    parser = argparse.ArgumentParser(description='Make lib file for sb toolkit')
    parser.add_argument("input", help="List of executables")
    parser.add_argument("lib_template", help="staphb toolkit lib template: \
    sb_lib_template.txt")

    return parser.parse_args()

args = get_args()

with open(args.lib_template, "r") as template_file:
    template = Template(template_file.read())

with open(args.input, 'r') as infile:
  for line in infile:
    sline = line.strip().split(",")
    variableMap = {}
    variableMap['exe'] = sline[1]
    variableMap['image'] = sline[0]
    variableMap['class'] = sline[0].capitalize()
    with open("sb_{0}.py".format(variableMap['image']), "w") as lib_file:
      out = template.substitute(variableMap)
      lib_file.write(out)
