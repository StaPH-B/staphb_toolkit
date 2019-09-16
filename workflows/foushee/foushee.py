#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import json
import sys
import csv
import datetime
import pathlib
import sys,os,re
import core.calldocker as cd
import core.sb_programs as sb_prog
import pprint
import itertools
from operator import itemgetter
from workflows.tredegar.tredegar import tredegar


sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

#load staphb libaries
from lib import sb_mash_species
from core import sb_programs
from core import fileparser
import xml.etree.ElementTree as ET

def group_by_emm(isolate_qual):
    emm_types=[]
    emm_groups={}
    for id in isolate_qual:
        if isolate_qual[id]["species_prediction"] != "Streptococcus_pyogenes":
            continue
        if isolate_qual[id]["subspecies_predictions"] not in emm_types:
            emm_groups.update({isolate_qual[id]["subspecies_predictions"]: [id]})
            emm_types.append((isolate_qual[id]["subspecies_predictions"]))
        else:
            emm_groups[isolate_qual[id]["subspecies_predictions"]].append(id)
    return emm_groups



################################
#main foushee function
################################
def foushee(memory,cpus,read_file_path,output_dir="",configuration=""):
    print("First run tredegar")
    isolate_qual=tredegar(memory=memory,cpus=cpus,read_file_path=read_file_path,output_dir=output_dir,configuration=configuration)

    emm_groups = group_by_emm(isolate_qual)

    for group in emm_groups.keys():
        if len(emm_groups[group]) > 1:
            print(f"{group}_assemblies.txt")
            for isolate in emm_groups[group]:
                assembly_path = os.path.join(*["shovill_output", isolate, "contigs.fa"])
                print(assembly_path)
        else:
            print("Only one isolate identified as "+ group)








