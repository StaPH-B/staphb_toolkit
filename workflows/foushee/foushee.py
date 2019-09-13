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
    print(emm_groups)
    return emm_groups



################################
#main foushee function
################################
def foushee(memory,cpus,read_file_path,output_dir="",configuration=""):
    print("First run tredegar")
    isolate_qual=tredegar(memory=memory,cpus=cpus,read_file_path=read_file_path,output_dir=output_dir,configuration=configuration)

    group_by_emm(isolate_qual)




