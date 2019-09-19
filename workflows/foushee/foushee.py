#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import json
import sys
import csv
import datetime
import pathlib
import sys,os
from workflows.tredegar.tredegar import tredegar


sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

#load staphb libaries
from lib import sb_mash_species
from core import sb_programs
from core import fileparser
import xml.etree.ElementTree as ET

def group_by_emm(isolate_qual):
    '''
    create data dictionary of emm groups in the format:
        {emm_typeX: [list of isolates ID'ed as emm_typeX],
        emm_typeY: [list of isolates ID'ed as emm_typeY],
        ...
        }
    '''
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

def ksnp3_input_file(emm_groups, out_dir):
    for group in emm_groups.keys():
        if len(emm_groups[group]) > 1:
            group_length=len(emm_groups[group])
            isolates_in_group=", ".join(emm_groups[group])
            ksnp3_in_file= os.path.join(out_dir, group, f"{group}_assemblies.txt")
            pathlib.Path(os.path.dirname(ksnp3_in_file)).mkdir(parents=True, exist_ok=True)
            print(f"{group_length} isolates identified as {group}: {isolates_in_group}")
            print(f"Creating ksnp3 input file: {ksnp3_in_file}")
            with open(ksnp3_in_file, "w") as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                for isolate in emm_groups[group]:
                    assembly_path = os.path.join(*["/datain", "shovill_output", isolate, "contigs.fa"])
                    writer.writerow([assembly_path, isolate])
        else:
            print("Only one isolate identified as "+ group)


################################
# main foushee function
################################
def foushee(memory,cpus,read_file_path,output_dir="",configuration=""):

    #if we don't have an output dir, use the cwd with a tredegar_output dir
    #if we do get the absolute path
    if not output_dir:
        output_dir = os.path.join(os.getcwd(),"foushee_output")
    else:
        output_dir = os.path.abspath(output_dir)
    print("First run tredegar")
    isolate_qual=tredegar(memory=memory,cpus=cpus,read_file_path=read_file_path,output_dir=output_dir,configuration=configuration)

    #get the configuration file
    if configuration:
        config_file_path = os.path.absolute(configuration)
    else:
        #use default
        config_file_path = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))),"foushee_config.json")

    with open(config_file_path) as config_file:
        ksnp3_config = json.load(config_file)
    ksnp3_params=ksnp3_config["parameters"]["ksnp3"]

    # create foushee output dir
    foushee_output = os.path.join(output_dir, "foushee_output")
    pathlib.Path(foushee_output).mkdir(parents=True, exist_ok=True)

    # create emm_group dict
    emm_groups = group_by_emm(isolate_qual)

    # create ksnp3 directory
    ksnp3_output = os.path.join(output_dir,"ksnp3_output")
    pathlib.Path(ksnp3_output).mkdir(parents=True, exist_ok=True)

    # create input file for kSNP3 analysis
    ksnp3_input_file(emm_groups, out_dir=ksnp3_output)

    # Run ksnp3 & snp_dists for each emm group:
    for group in emm_groups.keys():
        if len(emm_groups[group]) == 1:
            continue

        # setup mounting in docker container
        ksnp3_mounting = {os.path.abspath(output_dir): '/datain', os.path.join(os.path.abspath(output_dir), "ksnp3_output", group): '/dataout'}

        # generate command to run shovill on the id
        ksnp3_command = f"kSNP3 -in /datain/ksnp3_output/{group}/{group}_assemblies.txt -outdir /dataout/ {ksnp3_params}"

        # # generate ksnp3 object
        ksnp3_obj = sb_programs.Run(command=ksnp3_command, path=ksnp3_mounting, docker_image="ksnp3")

        # path for the ksnp3 result file
        ksnp3_result = os.path.join(*[output_dir, "ksnp3_output", group, "core_SNPs_matrix.fasta"])

        #run the SNP analysis process using ksnp3
        if not os.path.isfile(ksnp3_result):
            print(f"Performing SNP analysis for isolates identified as {group}")
            ksnp3_obj.run()
            ksnp3_tree = ksnp3_result.replace("core_SNPS_matrix.fasta", "tree.core.tre")
            ksnp3_tree_dest = os.path.join(*[output_dir, "foushee_output", f"{group}_core_SNPs_matrix.fasta"])
            os.link(ksnp3_tree, ksnp3_tree_dest)

        # create snp-dists output dir:
        snp_dists_output = os.path.join(os.path.abspath(output_dir), "snp_dists_output")
        pathlib.Path(os.path.join(snp_dists_output)).mkdir(parents=True, exist_ok=True)

        # setup mounting in docker container
        snp_dists_mounting = {os.path.join(os.path.abspath(output_dir), "ksnp3_output", group): '/datain', snp_dists_output: '/dataout'}

        # generate command to run shovill on the id
        snp_dists_command = f"bash -c 'snp-dists /datain/core_SNPs_matrix.fasta > /dataout/{group}_snp_distance_matrix.tsv'"

        # # generate snp_dists object
        snp_dists_obj = sb_programs.Run(command=snp_dists_command, path=snp_dists_mounting, docker_image="snp-dists")

        # path for the snp_dists result file
        snp_dists_result = os.path.join(*[snp_dists_output, f"{group}_snp_distance_matrix.tsv"])

        #run the SNP analysis process using snp_dists
        if not os.path.isfile(snp_dists_result):
            print(f"Performing SNP analysis for isolates identified as {group}")
            snp_dists_obj.run()
            os.link(snp_dists_result, snp_dists_result.replace("snp_dists_output", "foushee_output"))

        print(f"Done! Hard links for the ksnp3 core tree and snps-dists distance matrix made at {foushee_output}")










