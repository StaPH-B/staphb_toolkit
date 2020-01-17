#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import json
import csv
import pathlib
import sys,os
import logging
import getpass
import datetime

sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

# load staphb libaries
from workflows.tredegar.tredegar import tredegar
from core import sb_programs
from core import fileparser


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


def ksnp3_input_file(emm_groups, logger, out_dir):
    for group in emm_groups.keys():
        if len(emm_groups[group]) > 1:
            group_length=len(emm_groups[group])
            isolates_in_group=", ".join(emm_groups[group])
            ksnp3_in_file= os.path.join(out_dir, group, f"{group}_assemblies.txt")
            pathlib.Path(os.path.dirname(ksnp3_in_file)).mkdir(parents=True, exist_ok=True)
            logger.info(f"{group_length} isolates identified as {group}: {isolates_in_group}")
            logger.info(f"Creating ksnp3 input file: {ksnp3_in_file}")
            with open(ksnp3_in_file, "w") as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                for isolate in emm_groups[group]:
                    assembly_path = os.path.join(*["/datain", "shovill_output", isolate, "contigs.fa"])
                    writer.writerow([assembly_path, isolate])
        else:
            logger.info(f"Only one isolate identified as {group}: {''.join(emm_groups[group])}")


def ref_free_snp(output_dir, group, ksnp3_matrix, foushee_config, logger):

    # run the SNP analysis process using ksnp3 if output does not already exist
    if not os.path.isfile(ksnp3_matrix):
        logger.info(f"Performing SNP analysis for isolates identified as {group}")

        # setup mounting in docker container
        ksnp3_mounting = {os.path.abspath(output_dir): '/datain',
                          os.path.join(os.path.abspath(output_dir), "ksnp3_output", group): '/dataout'}

        # generate command to run shovill on the id
        ksnp3_configuration = foushee_config["parameters"]["ksnp3"]
        ksnp3_params = ksnp3_configuration["params"]
        ksnp3_command = f"kSNP3 -in /datain/ksnp3_output/{group}/{group}_assemblies.txt -outdir /dataout/ -k {ksnp3_params['kmer_length']} {ksnp3_params['core_snps_only']} "

        # # generate ksnp3 object and run it
        ksnp3_obj = sb_programs.Run(command=ksnp3_command, path=ksnp3_mounting, image=ksnp3_configuration["image"], tag=ksnp3_configuration["tag"])
        ksnp3_obj.run()


def snp_matrix(output_dir, group, snp_dists_output_dir, snp_dists_result, foushee_config, logger):
     #  run the SNP analysis process using snp_dists if output does not already exist
    if not os.path.isfile(snp_dists_result):
        logger.info(f"Performing SNP analysis for isolates identified as {group}")

        # create snp-dists output dir:
        pathlib.Path(os.path.join(snp_dists_output_dir)).mkdir(parents=True, exist_ok=True)

        # setup mounting in docker container
        snp_dists_mounting = {os.path.join(os.path.abspath(output_dir), "ksnp3_output", group): '/datain', snp_dists_output_dir: '/dataout'}

        # generate command to run shovill on the id
        snp_dists_configuration = foushee_config["parameters"]["snp-dists"]
        snp_dists_params = snp_dists_configuration["params"]
        snp_dists_command = f"bash -c 'snp-dists /datain/core_SNPs_matrix.fasta > /dataout/{group}_pairwise_snp_distance_matrix.tsv {snp_dists_params}'"

        # generate snp_dists object and then run it
        snp_dists_obj = sb_programs.Run(command=snp_dists_command, path=snp_dists_mounting, image=snp_dists_configuration["image"], tag=snp_dists_configuration["tag"])
        snp_dists_obj.run()

################################
# main foushee function
################################
def foushee(memory,cpus,read_file_path,output_dir="",configuration=""):

    #if we don't have an output dir, use the cwd with a tredegar_output dir
    #if we do get the absolute path
    if not output_dir:
        project = f"foushee_run_{datetime.datetime.today().strftime('%Y-%m-%d')}"
        output_dir = os.path.join(os.getcwd(), project)
    else:
        output_dir = os.path.abspath(output_dir)
        project = os.path.basename(output_dir)

    #get the configuration file
    if configuration:
        config_file_path = os.path.abspath(configuration)
    else:
        #use default
        config_file_path = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))),"foushee_config.json")

    with open(config_file_path, 'r') as config_file:
        foushee_config = json.load(config_file)

    # create foushee output dir
    foushee_output = os.path.join(output_dir, "foushee_output")
    pathlib.Path(foushee_output).mkdir(parents=True, exist_ok=True)

    # set logging file
    foushee_log_file = os.path.join(foushee_output, project + "_foushee.log")
    logFormatter = logging.Formatter("%(asctime)s: %(message)s")
    logger = logging.getLogger(__name__)

    fileHandler = logging.FileHandler(foushee_log_file)
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    logger.setLevel(logging.INFO)
    logger.info(
        f"{getpass.getuser()} ran Foushee as: foushee(memory={memory}, cpus={cpus}, read_file_path={read_file_path}, output_dir={output_dir}, configuration={configuration})")

    # process the raw reads
    fastq_files = fileparser.ProcessFastqs(read_file_path, output_dir=output_dir)

    # update configuration object with input files
    foushee_config = fastq_files.inputSubdomain(foushee_config)
    foushee_config["execution_info"]["run_id"] = project
    foushee_config["execution_info"]["user"] = getpass.getuser()
    foushee_config["execution_info"]["datetime"] = datetime.datetime.today().strftime('%Y-%m-%d')

    logger.info("First run tredegar")
    isolate_qual=tredegar(memory=memory,cpus=cpus,read_file_path=read_file_path,output_dir=output_dir,configuration=configuration)

    # create emm_group dict
    emm_groups = group_by_emm(isolate_qual)

    # create ksnp3 directory
    ksnp3_output = os.path.join(output_dir,"ksnp3_output")
    pathlib.Path(ksnp3_output).mkdir(parents=True, exist_ok=True)

    # create input file for kSNP3 analysis
    ksnp3_input_file(emm_groups, logger, out_dir=ksnp3_output)

    # Run ksnp3 & snp_dists for each emm group:
    for group in emm_groups.keys():
        if len(emm_groups[group]) == 1:
            continue

        # run reference free SNP calling with ksnp3
        ksnp3_matrix = os.path.join(*[output_dir, "ksnp3_output", group, "core_SNPs_matrix.fasta"])
        ref_free_snp(output_dir, group, ksnp3_matrix, foushee_config, logger)

        # hard link ksnp3 output into foushee directory
        ksnp3_tree = ksnp3_matrix.replace("core_SNPs_matrix.fasta", "tree.core.tre")
        ksnp3_tree_dest = os.path.join(foushee_output, f"{group}_tree.core.tre")
        if not os.path.isfile(ksnp3_tree_dest):
            os.link(ksnp3_tree, ksnp3_tree_dest)

        # generate snp distance matrix with snp-dists
        snp_dists_output_dir = os.path.join(os.path.abspath(output_dir), "snp_dists_output")
        snp_dists_result = os.path.join(*[snp_dists_output_dir, f"{group}_pairwise_snp_distance_matrix.tsv"])
        snp_matrix(output_dir, group, snp_dists_output_dir, snp_dists_result, foushee_config, logger)

        # hard link snp-dists output into foushee directory
        snp_dists_dest= os.path.join(foushee_output, f"{group}_pairwise_snp_distance_matrix.tsv")
        if not os.path.isfile(snp_dists_dest):
            os.link(snp_dists_result, snp_dists_dest)

        # update config file to include tredegar report
        foushee_config["file_io"]["output_files"]["ksnp_core_tree"] = os.path.join(foushee_output, f"{group}_tree.core.tre")
        foushee_config["file_io"]["output_files"]["pairwise_snps_dists_matrix"] = os.path.join(foushee_output, f"{group}_pairwise_snp_distance_matrix.tsv")
        foushee_config["file_io"]["output_files"]["log_file"] =foushee_log_file

        logger.info(f"Done! Hard links for the ksnp3 core tree and snps-dists distance matrix made at {foushee_output}")










