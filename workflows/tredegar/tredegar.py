#!/usr/bin/env python3

# author: Kevin Libuit
# email: kevin.libuit@dgs.virginia.gov

import os
import yaml
import sys
import csv
import datetime
import pathlib
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

# load staphb libaries
from lib import sb_mash_species
from core import sb_programs
from core import fileparser
import xml.etree.ElementTree as ET


# define function for running seqyclean
def clean_reads(id, output_dir, raw_read_file_path, fwd_read, rev_read, fwd_read_clean, config_file_path, tredegar_config):
    # path for the seqy clean result file
    seqy_clean_result = os.path.join(*[output_dir, "seqyclean_output", id, fwd_read_clean])

    # create and run seqyclean object if it results don't already exist
    if not os.path.isfile(seqy_clean_result):

        # create seqyclean output directory
        pathlib.Path(os.path.join(output_dir, "seqyclean_output")).mkdir(parents=True, exist_ok=True)

        # docker mounting dictionary
        seqyclean_mounting = {raw_read_file_path: '/datain', os.path.join(output_dir, "seqyclean_output"): '/dataout'}

        # command for creating the mash sketch
        seqyclean_params = tredegar_config["parameter_domain"]["seqyclean"]["params"]
        seqyclean_command = f"bash -c 'seqyclean -minlen 25 -qual -c /Adapters_plus_PhiX_174.fasta -1 /datain/{fwd_read} -2 /datain/{rev_read} -o /dataout/{id}/{id}_clean {seqyclean_params}'"

        # generate command to run seqyclean on the id
        seqyclean_obj = sb_programs.Run(command=seqyclean_command, path=seqyclean_mounting, docker_image="seqyclean", configuration=config_file_path)

        print(f"Cleaning read data with seqyclean. . .")
        seqyclean_obj.run()


# define function for running shovill
def assemble_contigs(id, output_dir, clean_read_file_path,fwd_read_clean, rev_read_clean, memory, cpus, assembly, config_file_path, tredegar_config):
    # create and run shovill object if results don't already exist
    if not os.path.isfile(assembly):

        # create shovill_output directory
        pathlib.Path(os.path.join(output_dir, "shovill_output")).mkdir(parents=True, exist_ok=True)

        # setup mounting in docker container
        shovill_mounting = {clean_read_file_path: '/datain', os.path.join(output_dir, "shovill_output"): '/dataout'}

        # generate command to run shovill on the id
        shovill_params = tredegar_config["parameter_domain"]["shovill"]["params"]
        shovill_command = f"shovill --outdir /dataout/{id}/ -R1 /datain/{fwd_read_clean} -R2 /datain/{rev_read_clean} --ram {memory} --cpus {cpus} --force {shovill_params}"

        # generate shovill object
        shovill_obj = sb_programs.Run(command=shovill_command, path=shovill_mounting, docker_image="shovill", configuration=config_file_path)


        print("Assemblying {id} with shovill. . .".format(id=id))
        shovill_obj.run()


def assembly_metrics(id, output_dir, assembly, quast_out_file, isolate_qual, config_file_path, tredegar_config):
    # create and run quast object if results don't already exist
    if not os.path.isfile(quast_out_file):

        # generate the path for quast output
        quast_output_path = os.path.join(output_dir, "quast_output")
        pathlib.Path(quast_output_path).mkdir(parents=True, exist_ok=True)

        # quast mounting dictionary paths
        quast_mounting = {os.path.dirname(assembly): '/datain', quast_output_path: '/dataout'}

        # create the quast command
        assembly_file_name = os.path.basename(assembly)
        quast_params = tredegar_config["parameter_domain"]["quast"]["params"]
        quast_command = f"bash -c 'quast.py /datain/{assembly_file_name} -o /dataout/{id} {quast_params}'"

        # create the quast object
        quast_obj = sb_programs.Run(command=quast_command, path=quast_mounting, docker_image="quast", configuration=config_file_path)

        print(f"Gathering {id} assembly quality metrics with Quast. . .")
        quast_obj.run()

    # open the quast results to capture relevant metrics
    with open(quast_out_file) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        for line in tsv_reader:
            if "Total length" in line[0]:
                genome_length=line[1]
                isolate_qual[id]["est_genome_length"] = genome_length
            if "# contigs" in line[0]:
                number_contigs=line[1]
                isolate_qual[id]["number_contigs"] = number_contigs
        if not genome_length:
            raise ValueError("Unable to predict genome length for isolate" + id)
        if not number_contigs:
            raise ValueError("Unable to predict number of contigs for isolate" + id)

def read_metrics(id, output_dir, raw_read_file_path, all_reads, isolate_qual, cgp_out, config_file_path, tredegar_config):
    # check for cg_pipeline output file if not exists run the cg_pipeline object
    if not os.path.isfile(cgp_out):
        # set  genome length
        genome_length = isolate_qual[id]["est_genome_length"]

        # create cg_pipeline output path
        cg_pipeline_output_path = os.path.join(output_dir, "cg_pipeline_output")
        pathlib.Path(cg_pipeline_output_path).mkdir(parents=True, exist_ok=True)

        # generate path mounting for container
        cg_mounting = {raw_read_file_path: '/datain', cg_pipeline_output_path: '/dataout'}

        # generate command for cg_pipeline
        cgp_params = tredegar_config["parameter_domain"]["cg_pipeline"]["params"]
        cgp_result_file = id + "_readMetrics.tsv"
        cg_command = f"bash -c \'run_assembly_readMetrics.pl {cgp_params} /datain/{all_reads} -e {genome_length} > /dataout/{cgp_result_file}\'"

        # generate the cg_pipeline object
        cg_obj = sb_programs.Run(command=cg_command, path=cg_mounting, docker_image="cg_pipeline", configuration=config_file_path)

        print(f"Getting {id} sequencing quality metrics with CG Pipeline")
        cg_obj.run()

    # open cg_pipeline results and capture relevant metrics
    with open(cgp_out) as tsv_file:
        tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

        for line in tsv_reader:
            if any(fwd_format in line["File"] for fwd_format in ["_1.fastq", "_R1.fastq", "_1P.fq.gz"]):
                isolate_qual[id]["r1_q"] = line["avgQuality"]
                isolate_qual[id]["est_cvg"] = float(line["coverage"])
            if any(rev_format in line["File"] for rev_format in ["_2.fastq", "_R2.fastq", "_2P.fq.gz"]):
                isolate_qual[id]["r2_q"] = line["avgQuality"]
                isolate_qual[id]["est_cvg"] += float(line["coverage"])

#define functions for determining ecoli serotype and sal serotype
def ecoli_serotype(output_dir,assembly,id,config_file_path, tredegar_config):
    # ambiguous allele calls
    matched_wzx = ["O2","O50","O17","O77","O118","O151","O169","O141ab","O141ac"]
    matched_wzy = ["O13","O135","O17","O44","O123","O186"]

    #path to serotypefinder results file, if it doesn't exist run the serotypefinder
    stf_out = f"{output_dir}/serotypefinder_output/{id}/results_tab.txt"
    if not os.path.isfile(stf_out):
        # output path for serotypefinder
        serotypefinder_output_path = os.path.join(output_dir, "serotypefinder_output")
        pathlib.Path(serotypefinder_output_path).mkdir(parents=True, exist_ok=True)

        # setup container mounting
        assembly_path = os.path.dirname(assembly)
        stf_mounting = {assembly_path: '/datain', serotypefinder_output_path: '/dataout'}

        # generate serotypefinder command
        assembly_name = os.path.basename(assembly)
        stf_params = tredegar_config["parameter_domain"]["serotypefinder"]["params"]
        stf_command = f"serotypefinder.pl -d /serotypefinder/database/ -i /datain/{assembly_name} -b /blast-2.2.26/ -o /dataout/{id} {stf_params}"

        # create serotypefinder object
        stf_obj = sb_programs.Run(command=stf_command, path=stf_mounting, docker_image="serotypefinder", configuration=config_file_path)
        print("Isolate {id} identified as E.coli. Running SerotypeFinder for serotype prediction".format(id=id))
        stf_obj.run()

    # process the results of serotypefinder as per literature guidelines (Joensen, et al. 2015, DOI: 10.1128/JCM.00008-15)
    with open(stf_out) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        wzx_allele = ""
        wzy_allele = ""
        wzm_allele = ""
        h_type = ""

        for line in tsv_reader:
            if "fl" in line [0]:
                h_type = line[5]

            if line[0] == "wzx":
                wzx_allele = line[5]
            if line[0] == "wzy":
                wzy_allele = line[5]
            if line[0] == "wzm":
                wzm_allele = line[5]

        o_type = wzx_allele
        if not wzx_allele:
            o_type = wzy_allele
        if not wzx_allele and not wzy_allele:
            o_type = wzm_allele

        if o_type in matched_wzx:
            o_type = wzy_allele
        if o_type in matched_wzy:
            o_type = wzx_allele
        serotype = f"{o_type}:{h_type}"

        # NA if no o-type or h-type identified
        if serotype == ":":
            serotype = "NA"

    return serotype

def salmonella_serotype(output_dir,raw_read_file_path,all_reads,id, config_file_path, tredegar_config):
    #path to seqsero results, if it doesn't exist run the seqsero object
    seqsero_out = f"{output_dir}/seqsero_output/{id}/Seqsero_result.txt"
    if not os.path.isfile(seqsero_out):
        # seqsero ouput path
        seqsero_output_path = os.path.join(output_dir, "seqsero_output")
        pathlib.Path(seqsero_output_path).mkdir(parents=True, exist_ok=True)

        # container mounting dictionary
        seqsero_mounting = {raw_read_file_path: '/datain', seqsero_output_path: '/dataout'}

        # container command
        seqsero_params = tredegar_config["parameter_domain"]["seqsero"]["params"]
        seqsero_command = f"bash -c 'SeqSero.py -m2 -i /datain/{all_reads} -d /dataout/{id} {seqsero_params}'"

        # generate seqsero object
        seqsero_obj = sb_programs.Run(command=seqsero_command, path=seqsero_mounting, docker_image="seqsero", configuration=config_file_path)

        print(f"Isolate {id} identified as identified as S.enterica. Running SeqSero for serotype prediction")
        seqsero_obj.run()

    #read the result file and return the serotype
    serotype = ""
    with open(seqsero_out) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        for line in tsv_reader:
            try:
                if "Predicted serotype" in line[0]:
                    serotype = line[1]
            except:
                pass
    return serotype

def gas_emmtype(output_dir,raw_read_file_path,id,fwd,rev, config_file_path, tredegar_config):
    #path to seqsero results, if it doesn't exist run the seqsero object
    emmtyper_out = f"{emmtyper_output_path}/{id}/{id}_1.results.xml"
    if not os.path.isfile(emmtyper_out):
        # seqsero ouput path
        emmtyper_output_path = os.path.join(output_dir, "emmtyper_output")
        pathlib.Path(emmtyper_output_path).mkdir(parents=True, exist_ok=True)
        # container mounting dictionary
        emmtyper_mounting = {raw_read_file_path: '/datain', emmtyper_output_path: '/dataout'}

        # container command
        emmtyper_params = tredegar_config["parameter_domain"]["emm-typing-tool"]["params"]
        emmtyper_command = f"emm_typing.py -1 /datain/{fwd} -2 /datain/{rev} -m /db/ -o /dataout/{id}/ {emmtyper_params}"

        # generate seqsero object
        emmtyper_obj = sb_programs.Run(command=emmtyper_command, path=emmtyper_mounting, docker_image="emm-typing-tool", configuration=config_file_path)

        print(f"Isolate {id} identified as identified as Streptococcus_pyogenes. Running emm-typing-tool for "
              "emm-type prediction")
        emmtyper_obj.run()

    #read the result file and return the serotype
    emm_type=""
    tree=ET.parse(emmtyper_out)
    root = tree.getroot()
    for result in root[1].findall("result"):
        if result.attrib['type'] == 'Final_EMM_type':

            emm_type=(result.attrib['value'])
    return emm_type.split(".")[0]

################################
#main tredegar function
################################
def tredegar(memory,cpus,read_file_path,output_dir="",configuration=""):
    # get the configuration file
    if configuration:
        config_file_path = os.path.abspath(configuration)
    else:
        # use default
        config_file_path = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))), "tredegar_config.yaml")
    # pull in configuration parameters
    with open(config_file_path, 'r') as config_file:
        tredegar_config = yaml.safe_load(config_file)

    # get the absolute path for the read_file_path
    read_file_path = os.path.abspath(read_file_path)

    # if we don't have an output dir, use the cwd with a tredegar_output dir
    if not output_dir:
        output_dir = os.path.join(os.getcwd(),"tredegar_output")
        project = datetime.datetime.today().strftime('%Y-%m-%d')
    else:
        # if we do, get the absolute path
        output_dir = os.path.abspath(output_dir)
        project = os.path.basename(output_dir)

    #process the raw reads
    fastq_files = fileparser.ProcessFastqs(read_file_path, output_dir=output_dir)

    #Set the read_file_path equal to the output_dir since reads have been copied/hard linked there
    if os.path.isdir(os.path.join(read_file_path,"AppResults")):
        read_file_path = output_dir
    else:
        fastq_files.link_reads(output_dir=output_dir)
        read_file_path=output_dir

    # Path to untrimmed reads
    read_file_path = os.path.join(read_file_path, "raw_reads")

    #Run MASH, CG_Pipeline, SeqSero, and SerotypeFinder results
    isolate_qual = {}
    genome_length = ""
    mash_species_obj = sb_mash_species.MashSpecies(path=read_file_path, output_dir=output_dir, configuration=config_file_path)

    #if we don't have mash species completed run it, otherwise parse the file and get the results
    mash_species_results = os.path.join(*[output_dir,'mash_output','mash_species.csv'])
    if not os.path.isfile(mash_species_results):
        mash_species = mash_species_obj.run()

    else:
        mash_species = {}
        with open(mash_species_results,'r') as csvin:
            reader = csv.reader(csvin,delimiter=',')
            for row in reader:
                mash_species[row[0]] = row[1]

    #dictionary of each set of reads found
    reads_dict = fastq_files.id_dict()

    #process each sample
    for id in reads_dict:
        if "seqsero_output" in os.path.abspath(reads_dict[id].fwd):
            continue
        print("Assessing read quality for isolate: " + id)

        #get read names
        fwd_read = os.path.basename(reads_dict[id].fwd)
        fwd_read_clean = fwd_read.replace("_1.fastq.gz", "_clean_PE1.fastq")
        rev_read = os.path.basename(reads_dict[id].rev)
        rev_read_clean = rev_read.replace("_2.fastq.gz", "_clean_PE2.fastq")
        all_reads = fwd_read.replace("_1.fastq", "*.fastq")

        # Change read dir since reads hardlinked/copied to an isolate sub dir
        raw_read_file_path = os.path.join(read_file_path, id)
        clean_read_file_path = os.path.join(read_file_path.replace("raw_reads","seqyclean_output"), id)

        #initalize result dictionary for this id
        isolate_qual[id] = {"r1_q": None, "r2_q": None, "est_genome_length": None,"est_cvg": None, "number_contigs": None, "species_prediction": None, "subspecies_predictions": "NA"}
        isolate_qual[id]["species_prediction"] = mash_species[id]

        # Clean read data with SeqyClean before assembling
        clean_reads(id, output_dir, raw_read_file_path, fwd_read, rev_read, fwd_read_clean, config_file_path, tredegar_config)

        # path for the assembly result file
        assembly = os.path.join(*[output_dir, "shovill_output", id, "contigs.fa"])

        # Assemble contigs using cleaned read data
        assemble_contigs(id, output_dir, clean_read_file_path, fwd_read_clean, rev_read_clean, memory, cpus, assembly, config_file_path, tredegar_config)

        # Get assembly metrics using quast
        quast_out_file = f"{output_dir}/quast_output/{id}/report.tsv"
        assembly_metrics(id, output_dir, assembly, quast_out_file, isolate_qual, config_file_path, tredegar_config)

        # Get read metrics using CG Pipeline
        cgp_out = f"{output_dir}/cg_pipeline_output/{id}_readMetrics.tsv"
        read_metrics(id, output_dir, raw_read_file_path, all_reads, isolate_qual, cgp_out, config_file_path, tredegar_config)

        # if the predicted species is ecoli run serotype finder
        if "Escherichia_coli" in isolate_qual[id]["species_prediction"]:
            isolate_qual[id]["subspecies_predictions"] = ecoli_serotype(output_dir,assembly,id, config_file_path, tredegar_config)

        # if the predicted species is salmonella enterica run seqsero
        if "Salmonella_enterica" in isolate_qual[id]["species_prediction"]:
            isolate_qual[id]["subspecies_predictions"] = salmonella_serotype(output_dir,raw_read_file_path,all_reads,id, config_file_path, tredegar_config)

        # if the predicted species is streptococcus pyogenes run seqsero
        if "Streptococcus_pyogenes" in isolate_qual[id]["species_prediction"]:
            isolate_qual[id]["subspecies_predictions"] = gas_emmtype(output_dir,raw_read_file_path,id,fwd_read,rev_read, tredegar_config)

    # Generate the Tredegar report
    tredegar_output = os.path.join(output_dir,"tredegar_output")
    pathlib.Path(tredegar_output).mkdir(parents=True, exist_ok=True)
    report_file = os.path.join(tredegar_output, project+"_tredegar_report.tsv")
    column_headers=["sample", "r1_q", "r2_q", "est_genome_length", "est_cvg", "number_contigs", "species_prediction", "subspecies_predictions"]

    # If we don't have a report, write one
    if not os.path.isfile(tredegar_output):
        with open(report_file, "w") as csvfile:
            w = csv.DictWriter(csvfile, column_headers, dialect=csv.excel_tab)
            w.writeheader()
            for key,val in sorted(isolate_qual.items()):
                row = {"sample":key}
                row.update(val)
                w.writerow(row)

    print(f"Tredegar is complete! Output saved to {tredegar_output}")
    return isolate_qual
