#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import json
import sys
import csv
import datetime
import pathlib
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

#load staphb libaries
from lib import sb_mash_species
from core import sb_programs
from core import fileparser

#define functions for determining ecoli serotype and sal serotype
def ecoli_serotype(output_dir,assembly,id,tredegar_config):
    #output path for serotypefinder
    serotypefinder_output_path = os.path.join(output_dir,"serotypefinder_output")
    pathlib.Path(serotypefinder_output_path).mkdir(parents=True, exist_ok=True)

    #get paratmeters
    stf_params=tredegar_config["parameters"]["serotypefinder"]

    #setup container mounting
    assembly_path = os.path.dirname(assembly)
    stf_mounting = {assembly_path: '/datain',
                   serotypefinder_output_path: '/dataout'}

    #generate serotypefinder command
    assembly_name = os.path.basename(assembly)
    stf_command = f"serotypefinder.pl -d /serotypefinder/database/ -i /datain/{assembly_name} -b /blast-2.2.26/ -o /dataout/{id} {stf_params}"

    #create serotypefinder object
    stf_obj = sb_programs.Run(command=stf_command, path=stf_mounting, docker_image="serotypefinder")

    #path to serotypefinder results file, if it doesn't exist run the serotypefinder
    stf_out = f"{output_dir}/serotypefinder_output/{id}/results_tab.txt"
    if not os.path.isfile(stf_out):
        print("Isolate {id} identified as E.coli. Running SerotypeFinder for serotype prediction".format(id=id))
        stf_obj.run()

    #process the results of the serotype finder file
    with open(stf_out) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        h_type = ''
        o_type = ''
        oh_dict = {}

        #parse each line of the result file and get the H and O type
        for line in tsv_reader:
            if 'Gene' not in line[0]:
                oh_dict[line[5]] = line[0]

        for oh_key in oh_dict:
            if 'H' in oh_key:
                h_type = h_type + oh_key + ' '
            if 'O' in oh_key:
                o_type = o_type + oh_key + ' '
        serotype_result = h_type + ' ' + o_type
        return serotype_result

def salmonella_serotype(output_dir,read_file_path,all_reads,id,tredegar_config):
    #seqsero ouput path
    seqsero2_output_path = os.path.join(output_dir,"seqsero2_output")
    pathlib.Path(seqsero2_output_path).mkdir(parents=True, exist_ok=True)

    #container mounting dictonary
    seqsero2_mounting = {read_file_path: '/datain',seqsero2_output_path: '/dataout'}

    #container command
    seqsero2_command = f"bash -c 'SeqSero2_package.py -i /datain/{all_reads} -t2 -d /dataout/{id}'"

    #generate seqsero object
    seqsero2_obj = sb_programs.Run(command=seqsero2_command, path=seqsero2_mounting, docker_image="seqsero2")

    #path to seqsero results, if it doesn't exist run the seqsero object
    seqsero2_out = f"{output_dir}/seqsero2_output/{id}/Seqsero_result.txt"
    if not os.path.isfile(seqsero2_out):
        print(f"Isolate {id} identified as identified as S.enterica. Running SeqSero for "
              "serotype prediction")
        seqsero2_obj.run()

    #read the result file and return the serotype
    serotype = ""
    with open(seqsero2_out) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        for line in tsv_reader:
            try:
                if "Predicted serotype" in line[0]:
                    serotype = line[1]
            except:
                pass
    return serotype

################################
#main tredegar function
################################
def tredegar(memory,cpus,read_file_path,output_dir="",configuration=""):
    #get the configuration file
    if configuration:
        config_file_path = os.path.absolute(configuration)
    else:
        #use default
        config_file_path = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))),"tredegar_config.json")
    #pull in configuration parameters
    with open(config_file_path) as config_file:
        tredegar_config = json.load(config_file)

    #get the absolute path for the read_file_path
    read_file_path = os.path.abspath(read_file_path)

    #if we don't have an output dir, use the cwd with a tredegar_output dir
    #if we do get the absolute path
    if not output_dir:
        output_dir = os.path.join(os.getcwd(),"tredegar_output")
        project = datetime.datetime.today().strftime('%Y-%m-%d')
    else:
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

    #Run MASH, CG_Pipeline, SeqSero, and SerotypeFinder results
    isolate_qual = {}
    genome_length = ""
    mash_species_obj = sb_mash_species.MashSpecies(path=read_file_path, output_dir=output_dir)

    #if we don't have mash species completed run it, otherwise parse the file and get the results
    mash_species_results = os.path.join(*[output_dir,'mash_output','mash_species.csv'])
    if not os.path.isfile(mash_species_results):
        mash_species = mash_species_obj.run(trim=True)
    else:
        mash_species = {}
        with open(mash_species_results,'r') as csvin:
            reader = csv.reader(csvin,delimiter=',')
            for row in reader:
                mash_species[row[0]] = row[1]

    #TODO add comment here on what this is/does
    matched_wzx = ["O2","O50","O17","O77","O118","O151","O169","O141ab","O141ac"]
    matched_wzy = ["O13","O135","O17","O44","O123","O186"]

    #dictonary of each set of reads found
    reads_dict = fastq_files.id_dict()

    #process each sample
    for id in reads_dict:
        if "seqsero_output" in os.path.abspath(reads_dict[id].fwd):
            continue
        print("Assessing read quality for isolate: " + id)

        #get read names
        # fwd_read = os.path.join("raw_reads", os.path.basename(reads_dict[id].fwd))
        # rev_read = os.path.join("raw_reads", os.path.basename(reads_dict[id].rev))
        # all_reads = os.path.join("raw_reads", reads_dict[id].fwd.replace("1.fastq", "*.fastq"))

        fwd_read = os.path.join("trimmomatic_output", id, f"{id}_1P.fq.gz")
        rev_read = os.path.join("trimmomatic_output", id, f"{id}_2P.fq.gz")
        all_reads = os.path.join("trimmomatic_output", id, f"{id}_*P.fq.gz")

        #initalize result dictonary for this id
        isolate_qual[id] = {"r1_q": None, "r2_q": None, "est_genome_length": None,"est_cvg": None, "predicted_species": None, "predicted_serotype": "NA"}
        isolate_qual[id]["predicted_species"] = mash_species[id]

        #create shovill_output directory
        pathlib.Path(os.path.join(output_dir, "shovill_output")).mkdir(parents=True, exist_ok=True)

        #setup mounting in docker container
        shovill_mounting = {read_file_path: '/datain', os.path.join(output_dir,"shovill_output"):'/dataout'}

        #generate command to run shovill on the id
        shovill_command = f"shovill --outdir /dataout/{id}/ -R1 /datain/{fwd_read} -R2 /datain/{rev_read} --ram {memory} --cpus {cpus} --force"

        #generate shovill object
        shovill_obj = sb_programs.Run(command=shovill_command, path=shovill_mounting, docker_image="shovill")

        #path for the assembly result file
        assembly_result_file_path = os.path.join(*[output_dir,"shovill_output",id,"contigs.fa"])

        #run the assembly process using shovill
        if not os.path.isfile(assembly_result_file_path):
            print("Assemblying {id} with shovill. . .".format(id=id))
            shovill_obj.run()

        #generate the path for quast output
        quast_output_path = os.path.join(output_dir,"quast_output")
        pathlib.Path(quast_output_path).mkdir(parents=True, exist_ok=True)

        #quast mounting dictonary paths
        quast_mounting = {output_dir: '/datain', quast_output_path: '/dataout'}

        #create the quast command
        assembly_file_name = os.path.basename(assembly_result_file_path)
        quast_command = f"bash -c 'quast.py /datain/shovill_output/{id}/{assembly_file_name} -o /dataout/{id}'"

        #create the quast object
        quast_obj = sb_programs.Run(command=quast_command, path=quast_mounting, docker_image="quast")

        #check for quast output, if it doesn't exist run the quast object
        quast_out_file = f"{output_dir}/quast_output/{id}/report.tsv"
        if not os.path.isfile(quast_out_file):
            print(f"Gathering {id} assembly quality metrics with Quast. . .")
            quast_obj.run()

        #open the output file and read in the genome length
        with open(quast_out_file) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for line in tsv_reader:
                if "Total length" in line[0]:
                    genome_length = line[1]
                    isolate_qual[id]["est_genome_length"] = genome_length
            if not genome_length:
                raise ValueError("Unable to predict genome length for isolate" + id)

        #create cg_pipeline output path
        cg_pipeline_output_path = os.path.join(output_dir,"cg_pipeline_output")
        pathlib.Path(cg_pipeline_output_path).mkdir(parents=True, exist_ok=True)
        #generate path mounting for container
        cg_mounting= {read_file_path: '/datain',cg_pipeline_output_path :'/dataout'}

        #generate command for cg_pipeline
        cg_params = tredegar_config["parameters"]["cg_pipeline"]
        cgp_result_file = id + "_readMetrics.tsv"
        cg_command = f"bash -c \'run_assembly_readMetrics.pl {cg_params} /datain/{all_reads} -e {genome_length} > /dataout/{cgp_result_file}\'"

        #generate the cg_pipeline object
        cg_obj= sb_programs.Run(command=cg_command, path=cg_mounting, docker_image="lyveset")

        #check for cg_pipeline output file if not exists run the cg_pipeline object
        cg_out = f"{output_dir}/cg_pipeline_output/{id}_readMetrics.tsv"
        if not os.path.isfile(cg_out):
            print(f"Getting {id} sequencing quality metrics with CG Pipeline")
            cg_obj.run()
            print("CG Pipeline complete.")
        
        #parse cg_pipeline results and store them
        with open(cg_out) as tsv_file:
            tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

            for line in tsv_reader:
                if any(fwd_format in line["File"] for fwd_format in ["_1.fastq", "_R1.fastq", "_1P.fq.gz"]):
                    isolate_qual[id]["r1_q"] = line["avgQuality"]
                    isolate_qual[id]["est_cvg"] = float(line["coverage"])
                if any(rev_format in line["File"] for rev_format in ["_2.fastq", "_R2.fastq", "_2P.fq.gz"]):
                    isolate_qual[id]["r2_q"] = line["avgQuality"]
                    isolate_qual[id]["est_cvg"] += float(line["coverage"])

        #if the predicted species i ecoli run serotype finder
        if "Escherichia_coli" in isolate_qual[id]["predicted_species"]:
            isolate_qual[id]["predicted_serotype"] = ecoli_serotype(output_dir,assembly_result_file_path,id,tredegar_config)

        #if the predicted species is salmonella run seqsero
        if "Salmonella_enterica" in isolate_qual[id]["predicted_species"]:
            isolate_qual[id]["predicted_serotype"] = salmonella_serotype(output_dir,read_file_path,all_reads,id,tredegar_config)

    # Generate the Tredegar report
    report_path = os.path.join(output_dir,"reports")
    pathlib.Path(report_path).mkdir(parents=True, exist_ok=True)
    report_file = os.path.join(report_path, project+"_tredegar_report.tsv")
    column_headers=["sample", "r1_q", "r2_q", "est_genome_length", "est_cvg", "predicted_species", "predicted_serotype"]

    #if we don't have a report, let's write one
    if not os.path.isfile(report_path):
        with open(report_file, "w") as csvfile:
            w = csv.DictWriter(csvfile, column_headers, dialect=csv.excel_tab)
            w.writeheader()
            for key,val in sorted(isolate_qual.items()):
                row = {"sample":key}
                row.update(val)
                w.writerow(row)

    print(f"Tredegar is complete! Output saved as {report_path}")
