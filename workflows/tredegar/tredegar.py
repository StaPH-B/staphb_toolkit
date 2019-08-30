#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import json
import sys
import argparse
import csv
import psutil
import datetime
import pathlib
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

from lib import sb_mash_species
from core import sb_programs
from core import fileparser

def main():
    parser = argparse.ArgumentParser(usage="tredegar.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path=args.input
    output_dir = args.o

    if not output_dir:
        output_dir = os.path.abspath("tredegar_output")

    if output_dir.endswith('/'):
        output_dir = output_dir[:-1]

    if '/' in output_dir and "tredegar_output" not in output_dir:
        project = output_dir.split('/')
        project = project[-1]
    elif "tredegar_output" in output_dir:
        project = datetime.datetime.today().strftime('%Y-%m-%d')
    else:
        project = output_dir

    runfiles = fileparser.RunFiles(path, output_dir=output_dir)
    if os.path.isdir(path + "/AppResults"):
        path = output_dir

    mem = psutil.virtual_memory()
    mem = str(mem.total >> 30)
    cpu = str(psutil.cpu_count())

    #Run MASH, CG_Pipeline, SeqSero, and SerotypeFinder results
    isolate_qual = {}
    genome_length = ""
    mash_species_obj = sb_mash_species.MashSpecies(path=os.path.abspath(path), output_dir=output_dir)
    mash_species = mash_species_obj.main()

    matched_wzx = ["O2","O50","O17","O77","O118","O151","O169","O141ab","O141ac"]
    matched_wzy = ["O13","O135","O17","O44","O123","O186"]

    config_file = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + "/tredegar_config.json"

    with open(config_file) as config_file:
        tredegar_config = json.load(config_file)


    for read in runfiles.reads:
        if "seqsero_output" in os.path.abspath(runfiles.reads[read].fwd):
            continue

        print("Analyzing isolate: " + read)
        id = runfiles.reads[read].id

        if os.path.isdir(path + "/AppResults"):
            path = output_dir

        fwd = os.path.abspath(runfiles.reads[read].fwd).replace(os.path.abspath(path), "")
        rev = os.path.abspath(runfiles.reads[read].rev).replace(os.path.abspath(path),"")
        if "R1.fastq" in fwd:
            all_reads = fwd.replace("R1.fastq", "*.fastq")
        else:
            all_reads = fwd.replace("_1.fastq", "*.fastq")

        out_dir = '/dataout'
        in_dir = '/datain'
        isolate_qual[id] = {"r1_q": None, "r2_q": None, "est_genome_length": None,
                            "est_cvg": None, "predicted_species": None, "predicted_serotype": "NA"}

        isolate_qual[id]["predicted_species"] = mash_species[id]

        pathlib.Path(output_dir + "/shovill_output/").mkdir(parents=True, exist_ok=True)
        shovill_mounting = {os.path.abspath(path): '/datain', os.path.abspath(output_dir) + "/shovill_output/":'/dataout'}
        shovill_command = "shovill --outdir {out_dir}/{id}/ -R1 {in_dir}/{fwd} -R2 {in_dir}/{rev} " \
                         "--ram {mem} --cpus {cpu} --force".format(
            in_dir=in_dir,out_dir=out_dir,fwd=fwd,rev=rev,mem=mem,cpu=cpu,id=id)

        shovill_obj = sb_programs.Run(command=shovill_command, path=shovill_mounting, docker_image="shovill")
        assembly = "/shovill_output/" + id + "/" + "contigs.fa"

        if not os.path.isfile(output_dir + assembly):
            print("Assemblying {id} with shovill. . .".format(id=id))
            shovill_obj.run()

        pathlib.Path(output_dir + "/quast_output/").mkdir(parents=True, exist_ok=True)
        quast_mounting = {os.path.abspath(output_dir): '/datain', os.path.abspath(output_dir) + "/quast_output/": '/dataout'}
        quast_command = "bash -c 'quast.py {in_dir}/{assembly} -o {out_dir}/{id}'".format(
            assembly=assembly, id=id, out_dir=out_dir, in_dir=in_dir)

        quast_obj = sb_programs.Run(command=quast_command, path=quast_mounting, docker_image="quast")
        quast_out = "%s/quast_output/%s/report.tsv" % (output_dir, id)

        if not os.path.isfile(quast_out):
            print("Gathering {id} assembly quality metrics with Quast. . .".format(id=id))
            quast_obj.run()

        with open(quast_out) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for line in tsv_reader:
                if "Total length" in line[0]:
                    isolate_qual[id]["est_genome_length"]= line[1]

        pathlib.Path(output_dir + "/cg_pipeline_output/").mkdir(parents=True, exist_ok=True)
        cg_mounting= {os.path.abspath(path): '/datain', os.path.abspath(output_dir) + "/cg_pipeline_output" :'/dataout'}
        with open("%s/quast_output/%s/report.tsv" % (output_dir, id)) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for line in tsv_reader:
                if "Total length" == line[0]:
                    genome_length = line[1]

        isolate_qual[id]["est_genome_length"] = genome_length

        if not genome_length:
            raise ValueError("Unable to predict genome length for isolate" + id)

        cg_command = "bash -c 'run_assembly_readMetrics.pl {cg_params} {in_dir}/{reads} -e {genome_length} > " \
                          "{out_dir}/{cgp_result}'".format(
            in_dir=in_dir,out_dir=out_dir,reads=all_reads,genome_length=genome_length,
            cgp_result=id + "_readMetrics.tsv", cg_params=tredegar_config["parameters"]["cg_pipeline"])

        cg_obj= sb_programs.Run(command=cg_command, path=cg_mounting, docker_image="lyveset")
        cg_out = "%s/cg_pipeline_output/%s_readMetrics.tsv"%(output_dir,id)

        if not os.path.isfile(cg_out):
            print("Getting {id} sequencing quality metrics with CG Pipeline".format(id=id))
            cg_obj.run()
            print("CG Pipeline complete.")

        with open(cg_out) as tsv_file:
            tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

            for line in tsv_reader:
                if any(fwd_format in line["File"] for fwd_format in ["_1.fastq", "_R1.fastq"]):
                    isolate_qual[id]["r1_q"] = line["avgQuality"]
                    isolate_qual[id]["est_cvg"] = float(line["coverage"])
                if any(rev_format in line["File"] for rev_format in ["_2.fastq", "_R2.fastq"]):
                    isolate_qual[id]["r2_q"] = line["avgQuality"]
                    isolate_qual[id]["est_cvg"] += float(line["coverage"])

        if "Escherichia_coli" in isolate_qual[id]["predicted_species"]:
            pathlib.Path(output_dir + "/serotypefinder_output/").mkdir(parents=True, exist_ok=True)
            stf_mounting = {os.path.abspath(output_dir): '/datain',
                           os.path.abspath(output_dir) + "/serotypefinder_output": '/dataout'}
            stf_command = "serotypefinder.pl -d /serotypefinder/database/ -i {in_dir}/{assembly} " \
                      "-b /blast-2.2.26/ -o {out_dir}/{id}/ {stf_params} ".format(
                in_dir=in_dir, out_dir=out_dir, assembly=assembly, id=id,
                stf_params=tredegar_config["parameters"]["serotypefinder"])

            stf_obj = sb_programs.Run(command=stf_command, path=stf_mounting, docker_image="serotypefinder")
            stf_out = "%s/serotypefinder_output/%s/results_tab.txt" % (output_dir, id)

            if not os.path.isfile(stf_out):
                print("Isolate {id} identified as E.coli. Running SerotypeFinder for serotype prediction".format(id=id))
                stf_obj.run()

            with open(stf_out) as tsv_file:
                tsv_reader = csv.reader(tsv_file, delimiter="\t")
                wzx_allele=""
                wzy_allele=""
                wzm_allele=""

                for line in tsv_reader:

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

                if o_type:
                    isolate_qual[id]["predicted_serotype"] = o_type

        if "Salmonella_enterica" in isolate_qual[id]["predicted_species"]:
            pathlib.Path(output_dir + "/seqsero_output/").mkdir(parents=True, exist_ok=True)
            seqsero_mounting = {os.path.abspath(path): '/datain',
                           os.path.abspath(output_dir) + "/seqsero_output": '/dataout'}
            seqsero_command = "bash -c 'SeqSero.py -m2 -i {in_dir}/{reads} -d {out_dir}/{id}'".format(
                in_dir=in_dir, out_dir=out_dir, assembly=assembly, id=id, reads=all_reads)
            seqsero_obj = sb_programs.Run(command=seqsero_command, path=seqsero_mounting, docker_image="seqsero")
            seqsero_out = "%s/seqsero_output/%s/Seqsero_result.txt" % (output_dir, id)

            if not os.path.isfile(seqsero_out):
                print("Isolate {id} identified as identified as S.enterica. Running SeqSero for "
                      "serotype prediction".format(id=id))
                seqsero_obj.run()

            with open(seqsero_out) as tsv_file:
                tsv_reader = csv.reader(tsv_file, delimiter="\t")
                for line in tsv_reader:
                    if "Predicted serotype(s)" in line[0]:
                        isolate_qual[id]["predicted_serotype"] = line[1]

    # Curate Tredegar report
    pathlib.Path(output_dir + "/reports/").mkdir(parents=True, exist_ok=True)
    tredegar_out = "{output_dir}/reports/{project}_tredegar_report.tsv".format(output_dir=output_dir,project=project)

    column_headers=["sample", "r1_q", "r2_q", "est_genome_length", "est_cvg", "predicted_species", "predicted_serotype"]

    if not os.path.isfile(tredegar_out):
        with open(tredegar_out, "w") as csvfile:
            w = csv.DictWriter(csvfile, column_headers, dialect=csv.excel_tab)
            w.writeheader()
            for key,val in sorted(isolate_qual.items()):
                row = {"sample":key}
                row.update(val)
                w.writerow(row)

    print("Tredegar is complete! Output saved as %s"%tredegar_out)


if __name__ == '__main__':
    main()
