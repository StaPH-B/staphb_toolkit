#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import argparse
import csv
import pandas
import datetime
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

from staphB_ToolKit.lib import mash, cg_pipeline

def main():
    parser = argparse.ArgumentParser(usage="belvidere.py <input> [options]")
    parser.add_argument("input", type=str, nargs='?', help="path to dir containing read files")
    parser.add_argument("-o", default="", nargs='?', type=str, help="Name of output_dir")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    path=args.input
    output_dir = args.o

    if not output_dir:
        output_dir = os.path.abspath("belvidere_output")

    if output_dir.endswith('/'):
        output_dir = output_dir[:-1]

    if '/' in output_dir and "belvidere_output" not in output_dir:
        project = output_dir.split('/')
        project = project[-1]
    elif "belvidere_output" in output_dir:
        project = datetime.datetime.today().strftime('%Y-%m-%d')
    else:
        project = output_dir

    #Gather MASH & CG_Pipeline results
    mash_obj = mash.Mash(path=path,output_dir=output_dir)
    mash_species = mash_obj.mash_species()

    isolate_qual = {}
    cgpipeline_obj = cg_pipeline.CGPipeline(path=path,output_dir=output_dir)
    cgpipeline_obj.read_metrics()

    for id in mash_species:
        isolate_qual[id] = {"species": None, "r1_q": None, "r1_avgReadLength": None, "r1_totalBases": None,
                            "r1_numReads": None, "r2_q": None, "r2_avgReadLength": None, "r2_totalBases": None,
                            "r2_numReads": None, "est_cvg": None}

        isolate_qual[id]["predicted_species"] = mash_species[id]

        with open("%s/cg_pipeline_output/%s_readMetrics.tsv"%(output_dir,id)) as tsv_file:
            tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

            for line in tsv_reader:
                if any(fwd_format in line["File"] for fwd_format in ["_1.fastq", "_R1.fastq"]):
                    isolate_qual[id]["r1_q"] = line["avgQuality"]
                    isolate_qual[id]["r1_avgReadLength"] = line["avgReadLength"]
                    isolate_qual[id]["r1_totalBases"] = line["totalBases"]
                    isolate_qual[id]["r1_numReads"] = line["numReads"]
                    isolate_qual[id]["est_cvg"] = float(line["coverage"])
                if any(rev_format in line["File"] for rev_format in ["_2.fastq", "_R2.fastq"]):
                    isolate_qual[id]["r2_q"] = line["avgQuality"]
                    isolate_qual[id]["r2_avgReadLength"] = line["avgReadLength"]
                    isolate_qual[id]["r2_totalBases"] = line["totalBases"]
                    isolate_qual[id]["r2_numReads"] = line["numReads"]
                    isolate_qual[id]["est_cvg"] += float(line["coverage"])

    # Curate Belvidere report

    reports_dir = output_dir + "/reports/"

    belvidere_out = "%s/%s_belvidere_report.csv"%(reports_dir,project)

    if not os.path.isdir(reports_dir):
        os.makedirs(reports_dir)
        print("Directory for WGS reports made:", reports_dir)

    # Change data dictionary to dataframe to csv
    df = pandas.DataFrame(isolate_qual).T[["predicted_species", "r1_q", "r1_avgReadLength","r1_totalBases",
                                           "r1_numReads", "r2_q", "r2_avgReadLength","r2_totalBases", "r2_numReads",
                                           "est_cvg"]]
    df.to_csv(belvidere_out)

    print("Belvidere is complete! Output saved as %s"%belvidere_out)


if __name__ == '__main__':
    main()
