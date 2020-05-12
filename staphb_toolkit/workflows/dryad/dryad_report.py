#!/usr/bin/env python3
#Author: Kelsey Florek and Abigail Shockey
#email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu
#description: A pipeline for constructing SNP based and  core gene set reference free phylogenies


import sys,os,re
import argparse
from shutil import which, copyfile
from datetime import date
import pexpect

import re, sys

def main():
    #get nextflow executable
    lib_path = os.path.abspath(os.path.dirname(__file__) + '/lib')
    dryad_path = os.path.abspath(os.path.dirname(__file__))
    nextflow_path = os.path.join(lib_path,'nextflow')

    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description='Rebuild a previously generated PDF report.')
    parser.add_argument('rmd',type=str,help="path to Rmarkdown file (.Rmd)",nargs='?', default=False)
    parser.add_argument('snp_matrix',type=str,help="path to snp matrix",nargs='?', default=False)
    parser.add_argument('cg_tree',type=str,help="path to core genome tree",nargs='?', default=False)
    parser.add_argument('--ar',type=str,help="path to ar TSV file")
    parser.add_argument('--profile', type=str,choices=["docker", "singularity"],help="specify nextflow profile, dryad_report will try to use docker first, then singularity")
    parser.add_argument('--get_config',action="store_true",help="get a Nextflow configuration template for dryad")
    parser.add_argument('--config','-c', type=str,help="Nextflow custom configuration")

    args = parser.parse_args()

    #give config to user if requested
    if args.get_config:
        config_path = os.path.join(dryad_path,"configs/dryad_config_template.config")
        dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_dryad.config")
        copyfile(config_path,dest_path)
        sys.exit()

    #check for reads_path
    if not args.rmd or not args.snp_matrix or not args.cg_tree:
        parser.print_help()
        sys.exit(1)


    #check if we are using docker or singularity
    if which('docker'):
        profile = '-profile docker'
    elif which('singularity'):
        profile = '-profile singularity'
    else:
        profile = ''

    #check for config or profile
    config = ""
    if args.config:
        config = "-C " + os.path.abspath(args.config)
        profile = ""
    elif args.profile:
        if which(args.profile):
            profile = '-profile ' + args.profile
        else:
            print(f"{args.profile} is not installed or found in PATH.")
    elif not profile:
        print('Singularity or Docker is not installed or not found in PATH.')
        sys.exit(1)

    #set work dir into local logs dir if profile not aws
    work = ""
    output_path = os.path.join(os.getcwd(),'rebuild_results')
    output_work = os.path.join(output_path,'report_work')
    if profile:
        work = f"-w {output_work}"

    rmd = os.path.abspath(args.rmd)
    logo_path = os.path.abspath(os.path.dirname(__file__) + 'assets/dryad_logo_250.png')
    snp_mat = "--snp_matrix " + os.path.abspath(args.snp_matrix)
    cg_tree = "--cg_tree " + os.path.abspath(args.cg_tree)
    if args.ar:
        ar_tsv = "--ar_tsv " + os.path.abspath(args.ar)
    else:
        ar_tsv = ""

    #build command
    command = nextflow_path
    command = command + f" {config} run {dryad_path}/rebuild_report.nf {profile} --logo {logo_path} --outdir {output_path} --rmd {rmd} {snp_mat} {cg_tree} {ar_tsv} {work}"

    #run command using nextflow in a subprocess
    print("Rebuilding Dryad Report:")
    child = pexpect.spawn(command)
    child.interact()
