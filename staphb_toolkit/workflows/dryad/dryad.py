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

    parser = MyParser(description='A comprehensive tree building program.')
    parser.add_argument('reads_path', type=str,help="path to the directory of raw reads in the fastq format",nargs='?', default=False)
    parser.add_argument('--output','-o',metavar="<output_path>",type=str,help="path to ouput directory, default \"dryad_results\"",default="dryad_results")
    parser.add_argument('--core-genome','-cg',default=False, action="store_true", help="construct a core-genome tree")
    parser.add_argument('--snp','-s',default=False, action="store_true", help="construct a SNP tree, requires a reference sequence in fasta format (-r)")
    parser.add_argument('-r',metavar='<path>', type=str,help="reference sequence for SNP pipeline")
    parser.add_argument('-ar',default=False, action="store_true", help="detect AR mechanisms")
    parser.add_argument('--sep',metavar="sep_chars",type=str,help="dryad identifies sample names from the name of the read file by splitting the name on the specified separating characters, default \"_\"",default="_")
    parser.add_argument('--profile', type=str,choices=["docker", "singularity"],help="specify nextflow profile, dryad will try to use docker first, then singularity")
    parser.add_argument('--config','-c', type=str,help="Nextflow custom configuration")
    parser.add_argument('--get_config',action="store_true",help="get a Nextflow configuration template for dryad")
    parser.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    parser.add_argument('--report',action="store_true",help="generte a pdf report")

    args = parser.parse_args()

    #give config to user if requested
    if args.get_config:
        config_path = os.path.join(dryad_path,"configs/dryad_config_template.config")
        dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_dryad.config")
        copyfile(config_path,dest_path)
        sys.exit()

    #check for reads_path
    if not args.reads_path:
        parser.print_help()
        print("Please specify a path to a directory containing the raw reads.")
        sys.exit(1)

    #check for reference sequence
    if args.snp and args.r == None:
        parser.print_help()
        print("Please specify a reference sequence for the SNP pipeline.")
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
    if profile:
        work = f"-w {args.output}/logs/work"

    #build nextflow command
    selections = ""
    if args.ar:
        selections += " --ar"
    if args.core_genome:
        selections += " --cg"
    if args.snp:
        selections += f" --snp --snp_reference {args.r}"
    if args.report and args.snp and args.core_genome:
        report_template_path = os.path.abspath(os.path.dirname(__file__) + '/report/report.Rmd')
        logo_path = os.path.abspath(os.path.dirname(__file__) + '/assets/dryad_logo_250.png')
        selections += f" --report {report_template_path} --logo {logo_path}"

    #path for multiqc config
    mqc_config_path = f"--multiqc_config " + os.path.join(dryad_path,"configs/multiqc_config.yaml")
    mqc_logo_path =  f"--multiqc_logo " + os.path.join(dryad_path,"assets/dryad_logo_250.png")

    #add other arguments
    other_args = f"--name_split_on {args.sep} --outdir {args.output}"

    #build command
    command = nextflow_path
    command = command + f" {config} run {dryad_path}/dryad.nf {profile} {args.resume} --reads {args.reads_path} {selections} {other_args} {mqc_config_path} {mqc_logo_path} -with-trace {args.output}/logs/dryad_trace.txt -with-report {args.output}/logs/dryad_execution_report.html {work}"

    #run command using nextflow in a subprocess
    print("Starting the Dryad pipeline:")
    child = pexpect.spawn(command)
    child.interact()
