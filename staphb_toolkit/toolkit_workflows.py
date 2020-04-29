#!/usr/bin/env python3

#authors:
# Kelsey Florek (kelsey.florek@slh.wisc.edu)
# Kevin Libuit (kevin.libuit@dgs.virginia.gov)

import sys,os,re
import argparse
from shutil import which, copyfile
from datetime import date
import pexpect
import staphb_toolkit.core.update_app as autoupdate

def main():
    #get nextflow executable
    lib_path = os.path.abspath(os.path.dirname(__file__) + '/' + 'lib')
    workflows_path = os.path.abspath(os.path.dirname(__file__) + '/' + 'workflows')
    nextflow_path = os.path.join(lib_path,'nextflow')

    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(usage="staphb-wf [optional arguments] <workflow> [workflow arguments]")
    parser.add_argument("--auto_update",default=False,action="store_true",help="Toggle automatic ToolKit updates. Default is off.")
    subparsers = parser.add_subparsers(title='workflows',metavar='',dest="subparser_name")

    #check if we are using docker or singularity
    if which('docker'):
        profile = '-profile docker'
    elif which('singularity'):
        profile = '-profile singularity'
    else:
        profile = ''

    #parser for workflows
    #tredegar-----------------------------------------
    parser_tredegar = subparsers.add_parser('tredegar', help='Quality control of WGS read data.', add_help=False)
    parser_tredegar.add_argument('reads_path', type=str,help="path to the location of the reads in a fastq format",nargs='?', default=False)
    parser_tredegar.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"tredegar_results\".",default="tredegar_results")
    parser_tredegar.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    parser_tredegar.add_argument('--config','-c', type=str,help="Nextflow custom configureation.")
    parser_tredegar.add_argument('--get_config',action="store_true",help="Get a Nextflow configuration template for tredegar.")
    parser_tredegar.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")

    #monroe-----------------------------------------
    parser_monroe = subparsers.add_parser('monroe', help='Consensus assembly for SARS-CoV-2 from ARTIC + Illumina protocols.', add_help=False)
    monroe_subparsers = parser_monroe.add_subparsers(title='monroe_commands',metavar='',dest='monroe_command')
    parser_monroe.add_argument('--get_config',type=str,choices=["pe_assembly", "ont_assembly", "cluster_analysis"],help="Get a Nextflow configuration template for the chosen workflow.")

    ##monroe_pe_assembly----------------------------
    subparser_monroe_pe_assembly = monroe_subparsers.add_parser('pe_assembly',help='Assembly SARS-CoV-2 genomes from paired-end read data generated from ARTIC amplicons', add_help=False)
    subparser_monroe_pe_assembly.add_argument('reads_path', type=str,help="path to the location of the reads in a fastq format")
    subparser_monroe_pe_assembly.add_argument('--primers', type=str,choices=["V1", "V2", "V3"], help="indicate which ARTIC primers were used (V1, V2, or V3)",required=True)
    subparser_monroe_pe_assembly.add_argument('--profile', type=str,choices=["docker", "singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    subparser_monroe_pe_assembly.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"monroe_results\".",default="monroe_results")
    subparser_monroe_pe_assembly.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    subparser_monroe_pe_assembly.add_argument('--config','-c', type=str,help="Nextflow custom configureation.")

    ##monroe_ont_assembly----------------------------
    subparser_monroe_ont_assembly = monroe_subparsers.add_parser('ont_assembly',help='Assembly SARS-CoV-2 genomes from ONT read data generated from ARTIC amplicons', add_help=False)
    subparser_monroe_ont_assembly.add_argument('reads_path', type=str,help="path to the location of the reads in a fastq or fast5 format")
    subparser_monroe_ont_assembly.add_argument('sequencing_summary', type=str,help="path to the location of the sequencing summary")
    subparser_monroe_ont_assembly.add_argument('--run_prefix', type=str,help="desired run prefix. Default \"artic_ncov19\"",default="artic_ncov19")
    subparser_monroe_ont_assembly.add_argument('--ont_basecalling', default=False, action="store_true",help="perform high accuracy basecalling using GPU (only use if you have setup a GPU compatable device)")
    subparser_monroe_ont_assembly.add_argument('--primers', type=str,choices=["V1", "V2", "V3"], help="indicate which ARTIC primers were used (V1, V2, or V3)",required=True)
    subparser_monroe_ont_assembly.add_argument('--config','-c', type=str,help="Nextflow custom configureation.")
    subparser_monroe_ont_assembly.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"monroe_results\".",default="monroe_results")
    subparser_monroe_ont_assembly.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    subparser_monroe_ont_assembly.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")

    ##monroe_cluster_analysis-----------------------
    subparser_monroe_cluster_analysis = monroe_subparsers.add_parser('cluster_analysis',help='Perform multiple sequence alinmment of SC2 assemblies to generate a SNP-distance matrix & ML phylogenetic tree', add_help=False)
    subparser_monroe_cluster_analysis.add_argument('assemblies_path', type=str,help="path to the location of the SC2 assemblies in a fasta format")
    subparser_monroe_cluster_analysis.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"monroe_results\".",default="monroe_results")
    subparser_monroe_cluster_analysis.add_argument('--report','-r', type=str, help="path to report rmarkdown", default=os.path.join(workflows_path,"monroe/report/report.Rmd"))
    subparser_monroe_cluster_analysis.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    subparser_monroe_cluster_analysis.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    subparser_monroe_cluster_analysis.add_argument('--config','-c', type=str,help="Nextflow custom configureation.")


    #foushee-----------------------------------------
    parser_foushee = subparsers.add_parser('foushee', help='Reference-free SNP calling for Streptococcus pyogenes isolates.', add_help=False)
    parser_foushee.add_argument('reads_path', type=str,help="path to the directory of raw reads in the fastq format",nargs='?', default=False)
    parser_foushee.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"tredegar_results\".",default="foushee_results")
    parser_foushee.add_argument('--report','-r', type=str, help="path to report rmarkdown", default=os.path.join(workflows_path,"foushee/report/report.Rmd"))
    parser_foushee.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    parser_foushee.add_argument('--config','-c', type=str,help="Nextflow custom configureation.")
    parser_foushee.add_argument('--get_config',action="store_true",help="Get a Nextflow configuration template for dryad.")
    parser_foushee.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    #dryad-----------------------------------------
    parser_dryad = subparsers.add_parser('dryad', help='A comprehensive tree building program.', add_help=False)
    parser_dryad.add_argument('reads_path', type=str,help="path to the directory of raw reads in the fastq format",nargs='?', default=False)
    parser_dryad.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"dryad_results\".",default="dryad_results")
    parser_dryad.add_argument('--core-genome','-cg',default=False, action="store_true", help="Construct a core-genome tree.")
    parser_dryad.add_argument('--snp','-s',default=False, action="store_true", help="Construct a SNP tree. Note: Requires a reference genome in fasta format (-r).")
    parser_dryad.add_argument('-r',metavar='<path>', type=str,help="Reference genome for SNP pipeline.")
    parser_dryad.add_argument('-ar',default=False, action="store_true", help="Detect AR mechanisms.")
    parser_dryad.add_argument('--sep',metavar="sep_chars",type=str,help="Dryad identifies sample names from the name of the read file by splitting the name on the specified separating characters, default \"_\".",default="_")
    parser_dryad.add_argument('--profile', type=str,choices=["docker", "singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    parser_dryad.add_argument('--config','-c', type=str,help="Nextflow custom configureation.")
    parser_dryad.add_argument('--get_config',action="store_true",help="Get a Nextflow configuration template for dryad.")
    parser_dryad.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")

    args = parser.parse_args()
    #check for updates
    if args.auto_update:
        #get current status
        update_status = autoupdate.check_update_status()
        if update_status:
            autoupdate.toggle_updater(False)
        else:
            autoupdate.toggle_updater(True)

    if autoupdate.check_update_status():
        autoupdate.check_for_updates()

    if args.subparser_name == None:
        parser.print_help()
        sys.exit(1)
    program = args.subparser_name

    #################################
    #Program specific execution code#
    #################################

    #tredegar--------------------------------

    if program == 'tredegar':
        #tredegar path
        tredegar_path = os.path.join(workflows_path,"tredegar/")

        #give config to user if requested
        if args.get_config:
            config_path = os.path.join(tredegar_path,"configs/tredegar_config_template.config")
            dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_tredegar.config")
            copyfile(config_path,dest_path)
            sys.exit()

        #check for reads_path
        if not args.reads_path:
            parser.print_help()
            print("Please specify a path to a directory containing the raw reads.")
            sys.exit(1)


        #check for config or profile
        config = ""
        if args.config:
            config = "-C " + os.path.abspath(args.config)
            profile = ""
        elif args.profile:
            profile = f"-profile {args.profile}"
        elif not profile:
            print('Singularity or Docker is not installed or not in found in PATH.')
            sys.exit(1)

        #set work dir into local logs dir if profile not aws
        work = ""
        if profile:
            work = f"-w {args.output}/logs/work"

        #build command
        command = nextflow_path
        command = command + f" {config} run {tredegar_path}/tredegar.nf {profile} {args.resume} --reads {args.reads_path} --outdir {args.output} -with-trace {args.output}/logs/Tredegar_trace.txt -with-report {args.output}/logs/Tredegar_execution_report.html {work}"
        print(command)

        #run command using nextflow in a subprocess
        print("Starting the Tredegar pipeline:")
        child = pexpect.spawn(command)
        child.interact()

    #monroe----------------------------------

    if program == 'monroe':
        #monroe path
        monroe_path = os.path.join(workflows_path,"monroe/")

        #give config to user if requested
        if args.get_config == "pe_assembly":
            config_path = os.path.join(monroe_path,"configs/pe_user_config.config")
            dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_pe_assembly.config")
            copyfile(config_path,dest_path)
            sys.exit()
        elif args.get_config == "ont_assembly":
            config_path = os.path.join(monroe_path,"configs/ont_user_config.config")
            dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_ont_assembly.config")
            copyfile(config_path,dest_path)
            sys.exit()
        elif args.get_config == "cluster_analysis":
            config_path = os.path.join(monroe_path,"configs/cluster_user_config.config")
            dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_cluster_analysis.config")
            copyfile(config_path,dest_path)
            sys.exit()

        if args.monroe_command == None:
            parser_monroe.print_help()
            sys.exit(1)

        #check for config or profile
        config = ""
        if args.config:
            config = "-C " + os.path.abspath(args.config)
            profile = ""
        elif args.profile:
            profile = f"-profile {args.profile}"
        elif not profile:
            print('Singularity or Docker is not installed or not in found in PATH.')
            sys.exit(1)

        #set work dir into local logs dir if profile not aws
        work = ""
        if profile:
            work = f"-w {args.output}/logs/work"

        if args.monroe_command == 'pe_assembly':
            #build command
            command = nextflow_path + f" {config} run {monroe_path}/monroe_pe_assembly.nf {profile} {args.resume} --pipe pe --reads {args.reads_path} --primers {args.primers} --outdir {args.output} -with-trace {args.output}/logs/Monroe_trace.txt -with-report {args.output}/logs/Monroe_execution_report.html {work}"
            #run command using nextflow in a subprocess
            print("Starting the Monroe paired-end assembly:")
            child = pexpect.spawn(command)
            child.interact()

        if args.monroe_command == 'cluster_analysis':
            #build command
            command = nextflow_path + f" {config} run {monroe_path}/monroe_cluster_analysis.nf {profile} {args.resume} --pipe cluster --assemblies {args.assemblies_path} --report {args.report} --outdir {args.output} -with-trace {args.output}/logs/Monroe_trace.txt -with-report {args.output}/logs/Monroe_execution_report.html {work}"
            #run command using nextflow in a subprocess
            print("Starting the Monroe cluster analysis:")
            child = pexpect.spawn(command)
            child.interact()

        if args.monroe_command == 'ont_assembly':
            #build command
            if(args.ont_basecalling):
                basecall = f"--basecalling --fast5_dir {args.reads_path}"
            else:
                basecall = f"--fastq_dir {args.reads_path}"

            command = nextflow_path + f" {config} run {monroe_path}/monroe_ont_assembly.nf {profile} {args.resume} {basecall} --pipe ont --sequencing_summary {args.sequencing_summary} --primers {args.primers} --outdir {args.output} --run_prefix {args.run_prefix} -with-trace {args.output}/logs/Monroe_trace.txt -with-report {args.output}/logs/Monroe_execution_report.html {work}"
            #run command using nextflow in a subprocess
            print("Starting the Monroe ONT assembly:")
            child = pexpect.spawn(command)
            child.interact()

    #foushee---------------------------------

    if program == 'foushee':
        #tredegar path
        foushee_path = os.path.join(workflows_path,"foushee/")

        #give config to user if requested
        if args.get_config:
            config_path = os.path.join(foushee_path,"configs/foushee_config_template.config")
            dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_foushee.config")
            copyfile(config_path,dest_path)
            sys.exit()

        #check for reads_path
        if not args.reads_path:
            parser.print_help()
            print("Please specify a path to a directory containing the raw reads.")
            sys.exit(1)


        #check for config or profile
        config = ""
        if args.config:
            config = "-C " + os.path.abspath(args.config)
            profile = ""
        elif args.profile:
            profile = f"-profile {args.profile}"
        elif not profile:
            print('Singularity or Docker is not installed or not in found in PATH.')
            sys.exit(1)

        #set work dir into local logs dir if profile not aws
        work = ""
        if profile:
            work = f"-w {args.output}/logs/work"

        #build command
        command = nextflow_path
        command = command + f" {config} run {foushee_path}/foushee.nf {profile} {args.resume} --reads {args.reads_path} --report {args.report} --outdir {args.output} -with-trace {args.output}/logs/Foushee_trace.txt -with-report {args.output}/logs/Foushee_execution_report.html {work}"
        print(command)

        #run command using nextflow in a subprocess
        print("Starting the Foushee pipeline:")
        child = pexpect.spawn(command)
        child.interact()

    #dryad-----------------------------------

    if program == 'dryad':
        #dryad path
        dryad_path = os.path.join(workflows_path,"dryad/")

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
            parser_dryad.print_help()
            print("Please specify a reference sequence for the SNP pipeline.")
            sys.exit(1)

        #check for config or profile
        config = ""
        if args.config:
            config = "-C " + os.path.abspath(args.config)
            profile = ""
        elif args.profile:
            profile = args.profile
        elif not profile:
            print('Singularity or Docker is not installed or not in found in PATH.')
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
        #add other arguments
        other_args = f"--name_split_on {args.sep} --outdir {args.output}"
        #build command
        command = nextflow_path
        command = command + f" {config} run {dryad_path}/dryad.nf {profile} {args.resume} --reads {args.reads_path} {selections} {other_args} {work}"

        #run command using nextflow in a subprocess
        print("Starting the Dryad pipeline:")
        child = pexpect.spawn(command)
        child.interact()
