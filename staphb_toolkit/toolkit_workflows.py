#!/usr/bin/env python3

#authors:
# Kelsey Florek (kelsey.florek@slh.wisc.edu)
# Kevin Libuit (kevin.libuit@dgs.virginia.gov)

import sys,os,re
import argparse
from shutil import which, copyfile
from datetime import date, datetime
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
    parser.add_argument("--update",default=False,action="store_true",help="Check for and install a ToolKit update.")
    parser.add_argument("--auto_update",default=False,action="store_true",help="Toggle automatic ToolKit updates. Default is off.")
    subparsers = parser.add_subparsers(title='workflows',metavar='',dest="subparser_name")

    #parser for workflows
    #tredegar-----------------------------------------
    parser_tredegar = subparsers.add_parser('tredegar', help='Quality control of WGS read data.', add_help=False)
    parser_tredegar.add_argument('reads_path', type=str,help="path to the location of the reads in a fastq format",nargs='?', default=False)
    parser_tredegar.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"tredegar_results\".",default="tredegar_results")
    parser_tredegar.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    parser_tredegar.add_argument('--config','-c', type=str,help="Nextflow custom configuration.")
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
    subparser_monroe_pe_assembly.add_argument('--config','-c', type=str,help="Nextflow custom configuration.")

    ##monroe_ont_assembly----------------------------
    ##medaka nanopolish
    subparser_monroe_ont_assembly = monroe_subparsers.add_parser('ont_assembly',help='Assembly SARS-CoV-2 genomes from ONT read data generated from ARTIC amplicons', add_help=False)
    subparser_monroe_ont_assembly.add_argument('fast5_path', type=str,help="path to the location of the reads in a fast5 format")
    subparser_monroe_ont_assembly.add_argument('--fastq_path', type=str,help="path to the location of the reads in a fastq format, needed if not performing bascalling")
    subparser_monroe_ont_assembly.add_argument('--polish_method',type=str,choices=["medaka","nanopolish"],help="polishing method, default: medaka",default="medaka")
    subparser_monroe_ont_assembly.add_argument('--summary', type=str,help="path to the location of the sequencing summary, only needed when using nanopolish from fastq data")
    subparser_monroe_ont_assembly.add_argument('--run_prefix', type=str,help="desired run prefix. Default \"artic_ncov19\"",default="artic_ncov19")
    subparser_monroe_ont_assembly.add_argument('--ont_basecalling', default=False, action="store_true",help="perform high accuracy basecalling using GPU (only use if you have setup a GPU compatable device)")
    subparser_monroe_ont_assembly.add_argument('--primers', type=str,choices=["V1", "V2", "V3"], help="indicate which ARTIC primers were used (V1, V2, or V3)",required=True)
    subparser_monroe_ont_assembly.add_argument('--config','-c', type=str,help="Nextflow custom configuration.")
    subparser_monroe_ont_assembly.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"monroe_results\".",default="monroe_results")
    subparser_monroe_ont_assembly.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    subparser_monroe_ont_assembly.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")

    ##monroe_cluster_analysis-----------------------
    subparser_monroe_cluster_analysis = monroe_subparsers.add_parser('cluster_analysis',help='Perform multiple sequence alinmment of SC2 assemblies to generate a SNP-distance matrix & ML phylogenetic tree', add_help=False)
    subparser_monroe_cluster_analysis.add_argument('assemblies_path', type=str,help="path to the location of the SC2 assemblies in a fasta format",nargs='?', default=False)
    subparser_monroe_cluster_analysis.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"monroe_results\".",default="monroe_results")
    subparser_monroe_cluster_analysis.add_argument('--rtemplate','-rt', type=str, help="path to report rmarkdown", default=os.path.join(workflows_path,"monroe/report/report.Rmd"))
    subparser_monroe_cluster_analysis.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    subparser_monroe_cluster_analysis.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    subparser_monroe_cluster_analysis.add_argument('--config','-c', type=str,help="Nextflow custom configuration.")
    subparser_monroe_cluster_analysis.add_argument('--get_rtemplate',action="store_true",help="Get a Rmd configuration template for report building.")

    #foushee-----------------------------------------
    parser_foushee = subparsers.add_parser('foushee', help='Reference-free SNP calling for Streptococcus pyogenes isolates.', add_help=False)
    parser_foushee.add_argument('reads_path', type=str,help="path to the directory of raw reads in the fastq format",nargs='?', default=False)
    parser_foushee.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"tredegar_results\".",default="foushee_results")
    parser_foushee.add_argument('--rtemplate','-rt', type=str, help="path to report rmarkdown", default=os.path.join(workflows_path,"foushee/report/report.Rmd"))
    parser_foushee.add_argument('--profile', type=str,choices=["docker","singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    parser_foushee.add_argument('--config','-c', type=str,help="Nextflow custom configuration.")
    parser_foushee.add_argument('--get_config',action="store_true",help="Get a Nextflow configuration template for foushee.")
    parser_foushee.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")

    #dryad-----------------------------------------
    parser_dryad = subparsers.add_parser('dryad', help='A comprehensive tree building program.', add_help=False)
    subparser_dryad = parser_dryad.add_subparsers(title='dryad apps',metavar='',dest='dryad_app')

    subparser_dryad_main = subparser_dryad.add_parser('main',help="dryad workflow")
    subparser_dryad_main.add_argument('reads_path', type=str,help="path to the directory of raw reads in the fastq format",nargs='?', default=False)
    subparser_dryad_main.add_argument('--output','-o',metavar="<output_path>",type=str,help="path to ouput directory, default \"dryad_results\"",default="dryad_results")
    subparser_dryad_main.add_argument('--core-genome','-cg',default=False, action="store_true", help="construct a core-genome tree")
    subparser_dryad_main.add_argument('--snp','-s',default=False, action="store_true", help="construct a SNP tree, requires a reference sequence in fasta format (-r)")
    subparser_dryad_main.add_argument('-r',metavar='<path>', type=str,help="reference sequence for SNP pipeline")
    subparser_dryad_main.add_argument('-ar',default=False, action="store_true", help="detect AR mechanisms")
    subparser_dryad_main.add_argument('--sep',metavar="sep_chars",type=str,help="dryad identifies sample names from the name of the read file by splitting the name on the specified separating characters, default \"_\"",default="_")
    subparser_dryad_main.add_argument('--profile', type=str,choices=["docker", "singularity"],help="specify nextflow profile, dryad will try to use docker first, then singularity")
    subparser_dryad_main.add_argument('--config','-c', type=str,help="Nextflow custom configuration")
    subparser_dryad_main.add_argument('--get_config',action="store_true",help="get a Nextflow configuration template for dryad")
    subparser_dryad_main.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")
    subparser_dryad_main.add_argument('--report',action="store_true",help="generte a pdf report")

    subparser_dryad_report = subparser_dryad.add_parser('rebuild_report',help='rebuild a previously generated PDF report')
    subparser_dryad_report.add_argument('rmd',type=str,help="path to Rmarkdown file (.Rmd)",nargs='?', default=False)
    subparser_dryad_report.add_argument('snp_matrix',type=str,help="path to snp matrix",nargs='?', default=False)
    subparser_dryad_report.add_argument('cg_tree',type=str,help="path to core genome tree",nargs='?', default=False)
    subparser_dryad_report.add_argument('--ar',type=str,help="path to ar TSV file")
    subparser_dryad_report.add_argument('--profile', type=str,choices=["docker", "singularity"],help="specify nextflow profile, dryad_report will try to use docker first, then singularity")
    subparser_dryad_report.add_argument('--get_config',action="store_true",help="get a Nextflow configuration template for dryad")
    subparser_dryad_report.add_argument('--config','-c', type=str,help="Nextflow custom configuration")

    #----------------------------------------------
    args = parser.parse_args()

    #check if we are using docker or singularity
    if which('docker'):
        profile = '-profile docker'
    elif which('singularity'):
        profile = '-profile singularity'
    else:
        profile = ''

    #check for updates
    if args.update:
        autoupdate.check_for_updates()

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

    #current time for log files
    exec_time = datetime.now().strftime("%y_%m_%d_%H_%M_%S_")

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
            parser_tredegar.print_help()
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
        if profile and not args.config:
            work = f"-w {args.output}/logs/work"

        #build command
        command = nextflow_path
        command = command + f" {config} run {tredegar_path}/tredegar.nf {profile} {args.resume} --reads {args.reads_path} --outdir {args.output} -with-trace {args.output}/logs/{exec_time}Tredegar_trace.txt -with-report {args.output}/logs/{exec_time}Tredegar_execution_report.html {work}"
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
        if profile and not args.config:
            work = f"-w {args.output}/logs/work"

        if args.monroe_command == 'pe_assembly':
            #build command
            command = nextflow_path + f" {config} run {monroe_path}/monroe_pe_assembly.nf {profile} {args.resume} --pipe pe --reads {args.reads_path} --primers {args.primers} --outdir {args.output} -with-trace {args.output}/logs/{exec_time}Monroe_trace.txt -with-report {args.output}/logs/{exec_time}Monroe_execution_report.html {work}"
            #run command using nextflow in a subprocess
            print("Starting the Monroe paired-end assembly:")
            child = pexpect.spawn(command)
            child.interact()

        if args.monroe_command == 'cluster_analysis':
            #give report template to user if requested
            if args.get_rtemplate:
                rtemplate_path = os.path.join(monroe_path,"report/report.Rmd")
                dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_cluster_analysis_report.Rmd")
                copyfile(rtemplate_path,dest_path)
                sys.exit()

            #check for assemblies path
            if not args.assemblies_path:
                subparser_monroe_cluster_analysis.print_help()
                print("Please specify a path to a directory containing the raw reads.")
                sys.exit(1)

            #build command
            command = nextflow_path + f" {config} run {monroe_path}/monroe_cluster_analysis.nf {profile} {args.resume} --pipe cluster --assemblies {args.assemblies_path} --report {args.rtemplate} --outdir {args.output} -with-trace {args.output}/logs/{exec_time}Monroe_trace.txt -with-report {args.output}/logs/{exec_time}Monroe_execution_report.html {work}"
            #run command using nextflow in a subprocess
            print("Starting the Monroe cluster analysis:")
            child = pexpect.spawn(command)
            child.interact()

        if args.monroe_command == 'ont_assembly':
            #build command
            if args.ont_basecalling:
                read_paths = f"--basecalling --fast5_dir {args.fast5_path}"
            elif args.fast5_path and args.fastq_path:
                read_paths = f"--fast5_dir {args.fast5_path} --fastq_dir {args.fastq_path}"
            else:
                subparser_monroe_ont_assembly.print_help()
                print("Please provide path to both fastq and fast5 or perform basecalling.")
                sys.exit(1)

            if args.polish_method == "nanopolish" and args.summary:
                seq_summary = f"--sequencing_summary {args.summary}"
            elif args.polish_method == "nanopolish" and not args.summary and not args.ont_basecalling:
                print("Nanopolish requires a sequencing summary file generated during basecalling.")
                sys.exit(1)
            else:
                seq_summary = ""

            command = nextflow_path + f" {config} run {monroe_path}/monroe_ont_assembly.nf {profile} {args.resume} {read_paths} --pipe ont {seq_summary} --polishing {args.polish_method} --primers {args.primers} --outdir {args.output} --run_prefix {args.run_prefix} -with-trace {args.output}/logs/{exec_time}Monroe_trace.txt -with-report {args.output}/logs/{exec_time}Monroe_execution_report.html {work}"
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
            parser_foushee.print_help()
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
        if profile and not args.config:
            work = f"-w {args.output}/logs/work"

        #build command
        command = nextflow_path
        command = command + f" {config} run {foushee_path}/foushee.nf {profile} {args.resume} --reads {args.reads_path} --report {args.rtemplate} --outdir {args.output} -with-trace {args.output}/logs/{exec_time}Foushee_trace.txt -with-report {args.output}/logs/{exec_time}Foushee_execution_report.html {work}"
        print(command)

        #run command using nextflow in a subprocess
        print("Starting the Foushee pipeline:")
        child = pexpect.spawn(command)
        child.interact()

    #dryad-----------------------------------

    if program == 'dryad':
        #dryad path
        dryad_path = os.path.join(workflows_path,"dryad/")

        if args.dryad_app == None:
            parser_dryad.print_help()
            sys.exit(1)

        #give config to user if requested
        if args.get_config:
            config_path = os.path.join(dryad_path,"configs/dryad_config_template.config")
            dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_dryad.config")
            copyfile(config_path,dest_path)
            sys.exit()

        #rebuild report subroutine
        if args.dryad_app == "rebuild_report":
            #check for rmd path
            if not args.rmd or not args.snp_matrix or not args.cg_tree:
                subparser_dryad_report.print_help()
                sys.exit(1)

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
            if profile and not args.config:
                work = f"-w {output_work}"

            rmd = os.path.abspath(args.rmd)
            logo_path = os.path.join(dryad_path, 'assets/dryad_logo_250.png')
            snp_mat = "--snp_matrix " + os.path.abspath(args.snp_matrix)
            cg_tree = "--cg_tree " + os.path.abspath(args.cg_tree)
            if args.ar:
                ar_tsv = "--ar_tsv " + os.path.abspath(args.ar)
            else:
                ar_tsv = ""

            #build command
            command = nextflow_path
            command = command + f" {config} run {dryad_path}/rebuild_report.nf {profile} --logo {logo_path} --outdir {output_path} --rmd {rmd} {snp_mat} {cg_tree} {ar_tsv} -with-trace {args.output}/logs/{exec_time}dryad_trace.txt -with-report {args.output}/logs/{exec_time}dryad_execution_report.html {work}"

            #run command using nextflow in a subprocess
            print("Rebuilding Dryad Report:")
            child = pexpect.spawn(command)
            child.interact()

        #dryad main application
        if args.dryad_app == "main":
            #check for reads_path
            if not args.reads_path:
                subparser_dryad_main.print_help()
                print("Please specify a path to a directory containing the raw reads.")
                sys.exit(1)

            #check for reference sequence
            if args.snp and args.r == None:
                subparser_dryad_main.print_help()
                print("Please specify a reference sequence for the SNP pipeline.")
                sys.exit(1)

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
            if profile and not args.config:
                work = f"-w {args.output}/logs/work"

            #build nextflow command
            logo_path = os.path.join(dryad_path,"assets/dryad_logo_250.png")
            selections = ""
            if args.ar:
                selections += " --ar"
            if args.core_genome:
                selections += " --cg"
            if args.snp:
                selections += f" --snp --snp_reference {args.r}"
            if args.report and args.snp and args.core_genome:
                report_template_path = os.path.join(dryad_path,"report/report.Rmd")
                selections += f" --report {report_template_path} --logo {logo_path}"

            #path for multiqc config
            mqc_config_path = f"--multiqc_config " + os.path.join(dryad_path,"configs/multiqc_config.yaml")
            mqc_logo_path =  f"--multiqc_logo " + logo_path

            #add other arguments
            other_args = f"--name_split_on {args.sep} --outdir {args.output}"

            #build command
            command = nextflow_path
            command = command + f" {config} run {dryad_path}/dryad.nf {profile} {args.resume} --reads {args.reads_path} {selections} {other_args} {mqc_config_path} {mqc_logo_path} -with-trace {args.output}/logs/{exec_time}dryad_trace.txt -with-report {args.output}/logs/{exec_time}dryad_execution_report.html {work}"

            #run command using nextflow in a subprocess
            print("Starting the Dryad pipeline:")
            child = pexpect.spawn(command)
            child.interact()
