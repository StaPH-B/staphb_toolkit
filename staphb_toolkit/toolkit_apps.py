#!/usr/bin/env python3

#authors:
# Kelsey Florek (kelsey.florek@slh.wisc.edu)
# Kevin Libuit (kevin.libuit@dgs.virginia.gov)

import sys,os,re
import argparse
from shutil import copy
import json
import staphb_toolkit.core.sb_programs as sb_prog
from staphb_toolkit.core.autopath import path_replacer
from staphb_toolkit.apps.sb_mash_species import MashSpecies
import staphb_toolkit.core.update_app as autoupdate
from datetime import date

#program dictionary
progs = {
'abricate':'Abricate - Mass screening of contigs for antimicrobial and virulence genes',
'augur': 'Pipeline components for real-time phylodynamic analysis',
'bbtools':'BBTools - Suite of fast, multithreaded bioinformatics tools for DNA and RNA sequence data',
'bwa':'BWA - mapping low-divergent sequences against a large reference genome',
'canu-racon':'Canu-Racon - Ultrafast consensus module for raw de novo assembly of long, uncorrected reads.',
'cfsan-snp':'CFSAN-SNP - SNP calling pipeline from the FDA CFSAN laboratory',
'circlator':'Circlator - A tool to circularize genome assemblies',
'clustalo':'ClustalO - A fast multiple sequence alignment program',
'emm-typing-tool':'Emm-typing-tool - Group A streptococci emm typing tool for NGS data',
'fastani':'FastANI - Fast whole-genome sequence average nucleotide identity (ANI) estimation',
'fastqc':'FastQC - A quality control tool for high throughput sequence data.',
'fasttree':'FastTree - Infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.',
'filtlong':'Filtlong - Quality filtering tool for long reads',
'flye':'Flye - De novo assembler for single molecule sequencing reads using repeat graphs',
'iqtree':'IQ-TREE - A fast and effective stochastic algorithm to infer phylogenetic trees by maximum likelihood.',
'kma':'KMA - Mapping method designed to map raw reads directly against redundant databases, in an ultra-fast manner using seed and extend.',
'ivar':'iVar - Computational package that contains functions broadly useful for viral amplicon-based sequencing.',
'ivar-SC2':'iVar - Computational package that contains functions broadly useful for viral amplicon-based sequencing: SARS-CoV-2 (SC2) reference sequence and ARTIC primers available in /reference ',
'kraken':'Kraken - Taxonomic sequence classification system ',
'kraken-build':'Kraken-build - Build a kraken database',
'kraken2':'Kraken2 - The second version of the Kraken taxonomic sequence classification system',
'kraken2-build':'Karken2-build - Build a kraken2 database',
'ksnp3':'kSNP3 - Identifies the pan-genome SNPs in a set of genome sequences, and estimates phylogenetic trees based upon those SNPs.',
'legsta':'Legsta - In silico Legionella pneumophila Sequence Based Typing',
'lyveset':'LYVE-SET - a method of using hqSNPs to create a phylogeny.',
'mash':'MASH - Fast genome and metagenome distance estimation using MinHash',
'mashtree':'MashTree - Create a tree using Mash distances',
'medaka':'Medaka - Sequence correction provided by ONT Research',
'minimap2':'Minimap2 - a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database.',
'mlst':'MLST - Scan contig files against PubMLST typing schemes',
'mugsy':'Mugsy - A multiple whole genome aligner.',
'multiqc':'MultiQC - Aggregate results from bioinformatics analyses across many samples into a single report.',
'nanoplot':'NanoPlot - Plotting scripts for long read sequencing data ',
'ncbi-amrfinder-plus':'NCBI AMRFinderPlus - Designed to find acquired antimicrobial resistance genes and some point mutations in protein or assembled nucleotide sequences.',
'orthofinder':'OrthoFinder - Phylogenetic orthology inference for comparative genomics',
'pilon':'Pilon - Automated genome assembly improvement and variant detection tool',
'plasmidseeker':'PlasmidSeeker - A k-mer based program for the identification of known plasmids from whole-genome sequencing reads',
'prokka':'Prokka - Rapid prokaryotic genome annotation',
'quast':'Quast - Genome assembly evaluation tool.',
'rasusa':'RASUA - Randomly subsample sequencing reads to a specified coverage',
'raxml':'RAxML -Maximum likelihood tree builder.',
'roary':'Roary - Rapid large-scale prokaryote pan genome analysis.',
'salmid':'SalmID - Rapid confirmation of Salmonella spp. and subspp. from sequence data',
'samtools':'Samtools - A suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories.',
'seqsero':'SeqSero - Salmonella serotyping from genome sequencing data.',
'seqsero2':'SeqSero2 - Salmonella serotype prediction from genome sequencing data.',
'seqyclean':'SeqyClean - Pre-process and clean NGS data in order to prepare for downstream analysis',
'seroba':'Seroba - k-mer based Pipeline to identify the Serotype from Illumina NGS reads ',
'serotypefinder':'SerotypeFinder - identifies the serotype in total or partial sequenced isolates of E. coli.',
'shovill':'Shovill - Faster SPAdes assembler',
'sistr':'SISTR - Salmonella in silico typing resource command-line tool',
'skesa':'SKESA - NCBI\'s de novo genome assemlber',
'snippy':'Snippy - Rapid haploid variant calling and core genome alignment',
'snp-dists':'SNP-dists - Pairwise SNP distance matrix from a FASTA sequence alignment',
'snp-sites': 'SNP-sites - Finds SNP sites from a multi-FASTA alignment file',
'spades':'SPAdes - St. Petersburg genome assembler',
'sra-toolkit':'SRA ToolKit - Collection of tools and libraries for using data in the INSDC Sequence Read Archives.',
'staramr':'StarAMR - Scans genome contigs against the ResFinder, PlasmidFinder, and PointFinder databases.',
'tiptoft':'TipToft - Predict plasmids from uncorrected long read data',
'trimmomatic':'Trimmoamtic - Flexible read trimming tool for Illumina NGS data',
'unicycler':'Unicycler - an assembly pipeline for bacterial genomes.',
'wtdbg2':'WTDBG2 - Fuzzy Bruijn graph approach to long noisy reads assembly'
}

def main():

    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            self.print_help()
            sys.stderr.write('\nerror: %s\n' % message)
            sys.exit(1)

    docker_config_path = os.path.abspath(os.path.dirname(__file__) + '/' + 'core/docker_config.json')

    parser = MyParser(usage="staphb-tk [optional arguments] <application> [application arguments]",add_help=True)
    subparsers = parser.add_subparsers(title='custom toolkit applications',metavar='',dest="subparser_name",parser_class=MyParser)
    parser.add_argument("--docker_config","-c", default=docker_config_path,metavar="<path>", help="Configuration file for container images and tags; if none provided, default container versions will be used.")
    parser.add_argument("--get_docker_config",default=False,action="store_true",help="Get the default docker container configureation file.")
    parser.add_argument("--list","-l",default=False,action="store_true",help="List all of the software available in the toolkit.")
    parser.add_argument("--update",default=False,action="store_true",help="Check for and install a ToolKit update.")
    parser.add_argument("--auto_update",default=False,action="store_true",help="Toggle automatic ToolKit updates. Default is off.")
    ###custom apps
    ## Mash Species
    parser_mash_species = subparsers.add_parser('mash_species',help='MASH_species uses a pre-sketched RefSeq database to identify the isolate species from paired-end read data.', usage="sb_mash_species <input> [options]")
    parser_mash_species.add_argument("input", type=str, nargs='?', help="path to dir containing paire-end read files")
    parser_mash_species.add_argument("-o",metavar='path', default="", type=str, help="Path for output directory",required=False)

    #parser for applications
    #-----------------------------------------
    parser_abricate = subparsers.add_parser('abricate', add_help=False)
    parser_augur = subparsers.add_parser('augur', add_help=False)
    parser_bbtols = subparsers.add_parser('bbtools', add_help=False)
    parser_bwa = subparsers.add_parser('bwa', add_help=False)
    parser_canuracon = subparsers.add_parser('canu-racon', add_help=False)
    parser_cfsansnp = subparsers.add_parser('cfsan-snp', add_help=False)
    parser_circlator = subparsers.add_parser('circlator', add_help=False)
    parser_clustalo = subparsers.add_parser('clustalo', add_help=False)
    parser_emmtypingtool = subparsers.add_parser('emm-typing-tool', add_help=False)
    parser_fastani = subparsers.add_parser('fastani', add_help=False)
    parser_fastqc = subparsers.add_parser('fastqc', add_help=False)
    parser_fasttree = subparsers.add_parser('fasttree', add_help=False)
    parser_filtong = subparsers.add_parser('filtlong', add_help=False)
    parser_flye = subparsers.add_parser('flye', add_help=False)
    parser_iqtree = subparsers.add_parser('iqtree', add_help=False)
    parser_ivar = subparsers.add_parser('ivar', add_help=False)
    parser_ivar_SC2 = subparsers.add_parser('ivar-SC2', add_help=False)
    parser_kma = subparsers.add_parser('kma', add_help=False)
    parser_kraken = subparsers.add_parser('kraken', add_help=False)
    parser_krakenbuild = subparsers.add_parser('kraken-build', add_help=False)
    parser_kraken2 = subparsers.add_parser('kraken2', add_help=False)
    parser_kraken2build = subparsers.add_parser('kraken2-build', add_help=False)
    parser_ksnp3 = subparsers.add_parser('ksnp3', add_help=False)
    parser_legsta = subparsers.add_parser('legsta', add_help=False)
    parser_lyveset = subparsers.add_parser('lyveset', add_help=False)
    parser_mash = subparsers.add_parser('mash', add_help=False)
    parser_mashtree = subparsers.add_parser('mashtree', add_help=False)
    parser_medaka = subparsers.add_parser('medaka', add_help=False)
    parser_minimap2 = subparsers.add_parser('minimap2', add_help=False)
    parser_mlst = subparsers.add_parser('mlst', add_help=False)
    parser_mugsy = subparsers.add_parser('mugsy', add_help=False)
    parser_multiqc = subparsers.add_parser('multiqc', add_help=False)
    parser_nanoplot = subparsers.add_parser('nanoplot', add_help=False)
    parser_ncbiamrfinder_plus = subparsers.add_parser('ncbi-amrfinder-plus', add_help=False)
    parser_orthofinder = subparsers.add_parser('orthofinder', add_help=False)
    parser_pilon = subparsers.add_parser('pilon', add_help=False)
    parser_plasmidseeker = subparsers.add_parser('plasmidseeker', add_help=False)
    parser_prokka = subparsers.add_parser('prokka', add_help=False)
    parser_quast = subparsers.add_parser('quast', add_help=False)
    parser_rasusa = subparsers.add_parser('rasusa', add_help=False)
    parser_raxml = subparsers.add_parser('raxml', add_help=False)
    parser_roary = subparsers.add_parser('roary', add_help=False)
    parser_salmid = subparsers.add_parser('salmid', add_help=False)
    parser_samtools = subparsers.add_parser('samtools', add_help=False)
    parser_seqsero = subparsers.add_parser('seqsero', add_help=False)
    parser_seqsero2 = subparsers.add_parser('seqsero2', add_help=False)
    parser_seqyclean = subparsers.add_parser('seqyclean', add_help=False)
    parser_seroba = subparsers.add_parser('seroba', add_help=False)
    parser_serotypefinder = subparsers.add_parser('serotypefinder', add_help=False)
    parser_shovill = subparsers.add_parser('shovill', add_help=False)
    parser_sistr = subparsers.add_parser('sistr', add_help=False)
    parser_skesa = subparsers.add_parser('skesa', add_help=False)
    parser_snippy = subparsers.add_parser('snippy', add_help=False)
    parser_snpdists = subparsers.add_parser('snp-dists', add_help=False)
    parser_snpsites = subparsers.add_parser('snp-sites', add_help=False)
    parser_spades = subparsers.add_parser('spades', add_help=False)
    parser_sratoolkit = subparsers.add_parser('sra-toolkit', add_help=False)
    parser_staramr = subparsers.add_parser('staramr', add_help=False)
    parser_tiptoft = subparsers.add_parser('tiptoft', add_help=False)
    parser_trimmomatic = subparsers.add_parser('trimmomatic', add_help=False)
    parser_unicycler = subparsers.add_parser('unicycler', add_help=False)
    parser_wtdbg2 = subparsers.add_parser('wtdbg2', add_help=False)

    #-----------------------------------------

    def print_prog_list():
        print("Available programs:")
        header = ["Command","Description","-------","-----------"]
        print(f"{header[0]:<25}{header[1]:^10}")
        print(f"{header[2]:<25}{header[3]:^10}")
        for key in progs:
            print(f"{key:<25}{progs[key]:^10}")
        return

    #handle the arguments and perform automatic path replacement
    parser_args = parser.parse_known_args()
    program = parser_args[0].subparser_name
    args = parser_args[1]

    #check for updates
    if parser_args[0].update:
        autoupdate.check_for_updates()
        sys.exit(0)

    if parser_args[0].auto_update:
        #get current status
        update_status = autoupdate.check_update_status()
        if update_status:
            autoupdate.toggle_updater(False)
        else:
            autoupdate.toggle_updater(True)

    if autoupdate.check_update_status():
        autoupdate.check_for_updates()

    #give user docker config if asked
    if parser_args[0].get_docker_config:
        cwd = os.getcwd()
        copy(docker_config_path,os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_docker_config.json"))
        sys.exit(0)

    #display list of programs if needed
    if parser_args[0].list:
        print_prog_list()
        sys.exit(0)

    if program == None:
        parser.print_help()
        sys.exit(1)

    #Run autopathing
    arg_string,path_map = path_replacer(args,os.getcwd())

    # set the configuration file
    if parser_args[0].docker_config == "/core/docker_config.json":
        # use default
        config_file_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + parser_args[0].docker_config
    else:
        config_file_path = os.path.abspath(parser_args[0].docker_config)

    with open(config_file_path, 'r') as config_file:
        config = json.load(config_file)

    #Custom program specific execution code
    #-----------------------------------------
    if program == 'mash_species':
        #get output dir if supplied, if not set it to cwd
        output_dir = None
        if not parser_args[0].o:
            output_dir = os.getcwd()
        else:
            try:
                output_dir = os.path.abspath(parser_args[0].o)
            except (AttributeError, TypeError):
                print("Please enter a valid output path.")
                sys.exit(1)

        #get input path, if not supplied print help
        try:
            path = os.path.abspath(parser_args[0].input)
        except (AttributeError, TypeError) as e:
            parser_mash_species.print_help()
            print("Please enter a valid input path.")
            sys.exit(1)

        #create and run the mash species object
        mash_species_obj = MashSpecies(path=path,output_dir=output_dir)
        mash_species_obj.run()


    #Program specific execution code
    #-----------------------------------------
    if program == 'ivar-SC2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = " "
        command = "ivar " + arg_string
        ivar_configuration = config["parameters"]["ivar-SC2"]
        ivar = sb_prog.Run(command=command, path=path_map, image=ivar_configuration["image"], tag=ivar_configuration["tag"])
        ivar.run()

    if program == 'ivar':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = " "
        command = "ivar " + arg_string
        ivar_configuration = config["parameters"]["ivar"]
        ivar = sb_prog.Run(command=command, path=path_map, image=ivar_configuration["image"], tag=ivar_configuration["tag"])
        ivar.run()

    if program == 'wtdbg2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "wtdbg2 " + arg_string
        wtdbg2_configuration = config["parameters"]["wtdbg2"]
        wtdbg2 = sb_prog.Run(command=command, path=path_map, image=wtdbg2_configuration["image"],
                             tag=wtdbg2_configuration["tag"])
        wtdbg2.run()

    if program == 'trimmomatic':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "trimmomatic " + arg_string
        trimmomatic_configuration = config["parameters"]["trimmomatic"]
        trimmomatic = sb_prog.Run(command=command, path=path_map, image=trimmomatic_configuration["image"], tag=trimmomatic_configuration["tag"])
        trimmomatic.run()

    if program == 'tiptoft':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "tiptoft " + arg_string
        tiptoft_configuration = config["parameters"]["tiptoft"]
        tiptoft = sb_prog.Run(command=command, path=path_map, image=tiptoft_configuration["image"], tag=tiptoft_configuration["tag"])
        tiptoft.run()

    if program == 'staramr':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "staramr " + arg_string
        staramr_configuration = config["parameters"]["staramr"]
        staramr = sb_prog.Run(command=command, path=path_map, image=staramr_configuration["image"], tag=staramr_configuration["tag"])
        staramr.run()

    if program == 'sra-toolkit':
        command = " " + arg_string
        sra_toolkit_configuration = config["parameters"]["sra-toolkit"]
        sra_toolkit = sb_prog.Run(command=command, path=path_map, image=sra_toolkit_configuration["image"], tag=sra_toolkit_configuration["tag"])
        if not re.search('[a-zA-Z]', arg_string):
            print("SRA toolkit tool must be specified, e.g. staphb-tk sra-toolkit fasterq-dump, staphb-tk sra-toolkit sra-pileup, etc. \n\nMore info on SRA Toolkit usage at: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc.")
        else:
            sra_toolkit.run()

    if program == 'snp-dists':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "snp-dists " + arg_string
        snp_dists_configuration = config["parameters"]["snp-dists"]
        snp_dists = sb_prog.Run(command=command, path=path_map, image=snp_dists_configuration["image"],
                                tag=snp_dists_configuration["tag"])
        snp_dists.run()

    if program == 'snp-sites':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "snp-sites " + arg_string
        snp_sites_configuration = config["parameters"]["snp-sites"]
        snp_sites = sb_prog.Run(command=command, path=path_map, image=snp_sites_configuration["image"],
                                tag=snp_sites_configuration["tag"])
        snp_sites.run()

    if program == 'snippy':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "snippy " + arg_string
        snippy_configuration = config["parameters"]["snippy"]
        snippy = sb_prog.Run(command=command, path=path_map, image=snippy_configuration["image"], tag=snippy_configuration["tag"])
        snippy.run()

    if program == 'skesa':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "skesa " + arg_string
        skesa_configuration = config["parameters"]["skesa"]
        skesa = sb_prog.Run(command=command, path=path_map, image=skesa_configuration["image"], tag=skesa_configuration["tag"])
        skesa.run()

    if program == 'sistr':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "sistr " + arg_string
        sistr_configuration = config["parameters"]["sistr"]
        sistr = sb_prog.Run(command=command, path=path_map, image=sistr_configuration["image"], tag=sistr_configuration["tag"])
        sistr.run()

    if program == 'seroba':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "seroba " + arg_string
        seroba_configuration = config["parameters"]["seroba"]
        seroba = sb_prog.Run(command=command, path=path_map, image=seroba_configuration["image"], tag=seroba_configuration["tag"])
        seroba.run()

    if program == 'seqsero2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "SeqSero2_package.py " + arg_string
        seqsero2_configuration = config["parameters"]["seqsero2"]
        seqsero2 = sb_prog.Run(command=command, path=path_map, image=seqsero2_configuration["image"],
                               tag=seqsero2_configuration["tag"])
        seqsero2.run()

    if program == 'salmid':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "SalmID.py " + arg_string
        salmid_configuration = config["parameters"]["salmid"]
        salmid = sb_prog.Run(command=command, path=path_map, image=salmid_configuration["image"], tag=salmid_configuration["tag"])
        salmid.run()

    if program == 'rasusa':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "rasusa " + arg_string
        rasusa_configuration = config["parameters"]["rasusa"]
        rasusa = sb_prog.Run(command=command, path=path_map, image=rasusa_configuration["image"], tag=rasusa_configuration["tag"])
        rasusa.run()

    if program == 'plasmidseeker':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "plasmidseeker.pl " + arg_string
        plasmidseeker_configuration = config["parameters"]["plasmidseeker"]
        plasmidseeker = sb_prog.Run(command=command, path=path_map, image=plasmidseeker_configuration["image"],
                                    tag=plasmidseeker_configuration["tag"])
        plasmidseeker.run()

    if program == 'pilon':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "--help"
        command = "pilon " + arg_string
        pilon_configuration = config["parameters"]["pilon"]
        pilon = sb_prog.Run(command=command, path=path_map, image=pilon_configuration["image"], tag=pilon_configuration["tag"])
        pilon.run()

    if program == 'orthofinder':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "orthofinder " + arg_string
        orthofinder_configuration = config["parameters"]["orthofinder"]
        orthofinder = sb_prog.Run(command=command, path=path_map, image=orthofinder_configuration["image"],
                                  tag=orthofinder_configuration["tag"])
        orthofinder.run()

    if program == 'ncbi-amrfinder-plus':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "--help"
        command = "amrfinder " + arg_string
        ncbi_amrfinderplus_configuration = config["parameters"]["ncbi-amrfinder-plus"]
        ncbi_amrfinderplus = sb_prog.Run(command=command, path=path_map, image=ncbi_amrfinderplus_configuration["image"], tag=ncbi_amrfinderplus_configuration["tag"])
        ncbi_amrfinderplus.run()

    if program == 'nanoplot':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "NanoPlot " + arg_string
        nanoplot_configuration = config["parameters"]["nanoplot"]
        nanoplot = sb_prog.Run(command=command, path=path_map, image=nanoplot_configuration["image"],
                               tag=nanoplot_configuration["tag"])
        nanoplot.run()

    if program == 'multiqc':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "multiqc " + arg_string
        multiqc_configuration = config["parameters"]["multiqc"]
        multiqc = sb_prog.Run(command=command, path=path_map, image=multiqc_configuration["image"],
                              tag=multiqc_configuration["tag"])
        multiqc.run()

    if program == 'mugsy':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mugsy " + arg_string
        mugsy_configuration = config["parameters"]["mugsy"]
        mugsy = sb_prog.Run(command=command, path=path_map, image=mugsy_configuration["image"],
                            tag=mugsy_configuration["tag"])
        mugsy.run()

    if program == 'mlst':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mlst " + arg_string
        mlst_configuration = config["parameters"]["mlst"]
        mlst = sb_prog.Run(command=command, path=path_map, image=mlst_configuration["image"], tag=mlst_configuration["tag"])
        mlst.run()

    if program == 'medaka':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "medaka " + arg_string
        medaka_configuration = config["parameters"]["medaka"]
        medaka = sb_prog.Run(command=command, path=path_map, image=medaka_configuration["image"],
                             tag=medaka_configuration["tag"])
        medaka.run()

    if program == 'mashtree':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mashtree " + arg_string
        mashtree_configuration = config["parameters"]["mashtree"]
        mashtree = sb_prog.Run(command=command, path=path_map, image=mashtree_configuration["image"],
                               tag=mashtree_configuration["tag"])
        mashtree.run()

    if program == 'legsta':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "legsta " + arg_string
        legsta_configuration = config["parameters"]["legsta"]
        legsta = sb_prog.Run(command=command, path=path_map, image=legsta_configuration["image"],
                             tag=legsta_configuration["tag"])
        legsta.run()

    if program == 'ksnp3':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = ""
        command = "kSNP3 " + arg_string
        ksnp3_configuration = config["parameters"]["ksnp3"]
        ksnp3 = sb_prog.Run(command=command, path=path_map, image=ksnp3_configuration["image"], tag=ksnp3_configuration["tag"])
        ksnp3.run()

    if program == 'kraken2-build':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken2-build " + arg_string
        kraken2_configuration = config["parameters"]["kraken2"]
        kraken2 = sb_prog.Run(command=command, path=path_map, image=kraken2_configuration["image"], tag=kraken2_configuration["tag"])
        kraken2.run()

    if program == 'kraken2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken2 " + arg_string
        kraken2_configuration = config["parameters"]["kraken2"]
        kraken2 = sb_prog.Run(command=command, path=path_map, image=kraken2_configuration["image"], tag=kraken2_configuration["tag"])
        kraken2.run()

    if program == 'kraken-build':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken-build " + arg_string
        kraken_configuration = config["parameters"]["kraken"]
        kraken = sb_prog.Run(command=command, path=path_map, image=kraken_configuration["image"], tag=kraken_configuration["tag"])
        kraken.run()

    if program == 'kraken':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken " + arg_string
        kraken_configuration = config["parameters"]["kraken"]
        kraken = sb_prog.Run(command=command, path=path_map, image=kraken_configuration["image"], tag=kraken_configuration["tag"])
        kraken.run()

    if program == 'kma':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kma " + arg_string
        kma_configuration = config["parameters"]["kma"]
        kma = sb_prog.Run(command=command, path=path_map, image=kma_configuration["image"], tag=kma_configuration["tag"])
        kma.run()

    if program == 'flye':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "flye " + arg_string
        flye_configuration = config["parameters"]["flye"]
        flye = sb_prog.Run(command=command, path=path_map, image=flye_configuration["image"], tag=flye_configuration["tag"])
        flye.run()

    if program == 'filtlong':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "filtlong " + arg_string
        filtlong_configuration = config["parameters"]["filtlong"]
        filtlong = sb_prog.Run(command=command, path=path_map, image=filtlong_configuration["image"], tag=filtlong_configuration["tag"])
        filtlong.run()

    if program == 'fastqc':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "fastqc " + arg_string
        fastqc_configuration = config["parameters"]["fastqc"]
        fastqc = sb_prog.Run(command=command, path=path_map, image=fastqc_configuration["image"],
                             tag=fastqc_configuration["tag"])
        fastqc.run()

    if program == 'fasttree':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "FastTree " + arg_string
        fasttree_configuration = config["parameters"]["fasttree"]
        fasttree = sb_prog.Run(command=command, path=path_map, image=fasttree_configuration["image"], tag=fasttree_configuration["tag"])
        fasttree.run()

    if program == 'fastani':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "fastANI " + arg_string
        fastani_configuration = config["parameters"]["fastani"]
        fastani = sb_prog.Run(command=command, path=path_map, image=fastani_configuration["image"], tag=fastani_configuration["tag"])
        fastani.run()

    if program == 'emm-typing-tool':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "emm_typing.py " + arg_string
        emm_typing_tool_configuration = config["parameters"]["emm-typing-tool"]
        emm_typing_tool = sb_prog.Run(command=command, path=path_map, image=emm_typing_tool_configuration["image"], tag=emm_typing_tool_configuration["tag"])
        emm_typing_tool.run()

    if program == 'circlator':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "circlator " + arg_string
        circlator_configuration = config["parameters"]["circlator"]
        circlator = sb_prog.Run(command=command, path=path_map, image=circlator_configuration["image"], tag=circlator_configuration["tag"])
        circlator.run()

    if program == 'cfsan-snp':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "cfsan_snp_pipeline " + arg_string
        cfsan_snp_configuration = config["parameters"]["cfsan-snp-pipeline"]
        cfsan_snp = sb_prog.Run(command=command, path=path_map, image=cfsan_snp_configuration["image"], tag=cfsan_snp_configuration["tag"])
        cfsan_snp.run()

    if program == 'canu-racon':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "racon " + arg_string
        canu_racon_configuration = config["parameters"]["canu-racon"]
        canu_racon = sb_prog.Run(command=command, path=path_map, image=canu_racon_configuration["image"], tag=canu_racon_configuration["tag"])
        canu_racon.run()

    if program == 'bbtools':
        command = " " + arg_string
        bbtools_configuration = config["parameters"]["bbtools"]
        bbtools = sb_prog.Run(command=command, path=path_map, image=bbtools_configuration["image"], tag=bbtools_configuration["tag"])
        if not re.search('[a-zA-Z]', arg_string):
            print("BBTools shell script must be specified, e.g. staphb-tk bbtools bbmap.sh, staphb-tk bbtools bbduk.sh, etc. \n\nMore info on BBTools at https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/.")
        else:
            bbtools.run()

    if program == 'raxml':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "raxmlHPC "+arg_string
        raxml_configuration = config["parameters"]["raxml"]
        raxml = sb_prog.Run(command=command, path=path_map, image=raxml_configuration["image"], tag =raxml_configuration["tag"])
        raxml.run()

    if program == 'spades':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "spades.py "+arg_string
        spades_configuration = config["parameters"]["spades"]
        spades = sb_prog.Run(command=command, path=path_map, image=spades_configuration["image"], tag =spades_configuration["tag"])
        spades.run()

    if program == 'mash':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mash "+arg_string
        mash_configuration = config["parameters"]["mash"]
        mash = sb_prog.Run(command=command, path=path_map, image=mash_configuration["image"], tag = mash_configuration["tag"])
        mash.run()

    if program == 'seqyclean':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "seqyclean "+arg_string
        seqyclean_configuration = config["parameters"]["seqyclean"]
        seqyclean = sb_prog.Run(command=command, path=path_map, image=seqyclean_configuration["image"], tag = seqyclean_configuration["tag"])
        seqyclean.run()

    if program == 'shovill':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "shovill " + arg_string
        shovill_configuration = config["parameters"]["shovill"]
        shovill = sb_prog.Run(command=command, path=path_map, image=shovill_configuration["image"], tag=shovill_configuration["tag"])
        shovill.run()

    if program == 'prokka':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "prokka " + arg_string
        prokka_configuration = config["parameters"]["prokka"]
        prokka = sb_prog.Run(command=command, path=path_map, image=prokka_configuration["image"], tag=prokka_configuration["tag"])
        prokka.run()

    if program == 'clustalo':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "clustalo " + arg_string
        clustalo_configuration = config["parameters"]["clustalo"]
        clustalo = sb_prog.Run(command=command, path=path_map, image=clustalo_configuration["image"], tag=clustalo_configuration["tag"])
        clustalo.run()

    if program == 'abricate':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "abricate " + arg_string
        abricate_configuration = config["parameters"]["abricate"]
        abricate = sb_prog.Run(command=command, path=path_map, image=abricate_configuration["image"], tag=abricate_configuration["tag"])
        abricate.run()

    if program == 'augur':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "augur " + arg_string
        augur_configuration = config["parameters"]["augur"]
        augur = sb_prog.Run(command=command, path=path_map, image=augur_configuration["image"], tag=augur_configuration["tag"])
        augur.run()

    if program == 'iqtree':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "iqtree " + arg_string
        iqtree_configuration = config["parameters"]["iqtree"]
        iqtree = sb_prog.Run(command=command, path=path_map, image=iqtree_configuration["image"], tag=iqtree_configuration["tag"])
        iqtree.run()

    if program == 'lyveset':
        command = "" + arg_string
        lyveset_configuration = config["parameters"]["lyveset"]
        lyveset = sb_prog.Run(command=command, path=path_map, image=lyveset_configuration["image"], tag=lyveset_configuration["tag"])
        if not re.search('[a-zA-Z]', arg_string):
            print("Lyev-SET perl script must be specified, e.g. staphb-tk lyveset launch_set.pl, staphb-tk lyveset set_manage.pl, staphb-tk lyveset run_assembly_readMeterics.pl. \n\nMore info on Lyve-SET usage at: github.com/lskatz/lyve-SET.")
        else:
            lyveset.run()

    if program == 'quast':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "quast.py " + arg_string
        quast_configuration = config["parameters"]["quast"]
        quast = sb_prog.Run(command=command, path=path_map, image=quast_configuration["image"], tag=quast_configuration["tag"])
        quast.run()

    if program == 'roary':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "roary " + arg_string
        roary_configuration = config["parameters"]["roary"]
        roary = sb_prog.Run(command=command, path=path_map, image=roary_configuration["image"], tag=roary_configuration["tag"])
        roary.run()

    if program == 'seqsero':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "SeqSero.py " + arg_string
        seqsero_configuration = config["parameters"]["seqsero"]
        seqsero = sb_prog.Run(command=command, path=path_map, image=seqsero_configuration["image"], tag=seqsero_configuration["tag"])
        seqsero.run()

    if program == 'samtools':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = ""
        command = "samtools " + arg_string
        samtools_configuration = config["parameters"]["samtools"]
        samtools = sb_prog.Run(command=command, path=path_map, image=samtools_configuration["image"], tag=samtools_configuration["tag"])
        samtools.run()

    if program == 'serotypefinder':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "serotypefinder.pl " + arg_string
        serotypefinder_configuration = config["parameters"]["serotypefinder"]
        serotypefinder = sb_prog.Run(command=command, path=path_map, image=serotypefinder_configuration["image"], tag=serotypefinder_configuration["tag"])
        serotypefinder.run()

    if program == 'bwa':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = ""
        command = "bwa " + arg_string
        bwa_configuration = config["parameters"]["bwa"]
        bwa = sb_prog.Run(command=command, path=path_map, image=bwa_configuration["image"], tag=bwa_configuration["tag"])
        bwa.run()

    if program == 'minimap2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "minimap2 " + arg_string
        minimap2_configuration = config["parameters"]["minimap2"]
        minimap2 = sb_prog.Run(command=command, path=path_map, image=minimap2_configuration["image"], tag=minimap2_configuration["tag"])
        minimap2.run()

    if program == 'unicycler':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "unicycler " + arg_string
        unicycler_configuration = config["parameters"]["unicycler"]
        unicycler = sb_prog.Run(command=command, path=path_map, image=unicycler_configuration["image"], tag=unicycler_configuration["tag"])
        unicycler.run()
