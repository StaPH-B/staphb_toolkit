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
'ariba':'Ariba - Antimicrobial Resistance Identification By Assembly',
'artic-ncov2019-medaka':'Artic conda environment for SC2 assembly with medaka',
'artic-ncov2019-nanopolish':'Artic conda environment for SC2 assembly with nanopolish',
'augur': 'Augur - Pipeline components for real-time phylodynamic analysis',
'bbtools':'BBTools - Suite of fast, multithreaded bioinformatics tools for DNA and RNA sequence data',
'bcftools':'BCFTools - Variant calling and manipulating files in the Variant Call Format (VCF) and its binary counterpart BCF',
'bwa':'BWA - mapping low-divergent sequences against a large reference genome',
'canu':'Canu = Long read assembly and polishing tools"',
'canu-racon':'Canu-Racon - Ultrafast consensus module for raw de novo assembly of long, uncorrected reads.',
'centroid' : 'centroid - a tool for determining an ideal reference genome from a set of fasta files' ,
'cfsan-snp':'CFSAN-SNP - SNP calling pipeline from the FDA CFSAN laboratory',
'circlator':'Circlator - A tool to circularize genome assemblies',
'clustalo':'ClustalO - A fast multiple sequence alignment program',
'colorid':'Colorid - Experiments with using BIGSI data structure for metagenomic and QC applications',
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
'kraken2-inspect': 'Kraken2-inspect - Inspect a kraken2 database',
'ksnp3':'kSNP3 - Identifies the pan-genome SNPs in a set of genome sequences, and estimates phylogenetic trees based upon those SNPs.',
'legsta':'Legsta - In silico Legionella pneumophila Sequence Based Typing',
'lyveset':'LYVE-SET - a method of using hqSNPs to create a phylogeny.',
'mafft':'MAFFT - multiple sequence alignment program for amino acid or nucleotide sequences',
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
'pangolin':'Pangolin - Phylogenetic Assignment of Named Global Outbreak LINeages',
'pilon':'Pilon - Automated genome assembly improvement and variant detection tool',
'plasmidseeker':'PlasmidSeeker - A k-mer based program for the identification of known plasmids from whole-genome sequencing reads',
'prokka':'Prokka - Rapid prokaryotic genome annotation',
'quast':'Quast - Genome assembly evaluation tool.',
'racon':'Racon - Long read assembly and polishing tools',
'rasusa':'RASUA - Randomly subsample sequencing reads to a specified coverage',
'raxml':'RAxML -Maximum likelihood tree builder.',
'roary':'Roary - Rapid large-scale prokaryote pan genome analysis.',
'salmid':'SalmID - Rapid confirmation of Salmonella spp. and subspp. from sequence data',
'samtools':'Samtools - A suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories.',
'seqsero':'SeqSero - Salmonella serotyping from genome sequencing data.',
'seqsero2':'SeqSero2 - Salmonella serotype prediction from genome sequencing data.',
'seqtk': 'SeqTK - Toolkit for processing sequences in FASTA/Q formats',
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
'trimmomatic':'Trimmomatic - Flexible read trimming tool for Illumina NGS data',
'unicycler':'Unicycler - an assembly pipeline for bacterial genomes.',
'vadr':'VADR - Software for viral annotations',
'vibrant': 'VIBRANT - a tool for automated recovery and annotation of bacterial and archaeal viruses, determination of genome completeness, and characterization of viral community function from metagenomic assemblies',
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

    parser = MyParser(description=f"StaPH-B ToolKit Programs v{autoupdate.version}",usage="staphb-tk [optional arguments] <application> [application arguments]",add_help=True)
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
    parser_ariba = subparsers.add_parser('ariba', add_help=False)
    parser_artic_ncov2019_medaka = subparsers.add_parser('artic-ncov2019-medaka', add_help=False)
    parser_artic_ncov2019_nanopolish = subparsers.add_parser('artic-ncov2019-nanopolish', add_help=False)
    parser_augur = subparsers.add_parser('augur', add_help=False)
    parser_bbtools = subparsers.add_parser('bbtools', add_help=False)
    parser_bcftools = subparsers.add_parser('bcftools', add_help=False)
    parser_bwa = subparsers.add_parser('bwa', add_help=False)
    parser_canu = subparsers.add_parser('canu', add_help=False)
    parser_canuracon = subparsers.add_parser('canu-racon', add_help=False)
    parser_centroid = subparsers.add_parser('centroid', add_help=False)
    parser_cfsansnp = subparsers.add_parser('cfsan-snp', add_help=False)
    parser_circlator = subparsers.add_parser('circlator', add_help=False)
    parser_clustalo = subparsers.add_parser('clustalo', add_help=False)
    parser_colorid = subparsers.add_parser('colorid', add_help=False)
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
    parser_kraken2inspect = subparsers.add_parser('kraken2-inspect', add_help=False)
    parser_ksnp3 = subparsers.add_parser('ksnp3', add_help=False)
    parser_legsta = subparsers.add_parser('legsta', add_help=False)
    parser_lyveset = subparsers.add_parser('lyveset', add_help=False)
    parser_mafft = subparsers.add_parser('mafft', add_help=False)
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
    parser_pangolin = subparsers.add_parser('pangolin', add_help=False)
    parser_pilon = subparsers.add_parser('pilon', add_help=False)
    parser_plasmidseeker = subparsers.add_parser('plasmidseeker', add_help=False)
    parser_prokka = subparsers.add_parser('prokka', add_help=False)
    parser_quast = subparsers.add_parser('quast', add_help=False)
    parser_racon = subparsers.add_parser('racon', add_help=False)
    parser_rasusa = subparsers.add_parser('rasusa', add_help=False)
    parser_raxml = subparsers.add_parser('raxml', add_help=False)
    parser_roary = subparsers.add_parser('roary', add_help=False)
    parser_salmid = subparsers.add_parser('salmid', add_help=False)
    parser_samtools = subparsers.add_parser('samtools', add_help=False)
    parser_seqsero = subparsers.add_parser('seqsero', add_help=False)
    parser_seqsero2 = subparsers.add_parser('seqsero2', add_help=False)
    parser_seqtk = subparsers.add_parser('seqtk', add_help=False)
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
    parser_vadr = subparsers.add_parser('vadr', add_help=False)
    parser_vibrant = subparsers.add_parser('vibrant', add_help=False)
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
        program_configuration = config["parameters"]["ivar-SC2"]

    if program == 'ivar':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = " "
        command = "ivar " + arg_string
        program_configuration = config["parameters"]["ivar"]

    if program == 'wtdbg2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "wtdbg2 " + arg_string
        program_configuration = config["parameters"]["wtdbg2"]

    if program == 'trimmomatic':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "trimmomatic " + arg_string
        program_configuration = config["parameters"]["trimmomatic"]

    if program == 'vibrant':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "python3 /VIBRANT/VIBRANT_run.py " + arg_string
        program_configuration = config["parameters"]["vibrant"]

    if program == 'vadr':
        if not re.search('[a-zA-Z]', arg_string):
            print("VADR perl script must be specified, e.g. staphb-tk vadr v-build.pl or staphb-tk vadr v-annotate.pl. \n\nMore info on VADR at https://github.com/ncbi/vadr.")
            sys.exit()
        command = " " + arg_string
        program_configuration = config["parameters"]["vadr"]

    if program == 'tiptoft':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "tiptoft " + arg_string
        program_configuration = config["parameters"]["tiptoft"]

    if program == 'staramr':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "staramr " + arg_string
        program_configuration = config["parameters"]["staramr"]

    if program == 'sra-toolkit':
        if not re.search('[a-zA-Z]', arg_string):
            print("SRA toolkit tool must be specified, e.g. staphb-tk sra-toolkit fasterq-dump, staphb-tk sra-toolkit sra-pileup, etc. \n\nMore info on SRA Toolkit usage at: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc.")
            sys.exit()
        command = " " + arg_string
        program_configuration = config["parameters"]["sra-toolkit"]

    if program == 'snp-dists':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "snp-dists " + arg_string
        program_configuration = config["parameters"]["snp-dists"]

    if program == 'snp-sites':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "snp-sites " + arg_string
        program_configuration = config["parameters"]["snp-sites"]

    if program == 'snippy':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "snippy " + arg_string
        program_configuration = config["parameters"]["snippy"]

    if program == 'skesa':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "skesa " + arg_string
        program_configuration = config["parameters"]["skesa"]

    if program == 'sistr':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "sistr " + arg_string
        program_configuration = config["parameters"]["sistr"]

    if program == 'seroba':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "seroba " + arg_string
        program_configuration = config["parameters"]["seroba"]

    if program == 'seqsero2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "SeqSero2_package.py " + arg_string
        program_configuration = config["parameters"]["seqsero2"]

    if program == 'seqtk':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = ""
        command = "seqtk " + arg_string
        program_configuration = config["parameters"]["seqtk"]

    if program == 'salmid':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "SalmID.py " + arg_string
        program_configuration = config["parameters"]["salmid"]

    if program == 'rasusa':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "rasusa " + arg_string
        program_configuration = config["parameters"]["rasusa"]

    if program == 'plasmidseeker':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "plasmidseeker.pl " + arg_string
        program_configuration = config["parameters"]["plasmidseeker"]

    if program == 'pangolin':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "--help"
        command = "pangolin " + arg_string
        program_configuration = config["parameters"]["pangolin"]

    if program == 'pilon':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "--help"
        command = "pilon " + arg_string
        program_configuration = config["parameters"]["pilon"]

    if program == 'orthofinder':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "orthofinder " + arg_string
        program_configuration = config["parameters"]["orthofinder"]

    if program == 'ncbi-amrfinder-plus':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "--help"
        command = "amrfinder " + arg_string
        program_configuration = config["parameters"]["ncbi-amrfinder-plus"]

    if program == 'nanoplot':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "NanoPlot " + arg_string
        program_configuration = config["parameters"]["nanoplot"]

    if program == 'multiqc':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "multiqc " + arg_string
        program_configuration = config["parameters"]["multiqc"]

    if program == 'mugsy':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mugsy " + arg_string
        program_configuration = config["parameters"]["mugsy"]

    if program == 'mlst':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mlst " + arg_string
        program_configuration = config["parameters"]["mlst"]

    if program == 'medaka':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "medaka " + arg_string
        program_configuration = config["parameters"]["medaka"]

    if program == 'mashtree':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mashtree " + arg_string
        program_configuration = config["parameters"]["mashtree"]

    if program == 'legsta':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "legsta " + arg_string
        program_configuration = config["parameters"]["legsta"]

    if program == 'ksnp3':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = ""
        command = "kSNP3 " + arg_string
        program_configuration = config["parameters"]["ksnp3"]

    if program == 'kraken2-inspect':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken2-inspect " + arg_string
        program_configuration = config["parameters"]["kraken2"]

    if program == 'kraken2-build':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken2-build " + arg_string
        program_configuration = config["parameters"]["kraken2"]

    if program == 'kraken2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken2 " + arg_string
        program_configuration = config["parameters"]["kraken2"]

    if program == 'kraken-build':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken-build " + arg_string
        program_configuration = config["parameters"]["kraken"]

    if program == 'kraken':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kraken " + arg_string
        program_configuration = config["parameters"]["kraken"]

    if program == 'kma':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "kma " + arg_string
        program_configuration = config["parameters"]["kma"]

    if program == 'flye':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "flye " + arg_string
        program_configuration = config["parameters"]["flye"]

    if program == 'filtlong':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "filtlong " + arg_string
        program_configuration = config["parameters"]["filtlong"]

    if program == 'fastqc':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "fastqc " + arg_string
        program_configuration = config["parameters"]["fastqc"]

    if program == 'fasttree':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "FastTree " + arg_string
        program_configuration = config["parameters"]["fasttree"]

    if program == 'fastani':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "fastANI " + arg_string
        program_configuration = config["parameters"]["fastani"]

    if program == 'emm-typing-tool':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "emm_typing.py " + arg_string
        program_configuration = config["parameters"]["emm-typing-tool"]

    if program == 'circlator':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "circlator " + arg_string
        program_configuration = config["parameters"]["circlator"]

    if program == 'cfsan-snp':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "cfsan_snp_pipeline " + arg_string
        program_configuration = config["parameters"]["cfsan-snp-pipeline"]

    if program == 'canu-racon':
        if not re.search('[a-zA-Z]', arg_string):
            print("This is a bundled application that requires a specific commands to be used (i.e. staphb-tk canu-racon canu -h) please see the documentation for Canu, Minimap2 and Racon to use.")
            sys.exit()
        command = " " + arg_string
        program_configuration = config["parameters"]["canu-racon"]

    if program == 'canu':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "canu " + arg_string
        program_configuration = config["parameters"]["canu"]

    if program == 'bbtools':
        if not re.search('[a-zA-Z]', arg_string):
            print("BBTools shell script must be specified, e.g. staphb-tk bbtools bbmap.sh, staphb-tk bbtools bbduk.sh, etc. \n\nMore info on BBTools at https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/.")
            sys.exit()
        command = " " + arg_string
        program_configuration = config["parameters"]["bbtools"]

    if program == 'bcftools':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "bcftools "+arg_string
        program_configuration = config["parameters"]["bcftools"]

    if program == 'raxml':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "raxmlHPC "+arg_string
        program_configuration = config["parameters"]["raxml"]

    if program == 'spades':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "spades.py "+arg_string
        program_configuration = config["parameters"]["spades"]

    if program == 'mash':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mash "+arg_string
        program_configuration = config["parameters"]["mash"]

    if program == 'seqyclean':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "seqyclean "+arg_string
        program_configuration = config["parameters"]["seqyclean"]

    if program == 'shovill':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "shovill " + arg_string
        program_configuration = config["parameters"]["shovill"]

    if program == 'prokka':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "prokka " + arg_string
        program_configuration = config["parameters"]["prokka"]

    if program == 'clustalo':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "clustalo " + arg_string
        program_configuration = config["parameters"]["clustalo"]

    if program == 'colorid':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "colorid " + arg_string
        program_configuration = config["parameters"]["colorid"]

    if program == 'abricate':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "abricate " + arg_string
        program_configuration = config["parameters"]["abricate"]

    if program == 'ariba':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "ariba " + arg_string
        program_configuration = config["parameters"]["ariba"]

    if program == 'artic-ncov2019-medaka':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "artic minion " + arg_string
        program_configuration = config["parameters"]["artic-ncov2019-medaka"]

    if program == 'artic-ncov2019-nanopolish':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "artic minion " + arg_string
        program_configuration = config["parameters"]["artic-ncov2019-medaka"]

    if program == 'augur':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "augur " + arg_string
        program_configuration = config["parameters"]["augur"]

    if program == 'iqtree':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "iqtree " + arg_string
        program_configuration = config["parameters"]["iqtree"]

    if program == 'lyveset':
        if not re.search('[a-zA-Z]', arg_string):
            print("Lyev-SET perl script must be specified, e.g. staphb-tk lyveset launch_set.pl, staphb-tk lyveset set_manage.pl, staphb-tk lyveset run_assembly_readMeterics.pl. \n\nMore info on Lyve-SET usage at: github.com/lskatz/lyve-SET.")
            sys.exit()
        command = "" + arg_string
        program_configuration = config["parameters"]["lyveset"]

    if program == 'quast':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "quast.py " + arg_string
        program_configuration = config["parameters"]["quast"]

    if program == 'racon':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "racon " + arg_string
        program_configuration = config["parameters"]["racon"]

    if program == 'roary':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "roary " + arg_string
        program_configuration = config["parameters"]["roary"]

    if program == 'seqsero':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "SeqSero.py " + arg_string
        program_configuration = config["parameters"]["seqsero"]

    if program == 'samtools':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = ""
        command = "samtools " + arg_string
        program_configuration = config["parameters"]["samtools"]

    if program == 'serotypefinder':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "serotypefinder.pl " + arg_string
        program_configuration = config["parameters"]["serotypefinder"]

    if program == 'bwa':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = ""
        command = "bwa " + arg_string
        program_configuration = config["parameters"]["bwa"]

    if program == 'minimap2':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "minimap2 " + arg_string
        program_configuration = config["parameters"]["minimap2"]

    if program == 'centroid':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "centroid.py " + arg_string
        program_configuration = config["parameters"]["centroid"]

    if program == 'unicycler':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "unicycler " + arg_string
        program_configuration = config["parameters"]["unicycler"]

    if program == 'mafft':
        if not re.search('[a-zA-Z]', arg_string):
            arg_string = "-h"
        command = "mafft " + arg_string
        program_configuration = config["parameters"]["mafft"]

    #Run the program
    #-----------------------------------------
    program_object = sb_prog.Run(command=command, path=path_map, image=program_configuration["image"], tag=program_configuration["tag"])
    program_object.run()
