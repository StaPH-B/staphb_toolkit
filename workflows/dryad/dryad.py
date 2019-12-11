#!/usr/bin/env python3
#Author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#description: A pipeline for constructing SNP based and  core gene set reference free phylogenies

import subprocess as sub
import sys, os
from shutil import copyfile

#local app libraries
from workflows.dryad.app.lib import cpu_count, checkexists
from workflows.dryad.app.lib import StatusTracker as ST
from workflows.dryad.app.trimming import q_trim
from workflows.dryad.app.core_genome import core_genome
from workflows.dryad.app.snp import snp

#define various pipelines
def dryad_cg(paired_reads,jobs,cpu_job,outdir,tracker):
    if not tracker.check_status('trimmed'):
        #start trimming
        q_trim(paired_reads,jobs,cpu_job,outdir,tracker)

    if tracker.check_status('cg'):
        #start core-genome pipeline
        core_genome(jobs,cpu_job,outdir,tracker)

def dryad_snp(paired_reads,jobs,cpu_job,outdir,tracker,reference):
    if not tracker.check_status('trimmed'):
        #start trimming
        q_trim(paired_reads,jobs,cpu_job,outdir,tracker)

    if tracker.check_status('snp'):
        #start snp pipeline
        snp(jobs,cpu_job,outdir,reference,tracker)

def dryad_all(paired_reads,jobs,cpu_job,outdir,tracker,reference):
    if not tracker.check_status('trimmed'):
        #start trimming
        q_trim(paired_reads,jobs,cpu_job,outdir,tracker)

    if tracker.check_status('cg'):
        #start core-genome pipeline
        core_genome(jobs,cpu_job,outdir,tracker)

    if tracker.check_status('snp'):
        #start snp pipeline
        snp(jobs,cpu_job,outdir,reference,tracker)

#main app function
def dryad(pipeline,threads,output_path,reads_path,ref_path=''):
    #get num of jobs and number of cpus per job
    jobs,cpu_job = cpu_count(threads)

    #get current working dir if output is empty
    try:
        out = os.path.abspath(output_path)
    except (AttributeError, TypeError) as err:
        out = os.getcwd()

    #open file and pull locations into a list
    with open(reads_path,'r') as f:
        r_list = []
        for line in f:
            if line.strip() != '':
                r_list.append(line.strip())

    #sort and join pairs
    r_list.sort()
    if len(r_list) % 2 != 0:
        print('There is an uneven number of read pairs in the read list. Exiting.')
        sys.exit()
    paired_reads = []
    [paired_reads.append([x,y]) for x,y in zip(r_list[0::2],r_list[1::2])]

    #initialize tracker
    tracker = ST()
    outdir = tracker.initialize(out,pipeline)

    #run pipeline
    if pipeline == 'cg':
        print(paired_reads,jobs,cpu_job,outdir,tracker)
        dryad_cg(paired_reads,jobs,cpu_job,outdir,tracker)
    elif pipeline == 'snp':
        print(paired_reads,jobs,cpu_job,outdir,tracker,ref_path)
        dryad_snp(paired_reads,jobs,cpu_job,outdir,tracker)
    elif pipeline == 'all':
        print(paired_reads,jobs,cpu_job,outdir,tracker,ref_path)
        dryad_all(paired_reads,jobs,cpu_job,outdir,tracker)
    else:
        print("not a valid pipeline")
        sys.exit(1)
