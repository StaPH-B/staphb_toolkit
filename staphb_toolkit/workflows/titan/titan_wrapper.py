#!/usr/bin/env python3

import re
import os,sys
import json
from pathlib import Path
from shutil import copy

#search a directory and return a dictonary of paired fastq files named by fastq name
def pe_search(path,bedfile_path):
    if not os.path.isdir(path):
        raise ValueError(path + " " + "not found.")
    path = os.path.abspath(path)
    fastqs = []
    fastq_dict = []
    for root,dirs,files in os.walk(path):
        for file in files:
            if '.fastq' in file or '.fastq.gz' in file or '.fq.gz' in file or '.fq' in file:
                fastqs.append(os.path.join(root,file))
    fastqs.sort()
    if not len(fastqs)%2 == 0:
        print("There is an uneven number of Fastq files in: '"+path+"'")
        sys.exit(1)

    for readA, readB in zip(fastqs[0::2], fastqs[1::2]):
        samplename = re.split("(_R1)|(_1)",os.path.basename(readA))[0]
        fastq_dict.append({"samplename":samplename,"read1_raw":readA,"read2_raw":readB,"primer_bed":bedfile_path})
    return(fastq_dict)

#search a directory and return a dictonary of fastq files (single end) named by fastq name
def se_search(path,bedfile_path):
    if not os.path.isdir(path):
        raise ValueError(path + " " + "not found.")
    path = os.path.abspath(path)
    fastqs = []
    fastq_dict = {}
    for root,dirs,files in os.walk(path):
        for file in files:
            if '.fastq' in file or '.fastq.gz' in file or '.fq.gz' in file or '.fq' in file:
                fastqs.append(os.path.join(root,file))
    fastqs.sort()
    for read in fastqs:
        samplename = re.split("(.fastq)|(.fq)",os.path.basename(read))[0]
        fastq_dict[samplename] = {"read":read,"primer_bed":bedfile_path}
    return(fastq_dict)

#search a path and jsonify the reads as input for wdl
def collect_input_data(reads_path,bedfile_path,pe = True):
    if pe == True:
        sample_dict = pe_search(reads_path,bedfile_path)
    else:
        sample_dict = se_search(reads_path,bedfile_path)

    input_data = {"cli_wrapper.inputSamples":sample_dict}
    return(input_data)

#grab and move data based on the WDL generated metadata.json file, this allows us to collect the output and organize it for the user
#NOTE: this will need to be updated if any additional steps are added to the workflow
def parseOutputMetadata(metaJSON,outPath):
    #prep output path
    Path(outPath).mkdir(parents=True,exist_ok=True)
    #load metadata json
    with open(metaJSON,'r') as jsonfile:
        data = json.load(jsonfile)

    #get sample ids
    sids = []
    for sample in data['inputs']['inputSamples']:
        sids.append(sample['samplename'])

    #copy final summary
    copy(data['outputs']['cli_wrapper.merged_metrics'],outPath)

    #copy consensus files
    p = os.path.join(outPath,"consensus/")
    Path(p).mkdir(parents=True,exist_ok=True)
    for f in data['outputs']['cli_wrapper.consensus_seq']:
        copy(f,p)

    #clean reads
    p = os.path.join(outPath,"cleaned_reads/")
    Path(p).mkdir(parents=True,exist_ok=True)
    for f in data['outputs']['cli_wrapper.read1_clean']:
        copy(f,p)
    for f in data['outputs']['cli_wrapper.read2_clean']:
        copy(f,p)

    #stats
    c = 0
    while c < len(sids):
        p = os.path.join(outPath,f"stats/{sids[c]}/")
        Path(p).mkdir(parents=True,exist_ok=True)
        f = data['outputs']['cli_wrapper.cov_hist'][c]
        copy(f,p)
        f = data['outputs']['cli_wrapper.amp_coverage'][c]
        copy(f,p)
        f = data['outputs']['cli_wrapper.samtools_stats'][c]
        copy(f,p)
        f = data['outputs']['cli_wrapper.cov_stats'][c]
        copy(f,p)
        f = data['outputs']['cli_wrapper.samtools_flagstat'][c]
        copy(f,p)
        f = data['outputs']['cli_wrapper.kraken_report'][c]
        copy(f,p)
        f = data['outputs']['cli_wrapper.pango_lineage_report'][c]
        copy(f,p)
        c += 1

    #untrimmed bam
    p = os.path.join(outPath,"untrimmed_bam/")
    Path(p).mkdir(parents=True,exist_ok=True)
    for f in data['outputs']['cli_wrapper.sorted_bam']:
        copy(f,p)
    for f in data['outputs']['cli_wrapper.sorted_bai']:
        copy(f,p)

    #trimmed bam
    p = os.path.join(outPath,"trimmed_bam/")
    Path(p).mkdir(parents=True,exist_ok=True)
    for f in data['outputs']['cli_wrapper.trim_sorted_bam']:
        copy(f,p)
    for f in data['outputs']['cli_wrapper.trim_sorted_bai']:
        copy(f,p)

#merge optional inputs with inputs and output json string
def mergeOptionalInputs(inputDATA,optionsDATA):
    merged_data = {**inputDATA, **optionsDATA}
    return(json.dumps(merged_data,sort_keys=True, indent=4))
