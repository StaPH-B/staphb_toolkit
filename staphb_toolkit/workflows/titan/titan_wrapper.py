#!/usr/bin/env python3

import re
import os,sys
import json

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

#search a directory and return a dictonary of fastq files named by fastq name
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

def create_input_json(reads_path,bedfile_path,pe = True):
    if pe == True:
        sample_dict = pe_search(reads_path,bedfile_path)
    else:
        sample_dict = se_search(reads_path,bedfile_path)

    input_json = {"cli_wrapper.inputSamples":sample_dict}
    return(json.dumps(input_json,sort_keys=True, indent=4))

if __name__ == '__main__' :
    with open("titan_input.json",'w') as outfile:
        outfile.write(create_input_json(sys.argv[1],"/test/path/bed.bed"))
