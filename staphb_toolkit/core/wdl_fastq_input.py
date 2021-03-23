#!/usr/bin/env python3

import re
import os,sys

#search a directory and return a dictonary of paired fastq files named by fastq name
def pe_search(path):
    if not os.path.isdir(path):
        raise ValueError(path + " " + "not found.")
    path = os.path.abspath(path)
    fastqs = []
    fastq_dict = {}
    for root,dirs,files in os.walk(path):
        for file in files:
            if '.fastq' in file or '.fastq.gz' in file or '.fq.gz' in file or '.fq' in file:
                fastqs.append(file)
    fastqs.sort()
    if not len(fastqs)%2 == 0:
        print("There is an uneven number of Fastq files in: '"+path+"'")
        sys.exit(1)

    for readA, readB in zip(fastqs[0::2], fastqs[1::2]):
        samplename = re.split("(_R1)|(_1)",readA)[0]
        fastq_dict[samplename] = {"read_1":readA,"read_2":readB}
    return(fastq_dict)

#search a directory and return a dictonary of fastq files named by fastq name
def se_search(path):
    if not os.path.isdir(path):
        raise ValueError(path + " " + "not found.")
    path = os.path.abspath(path)
    fastqs = []
    fastq_dict = {}
    for root,dirs,files in os.walk(path):
        for file in files:
            if '.fastq' in file or '.fastq.gz' in file or '.fq.gz' in file or '.fq' in file:
                fastqs.append(file)
    fastqs.sort()
    for read in fastqs:
        samplename = re.split("(.fastq)|(.fq)",read)[0]
        fastq_dict[samplename] = read
    return(fastq_dict)
