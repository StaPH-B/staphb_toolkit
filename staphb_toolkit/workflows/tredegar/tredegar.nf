#!/usr/bin/env nextflow

//Description: Workflow for quality control of raw illumina reads
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs(  "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nEnsure read data is compressed (gzip) before rerunning." }
    .set { raw_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(reads) into raw_reads_mash, raw_reads_trim, raw_reads_qc, raw_reads_gas, raw_reads_salmonella

  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
  """
  """
  }
}
//Create Mash Sketches for all isolates
process mash_dist{
  tag "$name"

  input:
  set val(name), file(reads) from raw_reads_mash

  output:
  tuple val("${name}"), file("${name}_top_hits.tab") into mash_dist

  script:
  """
  mash sketch -r -m 2 -o ${name}_sketch.msh ${reads}
  mash dist ${params.mash_db} ${name}_sketch.msh > ${name}_distance.tab && sort -gk3 ${name}_distance.tab | head > ${name}_top_hits.tab
  """
}

//Curate top Mash Tredegar_results
process mash_species {
  publishDir "${params.outdir}/logs/mash/", mode: 'copy'

  input:
  file(mash_hits) from mash_dist.collect()

  output:
  file "mash_species.tsv" into mash_species_ecoli, mash_species_GAS, mash_species_salmonella, mash_species_report

  script:
  """
#!/usr/bin/env python3
import re
import glob

mash_list = glob.glob("*_top_hits.tab")

mash_species = {}

for file in mash_list:
  with open(file, 'r') as f:
    for line in f:
      top_hit = line

      # capture only the genus and species of top hit
      top_hit = re.sub(r'.*-\\.-', '', top_hit)
      top_hit=top_hit.split()
      top_hit=top_hit[0]
      top_hit=re.match(r'^[^_]*_[^_]*', top_hit).group(0)
      top_hit=re.sub(r'.fna', '', top_hit)

      # Ensure top hit has a definitive species assignment
      if "_sp." not in top_hit:
        break
  # specify the top hit as the species for this id
  id = file.split('_top_hit')[0]
  mash_species[id] = top_hit

with open("mash_species.tsv", 'w') as f:
  f.write("Isolate,Predicted Species\\n")
  for key in mash_species.keys():
    f.write("%s,%s\\n"%(key,mash_species[key]))

  """

}

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"

  input:
  set val(name), file(reads) from raw_reads_trim

  output:
  tuple name, file("${name}_trimmed{_1,_2}.fastq.gz") into trimmed_reads

  script:
  """
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${task.cpus} ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} > ${name}.trim.stats.txt
  mv ${name}*1P.fastq.gz ${name}_trimmed_1.fastq.gz
  mv ${name}*2P.fastq.gz ${name}_trimmed_2.fastq.gz
  """
}
//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads
  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}

//Assemble cleaned reads with Shovill
process shovill {
  tag "$name"
  publishDir "${params.outdir}/logs/shovill", mode: 'copy'

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_serotypefinder

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory}  --outdir . --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv contigs.fa ${name}.contigs.fa
  """
}

//Assembly Quality Report
process quast {
  publishDir "${params.outdir}/logs/quast/",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file "${name}_report.tsv" into quast_results, quast_results_report

  script:
  """
  quast.py ${assembly} -o .
  mv report.tsv ${name}_report.tsv
  """


}

//QC of read data
process cg_pipeline {
  publishDir "${params.outdir}/logs/cg_pipeline",mode:'copy'

  input:
  set val(name), file(reads) from raw_reads_qc
  file(quast_report) from quast_results

  output:
  file "${name}_readMeterics.tsv" into cg_pipeline_results

  script:

  """
#!/usr/bin/env python3
import os
import csv
import glob

quast_report = "${quast_report}"
name = "${name}"
subsample = "${params.subsample}"
reads  = "${reads}"
genome_length = ""

# Set genome length from quast output
with open(quast_report) as tsv:
  tsv_reader = csv.reader(tsv, delimiter="\t")
  for line in tsv_reader:
    if "Total length" == line[0]:
      genome_length=line[1]
  if not genome_length:
    raise ValueError("Unable to predict genome length for isolate".format({name}))

#Run CG Pipeline
os.system("run_assembly_readMetrics.pl {} {} -e {} > {}_readMeterics.tsv".format(subsample,reads,genome_length,name))
"""
}


//GAS serotyping
process emmtype_finder {
  publishDir "${params.outdir}/logs/emmtyper/",mode:'copy'

  input:
  set val(name), file(reads) from raw_reads_gas
  file(mash_species) from mash_species_GAS

  output:
  file "${name}*.results.xml" into emmtyper_results optional true

  script:
  """
#!/usr/bin/env python
import os
import csv
import glob

reads =  glob.glob("*fastq*")
mash_species = "${mash_species}"
name = "${name}"
db = "${params.emmtyper_db}"
# Run Emmtyper if isolate predicted as Streptococcus_pyogenes
with open(mash_species) as tsv:
  tsv_reader = csv.reader(tsv, delimiter=",")
  for line in tsv_reader:
    if line[0] == name and line[1] == "Streptococcus_pyogenes":
        os.system("emm_typing.py -1 {} -2 {} -o . -m {}".format(reads[0],reads[1],db))
"""
}

//Salmonella serotyping
process seqsero {
  publishDir "${params.outdir}/logs/seqsero/",mode:'copy'

  input:
  set val(name), file(reads) from raw_reads_salmonella
  file(mash_species) from mash_species_salmonella

  output:
  file "${name}_seqsero.txt" into seqsero_results optional true

  script:
  """
#!/usr/bin/env python
import os
import csv
import glob

mash_species = "${mash_species}"
name = "${name}"
reads = " ".join(glob.glob("*.fastq.gz"))
db = "${params.emmtyper_db}"
# Run SeqSero if isolate predicted as Salmonella enterica
with open(mash_species) as tsv:
  tsv_reader = csv.reader(tsv, delimiter=",")
  for line in tsv_reader:
    if line[0] == name and line[1] == "Salmonella_enterica":
      os.system("SeqSero.py -m2 -i {} -d ./{}".format(reads, name))
      os.rename("./{}/Seqsero_result.txt".format(name), "./{}_seqsero.txt".format(name))

"""
}

//Ecoli serotyping
process serotypefinder {
  publishDir "${params.outdir}/logs/serotypefinder/",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_serotypefinder
  file(mash_species) from mash_species_salmonella

  output:
  file "${name}_serotypefinder.txt" into serotypefinder_results optional true

  script:
  """
#!/usr/bin/env python
import os
import csv

mash_species = "${mash_species}"
name = "${name}"
assembly = "${assembly}"
db = "${params.serotypefinder_db}"
agreement = "${params.serotypefinder_agreement}"
coverage = "${params.serotypefinder_coverage}"

# Run SerotypeFinder if isolate predicted as E.coli
with open(mash_species) as tsv:
  tsv_reader = csv.reader(tsv, delimiter=",")
  for line in tsv_reader:
    print(line[0])
    print(name)
    if line[0] == name and line[1] == "Escherichia_coli":
      print("serotypefinder.pl -d {} -i {} -b /blast-2.2.26/ -o ./{} -s ecoli -k {} -l {}".format(db, assembly, name, agreement, coverage))
      os.system("serotypefinder.pl -d {} -i {} -b /blast-2.2.26/ -o ./{} -s ecoli -k {} -l {}".format(db, assembly, name, agreement, coverage))
      os.rename("./{}/results_table.txt".format(name), "{}_serotypefinder.txt".format(name))

"""
}

//Collect and format all output

// First set falg files for optional-output processes
STF_EMPTY = file("${params.outdir}/logs/Tredegar_trace.txt")
SS_EMPTY = file("${params.outdir}/logs/execution_report.html")
ET_EMPTY = file("${params.outdir}/shovill/*contigs.fa")
process results{
  publishDir "${params.outdir}", mode: 'copy'
  echo true


  input:
  file(cg_pipeline_results) from cg_pipeline_results.collect()
  file(quast_report) from quast_results_report.collect()
  file(mash_species) from mash_species_report
  file(seortypefinder_result) from serotypefinder_results.collect().ifEmpty(STF_EMPTY)
  file(seqsero_results) from seqsero_results.collect().ifEmpty(SS_EMPTY)
  file(emmtyper_results) from emmtyper_results.collect().ifEmpty(ET_EMPTY)

  output:
  file "Tredegar_results.tsv"

  script:
  """
#!/usr/bin/env python3
import os, sys
import glob, csv
import xml.etree.ElementTree as ET
class result_values:
    def __init__(self,id):
        self.id = id
        self.r1_q = "NA"
        self.r2_q = "NA"
        self.est_genome_length = "NA"
        self.est_cvg = "NA"
        self.number_contigs = "NA"
        self.species_prediction = "NA"
        self.subspecies_prediction = "NA"

#get list of result files
cg_results = glob.glob("*_readMeterics.tsv")
quast_results = glob.glob("*_report.tsv")
mash_species = "mash_species.tsv"
serotype_finder_results = glob.glob("*_serotypefinder.txt")
seqsero_results = glob.glob("*_seqsero.txt")
emmtyper_results = glob.glob("*_R1.results.xml")

results = {}

# collect cg_pipeline results
for file in cg_results:
    id = file.split("_readMeterics.tsv")[0]
    result = result_values(id)
    with open(file,'r') as tsv_file:
        tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
        for line in tsv_reader:
            if "_R1" in line["File"]:
                result.r1_q = line["avgQuality"]
                result.est_cvg = float(line["coverage"])
            if "_R2" in line["File"]:
                result.r2_q = line["avgQuality"]
                result.est_cvg += float(line["coverage"])

    # collect quast results
    file = "{}_report.tsv".format(id)
    if not os.path.isfile(file):
        result.est_genome_length = "ASSEMBLY_FAILED"
        result.number_contigs = "ASSEMBLY_FAILED"
    else:
        with open(file, 'r') as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for line in tsv_reader:
                if "Total length" in line[0]:
                    result.est_genome_length = line[1]
                if "# contigs" in line[0]:
                    result.number_contigs = line[1]

    # collect mash_species result
    file = "mash_species.tsv"
    with open(file, 'r') as tsv_file:
      tsv_reader = csv.reader(tsv_file, delimiter=",")
      for line in tsv_reader:
        if line[0] == id:
            result.species_prediction = line[1]

    # collect emmtyper results
    if not glob.glob("{}*.results.xml".format(id)):
        pass
    else:
        file = glob.glob("{}*.results.xml".format(id))[0]
        tree=ET.parse(file)
        root = tree.getroot()
        for emm in root[1].findall("result"):
            if emm.attrib['type'] == 'Final_EMM_type':
                emm_type=(emm.attrib['value'])
                result.subspecies_prediction = emm_type.split(".")[0]

    # collect serotypefinder results
    file = "{}_serotypefinder.txt".format(id)
    if not os.path.isfile(file):
        pass
    else:
    # process the results of serotypefinder as per literature guidelines (Joensen, et al. 2015, DOI: 10.1128/JCM.00008-15)
        with open(file) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            wzx_allele = ""
            wzy_allele = ""
            wzm_allele = ""
            h_type = ""
            o_type = ""
            matched_wzx = ["O2", "O50", "O17", "O77", "O118", "O151", "O169", "O141ab", "O141ac"]
            matched_wzy = ["O13", "O135", "O17", "O44", "O123", "O186"]

            for line in tsv_reader:
                if len(line) == 0:
                    pass
                else:
                    if "fl" in line[0]:
                        h_type = line[5]

                    if line[0] == "wzx":
                        wzx_allele = line[5]
                    if line[0] == "wzy":
                        wzy_allele = line[5]
                    if line[0] == "wzm":
                        wzm_allele = line[5]

                    o_type = wzx_allele
                    if not wzx_allele:
                        o_type = wzy_allele
                    if not wzx_allele and not wzy_allele:
                        o_type = wzm_allele

                    if o_type in matched_wzx:
                        o_type = wzy_allele
                    if o_type in matched_wzy:
                        o_type = wzx_allele
                    serotype = "{}:{}".format(o_type,h_type)

                    # NA if no o-type or h-type identified
                    if serotype == ":":
                        serotype = "NA"
                    result.subspecies_prediction = serotype

    # collect serotypefinder results
    file = "{}_seqsero.txt".format(id)
    if not os.path.isfile(file):
        pass
    else:
        with open(file) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter="\t")
            for line in tsv_reader:
                try:
                    if "Predicted serotype" in line[0]:
                        result.subspecies_prediction = line[1]
                except:
                    pass

    results[id] = result

#create output file
with open("Tredegar_results.tsv",'w') as csvout:
    writer = csv.writer(csvout,delimiter='\t')
    writer.writerow(["sample","rq_1", "r2_q", "est_genome_length", "est_cvg", "number_contigs", "species_prediction", "subspecies_prediction"])
    for id in results:
        result = results[id]
        writer.writerow([result.id,result.r1_q,result.r2_q,result.est_genome_length,result.est_cvg,result.number_contigs,result.species_prediction,result.subspecies_prediction])

"""

}
