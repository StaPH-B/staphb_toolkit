#!/usr/bin/env nextflow

//Description: Workflow for quality control of raw illumina reads
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs(  "${params.reads}/*{R1,R2,_1,_2}*.fastq.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file("*{R1,R2,_1,_2}.fastq.gz") into raw_reads_clean, raw_reads_mash

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
//Trim reads and remove adapters
process seqyclean {
  tag "$name"
  publishDir "${params.outdir}/seqyclean/", mode: 'copy'

  input:
  set val(name), file(reads) from raw_reads_clean

  output:
  tuple val("${name}"), "${name}_clean_PE{1,2}.fastq" into cleaned_reads, cleaned_reads_mash, cleaned_reads_gas

  script:
  """
  seqyclean -1 ${name}_R1.fastq.gz -2 ${name}_R2.fastq.gz -minlen ${params.min_read_length} -o ${name}_clean -c ${params.contaminants} ${params.quality_trimming}
  """
}

//Create Mash Sketches for all isolates
process mash_dist{
  tag "$name"

  input:
  set val(name), file(reads) from cleaned_reads_mash

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
  publishDir "${params.outdir}/mash/", mode: 'copy'

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
  id = file.split('_')[0]
  mash_species[id] = top_hit

with open("mash_species.tsv", 'w') as f:
  f.write("Isolate,Predicted Species\\n")
  for key in mash_species.keys():
    f.write("%s,%s\\n"%(key,mash_species[key]))

  """

}

//Assemble cleaned reads with Shovill
process shovill {
  tag "$name"
  publishDir "${params.outdir}/shovill", mode: 'copy'

  memory '8 GB'
  ram=6

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}_contigs.fa") into assembled_genomes_quality, assembled_genomes_ksnp

  shell:
  '''
  ram=`awk '/MemTotal/ { printf "%.0f \\n", $2/1024/1024 - 1 }' /proc/meminfo`
  shovill --cpus 0 --ram $ram  --outdir . --R1 !{reads[0]} --R2 !{reads[1]} --force
  mv contigs.fa !{name}_contigs.fa
  '''
}

//Assembly Quality Report
process quast {
  publishDir "${params.outdir}/quast/",mode:'copy'

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

//GAS serotyping
process emmtype_finder {
  publishDir "${params.outdir}/emmtyper/",mode:'copy'

  input:
  set val(name), file(reads) from cleaned_reads_gas
  file(mash_species) from mash_species_GAS

  output:
  file "${name}*.results.xml" into emmtyper_results, emmtyper_results_ksnp

  script:
  """
#!/usr/bin/env python
import os
import csv
import glob

reads =  glob.glob("*fastq")
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

//Collect and qc metrics
EMPTY = file('empty')
process results{
  publishDir "${params.outdir}", mode: 'copy'
  echo true


  input:
  file(quast_report) from quast_results_report.collect()
  file(mash_species) from mash_species_report
  file(emmtyper_results) from emmtyper_results.collect().ifEmpty(EMPTY)

  output:
  file "assembly_metrics.csv" into qc_metrics

  script:
  """
#!/usr/bin/env python3
import os, sys
import glob, csv
import xml.etree.ElementTree as ET
class result_values:
    def __init__(self,id):
        self.id = id
        self.est_genome_length = "NA"
        self.est_cvg = "NA"
        self.number_contigs = "NA"
        self.species_prediction = "NA"
        self.subspecies_prediction = "NA"


#get list of result files
quast_results = glob.glob("*_report.tsv")
mash_species = "mash_species.tsv"
emmtyper_results = glob.glob("*.results.xml")


results = {}

# collect cg_pipeline results
for file in quast_results:
    id = file.split("_report.tsv")[0]
    result = result_values(id)
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
    file = glob.glob("{}*.results.xml".format(id))[0]
    if not os.path.isfile(file):
        pass
    else:
        tree=ET.parse(file)
        root = tree.getroot()
        for emm in root[1].findall("result"):
            if emm.attrib['type'] == 'Final_EMM_type':
                emm_type=(emm.attrib['value'])
                result.subspecies_prediction = emm_type.split(".")[0]

    results[id] = result

#create output file
with open("assembly_metrics.csv",'w') as csvout:
    writer = csv.writer(csvout,delimiter=',')
    writer.writerow(["sample", "est_genome_length", "est_cvg", "number_contigs", "species_prediction", "subspecies_prediction"])
    for id in results:
        result = results[id]
        writer.writerow([result.id,result.est_genome_length,result.est_cvg,result.number_contigs,result.species_prediction,result.subspecies_prediction])

"""
}

process ksnp3 {
  publishDir "${params.outdir}/ksnp3", mode: 'copy'

  input:
  file(qc_metrics) from qc_metrics
  file(assembly) from assembled_genomes_ksnp.collect()

  output:
  file("*core_SNPs_matrix.fasta") into ksnp_matrix
  file("*_tree.core.tre")

  script:
  """
#!/usr/bin/env python
import os
import glob, csv
import re
import xml.etree.ElementTree as ET
class result_values:
    def __init__(self,id):
        self.subspecies_prediction = "NA"

qc_metrics = "assembly_metrics.csv"
assemblies = glob.glob("*_contigs.fa")
kmer_length = "${params.kmer_length}"
kmer_snps = "${params.kmer_snps}"

with open(qc_metrics,'r') as csv_file:
    csv_reader = list(csv.DictReader(csv_file, delimiter=","))
    for line in csv_reader:
        if re.search(r"emm\\d", line["subspecies_prediction"]) is not None:
            emm_type = line["subspecies_prediction"]
            isolate = line["sample"]
            assembly = os.path.abspath("{}_contigs.fa".format(isolate))

            with open("{}_assemblies.txt".format(emm_type), 'a') as file:
                file.write("{}\\t{}\\n".format(assembly, isolate))
        else:
            print(str(re.search(r"emm/d", line["subspecies_prediction"])))
            print(line["subspecies_prediction"])

assembly_paths = glob.glob("*_assemblies.txt")

for file in assembly_paths:
    if len(open(file).readlines(  )) < 2:
        pass
    else:
        group = file.split("_assemblies.txt")[0]
        file = os.path.abspath(file)
        if not os.path.exists("ksnp3"):
            os.makedirs("ksnp3")
        os.system("kSNP3 -in {}_assemblies.txt -outdir ksnp3/{} -k {} {}".format(group,group,kmer_length,kmer_snps))
        os.rename("ksnp3/{}/core_SNPs_matrix.fasta".format(group), "{}_core_SNPs_matrix.fasta".format(group))
        os.rename("ksnp3/{}/tree.core.tre".format(group), "{}_tree.core.tre".format(group))
  """
}

// Generate SNP matrix from ksnp3 matrix
process snp_dists{
  publishDir "${params.outdir}", mode: 'copy'
  echo true

  input:
  file(alignment) from ksnp_matrix.collect()

  output:
  file "*pairwise_snp_distance_matrix.tsv"

  shell:
  """
  for i in *_core_SNPs_matrix.fasta
  do
      group=\$(echo \${i} | sed 's/_core_SNPs_matrix\\.fasta//')
      snp-dists \${group}_core_SNPs_matrix.fasta > \${group}_pairwise_snp_distance_matrix.tsv
  done
  """
}