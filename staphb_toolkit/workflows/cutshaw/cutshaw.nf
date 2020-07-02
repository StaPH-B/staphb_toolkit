#!/usr/bin/env nextflow

//Description: Workflow for quality control of raw illumina reads
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.reads = ""
params.isolate_key = ""
params.pt_genomes =  ""
params.outdir = ""

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs(  "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

Channel
    .fromPath(params.isolate_key)
    .set { isolate_key }

Channel
    .fromPath(params.pt_genomes)
    .set { pt_genomes }

Channel
    .fromPath(params.pt_genomes)
    .set { pt_genomes_snp }

// process test {
//   input:
//   file(isolate_key) from isolate_key
//   file(genomes) from pt_genomes.collect()
//
//   script:
//   """
//   echo "${isolate_key}"
//   ls
//   asdfasdfsdf
//   """
// }


//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(reads) into raw_reads_assemble, raw_reads_snp, raw_reads_qc

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

//Assemble reads with Spades
process spades {
  tag "$name"
  publishDir "${params.outdir}/shovill", mode: 'copy'

  input:
  set val(name), file(reads) from raw_reads_assemble

  output:
  tuple name, file("${name}_contigs.fasta") into assembled_genomes_quality, assembled_genomes_fastani

  shell:
  '''
  ram=`awk '/MemTotal/ { printf "%.0f \\n", $2/1024/1024 - 1 }' /proc/meminfo`
  spades.py --memory $ram -1 !{reads[0]} -2 !{reads[1]} -o ./spades_out
  mv ./spades_out/contigs.fasta !{name}_contigs.fasta
  '''
}

//Assembly Quality Report
process quast {

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


// Run FastANI
process fastani{

  input:
  set val(name), file(assembly) from assembled_genomes_fastani
  file(genomes) from pt_genomes.collect()

  output:
  file "fastani_*.out" into fastani_results

  shell:
  """
  fastANI -q ${assembly} --rl ./reference_list.txt -o ./fastani_${name}.out
  """
}

process cfsan_snp{

  input:
  set val(name), file(reads) from raw_reads_snp
  file(isolate_key) from isolate_key.collect()
  file(genomes) from pt_genomes_snp.collect()

  output:
  file("*metrics.tsv") into cfsan_metrics

  script:
  """
#!/usr/bin/env python
import subprocess
import shutil
import os
import glob
import csv

isolate_key = "${isolate_key}"
reads = glob.glob("*.gz")
sample = "${name}"

with open(isolate_key,'r') as csv_file:
    csv_reader = list(csv.DictReader(csv_file, delimiter=","))
    for line in csv_reader:
        if line["sample_id"] in sample:
            reference_assembly = line["isolate_id"] + ".fasta"
            sample = line["sample_id"]

    if not reference_assembly:
        print("no PT reference found; please ensure isolate_key file is properly formatted")
        exit()

os.makedirs(os.path.join("input_reads", sample))
for file in reads:
    shutil.move(file, "input_reads/%s/"%sample)

os.system("run_snp_pipeline.sh -m soft -o . -s ./input_reads {}".format(reference_assembly))
os.rename("metrics.tsv", "%s_metrics.tsv"%sample)
"""
}

// //Collect and qc metrics
// EMPTY = file('empty')
// process assembly_results{
//   publishDir "${params.outdir}", mode: 'copy'
//
//
//   input:
//   file(quast_report) from quast_results_report.collect()
//   file(mash_species) from mash_species_report
//   file(emmtyper_results) from emmtyper_results.collect().ifEmpty(EMPTY)
//
//   output:
//   file "assembly_metrics.csv" into qc_metrics
//
//   script:
//   """
// #!/usr/bin/env python3
// import os, sys
// import glob, csv
// import xml.etree.ElementTree as ET
// class result_values:
//     def __init__(self,id):
//         self.id = id
//         self.est_genome_length = "NA"
//         self.number_contigs = "NA"
//         self.species_prediction = "NA"
//         self.subspecies_prediction = "NA"
//
//
// #get list of result files
// quast_results = glob.glob("*_report.tsv")
// mash_species = "mash_species.tsv"
// emmtyper_results = glob.glob("*.results.xml")
//
//
// results = {}
//
// # collect cg_pipeline results
// for file in quast_results:
//     id = file.split("_report.tsv")[0]
//     result = result_values(id)
//     if not os.path.isfile(file):
//         result.est_genome_length = "ASSEMBLY_FAILED"
//         result.number_contigs = "ASSEMBLY_FAILED"
//     else:
//         with open(file, 'r') as tsv_file:
//             tsv_reader = csv.reader(tsv_file, delimiter="\t")
//             for line in tsv_reader:
//                 if "Total length" in line[0]:
//                     result.est_genome_length = line[1]
//                 if "# contigs" in line[0]:
//                     result.number_contigs = line[1]
//spades
//     # collect mash_species result
//     file = "mash_species.tsv"
//     with open(file, 'r') as tsv_file:
//       tsv_reader = csv.reader(tsv_file, delimiter=",")
//       for line in tsv_reader:
//         if line[0] == id:
//             result.species_prediction = line[1]
//
//     # collect emmtyper results
//     file = glob.glob("{}*.results.xml".format(id))[0]
//     if not os.path.isfile(file):
//         pass
//     else:
//         tree=ET.parse(file)
//         root = tree.getroot()
//         for emm in root[1].findall("result"):
//             if emm.attrib['type'] == 'Final_EMM_type':
//                 emm_type=(emm.attrib['value'])
//                 result.subspecies_prediction = emm_type.split(".")[0]
//
//     results[id] = result
//
// #create output file
// with open("assembly_metrics.csv",'w') as csvout:
//     writer = csv.writer(csvout,delimiter=',')
//     writer.writerow(["sample", "est_genome_length", "number_contigs", "species_prediction", "subspecies_prediction"])
//     for id in results:
//         result = results[id]
//         writer.writerow([result.id,result.est_genome_length,result.number_contigs,result.species_prediction,result.subspecies_prediction])
//
// """
// }
//
//
// process render{
//   publishDir "${params.outdir}/cluster_analysis/", mode: 'copy', pattern: "*.pdf"
//   publishDir "${params.outdir}/cluster_analysis/images", mode: 'copy', pattern: "*.png"
//   publishDir "${params.outdir}/cluster_analysis/snps/", mode: 'copy', pattern: "*ordered_snp_distance_matrix.tsv", overwrite: false
//
//   input:
//   file(pairwise_snp_distance_matrix) from matrix.collect()
//   file(tree_file) from ksnp_tree.collect()
//   file(rmd) from report
//
//   output:
//   file "*foushee_cluster_report.pdf"
//   file "*_core_parsimony_tree.png"
//   file "*SNP_heatmap.png"
//   file "*snp_distance_matrix.tsv"
//   shell:
// """
// for i in *pairwise_snp_distance_matrix.tsv
// do
//   emm_type=\$(echo \$i | cut -d _ -f 1)
//   if [ -s \${emm_type}_tree.core.tre ]
//   then
//     cp ${rmd} ./\${emm_type}_report_template.Rmd
//     Rscript /reports/render.R \${emm_type}_pairwise_snp_distance_matrix.tsv \${emm_type}_tree.core.tre ./\${emm_type}_report_template.Rmd
//     mv report.pdf \${emm_type}_foushee_cluster_report.pdf
//     mv core_parsimony_tree.png \${emm_type}_core_parsimony_tree.png
//     mv SNP_heatmap.png \${emm_type}_SNP_heatmap.png
//     mv snp_distance_matrix.tsv \${emm_type}_ordered_snp_distance_matrix.tsv
//   else
//   # also publish matrices for emmtypes with less than 3 isolates & remove empty tree file
//     mv \${emm_type}_pairwise_snp_distance_matrix.tsv \${emm_type}_ordered_snp_distance_matrix.tsv
//     rm \${emm_type}_tree.core.tre
//   fi
// done
// """
// }
