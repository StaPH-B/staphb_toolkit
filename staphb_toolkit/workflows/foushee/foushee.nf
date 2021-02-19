#!/usr/bin/env nextflow

//Description: Workflow for quality control of raw illumina reads
//Author: Kevin Libuit
//eMail: kevin.libuit@dgs.virginia.gov

//starting parameters
params.reads = ""
params.outdir = ""
params.report = ""


//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs(  "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

Channel
    .fromPath(params.report)
    .set { report }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(reads) into raw_reads_trim, raw_reads_mash

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
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads, cleaned_reads_mash, cleaned_reads_gas
  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
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
  id = file.split('_top_hits')[0]
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

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}_contigs.fa") into assembled_genomes_quality, assembled_genomes_ksnp

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory}  --outdir . --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv contigs.fa ${name}_contigs.fa
  """
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

reads =  glob.glob("*fastq*")
mash_species = "${mash_species}"
print(reads)
name = "${name}"
db = "${params.emmtyper_db}"
# Run Emmtyper if isolate predicted as Streptococcus_pyogenes
with open(mash_species) as tsv:
  tsv_reader = csv.reader(tsv, delimiter=",")
  for line in tsv_reader:
    if line[0] == name and line[1] == "Streptococcus_pyogenes":
        os.system("emm_typing.py -1 {} -2 {} -o . -m {}".format(reads[0],reads[1],db))
    else:
        print(name + " not identified as Streptococcus pyogenes")
"""
}

//Collect and qc metrics
EMPTY = file('empty')
process assembly_results{
  publishDir "${params.outdir}", mode: 'copy'


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
    writer.writerow(["sample", "est_genome_length", "number_contigs", "species_prediction", "subspecies_prediction"])
    for id in results:
        result = results[id]
        writer.writerow([result.id,result.est_genome_length,result.number_contigs,result.species_prediction,result.subspecies_prediction])

"""
}

process ksnp3 {
  publishDir "${params.outdir}/cluster_analysis/snps", mode: 'copy'

  input:
  file(qc_metrics) from qc_metrics
  file(assembly) from assembled_genomes_ksnp.collect()

  output:
  file("*core_SNPs_matrix.fasta") into ksnp_matrix
  file("*_tree.core.tre") into ksnp_tree

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

  input:
  file(alignment) from ksnp_matrix.collect()

  output:
  file "*pairwise_snp_distance_matrix.tsv" into matrix

  shell:
  """
  for i in *_core_SNPs_matrix.fasta
  do
      group=\$(echo \${i} | sed 's/_core_SNPs_matrix\\.fasta//')
      snp-dists \${group}_core_SNPs_matrix.fasta > \${group}_pairwise_snp_distance_matrix.tsv
  done
  """
}

process render{
  publishDir "${params.outdir}/cluster_analysis/", mode: 'copy', pattern: "*.pdf"
  publishDir "${params.outdir}/cluster_analysis/images", mode: 'copy', pattern: "*.png"
  publishDir "${params.outdir}/cluster_analysis/snps/", mode: 'copy', pattern: "*ordered_snp_distance_matrix.tsv", overwrite: false

  input:
  file(pairwise_snp_distance_matrix) from matrix.collect()
  file(tree_file) from ksnp_tree.collect()
  file(rmd) from report

  output:
  file "*foushee_cluster_report.pdf"
  file "*_core_parsimony_tree.png"
  file "*SNP_heatmap.png"
  file "*snp_distance_matrix.tsv"
  shell:
"""
for i in *pairwise_snp_distance_matrix.tsv
do
  emm_type=\$(echo \$i | cut -d _ -f 1)
  if [ -s \${emm_type}_tree.core.tre ]
  then
    cp ${rmd} ./\${emm_type}_report_template.Rmd
    Rscript /reports/render.R \${emm_type}_pairwise_snp_distance_matrix.tsv \${emm_type}_tree.core.tre ./\${emm_type}_report_template.Rmd
    mv report.pdf \${emm_type}_foushee_cluster_report.pdf
    mv core_parsimony_tree.png \${emm_type}_core_parsimony_tree.png
    mv SNP_heatmap.png \${emm_type}_SNP_heatmap.png
    mv snp_distance_matrix.tsv \${emm_type}_ordered_snp_distance_matrix.tsv
  else
  # also publish matrices for emmtypes with less than 3 isolates & remove empty tree file
    mv \${emm_type}_pairwise_snp_distance_matrix.tsv \${emm_type}_ordered_snp_distance_matrix.tsv
    rm \${emm_type}_tree.core.tre
  fi
done
"""
}
