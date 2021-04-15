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
    .fromPath(params.isolate_key)
    .set { isolate_key_quast }

Channel
    .fromPath(params.isolate_key)
    .set { isolate_key_qc }

Channel
    .fromPath(params.pt_genomes)
    .set { pt_genomes }

Channel
    .fromPath(params.pt_genomes)
    .set { pt_genomes_quast }

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
  publishDir "${params.outdir}/logs/spades", mode: 'copy'

  input:
  set val(name), file(reads) from raw_reads_assemble

  output:
  tuple name, file("${name}_contigs.fasta") into assembled_genomes_quality, assembled_genomes_fastani

  script:
  """
  spades.py --memory ${task.memory.toGiga()} -1 ${reads[0]} -2 ${reads[1]} -o ./spades_out
  mv ./spades_out/contigs.fasta ${name}_contigs.fasta
  """
}

//Assembly Quality Report
process quast {

  input:
  set val(name), file(assembly) from assembled_genomes_quality
  file(isolate_key) from isolate_key_quast.collect()
  file(genomes) from pt_genomes_quast.collect()


  output:
  file "*quast_report.tsv" into quast_results, quast_results_report

  script:
  """
#!/usr/bin/env python
import subprocess
import shutil
import os
import glob
import csv

isolate_key = "${isolate_key}"
assembly = "${assembly}"
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

os.system("quast.py {} -r {} -o .".format(assembly, reference_assembly))
os.rename("report.tsv", "%s_quast_report.tsv"%sample)
"""
}

//QC of read data
process cg_pipeline{

  input:
  set val(name), file(reads) from raw_reads_qc
  file(quast_report) from quast_results

  output:
  file "${name}_readMetrics.tsv" into cg_pipeline_results

  script:

  """
#!/usr/bin/env python3
import os
import csv
import glob

quast_report = "${quast_report}"
name = "${name}"
reads  = "${reads}"
subsample = "${params.subsample}"
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
os.system("run_assembly_readMetrics.pl {} {} -e {} > {}_readMetrics.tsv".format(subsample,reads,genome_length,name))
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

//Collect and qc metrics
process qc_results{
  publishDir "${params.outdir}/logs", mode: 'copy'

  input:
  file(quast_report) from quast_results_report.collect()
  file(cg_results) from cg_pipeline_results.collect()
  file(fastani_results) from fastani_results.collect()
  file(cfsan_results) from cfsan_metrics.collect()
  file(isolate_key) from isolate_key_qc.collect()

  output:
  file "seq_results.tsv" into qc_metrics

  script:
  """
#!/usr/bin/env python3
import os, sys
import glob, csv
import xml.etree.ElementTree as ET
class result_values:
    def __init__(self,id):
        self.SampleID = id
        self.IsolateID = "NA"
        self.Organism = "NA"
        self.Genus = "NA"
        self.Sequencer = "Illumina MiSeq"
        self.Machine = "NA"
        self.FlowCell = "NA"
        self.LibKit = "Netera XT"
        self.Chemistry = "NA"
        self.RunDate = "NA"
        self.SequencedBy = "NA"
        self.SamplesPerRun = 16
        self.SeqLength = "NA"
        self.Reads = "NA"
        self.MeanR1Qual = "NA"
        self.MeanR2Qual = "NA"
        self.PercMapped = "NA"
        self.MeanDepth = "NA"
        self.CovLT10 = 0
        self.SNPs = "NA"
        self.MeanInsert = "NA"
        self.NG50 = "NA"
        self.GenomeFraction = "NA"
        self.Contigs = "NA"
        self.LengthDelta = "NA"
        self.UnalignedLength = "NA"
        self.MostAbundantOrganism = "NA"
        self.Misannotated = "NA"
        self.Coverage = "NA"

#get list of result files
quast_results = glob.glob("*_report.tsv")
cg_results =  glob.glob("*_readMetrics.tsv")
fastani_results = glob.glob("fastani_*.out")
cfsan_results = glob.glob("*_metrics.tsv")
isolate_key = "${isolate_key}"


results = {}

# collect cfsan results
for file in cfsan_results:
    id = file.split("_metrics.tsv")[0]
    result = result_values(id)
    with open(file, 'r') as tsv_file:
        tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

        result.MeanDepth = tsv_reader[0]["Average_Pileup_Depth"]
        result.MeanInsert = int(round(float(tsv_reader[0]["Average_Insert_Size"])))
        result.PercMapped = tsv_reader[0]["Percent_of_Reads_Mapped"]
        result.SNPs = tsv_reader[0]["Phase2_Preserved_SNPs"]
        result.FlowCell = tsv_reader[0]["Flowcell"]
        result.Machine = tsv_reader[0]["Machine"]

    # collect cg_pipeline results
    for cg_result in cg_results:
        if id in cg_result:
            file = cg_result
            with open(file, 'r') as tsv_file:
                tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))

                for line in tsv_reader:
                    if any(fwd_format in line["File"] for fwd_format in ["_1.fastq", "_R1.fastq", "_R1_001.fastq.gz"]):
                        result.MeanR1Qual = line["avgQuality"]
                        result.Reads = int(line["numReads"])
                        result.Coverage = float(line["coverage"])
                    if any(rev_format in line["File"] for rev_format in ["_2.fastq", "_R2.fastq", "_R2_001.fastq.gz"]):
                        result.MeanR2Qual = line["avgQuality"]
                        result.r2_totalBases = line["totalBases"]
                        result.Reads += int(line["numReads"])
                        result.Coverage += float(line["coverage"])

    # collect cg_pipeline results
    for quast_result in quast_results:
        if id in quast_result:
            file = quast_result

            with open(file, 'r') as tsv_file:
                tsv_reader = csv.reader(tsv_file, delimiter="\t")
                for line in tsv_reader:
                    if "# contigs" in line[0]:
                        result.Contigs = line[1]
                    if "Genome fraction (%)" in line[0]:
                        result.GenomeFraction = line[1]
                    if "NG50" in line[0]:
                        result.NG50 = line[1]
                    if "Unaligned length" in line[0]:
                        result.UnalignedLength = line[1]
                    if "Total length" in line[0]:
                        total_length = line[1]
                    if "Reference length" in line[0]:
                        reference_length = line[1]
                result.LengthDelta = int(abs(int(total_length) - int(reference_length)))

    # collect fastani results
    for fastani_result in fastani_results:
        if id in fastani_result:
            file = fastani_result

            with open(file, 'r') as file:
                tsv_reader = csv.reader(file, delimiter="\t", quotechar='"')
                result.Organism = ""
                predicted_reference = str(os.path.basename(next(tsv_reader)[1]))
                if "SAP18-0432" in predicted_reference:
                    result.Organism = "Salmonella enterica subsp. enterica serover Enteritidis"

                elif "SAP18-H9654" in predicted_reference:
                    result.Organism = "Salmonella enterica subsp. enterica serover Enteritidis"
                elif "SAP18-6199" in predicted_reference:
                    result.Organism = "Salmonella enterica subsp. enterica serover Typhimurium"
                elif "SAP18-8729" in predicted_reference:
                    result.Organism = "Salmonella enterica subsp. enterica serover Newport"
                elif "LMP18-H2446" in predicted_reference:
                    result.Organism = "Listeria monocytogenes"
                elif "LMP18-H8393" in predicted_reference:
                    result.Organism = "Listeria monocytogenes"
                else:
                    raise ValueError("Sample %s not identified as a 2018 PT isolate"%id)
                result.Genus = result.Organism.split()[0]
            # set indicated isolate ID
            with open(isolate_key,'r') as csv_file:
                csv_reader = list(csv.DictReader(csv_file, delimiter=","))
                for line in csv_reader:
                    if line["sample_id"] in result.SampleID:
                        result.IsolateID = line["isolate_id"]

            if result.IsolateID in predicted_reference:
                result.Misannotated = "FALSE"
            else:
                print("Isolate: " + result.IsolateID + " Predicted: " + predicted_reference)
                result.Misannotated = "TRUE"

    results[id] = result

#create output file
with open("seq_results.tsv",'w') as csvout:
    writer = csv.writer(csvout,delimiter='\t')
    writer.writerow(["SampleID", "IsolateID", "Organism", "Genus","Sequencer","Machine", "FlowCell", "LibKit", "Chemistry", "RunDate", "SequencedBy", "SamplesPerRun", "Reads", "SeqLength", "MeanR1Qual", "MeanR2Qual", "PercMapped", "MeanDepth", "CovLT10", "SNPs", "MeanInsert", "NG50", "GenomeFraction", "Contigs", "LengthDelta","UnalignedLength", "MostAbundantOrganism", "Misannotated", "Coverage"])
    for id in results:
        result = results[id]
        writer.writerow([result.SampleID,result.IsolateID,result.Organism,result.Genus,result.Sequencer,result.Machine,result.FlowCell,result.LibKit,result.Chemistry,result.RunDate,result.SequencedBy,result.SamplesPerRun,result.Reads,result.SeqLength,result.MeanR1Qual,result.MeanR2Qual,result.PercMapped,result.MeanDepth,result.CovLT10,result.SNPs,result.MeanInsert,result.NG50,result.GenomeFraction,result.Contigs,result.LengthDelta,result.UnalignedLength,result.MostAbundantOrganism,result.Misannotated,result.Coverage])

"""
}


process render{
  publishDir "${params.outdir}/", mode: 'copy'


  input:
  file(seq_results) from qc_metrics

  output:
  file "Cutshaw-report.pdf"
  shell:
"""
create_competency_report.sh -p "Cutshaw-report" -t "${params.title}" -T /wgs_competency/report_template.Rnw -o . -s ${seq_results} -S /wgs_competency/seq_stats.tsv
"""
}
