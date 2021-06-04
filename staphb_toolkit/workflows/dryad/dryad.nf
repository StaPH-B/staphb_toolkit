#!/usr/bin/env nextflow

//Description: Workflow for building various trees from raw illumina reads
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
    .set { raw_reads }

if (params.snp) {
    Channel
        .fromPath(params.snp_reference)
        .set { snp_reference }
}

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(outfiles) into read_files_fastqc, read_files_trimming

  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    outfiles = ["${name}_R1.fastq.gz","${name}_R2.fastq.gz"]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
    outfiles = reads
    """
    """
  }
}

//Step1a: FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/logs/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from read_files_fastqc

  output:
  file("*_fastqc.{zip,html}") into fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"
  if(params.savetrimmedreads){
    publishDir "${params.outdir}/trimmed", mode: 'copy'
  }
  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_trimmed{_1,_2}.fastq.gz") into trimmed_reads
  file("${name}.trim.stats.txt") into trimmomatic_stats

  script:
  """
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${task.cpus} ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} 2> ${name}.trim.stats.txt
  mv ${name}*1P.fastq.gz ${name}_trimmed_1.fastq.gz
  mv ${name}*2P.fastq.gz ${name}_trimmed_2.fastq.gz
  """
}
//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"
  publishDir "${params.outdir}/logs/cleanedreads", mode: 'copy',pattern:"*.stats.txt"

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads_cg
  file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads_snp
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats

  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}

if (params.snp) {
    //SNP Step1: Run CFSAN-SNP Pipeline
    process cfsan {
      publishDir "${params.outdir}/results", mode: 'copy'

      input:
      file(reads) from cleaned_reads_snp.collect()
      file(reference) from snp_reference

      output:
      file("snp_distance_matrix.tsv") into snp_mat
      file("snpma.fasta") into snp_alignment

      when:
      params.snp == true

      script:
      """
      #!/usr/bin/env python
      import subprocess
      import glob
      import os

      fwd_reads = glob.glob("*_1.clean.fastq.gz")
      fwd_reads.sort()

      readDict = {}
      for file in fwd_reads:
        sid = os.path.basename(file).split('_')
        sid = sid[0:(len(sid)-1)]
        sid = '_'.join(sid)
        fwd_read = glob.glob(sid+"_1.clean.fastq.gz")[0]
        rev_read = glob.glob(sid+"_2.clean.fastq.gz")[0]
        readDict[sid] = [fwd_read,rev_read]

      os.mkdir("input_reads")
      for key in readDict:
        print key
        os.mkdir(os.path.join("input_reads",key))
        os.rename(readDict[key][0],os.path.join(*["input_reads",key,readDict[key][0]]))
        os.rename(readDict[key][1],os.path.join(*["input_reads",key,readDict[key][1]]))

      command = "cfsan_snp_pipeline run ${reference} -o . -s input_reads"
      process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
      output, error = process.communicate()
      print output
      print error
      """
    }

    //SNP Step2: Run IQTREE on snp alignment
    process snp_tree {
      publishDir "${params.outdir}/results", mode: 'copy'

      input:
      file(snp_fasta) from snp_alignment

      output:
      file("snp.tree") optional true

      script:
        """
        numGenomes=`grep -o '>' snpma.fasta | wc -l`
        if [ \$numGenomes -gt 3 ]
        then
          iqtree -nt AUTO -s snpma.fasta -m ${params.cg_tree_model} -bb 1000
          mv snpma.fasta.contree snp.tree
        fi
        """
    }
}

//CG Step1: Assemble trimmed reads with Shovill and map reads back to assembly
process shovill {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/results/assembled", mode: 'copy',pattern:"*.fa"
  publishDir "${params.outdir}/results/alignments", mode: 'copy',pattern:"*.sam"

  input:
  set val(name), file(reads) from cleaned_reads_cg

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_annotation, assembled_genomes_ar, assembled_genomes_mash, assembled_genomes_mlst
  tuple name, file("${name}.sam") into sam_files

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory} --outdir ./output --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv ./output/contigs.fa ${name}.contigs.fa
  bwa index ${name}.contigs.fa
  bwa mem ${name}.contigs.fa ${reads[0]} ${reads[1]} > ${name}.sam
  """
}

//Index and sort bam file then calculate coverage
process samtools {
  tag "$name"

  publishDir "${params.outdir}/results/alignments", mode: 'copy', pattern:"*.bam"
  publishDir "${params.outdir}/results/coverage", mode: 'copy', pattern:"*_depth.tsv*"

  input:
  set val(name), file(sam) from sam_files

  output:
  file("${name}_depth.tsv") into cov_files

  shell:
  """
  samtools view -S -b ${name}.sam > ${name}.bam
  samtools sort ${name}.bam > ${name}.sorted.bam
  samtools index ${name}.sorted.bam
  samtools depth -a ${name}.sorted.bam > ${name}_depth.tsv
  """
}

//Calculate median coverage
process coverage_stats {
  publishDir "${params.outdir}/results/coverage", mode: 'copy'

  input:
  file(cov) from cov_files.collect()

  output:
  file('coverage_stats.txt')

  script:
  '''
  #!/usr/bin/env python3
  import glob
  import os
  from numpy import median
  from numpy import average

  results = []

  files = glob.glob("*_depth.tsv*")
  for file in files:
    nums = []
    sid = os.path.basename(file).split('_')[0]
    with open(file,'r') as inFile:
      for line in inFile:
        nums.append(int(line.strip().split()[2]))
      med = median(nums)
      avg = average(nums)
      results.append(f"{sid}\\t{med}\\t{avg}\\n")

  with open('coverage_stats.txt', 'w') as outFile:
    outFile.write("Sample\\tMedian Coverage\\tAverage Coverage\\n")
    for result in results:
      outFile.write(result)
  '''
}

process mash {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/results/mash",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_mash

  output:
  file("${name}.mash.txt") into mash_result

  script:
  """
  mash dist /db/RefSeqSketchesDefaults.msh ${assembly} > ${name}.txt
  sort -gk3 ${name}.txt | head > ${name}.mash.txt
  """
}

//CG Step2a: Assembly Quality Report
process quast {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/logs/quast",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_report

  script:
  """
  quast.py ${assembly} -o .
  mv report.txt ${name}.quast.tsv
  """
}


//CG Step2b: Annotate with prokka
process prokka {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/results/annotated",mode:'copy'

  input:
  file(mash) from mash_result
  set val(name), file(assembly) from assembled_genomes_annotation


  output:
  file("${name}.gff") into annotated_genomes
  file("${name}.prokka.stats.txt") into prokka_stats

  when:
  params.cg == true

  script:
  """
  #genus=`awk 'FNR == 1 {split(\$1,a,"-\\.-|\\.fna");print a[2]}'  ${name}.mash.txt | awk '{split(\$1,a,"_");print a[1]}'`
  #species=`awk 'FNR == 1 {split(\$1,a,"-\\.-|\\.fna");print a[2]}'  ${name}.mash.txt | awk '{split(\$1,a,"_");print a[2]}'`
  prokka --cpu ${task.cpus} --force --compliant --prefix ${name} --genus "" --species ""  --mincontiglen 500 --outdir . ${assembly} > ${name}.log
  mv ${name}.txt ${name}.prokka.stats.txt
  """
}

//CG Step3: Align with Roary
process roary {
  publishDir "${params.outdir}/results",mode:'copy'

  numGenomes = 0
  input:
  file(genomes) from annotated_genomes.collect()

  output:
  file("core_gene_alignment.aln") into core_aligned_genomes
  file("core_genome_statistics.txt") into core_aligned_stats

  script:
  if(params.roary_mafft == true){
    mafft="-n"
  }else{mafft=""}
  """
  roary -e ${mafft} -p ${task.cpus} ${genomes}
  mv summary_statistics.txt core_genome_statistics.txt
  """
}

//CG Step4: IQTree for core-genome
process cg_tree {
  publishDir "${params.outdir}/results",mode:'copy'

  input:
  file(alignedGenomes) from core_aligned_genomes

  output:
  file("core_genome.tree") optional true into cgtree

  script:
    """
    numGenomes=`grep -o '>' core_gene_alignment.aln | wc -l`
    if [ \$numGenomes -gt 3 ]
    then
      iqtree -nt AUTO -s core_gene_alignment.aln -m ${params.cg_tree_model} -bb 1000
      mv core_gene_alignment.aln.contree core_genome.tree
    fi
    """
}

//AR Step1: Find AR genes with amrfinder+
process amrfinder {
  tag "$name"
  publishDir "${params.outdir}/results/amrfinder",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_ar

  output:
  file("${name}.tsv") into ar_predictions

  when:
  params.ar == true

  script:
  """
  amrfinder -n ${assembly} -o ${name}.tsv
  """
}

//AR Step 2: Summarize amrfinder+ results as a binary presence/absence matrix
process amrfinder_summary {
  tag "$name"
  publishDir "${params.outdir}/results",mode:'copy'

  input:
  file(predictions) from ar_predictions.collect()

  output:
  file("ar_predictions.tsv") into ar_tsv

  when:
  params.ar == true

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd

  files = glob.glob("*.tsv")
  dfs = []
  for file in files:
      sample_id = os.path.basename(file).split(".")[0]
      print(sample_id)
      df = pd.read_csv(file, header=0, delimiter="\\t")
      df.columns=df.columns.str.replace(' ', '_')
      print(df)
      df = df.assign(Sample=sample_id)
      df = df[['Sample','Gene_symbol','%_Coverage_of_reference_sequence','%_Identity_to_reference_sequence']]
      df = df.rename(columns={'%_Coverage_of_reference_sequence':'Coverage','%_Identity_to_reference_sequence':'Identity','Gene_symbol':'Gene'})
      dfs.append(df)

  concat = pd.concat(dfs)
  concat.to_csv('ar_predictions.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

Channel
  .fromPath(params.multiqc_config)
  .set { mqc_config }

Channel
  .fromPath(params.multiqc_logo)
  .set { mqc_logo }

//Collect Results
process multiqc {
  publishDir "${params.outdir}/results", mode: 'copy'
  echo true

  input:
  //file multiqc_config
  file logo from mqc_logo
  file config from mqc_config
  file fastq_results from fastqc_results.collect()
  path trimmomatic from trimmomatic_stats.collect()
  path q_report from quast_report
  file cg_stats from core_aligned_stats
  file phix_removal from phix_cleanning_stats.collect()
  file adapter_removal from adapter_cleanning_stats.collect()
  file prokka_logs from prokka_stats.collect()


  output:
  file("*multiqc_report.html") into multiqc_report
  file("*_data")

  script:
  """
  multiqc . >/dev/null 2>&1
  """
}

process mlst {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/results",mode:'copy'

  input:
  file(assemblies) from assembled_genomes_mlst.collect()

  output:
  file("mlst.tsv") into mlst_results

  script:
  """
  mlst --nopath *.fa > mlst.tsv
  """
}

process mlst_formatting {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/results",mode:'copy'

  input:
  file(mlst) from mlst_results

  output:
  file("mlst_formatted.tsv")

  script:
  """
  #!/usr/bin/env python3
  import csv

  string_map = {}

  with open('mlst.tsv','r') as csvfile:
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
    csvfile.seek(0)
    reader = csv.reader(csvfile,dialect,delimiter='\t')
    for row in reader:
      id_string = row[0]
      sp_string = row[1]
      st_string = row[2]
      string_map[id_string] = [sp_string,st_string]

  mlst = []
  for key in string_map:
    id = key
    id = id.replace('.contigs.fa','')
    species = string_map[key][0]
    st = string_map[key][1]
    if species == 'abaumannii':
        st = 'PubMLST ST' + str(st) + ' (Oxford)'
    if species == 'abaumannii_2':
        st = 'PubMLST ST' + str(st) + ' (Pasteur)'
    else:
        st = 'PubMLST ST' + str(st)
    if '-' in st:
        st = 'NA'
    mlst.append(f'{id}\\t{st}\\n')

  with open('mlst_formatted.tsv','w') as outFile:
    outFile.write('Sample\\tMLST Scheme\\n')
    for scheme in mlst:
      outFile.write(scheme)

  """
}

if (params.report && !params.ar) {

  report = file(params.report)
  logo = file(params.logo)

  process render{
    publishDir "${params.outdir}/results", mode: 'copy'
    stageInMode = "copy"

    input:
    file snp from snp_mat
    file tree from cgtree
    file rmd from report
    file dryad_logo from logo

    output:
    file("cluster_report.pdf")
    file("report_template.Rmd")

    shell:
    """
    Rscript /reports/render.R ${snp} ${tree} ${rmd}
    mv report.pdf cluster_report.pdf
    mv ${rmd} report_template.Rmd
    """
  }
}

if (params.report && params.ar) {

  report = file(params.report)
  logo = file(params.logo)

  process renderWithAR{
    publishDir "${params.outdir}/results", mode: 'copy'
    stageInMode = "copy"

    input:
    file snp from snp_mat
    file tree from cgtree
    file ar from ar_tsv
    file rmd from report
    file dryad_logo from logo

    output:
    file("cluster_report.pdf")
    file("report_template.Rmd")

    shell:
    """
    Rscript /reports/render.R ${snp} ${tree} ${rmd} ${ar}
    mv report.pdf cluster_report.pdf
    mv ${rmd} report_template.Rmd
    """
  }
}
