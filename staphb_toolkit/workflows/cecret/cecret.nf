#!/usr/bin/env nextflow

println("Currently using the Cecret workflow for use with amplicon-based Illumina hybrid library prep on MiSeq\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.2.1.2021117")
println("")

params.reads = workflow.launchDir + '/reads'
params.single_reads = workflow.launchDir + '/single_reads'
params.fastas = workflow.launchDir + '/fastas'
if ( params.reads == params.single_reads ) {
  println("'params.reads' and 'params.single_reads' cannot point to the same directory!")
  println("'params.reads' is set to " + params.reads)
  println("'params.single_reads' is set to " + params.single_reads)
  exit 1
}
params.outdir = workflow.launchDir + '/cecret'

Channel
  .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fastq.gz,fq,fq.gz}", size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1], "paired" ) }
  .view { "Fastq file found : ${it[0]}" }
  .into { paried_reads_check ; paired_reads }

Channel
  .fromPath("${params.single_reads}/*.{fastq,fastq.gz,fq,fq.gz}")
  .map { reads -> tuple(reads.simpleName, reads, "single" ) }
  .view { "Fastq file found : ${it[0]}" }
  .into { single_reads_check ; single_reads }

Channel
  .fromPath("${params.fastas}/*{.fa,.fasta,.fna}", type:'file')
  .map { fasta -> tuple(fasta.baseName, fasta ) }
  .view { "Fasta file found : ${it[0]}" }
  .into { fastas_check ; fastas }

paried_reads_check
  .concat(single_reads_check)
  .concat(fastas_check)
  .ifEmpty{
    println("FATAL : No input files were found!")
    println("No paired-end fastq files were found at ${params.reads}." )
    println("Set 'params.reads' to directory with paired-end reads")
    println("No single-end fastq files were found at ${params.single_reads}." )
    println("Set 'params.single_reads' to directory with single-end reads")
    println("No fasta files were found at ${params.fastas}." )
    println("Set 'params.fastas' to directory with fastas.")
    exit 1
  }

// reference files for SARS-CoV-2 (part of the github repository)
params.reference_genome = workflow.projectDir + "/configs/MN908947.3.fasta"
params.gff_file = workflow.projectDir + "/configs/MN908947.3.gff"
params.primer_bed = workflow.projectDir + "/configs/artic_V3_nCoV-2019.bed"
params.amplicon_bed = workflow.projectDir + "/configs/nCoV-2019.insert.bed"

params.trimmer = 'ivar'
params.cleaner = 'seqyclean'
params.aligner = 'bwa'
params.msa = 'mafft'

// to toggle off processes
params.bcftools_variants = false // fails to download a lot
params.fastqc = true
params.ivar_variants = true
params.samtools_stats = true
params.samtools_coverage = true
params.samtools_depth = true
params.samtools_flagstat = true
params.samtools_ampliconstats = true
params.samtools_plot_ampliconstats = true
params.bedtools_multicov = true
params.nextclade = true
params.pangolin = true
params.bamsnap = false // can be really slow
params.rename = false
params.filter = false
params.vadr = true

// for optional contamination determination
params.kraken2 = false
params.kraken2_db = ''
params.kraken2_organism = "Severe acute respiratory syndrome-related coronavirus"

// for optional route of tree generation and counting snps between samples
params.relatedness = false
params.snpdists = true
params.iqtree2 = true

// for optional renaming of files for GISAID and GenBank submissions
params.sample_file = workflow.launchDir + '/covid_samples.csv'
params.gisaid_threshold = '25000'
params.genbank_threshold = '15000'

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

// This is where the results will be
println("The files and directory for results is " + params.outdir)
println("A table summarizing results will be created: ${params.outdir}/summary.txt and ${workflow.launchDir}/cecret_run_results.txt\n")

Channel
  .fromPath(params.reference_genome, type:'file')
  .ifEmpty{
    println("No reference genome was selected. Set with 'params.reference_genome'")
    exit 1
  }
  .view { "Reference Genome : $it"}
  .into { reference_genome ; reference_genome2 ; reference_genome_msa ; reference_genome_bamsnap ; reference_genome3 }

Channel
  .fromPath(params.gff_file, type:'file')
  .view { "GFF file for Reference Genome : $it"}
  .set { gff_file }

Channel
  .fromPath(params.primer_bed, type:'file')
  .ifEmpty{
    println("A bedfile for primers is required. Set with 'params.primer_bed'.")
    exit 1
  }
  .view { "Primer BedFile : $it"}
  .into { primer_bed ; primer_bed_bedtools ; primer_bed_ampliconstats }

amplicon_bed = params.bedtools_multicov
  ? Channel.fromPath(params.amplicon_bed, type:'file').view { "Amplicon BedFile : $it"}
  : Channel.empty()

kraken2_db = params.kraken2
  ? Channel.fromPath(params.kraken2_db, type:'dir').view { "Kraken2 database : $it" }
  : Channel.empty()

paired_reads
  .concat(single_reads)
  .ifEmpty{ println("No fastq or fastq.gz files were found at ${params.reads} or ${params.single_reads}") }
  .into { fastq_reads_seqyclean ; fastq_reads_fastp ; fastq_reads_fastqc ; fastq_reads_rename }

println("") // just for aesthetics

params.fastqc_options = ''
process fastqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag "$sample"
  cpus 1
  container 'staphb/fastqc:latest'

  when:
  params.fastqc && sample != null

  input:
  tuple val(sample), file(raw), val(type) from fastq_reads_fastqc

  output:
  file("${task.process}/*.{html,zip}")
  tuple sample, env(raw_1) into fastqc_1_results
  tuple sample, env(raw_2) into fastqc_2_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastqc --version >> $log_file

    fastqc !{params.fastqc_options} \
      --outdir !{task.process} \
      --threads !{task.cpus} \
      !{raw} \
      2>> $err_file >> $log_file

    zipped_fastq=($(ls !{task.process}/*fastqc.zip) "")

    raw_1=$(unzip -p ${zipped_fastq[0]} */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' )
    raw_2=NA
    if [ -f "${zipped_fastq[1]}" ] ; then raw_2=$(unzip -p !{task.process}/*fastqc.zip */fastqc_data.txt | grep "Total Sequences" | awk '{ print $3 }' ) ; fi

    if [ -z "$raw_1" ] ; then raw_1="0" ; fi
    if [ -z "$raw_2" ] ; then raw_2="0" ; fi
  '''
}

if ( params.cleaner == 'seqyclean' ) {
  params.seqyclean_contaminant_file="/Adapters_plus_PhiX_174.fasta"
  params.seqyclean_options = '-minlen 25 -qual'
  process seqyclean {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 1
    container 'staphb/seqyclean:latest'

    when:
    params.cleaner == 'seqyclean' && sample != null

    input:
    tuple val(sample), file(reads), val(paired_single) from fastq_reads_seqyclean

    output:
    tuple sample, file("${task.process}/${sample}_clean_PE{1,2}.fastq.gz") optional true into paired_files
    tuple sample, file("${task.process}/${sample}_cln_SE.fastq.gz") optional true into single_files
    tuple sample, file("${task.process}/${sample}_clean_PE{1,2}.fastq.gz"), val(paired_single) optional true into paired_files_kraken2
    tuple sample, file("${task.process}/${sample}_cln_SE.fastq.gz"), val(paired_single) optional true into single_files_kraken2
    file("${task.process}/${sample}_cl*n_SummaryStatistics.tsv") into seqyclean_files
    file("${task.process}/${sample}_cl*n_SummaryStatistics.txt")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
    tuple sample, env(perc_kept) into seqyclean_perc_kept_results
    tuple sample, env(kept) into seqyclean_pairskept_results
    tuple sample, env(cleaner_version) into cleaner_version

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      echo "seqyclean version: $(seqyclean -h | grep Version)" >> $log_file
      cleaner_version="seqyclean : $(seqyclean -h | grep Version)"

      kept=''
      perc_kept=''

      if [ "!{paired_single}" == "single" ]
      then
        seqyclean !{params.seqyclean_options} \
          -c !{params.seqyclean_contaminant_file} \
          -U !{reads} \
          -o !{task.process}/!{sample}_cln \
          -gz \
          2>> $err_file >> $log_file
        kept=$(cut -f 36 !{task.process}/!{sample}_cln_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
        perc_kept=$(cut -f 37 !{task.process}/!{sample}_cln_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
      else
        seqyclean !{params.seqyclean_options} \
          -c !{params.seqyclean_contaminant_file} \
          -1 !{reads[0]} -2 !{reads[1]} \
          -o !{task.process}/!{sample}_clean \
          -gz \
          2>> $err_file >> $log_file
        kept=$(cut -f 58 !{task.process}/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
        perc_kept=$(cut -f 59 !{task.process}/!{sample}_clean_SummaryStatistics.tsv | grep -v "Kept" | head -n 1)
      fi

      if [ -z "$kept" ] ; then kept="0" ; fi
      if [ -z "$perc_kept" ] ; then perc_kept="0" ; fi
    '''
  }

  fastp_results=Channel.empty()

  seqyclean_files
    .collectFile(name: "Combined_SummaryStatistics.tsv",
      keepHeader: true,
      sort: true,
      storeDir: "${params.outdir}/seqyclean")

} else if ( params.cleaner == 'fastp' ) {
  params.fastp_options = ''
  process fastp {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 1
    container 'bromberglab/fastp:latest'

    input:
    tuple val(sample), file(reads), val(paired_single) from fastq_reads_fastp

    output:
    tuple sample, file("${task.process}/${sample}_clean_PE{1,2}.fastq.gz") optional true into paired_files
    tuple sample, file("${task.process}/${sample}_cln.fastq.gz") optional true into single_files
    tuple sample, file("${task.process}/${sample}_clean_PE{1,2}.fastq.gz"), val(paired_single) optional true into paired_files_kraken2
    tuple sample, file("${task.process}/${sample}_cln.fastq.gz"), val(paired_single) optional true into single_files_kraken2
    file("${task.process}/${sample}_fastp.{html,json}")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
    tuple sample, env(passed_reads) into fastp_results
    tuple sample, env(cleaner_version) into cleaner_version

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      fastp --version >> $log_file 2>> $err_file
      cleaner_version="fastp : $(fastp --version | head -n 1)"

      if [ "!{paired_single}" == "single" ]
      then
        fastp !{params.fastp_options} \
          -i !{reads} \
          -o !{task.process}/!{sample}_cln.fastq.gz \
          -h !{task.process}/!{sample}_fastp.html \
          -j !{task.process}/!{sample}_fastp.json \
          2>> $err_file >> $log_file
      else
        fastp !{params.fastp_options} \
          -i !{reads[0]} \
          -I !{reads[1]} \
          -o !{task.process}/!{sample}_clean_PE1.fastq.gz \
          -O !{task.process}/!{sample}_clean_PE2.fastq.gz \
          -h !{task.process}/!{sample}_fastp.html \
          -j !{task.process}/!{sample}_fastp.json \
          2>> $err_file >> $log_file
      fi

      passed_reads=$(grep "reads passed filter" $err_file | tail -n 1 | cut -f 2 -d ":" | sed 's/ //g' )
      if [ -z "$passed_reads" ] ; then passed_reads="0" ; fi
    '''
  }
  seqyclean_perc_kept_results=Channel.empty()
  seqyclean_pairskept_results=Channel.empty()
}

if (params.aligner == 'bwa') {
  process bwa {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/bwa/*.{log,err}"
    tag "${sample}"
    cpus params.maxcpus
    container 'staphb/bwa:latest'

    input:
    tuple val(sample), file(reads), file(reference_genome) from paired_files.concat(single_files).combine(reference_genome)

    output:
    tuple sample, file("aligned/${sample}.sam") into sams, sams_filter
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
    tuple sample, env(bwa_version) into aligner_version

    shell:
    '''
      mkdir -p aligned logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
      bwa_version="bwa : "$(bwa 2>&1 | grep Version)

      # index the reference fasta file
      bwa index !{reference_genome}

      # bwa mem command
      bwa mem -t !{task.cpus} !{reference_genome} !{reads} 2>> $err_file > aligned/!{sample}.sam
    '''
  }
} else if (params.aligner == 'minimap2') {
  params.minimap2_options = '-K 20M'
  process minimap2 {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/minimap2/*.{log,err}"
    tag "${sample}"
    cpus params.maxcpus
    container 'staphb/minimap2:latest'

    input:
    tuple val(sample), file(reads), file(reference_genome) from paired_files.concat(single_files).combine(reference_genome)

    output:
    tuple sample, file("aligned/${sample}.sam") into sams, sams_filter
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
    tuple sample, env(minimap2_version) into aligner_version

    shell:
    '''
      mkdir -p aligned logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      minimap2 --version >> $log_file
      minimap2_version=$(echo "minimap2 : "$(minimap2 --version))

      minimap2 !{params.minimap2_options} \
        -ax sr -t !{task.cpus} \
        -o aligned/!{sample}.sam \
        !{reference_genome} !{reads} 2>> $err_file >> $log_file
    '''
  }
}

process sort {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.maxcpus
  container 'staphb/samtools:latest'

  input:
  tuple val(sample), file(sam) from sams

  output:
  tuple sample, file("aligned/${sample}.sorted.bam") into pre_trim_bams, pre_trim_bams2
  tuple sample, file("aligned/${sample}.sorted.bam"), file("aligned/${sample}.sorted.bam.bai") into pre_trim_bams_bamsnap
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p aligned logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools sort !{sam} 2>> $err_file | \
      samtools view -F 4 -o aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file

    # indexing the bams
    samtools index aligned/!{sample}.sorted.bam 2>> $err_file >> $log_file
  '''
}

params.filter_options = ''
process filter {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/samtools:latest'

  when:
  params.filter

  input:
  tuple val(sample), file(sam) from sams_filter

  output:
  tuple sample, file("${task.process}/${sample}_filtered_{R1,R2}.fastq.gz") optional true into filtered_reads
  file("${task.process}/${sample}_filtered_unpaired.fastq.gz") optional true
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools sort -n !{sam} 2>> $err_file | \
      samtools fastq -F 4 !{params.filter_options} \
      -s !{task.process}/!{sample}_filtered_unpaired.fastq.gz \
      -1 !{task.process}/!{sample}_filtered_R1.fastq.gz \
      -2 !{task.process}/!{sample}_filtered_R2.fastq.gz \
      2>> $err_file >> $log_file
  '''
}

if (params.trimmer == 'ivar' ) {
  params.ivar_trim_options = ''
  process ivar_trim {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 1
    container 'staphb/ivar:latest'

    input:
    tuple val(sample), file(bam), file(primer_bed) from pre_trim_bams.combine(primer_bed)

    output:
    tuple sample, file("${task.process}/${sample}.primertrim.sorted.bam") into trimmed_bams, trimmed_bams4, trimmed_bams5
    tuple sample, file("${task.process}/${sample}.primertrim.sorted.bam"), file("ivar_trim/${sample}.primertrim.sorted.bam.bai") into bam_bai
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
    tuple sample, env(trimmer_version) into trimmer_version

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      ivar version >> $log_file
      trimmer_version="ivar : $(ivar version)"

      # trimming the reads
      ivar trim !{params.ivar_trim_options} -e -i !{bam} -b !{primer_bed} -p !{task.process}/!{sample}.primertrim 2>> $err_file >> $log_file

      # sorting and indexing the trimmed bams
      samtools sort !{task.process}/!{sample}.primertrim.bam -o !{task.process}/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
      samtools index !{task.process}/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
    '''
  }

} else if (params.trimmer == 'samtools') {
  params.samtools_ampliconclip_options = ''
  process samtools_ampliconclip {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 1
    container 'staphb/samtools:latest'

    input:
    tuple val(sample), file(bam), file(primer_bed) from pre_trim_bams.combine(primer_bed)

    output:
    tuple sample, file("${task.process}/${sample}.primertrim.sorted.bam") into trimmed_bams, trimmed_bams4, trimmed_bams5
    tuple sample, file("${task.process}/${sample}.primertrim.sorted.bam"), file("${task.process}/${sample}.primertrim.sorted.bam.bai") into bam_bai
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
    tuple sample, env(trimmer_version) into trimmer_version

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      samtools --version >> $log_file
      trimmer_version="samtools ampliconclip : $(samtools --version | head -n 1)"

      # trimming the reads
      samtools ampliconclip !{params.samtools_ampliconclip_options} -b !{primer_bed} !{bam} 2>> $err_file | \
        samtools sort 2>> $err_file |  \
        samtools view -F 4 -o !{task.process}/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file

      samtools index !{task.process}/!{sample}.primertrim.sorted.bam 2>> $err_file >> $log_file
    '''
  }
} else if (params.trimmer == 'none') {
  pre_trim_bams
    .into { trimmed_bams ; trimmed_bams4 ; trimmed_bams5 }
}

trimmed_bams
 .combine(reference_genome2)
 .into { trimmed_bams_genome ; trimmed_bams_ivar_consensus ; trimmed_bams_bcftools_variants }

params.minimum_depth = 100
params.mpileup_depth = 8000
params.ivar_variants_options = '-q 20 -t 0.6'
process ivar_variants {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/ivar:latest'
  memory {2.GB * task.attempt}
  errorStrategy 'retry'
  maxRetries 2

  when:
  params.ivar_variants

  input:
  tuple val(sample), file(bam), file(reference_genome), file(gff_file) from trimmed_bams_genome.combine(gff_file)

  output:
  tuple sample, file("${task.process}/${sample}.variants.tsv")
  tuple sample, file("${task.process}/${sample}.ivar_variants.vcf") into ivar_variant_file
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(variants_num) into ivar_variants_results

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{reference_genome} !{bam} 2>> $err_file | \
      ivar variants -p !{task.process}/!{sample}.variants !{params.ivar_variants_options} -m !{params.minimum_depth} -r !{reference_genome} -g !{gff_file} 2>> $err_file >> $log_file

    variants_num=$(grep "TRUE" !{task.process}/!{sample}.variants.tsv | wc -l)

    if [ -z "$variants_num" ] ; then variants_num="0" ; fi

    echo '##fileformat=VCFv4.2' > !{task.process}/!{sample}.ivar_variants.vcf
    echo '##source=iVar' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Deapth of alternate base on reverse reads">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo '##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">' >> !{task.process}/!{sample}.ivar_variants.vcf
    echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t!{sample}' >> !{task.process}/!{sample}.ivar_variants.vcf
    tail -n+2 !{task.process}/!{sample}.variants.tsv | \
      awk '{print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\t.\tREF_DP=" $5 ";REF_RV=" $6 ";REF_QUAL=" $7 ";ALT_DP=" $8 ";ALT_RV=" $9 ";ALT_QUAL=" $10 "\tGT:PL\t1/1:" $12 "," $12-$8 "," $8 }' \
      >> !{task.process}/!{sample}.ivar_variants.vcf
  '''
}

params.ivar_consensus_options = '-q 20 -t 0.6 -n N'
process ivar_consensus {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/ivar_consensus/*.{log,err}"
  publishDir "${params.outdir}", mode: 'copy', pattern: "consensus/*.consensus.fa"
  tag "${sample}"
  cpus 1
  container 'staphb/ivar:latest'
  memory {2.GB * task.attempt}
  errorStrategy {'retry'}
  maxRetries 2

  input:
  tuple val(sample), file(bam), file(reference_genome) from trimmed_bams_ivar_consensus

  output:
  tuple sample, file("consensus/${sample}.consensus.fa") into consensus_rename, consensus_pangolin, consensus_vadr, consensus_nextclade
  file("consensus/${sample}.consensus.fa") into consensus_msa
  file("consensus/${sample}.consensus.qual.txt")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results
  tuple sample, env(ivar_version) into ivar_version

  shell:
  '''
    mkdir -p consensus logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file
    ivar version >> $log_file
    ivar_version=$(ivar version | grep "version")

    samtools mpileup -A -d !{params.mpileup_depth} -B -Q 0 --reference !{reference_genome} !{bam} 2>> $err_file | \
      ivar consensus !{params.ivar_consensus_options} -m !{params.minimum_depth} -p consensus/!{sample}.consensus 2>> $err_file >> $log_file

    if [ -f "consensus/!{sample}.consensus.fa" ]
    then
      num_N=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o 'N' | wc -l )
      num_ACTG=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o -E "C|A|T|G" | wc -l )
      num_degenerate=$(grep -v ">" consensus/!{sample}.consensus.fa | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )

      if [ -z "$num_N" ] ; then num_N="0" ; fi
      if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
      if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
    else
      num_N="0"
      num_ACTG="0"
      num_degenerate="0"
    fi

    if [ -z "$num_N" ] ; then num_N="0" ; fi
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
    num_total=$(( $num_N + $num_degenerate + $num_ACTG ))
  '''
}

process fasta_prep {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${fasta}"
  cpus 1
  container 'staphb/parallel-perl:latest'

  when:
  sample != null && sample != 'input.1'

  input:
  tuple val(sample), file(fasta) from fastas

  output:
  tuple sample, file("${task.process}/${fasta}") optional true into fastas_rename, fastas_pangolin, fastas_vadr, fastas_nextclade
  file("${task.process}/${fasta}") into fastas_msa

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    if [ !{fasta} != 'null' ]
    then
      echo ">!{sample}" > !{task.process}/!{fasta}
      grep -v ">" !{fasta} | fold -w 75 >> !{task.process}/!{fasta}
    fi
  '''
}

process bcftools_variants {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/bcftools:latest'

  when:
  params.bcftools_variants

  input:
  tuple val(sample), file(bam), file(reference_genome) from trimmed_bams_bcftools_variants

  output:
  tuple sample, file("${task.process}/${sample}.vcf") into bcftools_variants_file
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(variants_num) into bcftools_variants_results

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bcftools --version >> $log_file

    bcftools mpileup -A -d !{params.mpileup_depth} -B -Q 0 -f !{reference_genome} !{bam} 2>> $err_file | \
      bcftools call -mv -Ov -o !{task.process}/!{sample}.vcf 2>> $err_file >> $log_file

    variants_num=$(grep -v "#" bcftools_variants/!{sample}.vcf | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
  '''
}

pre_trim_bams_bamsnap
  .join(ivar_variant_file, remainder: true, by:0)
  .join(bcftools_variants_file, remainder: true, by:0)
  .combine(reference_genome_bamsnap)
  .set { bamsnap_files }

params.bamsnap_options = ''
process bamsnap {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  errorStrategy 'ignore'
  container 'danielmsk/bamsnap:latest'
  time '1h'

  when:
  params.bamsnap

  input:
  tuple val(sample), file(bam), file(bai), file(ivar), file(bcftools), file(reference_genome) from bamsnap_files

  output:
  file("${task.process}/${sample}/{ivar,bcftools}/*.{png,log}") optional true
  file("${task.process}/${sample}/*.{png,log}") optional true
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bamsnap --version >> $log_file

    if [[ "!{ivar}" != *"input"* ]]
    then
      mkdir -p bamsnap/!{sample}
      bamsnap \
        -draw coordinates bamplot coverage base \
        !{params.bamsnap_options} \
        -process !{task.cpus} \
        -ref !{reference_genome} \
        -bam !{bam} \
        -vcf !{ivar} \
        -out !{task.process}/!{sample}/ivar \
        -imagetype png \
        -save_image_only 2>> $err_file >> $log_file
    fi

    if [[ "!{bcftools}" != *"input"* ]]
    then
      mkdir -p bamsnap/!{sample}
      bamsnap \
        -draw coordinates bamplot coverage base \
        !{params.bamsnap_options} \
        -process !{task.cpus} \
        -ref !{reference_genome} \
        -bam !{bam} \
        -vcf !{bcftools} \
        -out !{task.process}/!{sample}/bcftools \
        -imagetype png \
        -save_image_only 2>> $err_file >> $log_file
    fi
  '''
}

pre_trim_bams2
   .combine(trimmed_bams4, by: 0)
   .into { pre_post_bams ; pre_post_bams2 ; pre_post_bams3 ; pre_post_bams4 }

params.samtools_stats_options = ''
process samtools_stats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/samtools:latest'

  when:
  params.samtools_stats

  input:
  tuple val(sample), file(aligned), file(trimmed) from pre_post_bams

  output:
  file("${task.process}/aligned/${sample}.stats.txt")
  file("${task.process}/trimmed/${sample}.stats.trim.txt")
  tuple sample, env(insert_size_before_trimming) into samtools_stats_before_size_results
  tuple sample, env(insert_size_after_trimming) into samtools_stats_after_size_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/aligned !{task.process}/trimmed logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools stats !{params.samtools_stats_options} !{aligned} 2>> $err_file > !{task.process}/aligned/!{sample}.stats.txt
    samtools stats !{params.samtools_stats_options} !{trimmed} 2>> $err_file > !{task.process}/trimmed/!{sample}.stats.trim.txt

    insert_size_before_trimming=$(grep "insert size average" !{task.process}/aligned/!{sample}.stats.txt | cut -f 3)
    insert_size_after_trimming=$(grep "insert size average" !{task.process}/trimmed/!{sample}.stats.trim.txt | cut -f 3)
  '''
}

params.samtools_coverage_options = ''
process samtools_coverage {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/samtools:latest'

  when:
  params.samtools_coverage

  input:
  tuple val(sample), file(aligned), file(trimmed) from pre_post_bams2

  output:
  file("${task.process}/{aligned,trimmed}/${sample}.cov.{txt,hist}")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(coverage) into samtools_coverage_results
  tuple sample, env(depth) into samtools_covdepth_results

  shell:
  '''
    mkdir -p !{task.process}/aligned !{task.process}/trimmed logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools coverage !{params.samtools_coverage_options} !{aligned} -m -o !{task.process}/aligned/!{sample}.cov.hist 2>> $err_file >> $log_file
    samtools coverage !{params.samtools_coverage_options} !{aligned}    -o !{task.process}/aligned/!{sample}.cov.txt 2>> $err_file >> $log_file
    samtools coverage !{params.samtools_coverage_options} !{trimmed} -m -o !{task.process}/trimmed/!{sample}.cov.trim.hist 2>> $err_file >> $log_file
    samtools coverage !{params.samtools_coverage_options} !{trimmed}    -o !{task.process}/trimmed/!{sample}.cov.trim.txt 2>> $err_file >> $log_file

    coverage=$(cut -f 6 !{task.process}/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    depth=$(cut -f 7 !{task.process}/trimmed/!{sample}.cov.trim.txt | tail -n 1)
    if [ -z "$coverage" ] ; then coverage_trim="0" ; fi
    if [ -z "$depth" ] ; then depth_trim="0" ; fi
  '''
}

params.samtools_flagstat_options = ''
process samtools_flagstat {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/samtools:latest'

  input:
  tuple val(sample), file(aligned), file(trimmed) from pre_post_bams3

  when:
  params.samtools_flagstat

  output:
  file("${task.process}/{aligned,trimmed}/${sample}.flagstat.txt")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/aligned !{task.process}/trimmed logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools flagstat !{params.samtools_flagstat_options} \
      !{aligned} \
      2>> $err_file > !{task.process}/aligned/!{sample}.flagstat.txt

    samtools flagstat !{params.samtools_flagstat_options} \
      !{trimmed} \
      2>> $err_file > !{task.process}/trimmed/!{sample}.flagstat.txt
  '''
}

params.samtools_depth_options = ''
process samtools_depth {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/samtools:latest'

  input:
  tuple val(sample), file(aligned), file(trimmed) from pre_post_bams4

  when:
  params.samtools_depth

  output:
  file("${task.process}/{aligned,trimmed}/${sample}.depth.txt")
  tuple sample, env(depth) into samtools_depth_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/aligned !{task.process}/trimmed logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools depth !{params.samtools_depth_options} \
      !{aligned} \
      2>> $err_file > !{task.process}/aligned/!{sample}.depth.txt

    samtools depth !{params.samtools_depth_options} \
      !{trimmed} \
      2>> $err_file > !{task.process}/trimmed/!{sample}.depth.txt

    depth=$(awk '{ if ($3 > !{params.minimum_depth} ) print $0 }' !{task.process}/trimmed/!{sample}.depth.txt | wc -l )
    if [ -z "$depth" ] ; then depth="0" ; fi
  '''
}

params.kraken2_options = ''
process kraken2 {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.maxcpus
  container 'staphb/kraken2:latest'

  when:
  params.kraken2

  input:
  tuple val(sample), file(clean), val(paired_single), path(kraken2_db) from paired_files_kraken2.concat(single_files_kraken2).combine(kraken2_db)

  output:
  file("${task.process}/${sample}_kraken2_report.txt")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(percentage_cov) into kraken2_sars_results
  tuple sample, env(percentage_human) into kraken2_human_results

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    kraken2 --version >> $log_file

    if [ ! -d !{kraken2_db} ]
    then
      echo "Kraken2 database could not be found. Please specify with params.kraken2_db" | tee -a $err_file
    fi

    if [ "!{paired_single}" == "single" ]
    then
      kraken2 !{params.kraken2_options} \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{kraken2_db} \
        !{clean} \
        --report !{task.process}/!{sample}_kraken2_report.txt \
        2>> $err_file >> $log_file
    else
      kraken2 !{params.kraken2_options} \
        --paired \
        --classified-out cseqs#.fq \
        --threads !{task.cpus} \
        --db !{kraken2_db} \
        !{clean} \
        --report !{task.process}/!{sample}_kraken2_report.txt \
        2>> $err_file >> $log_file
    fi

    percentage_human=$(grep "Homo sapiens" !{task.process}/!{sample}_kraken2_report.txt | awk '{print $1}')
    percentage_cov=$(grep "!{params.kraken2_organism}" !{task.process}/!{sample}_kraken2_report.txt | awk '{print $1}')

    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_cov" ] ; then percentage_cov="0" ; fi
  '''
}

params.bedtools_multicov_options = '-f .1'
process bedtools_multicov {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/bedtools:latest'

  when:
  params.bedtools_multicov

  input:
  tuple val(sample), file(bam), file(bai), file(amplicon_bed) from bam_bai.combine(amplicon_bed)

  output:
  file("${task.process}/${sample}.multicov.txt")
  tuple sample, env(num_failed_amplicons) into bedtools_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    bedtools --version >> $log_file

    bedtools multicov !{params.bedtools_multicov_options} \
      -bams !{bam} \
      -bed !{amplicon_bed} \
      2>> $err_file >> !{task.process}/!{sample}.multicov.txt

    result_column=$(head -n 1 !{task.process}/!{sample}.multicov.txt | awk '{print NF}' )
    num_failed_amplicons=$(cat !{task.process}/!{sample}.multicov.txt | tr ' ' '\t' | cut -f $result_column | awk '{ if ( $1 < 20 ) print $0 }' | wc -l )
    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons="NA" ; fi
  '''
}

params.samtools_ampliconstats_options = ''
process samtools_ampliconstats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/samtools:latest'

  when:
  params.samtools_ampliconstats

  input:
  tuple val(sample), file(bam), file(primer_bed) from trimmed_bams5.combine(primer_bed_ampliconstats)

  output:
  tuple sample, file("${task.process}/${sample}_ampliconstats.txt") into samtools_ampliconstats_files
  tuple sample, env(num_failed_amplicons) into samtools_ampliconstats_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    samtools ampliconstats !{params.samtools_ampliconstats_options} \
      !{primer_bed} \
      !{bam} \
      2>> $err_file > !{task.process}/!{sample}_ampliconstats.txt

    num_failed_amplicons=$(grep ^FREADS !{task.process}/!{sample}_ampliconstats.txt | cut -f 2- | tr '\t' '\n' | awk '{ if ($1 < 20) print $0 }' | wc -l)
    if [ -z "$num_failed_amplicons" ] ; then num_failed_amplicons=0 ; fi
  '''
}

params.samtools_plot_ampliconstats_options = '-size 1200,900 -size2 1200,900 -size3 1200,900'
process samtools_plot_ampliconstats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/samtools:latest'
  errorStrategy 'ignore'

  when:
  params.samtools_plot_ampliconstats

  input:
  tuple val(sample), file(ampliconstats) from samtools_ampliconstats_files

  output:
  file("${task.process}/${sample}*")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    samtools --version >> $log_file

    plot-ampliconstats !{params.samtools_plot_ampliconstats_options} \
      !{task.process}/!{sample} \
      !{ampliconstats}
  '''
}

params.pangolin_options = ''
process pangolin {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/pangolin:latest'

  when:
  params.pangolin

  input:
  tuple val(sample), file(fasta) from consensus_pangolin.concat(fastas_pangolin)

  output:
  file("${task.process}/${sample}/lineage_report.csv") into pangolin_files
  tuple sample, env(pangolin_lineage) into pangolin_lineage
  tuple sample, env(pangolin_status) into pangolin_status
  tuple sample, env(pangolin_scorpio) into pangolin_scorpio
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(pangolin_version) into pangolin_version
  tuple sample, env(pangolearn_version) into pangolearn_version
  tuple sample, env(constellations_version) into constellations_version
  tuple sample, env(scorpio_version) into scorpio_version

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    pangolin --version >> $log_file
    pangolin --pangoLEARN-version >> $log_file

    pangolin_version=$(pangolin --all-versions | grep pangolin | cut -d ' ' -f 2)
    pangolearn_version=$(pangolin --all-versions | grep pangolearn | cut -d ' ' -f 2)
    constellations_version=$(pangolin --all-versions | grep constellations | cut -d ' ' -f 2)
    scorpio_version=$(pangolin --all-versions | grep scorpio | cut -d ' ' -f 2)

    pangolin !{params.pangolin_options} \
      --outdir !{task.process}/!{sample}   \
      !{fasta} \
      2>> $err_file >> $log_file

    lineage_column=$(head -n 1 !{task.process}/!{sample}/lineage_report.csv | tr ',' '\\n' | grep -n "lineage"      | cut -f 1 -d ":" )
    status_column=$(head  -n 1 !{task.process}/!{sample}/lineage_report.csv | tr ',' '\\n' | grep -n "status"       | cut -f 1 -d ":" )
    scorpio_column=$(head -n 1 !{task.process}/!{sample}/lineage_report.csv | tr ',' '\\n' | grep -n "scorpio_call" | cut -f 1 -d ":" )

    if [ -n "$lineage_column" ]
    then
      pangolin_lineage=$(grep "Consensus_!{sample}.consensus_threshold" !{task.process}/!{sample}/lineage_report.csv | cut -f $lineage_column -d ",")
    else
      pangolin_lineage="Not Found"
    fi

    if [ -n "$status_column" ]
    then
      pangolin_status=$(grep "Consensus_!{sample}.consensus_threshold" !{task.process}/!{sample}/lineage_report.csv | cut -f $status_column -d ",")
    else
      pangolin_lineage="Not Found"
    fi

    if [ -n "$scorpio_column" ]
    then
      pangolin_scorpio=$(grep "Consensus_!{sample}.consensus_threshold" !{task.process}/!{sample}/lineage_report.csv | cut -f $scorpio_column -d ",")
    else
      pangolin_lineage="Not Found"
    fi

    if [ -z "$pangolin_lineage" ] ; then pangolin_lineage="NA" ; fi
    if [ -z "$pangolin_status" ]  ; then pangolin_status="NA"  ; fi
    if [ -z "$pangolin_scorpio" ] ; then pangolin_scorpio="NA" ; fi
  '''
}

pangolin_files
  .collectFile(name: "combined_lineage_report.csv",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/pangolin")

params.nextclade_prep_options = '--name sars-cov-2'
process nextclade_prep {
  tag "Downloading SARS-CoV-2 dataset"
  cpus 1
  container 'nextstrain/nextclade:latest'

  when:
  params.nextclade || params.msa == 'nextalign'

  output:
  path("${task.process}") into prepped_nextclade, prepped_nextalign
  file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    nextclade --version >> $log_file

    nextclade dataset get !{params.nextclade_prep_options} --output-dir !{task.process}
  '''
}

params.nextclade_options = ''
process nextclade {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'nextstrain/nextclade:latest'

  when:
  params.nextclade

  input:
  tuple val(sample), file(fasta), path(dataset) from consensus_nextclade.concat(fastas_nextclade).combine(prepped_nextclade)

  output:
  file("${task.process}/${sample}/${sample}_nextclade.csv") into nextclade_files
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
  tuple sample, env(nextclade_clade) into nextclade_clade_results
  tuple sample, env(nextclade_version) into nextclade_version

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    nextclade --version >> $log_file
    nextclade_version=$(nextclade --version)

    nextclade !{params.nextclade_options} \
      --input-fasta=!{fasta} \
      --input-dataset !{dataset} \
      --output-json=!{task.process}/!{sample}/!{sample}_nextclade.json \
      --output-csv=!{task.process}/!{sample}/!{sample}_nextclade.csv \
      --output-tsv=!{task.process}/!{sample}/!{sample}_nextclade.tsv \
      --output-tree=!{task.process}/!{sample}/!{sample}_nextclade.auspice.json \
      --output-dir=!{task.process}/!{sample} \
      --output-basename=!{sample} \
      2>> $err_file >> $log_file

    nextclade_column=$(head -n 1 !{task.process}/!{sample}/!{sample}_nextclade.csv | tr ';' '\\n' | grep -wn "clade" | cut -f 1 -d ":" )
    if [ -n "$nextclade_column" ]
    then
      nextclade_clade=$(cat !{task.process}/!{sample}/!{sample}_nextclade.csv | grep !{sample} | cut -f $nextclade_column -d ";" | sed 's/,/;/g' | sed 's/"//g' | head -n 1 )
    else
      nextclade_clade="Not Found"
    fi

    if [ -z "$nextclade_clade" ] ; then nextclade_clade="clade" ; fi
  '''
}

nextclade_files
  .collectFile(name: "combined_nextclade_report.txt",
    keepHeader: true,
    sort: true,
    storeDir: "${params.outdir}/nextclade")

params.vadr_options = '--split --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn'
params.vadr_reference = 'sarscov2'
params.vadr_mdir = '/opt/vadr/vadr-models'
params.vadr_trim_options = '--minlen 50 --maxlen 30000'
process vadr {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/vadr:latest'

  when:
  params.vadr

  input:
  tuple val(sample), file(fasta) from consensus_vadr.concat(fastas_vadr)

  output:
  file("${task.process}/${sample}/*") optional true
  file("${task.process}/${sample}/${sample}.vadr.fail.fa") optional true into vadr_files_fail_fasta
  file("${task.process}/${sample}/${sample}.vadr.fail.list") optional true into vadr_files_fail_list
  file("${task.process}/${sample}/${sample}.vadr.pass.fa") optional true into vadr_files_pass_fasta
  file("${task.process}/${sample}/${sample}.vadr.pass.list") optional true into vadr_files_pass_list
  file("${task.process}/${sample}/${sample}.vadr.sqc") optional true into vadr_files_sqc
  tuple sample, env(pass_fail) into vadr_results
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process} !{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null
    echo "no version" >> $log_file
    v-annotate.pl -h >> $log_file

    fasta-trim-terminal-ambigs.pl !{params.vadr_trim_options} \
      !{fasta} > trimmed_!{fasta}

    if [ -s "trimmed_!{fasta}" ]
    then
      v-annotate.pl !{params.vadr_options} \
        --cpu !{task.cpus} \
        --noseqnamemax \
        --mkey !{params.vadr_reference} \
        --mdir !{params.vadr_mdir} \
        trimmed_!{fasta} \
        !{task.process}/!{sample} \
        2>> $err_file >> $log_file

      pass_fail=$(grep "Consensus_!{sample}.consensus_threshold" !{task.process}/!{sample}/!{sample}.vadr.sqc | awk '{print $4}')
    else
      pass_fail="FAIL"
    fi
    if [ -z "$pass_fail" ] ; then pass_fail="NA" ; fi
  '''
}

vadr_files_fail_fasta
  .collectFile(name: "combined_vadr.fail.fasta",
    keepHeader: false,
    sort: false,
    storeDir: "${params.outdir}/vadr")

vadr_files_fail_list
  .collectFile(name: "combined_vadr.fail.list",
    keepHeader: false,
    sort: true,
    storeDir: "${params.outdir}/vadr")

vadr_files_pass_fasta
  .collectFile(name: "combined_vadr.pass.fasta",
    keepHeader: false,
    sort: false,
    storeDir: "${params.outdir}/vadr")

vadr_files_pass_list
  .collectFile(name: "combined_vadr.pass.list",
    keepHeader: false,
    sort: true,
    storeDir: "${params.outdir}/vadr")

vadr_files_sqc
  .collectFile(name: "combined_vadr.sqc",
    keepHeader: true,
    sort: true,
    skip: 3,
    storeDir: "${params.outdir}/vadr")

consensus_results
//tuple sample, env(num_N), env(num_ACTG), env(num_degenerate), env(num_total) into consensus_results
  .join(fastqc_1_results, remainder: true, by: 0)
  .join(fastqc_2_results, remainder: true, by: 0)
  .join(seqyclean_pairskept_results, remainder: true, by: 0)
  .join(seqyclean_perc_kept_results, remainder: true, by: 0)
  .join(fastp_results, remainder: true, by: 0)
  .join(ivar_variants_results, remainder: true, by: 0)
  .join(bcftools_variants_results, remainder: true, by:0)
  .join(samtools_coverage_results, remainder: true, by: 0)
  .join(samtools_covdepth_results, remainder: true, by: 0)
  .join(samtools_depth_results, remainder: true, by: 0)
  .join(samtools_stats_before_size_results, remainder: true, by: 0)
  .join(samtools_stats_after_size_results, remainder: true, by: 0)
  .join(kraken2_human_results, remainder: true, by: 0)
  .join(kraken2_sars_results, remainder: true, by: 0)
  .join(nextclade_clade_results, remainder: true, by: 0)
  .join(bedtools_results, remainder: true, by: 0)
  .join(samtools_ampliconstats_results, remainder: true, by: 0)
  .join(aligner_version, remainder: true, by: 0)
  .join(trimmer_version, remainder: true, by: 0)
  .join(cleaner_version, remainder: true, by: 0)
  .join(ivar_version, remainder: true, by: 0)
  .join(pangolin_lineage, remainder: true, by: 0)
  .join(pangolin_status, remainder: true, by: 0)
  .join(pangolin_scorpio, remainder: true, by: 0)
  .join(vadr_results, remainder: true, by: 0)
  .join(pangolin_version, remainder: true, by: 0)
  .join(pangolearn_version, remainder: true, by: 0)
  .join(constellations_version, remainder: true, by: 0)
  .join(scorpio_version, remainder: true, by: 0)
  .join(nextclade_version, remainder: true, by: 0)
  .set { results }

process summary {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  tag "${sample}"
  cpus 1
  container 'staphb/parallel-perl:latest'

  input:
  tuple val(sample), val(num_N), val(num_ACTG), val(num_degenerate), val(num_total),
    val(raw_1),
    val(raw_2),
    val(pairskept),
    val(perc_kept),
    val(reads_passed),
    val(ivar_variants),
    val(bcftools_variants),
    val(coverage),
    val(covdepth),
    val(depth),
    val(samtools_stats_before_size_results),
    val(samtools_stats_after_size_results),
    val(percentage_human),
    val(percentage_cov),
    val(nextclade_clade),
    val(bedtools_num_failed_amplicons),
    val(samtools_num_failed_amplicons),
    val(aligner_version),
    val(trimmer_version),
    val(cleaner_version),
    val(ivar_version),
    val(pangolin_lineage),
    val(pangolin_status),
    val(pangolin_scorpio),
    val(vadr_results),
    val(pangolin_version),
    val(pangolearn_version),
    val(constellations_version),
    val(scorpio_version),
    val(nextclade_version) from results

  output:
  file("${task.process}/${sample}.summary.csv") into summary
  file("${task.process}/${sample}.summary.txt") into summary2
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    date | tee -a $log_file $err_file > /dev/null

    sample_id=($(echo !{sample} | cut -f 1 -d "_" ))

    header="sample_id,sample"
    result="${sample_id},!{sample}"

    header="$header,pangolin_lineage,pangolin_status,pangolin_scorpio_call"
    result="$result,!{pangolin_lineage},!{pangolin_status},!{pangolin_scorpio}"

    header="$header,nextclade_clade"
    result="$result,!{nextclade_clade}"

    header="$header,fastqc_raw_reads_1,fastqc_raw_reads_2"
    result="$result,!{raw_1},!{raw_2}"

    header="$header,seqyclean_pairs_kept_after_cleaning,seqyclean_percent_kept_after_cleaning"
    result="$result,!{pairskept},!{perc_kept}"

    header="$header,fastp_reads_passed"
    result="$result,!{reads_passed}"

    header="$header,depth_after_trimming,1X_coverage_after_trimming"
    result="$result,!{covdepth},!{coverage}"

    header="$header,num_pos_!{params.minimum_depth}X"
    result="$result,!{depth}"

    header="$header,insert_size_before_trimming,insert_size_after_trimming"
    result="$result,!{samtools_stats_before_size_results},!{samtools_stats_after_size_results}"

    organism=$(echo "!{params.kraken2_organism}" | sed 's/ /_/g')
    header="$header,%_human_reads,percent_${organism}_reads"
    result="$result,!{percentage_human},!{percentage_cov}"

    header="$header,ivar_num_variants_identified,bcftools_variants_identified"
    result="$result,!{ivar_variants},!{bcftools_variants}"

    header="$header,bedtools_num_failed_amplicons,samtools_num_failed_amplicons"
    result="$result,!{bedtools_num_failed_amplicons},!{samtools_num_failed_amplicons}"

    header="$header,vadr_conclusion"
    result="$result,!{vadr_results}"

    header="$header,num_N,num_degenerage,num_non-ambiguous,num_total"
    result="$result,!{num_N},!{num_degenerate},!{num_ACTG},!{num_total}"

    header="$header,pangolin_version,pangolearn_version,constellations_version,scorpio_version"
    result="$result,!{pangolin_version},!{pangolearn_version},!{constellations_version},!{scorpio_version}"

    header="$header,nextclade_version"
    result="$result,!{nextclade_version}"

    header="$header,cleaner_version,aligner_version,trimmer_version,ivar_version"
    result="$result,!{cleaner_version},!{aligner_version},!{trimmer_version},!{ivar_version}"

    echo $header > !{task.process}/!{sample}.summary.csv
    echo $result >> !{task.process}/!{sample}.summary.csv

    cat !{task.process}/!{sample}.summary.csv | tr ',' '\t' > !{task.process}/!{sample}.summary.txt
  '''
}

summary
  .collectFile(name: "summary.csv",
      keepHeader: true,
      sort: true,
      skip: 1,
      storeDir: "${params.outdir}")

summary2
  .collectFile(name: "cecret_run_results.txt",
    keepHeader: true,
    sort: true,
    skip: 1,
    storeDir: "${workflow.launchDir}")

if (params.relatedness) {
  if ( params.msa == 'mafft' ) {
    params.mafft_options = '--maxambiguous 0.5'
    process mafft {
      publishDir "${params.outdir}", mode: 'copy'
      tag "Multiple Sequence Alignment"
      cpus params.maxcpus
      container 'staphb/mafft:latest'
      errorStrategy 'retry'
      maxRetries 2

      input:
      file(consensus) from consensus_msa.concat(fastas_msa).collectFile(name:"ultimate.fasta")
      file(reference_genome) from reference_genome_msa

      output:
      file("${task.process}/mafft_aligned.fasta") into msa_file, msa_file2
      file("logs/${task.process}/mafft.${workflow.sessionId}.{log,err}")

      shell:
      '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/mafft.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/mafft.!{workflow.sessionId}.err

        date | tee -a $log_file $err_file > /dev/null
        echo "mafft version:" >> $log_file
        mafft --version 2>&1 >> $log_file

        mafft --auto \
          !{params.mafft_options} \
          --thread !{task.cpus} \
          --addfragments !{consensus} \
          !{reference_genome} \
          > !{task.process}/mafft_aligned.fasta \
          2>> $err_file
      '''
    }
  } else if ( params.msa == 'nextalign' ) {
    params.nextalign_options = '--genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --include-reference'
    process nextalign {
      publishDir "${params.outdir}", mode: 'copy'
      tag "Multiple Sequence Alignment"
      cpus params.maxcpus
      container 'nextstrain/nextalign:latest'

      input:
      file(consensus) from consensus_msa.concat(fastas_msa).collectFile(name:"ultimate.fasta")
      path(dataset) from prepped_nextalign

      output:
      file("${task.process}/nextalign.aligned.fasta") into msa_file, msa_file2
      file("${task.process}/{*.fasta,nextalign.*.csv}")
      file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

      shell:
      '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

        date | tee -a $log_file $err_file > /dev/null
        echo "nextalign version:" >> $log_file
        nextalign --version-detailed 2>&1 >> $log_file

        nextalign !{params.nextalign_options} \
          --sequences !{consensus} \
          --reference !{dataset}/reference.fasta \
          --genemap !{dataset}/genemap.gff \
          --jobs !{task.cpus} \
          --output-dir !{task.process} \
          --output-basename nextalign \
          >> $log_file 2>> $err_file
      '''
    }
  }

  params.snpdists_options = ''
  process snpdists {
    publishDir "${params.outdir}", mode: 'copy'
    tag "createing snp matrix with snp-dists"
    cpus 1
    container 'staphb/snp-dists:latest'

    when:
    params.snpdists

    input:
    file(msa) from msa_file

    output:
    file("snp-dists/snp-dists.txt")
    file("logs/${task.process}/snp-dists.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p snp-dists logs/!{task.process}
      log_file=logs/!{task.process}/snp-dists.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/snp-dists.!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      snp-dists -v >> $log_file

      snp-dists !{params.snpdists_options} !{msa} > snp-dists/snp-dists.txt 2> $err_file
    '''
  }
  params.iqtree2_outgroup = 'MN908947'
  params.iqtree2_options = '-ninit 2 -n 2 -me 0.05 -m GTR'
  process iqtree2 {
    publishDir "${params.outdir}", mode: 'copy'
    tag "Creating phylogenetic tree with iqtree"
    cpus params.maxcpus
    container 'staphb/iqtree2:latest'

    when:
    params.iqtree2

    input:
    file(msa) from msa_file2

    output:
    file("${task.process}/${task.process}.{iqtree,treefile,mldist,log}")
    file("logs/${task.process}/${task.process}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{task.process}.!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null
      iqtree2 --version >> $log_file

      if [ -n "!{params.iqtree2_outgroup}" ] && [ "!{params.iqtree2_outgroup}" != "null" ]
      then
        outgroup="-o !{params.iqtree2_outgroup}"
        cat !{msa} | sed 's/!{params.iqtree2_outgroup}.*/!{params.iqtree2_outgroup}/g' > !{msa}.renamed
      else
        outgroup=""
        mv !{msa} !{msa}.renamed
      fi

      # creating a tree
    	iqtree2 !{params.iqtree2_options} \
        -nt AUTO \
        -ntmax !{task.cpus} \
        -s !{msa}.renamed \
        -pre !{task.process}/iqtree2 \
        $outgroup \
        >> $log_file 2>> $err_file
    '''
  }
}

if (params.rename) {
  Channel
    .fromPath(params.sample_file, type:'file')
    .ifEmpty{
      println("No sample file was found. Set with 'params.sample_file'")
      exit 1
    }
    .view { "Sample File : $it" }
    .set{ sample_file }

  fastq_reads_rename
    .join(consensus_rename, by:0)
    .join(fastas_rename, by:0)
    .join(filtered_reads, remainder: true, by: 0)
    .combine(sample_file)
    .set { rename }

  process rename {
    publishDir "${params.outdir}", mode: 'copy'
    tag "Renaming files for ${sample}"
    cpus 1
    container 'staphb/parallel-perl:latest'
    stageInMode 'copy'

    input:
    tuple val(sample), file(reads), val(paired_single), file(consensus), file(filtered_reads), file(sample_file) from rename

    output:
    file("submission_files/*{genbank,gisaid}.fa") optional true
    file("submission_files/*.fastq.gz")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{err,log}")

    shell:
    '''
      mkdir -p submission_files logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      date | tee -a $log_file $err_file > /dev/null

      !{workflow.projectDir}/bin/genbank_submission.sh \
        -f !{sample_file} \
        -c . \
        -d . \
        -s !{params.gisaid_threshold} \
        -g !{params.genbank_threshold} \
        -o submission_files \
        2>> $err_file >> $log_file
    '''
  }
}

workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("A summary of results can be found in a tab-delimited file: ${workflow.launchDir}/run_results.txt")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
