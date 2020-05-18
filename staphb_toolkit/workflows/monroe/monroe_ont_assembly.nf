#!/usr/bin/env nextflow

//Description: Adaptation of ARTIC Network nCoV-2019 Bioinformatics SOP
//Available: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Authors of this Nextflow: Kelsey Florek and Abigail Shockey
//Email: kelsey.florek@slh.wisc.edu

//starting parameters
params.fast5_dir = ""
params.fastq_dir = ""
params.sequencing_summary = ""
params.outdir = ""
params.primers = ""
params.run_prefix = "artic_ncov19"
params.pipe = ""

// Input channels
Channel
    .value( "${params.primers}")
    .ifEmpty { exit 1, "Primers used must be included." }
    .set { polish_primers }

// If we have fast5 files then start with basecalling
if(params.basecalling){
  Channel
      .fromPath( "${params.fast5_dir}")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .into { raw_fast5; polish_fast5 }

  process guppy_basecalling {
    input:
      file(fast5s) from raw_fast5.collect()

    output:
      file "fastq/*.fastq" into fastq_reads
      file "fastq/sequencing_summary.txt" into sequencing_summary

    script:
      if(params.basecalling_mode == "fast"){
        """
        guppy_basecaller --chunk_size ${params.chunk_size} --chunks_per_runner ${params.chunks_per_runner} --gpu_runners_per_device ${params.gpu_runners_per_device} -c /opt/ont/guppy/data/dna_r9.4.1_450bps_fast.cfg -i . -s fastq -x auto -r
        """
      }else{
        """
        guppy_basecaller --chunk_size ${params.chunk_size} --chunks_per_runner ${params.chunks_per_runner} --gpu_runners_per_device ${params.gpu_runners_per_device} -c /opt/ont/guppy/data/dna_r9.4.1_450bps_hac.cfg -i . -s fastq -x auto -r
        """
      }
  }
}

// If we have already basecalled get fastqs and fast5s for polishing
else {
  Channel
      .fromPath( "${params.fastq_dir}/*.fastq*")
      .ifEmpty { exit 1, "Cannot find any fastq files in: ${params.fastq_dir} Path must not end with /" }
      .set { fastq_reads }

  Channel
      .fromPath( "${params.fast5_dir}")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .set { polish_fast5 }

if(params.polishing == "nanopolish"){
  Channel
      .fromPath( "${params.sequencing_summary}")
      .ifEmpty { exit 1, "Cannot find sequencing summary in: ${params.sequencing_summary}" }
      .set { sequencing_summary }
  }
}

// Demultiplex fastqs
process guppy_demultiplexing {
  publishDir "${params.outdir}/demultiplexing", mode: 'copy'

  input:
    file(fastqs) from fastq_reads.collect()

  output:
    path("output_directory/barcode*",type:'dir',maxDepth:1) into demultiplexed_reads

  script:
    """
      guppy_barcoder --require_barcodes_both_ends -i . -s output_directory --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
    """
}

// Run artic gupplyplex
process artic_guppyplex {
  publishDir "${params.outdir}/guppyplex", mode: 'copy'

  //guppy plex likes to occasionally fail when running in the cloud, set to retry here once
  errorStrategy { task.attempt > 2 ? 'ignore' : 'retry' }
  maxRetries 1

  input:
    file(reads) from demultiplexed_reads.flatten()

  output:
    file "${params.run_prefix}_barcode*.fastq" into polish_files

  script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --directory barcode* \
    --prefix ${params.run_prefix}
    """
}

// Run artic pipeline using nanopolish
if(params.polishing=="nanopolish"){
  process artic_nanopolish_pipeline {
    publishDir "${params.outdir}/pipeline_nanopolish", mode: 'copy'

    input:
      val primers from polish_primers
      tuple file(fastq), path(fast5path), file(sequencing_summary) from polish_files .combine(polish_fast5) .combine(sequencing_summary)

    output:
      file "*{.primertrimmed.bam,.vcf,.variants.tab}"
      file "*.consensus.fasta" into consensus_fasta

    script:
      """
      # get samplename by dropping file extension
      filename=${fastq}
      samplename=\${filename%.*}

      cpus=`grep -c ^processor /proc/cpuinfo`

      artic minion --normalise ${params.normalise} --threads \$cpus --scheme-directory /artic-ncov2019/primer_schemes --fast5-directory ${fast5path}  --sequencing-summary ${sequencing_summary} --read-file ${fastq} nCoV-2019/${primers} \$samplename
      """
  }
}

// Run artic pipeline using medaka
else {
  process artic_medaka_pipeline {
    publishDir "${params.outdir}/pipeline_medaka", mode: 'copy'

    input:
      val primers from polish_primers
      file(fastq) from polish_files

    output:
      file "*{.primertrimmed.rg,.primers.vcf,.vcf.gz,.trimmed.rg,.fail.vcf}*"
      file "*.consensus.fasta" into consensus_fasta

    script:
      """
      # get samplename by dropping file extension
      filename=${fastq}
      samplename=\${filename%.*}
      cpus=`grep -c ^processor /proc/cpuinfo`

      artic minion --medaka --normalise ${params.normalise} --threads \$cpus --scheme-directory /artic-ncov2019/primer_schemes --read-file ${fastq} nCoV-2019/${primers} \$samplename
      """
  }
}
