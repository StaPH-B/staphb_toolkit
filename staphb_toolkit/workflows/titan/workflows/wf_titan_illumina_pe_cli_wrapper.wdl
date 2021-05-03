version 1.0

import "wf_titan_illumina_pe.wdl" as titan_illumina_pe
import "../tasks/task_titan_illumina_pe_output_summary.wdl" as summary

struct InputJSON {
  File read1_raw
  File read2_raw
  String samplename
  File primer_bed
}

workflow cli_wrapper {
  input {
    Array[InputJSON] inputSamples
  }

  scatter (sample in inputSamples){
    call titan_illumina_pe.titan_illumina_pe{
      input:
        samplename = sample.samplename,
        seq_method = "Illumina paired-end",
        read1_raw = sample.read1_raw,
        read2_raw = sample.read2_raw,
        primer_bed = sample.primer_bed,
        pangolin_docker_image = "staphb/pangolin:2.3.2-pangolearn-2021-02-21"
    }

    call summary.sample_metrics {
      input:
        samplename        = sample.samplename,
        pangolin_lineage  = titan_illumina_pe.pangolin_lineage,
        pangolin_aLRT     = titan_illumina_pe.pangolin_aLRT,
        nextclade_clade   = titan_illumina_pe.nextclade_clade,
        nextclade_aa_subs = titan_illumina_pe.nextclade_aa_subs,
        nextclade_aa_dels = titan_illumina_pe.nextclade_aa_dels,
        fastqc_raw_pairs  = titan_illumina_pe.fastqc_raw_pairs,
        seqy_pairs        = titan_illumina_pe.seqy_pairs,
        seqy_percent      = titan_illumina_pe.seqy_percent,
        kraken_human      = titan_illumina_pe.kraken_human,
        kraken_sc2        = titan_illumina_pe.kraken_sc2,
        number_N          = titan_illumina_pe.number_N,
        number_ATCG       = titan_illumina_pe.number_ATCG,
        number_Degenerate = titan_illumina_pe.number_Degenerate,
        number_Total      = titan_illumina_pe.number_Total,
        coverage          = titan_illumina_pe.coverage,
        depth             = titan_illumina_pe.depth,
        meanbaseq_trim    = titan_illumina_pe.meanbaseq_trim,
        meanmapq_trim     = titan_illumina_pe.meanmapq_trim,
        coverage_trim     = titan_illumina_pe.coverage_trim,
        depth_trim        = titan_illumina_pe.depth_trim,
        amp_fail          = titan_illumina_pe.amp_fail
    }
  }

  call summary.merge_metrics {
    input:
      all_metrics = sample_metrics.single_metrics
  }

  output {
    Array[File]    read1_clean          = titan_illumina_pe.read1_clean
    Array[File]    read2_clean          = titan_illumina_pe.read2_clean
    Array[File]    kraken_report        = titan_illumina_pe.kraken_report
    Array[File]    sorted_bam           = titan_illumina_pe.sorted_bam
    Array[File]    sorted_bai           = titan_illumina_pe.sorted_bai
    Array[File]    trim_sorted_bam      = titan_illumina_pe.trim_sorted_bam
    Array[File]    trim_sorted_bai      = titan_illumina_pe.trim_sorted_bai
    Array[File]    consensus_seq        = titan_illumina_pe.consensus_seq
    Array[File]    samtools_stats       = titan_illumina_pe.consensus_stats
    Array[File]    cov_hist             = titan_illumina_pe.cov_hist
    Array[File]    cov_stats            = titan_illumina_pe.cov_stats
    Array[File]    samtools_flagstat    = titan_illumina_pe.consensus_flagstat
    Array[File]    pango_lineage_report = titan_illumina_pe.pango_lineage_report
    Array[File]    amp_coverage         = titan_illumina_pe.amp_coverage
    File           merged_metrics       = merge_metrics.run_results
  }

}
