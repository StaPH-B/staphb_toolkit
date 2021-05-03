version 1.0

task sample_metrics {

  input {
    String    samplename
    String    pangolin_lineage
    Float     pangolin_aLRT
    String    nextclade_clade
    String    nextclade_aa_subs
    String    nextclade_aa_dels
    Int?       fastqc_raw_pairs
    Int?       seqy_pairs
    Float?     seqy_percent
    Float?     kraken_human
    Float?     kraken_sc2
    Int       number_N
    Int       number_ATCG
    Int       number_Degenerate
    Int       number_Total
    Float     coverage
    Float     depth
    Float     meanbaseq_trim
    Float     meanmapq_trim
    Float     coverage_trim
    Float     depth_trim
    Int       amp_fail
    Float?    coverage_threshold = 95.00
    Float?    meanbaseq_threshold = 30.00
    Float?    meanmapq_threshold = 30.00
  }

  command <<<
  python3<<CODE

  if ~{coverage_trim} >= ~{coverage_threshold} and ~{meanbaseq_trim} >= ~{meanbaseq_threshold} and ~{meanmapq_trim} >= ~{meanmapq_threshold}:
    assembly_status = "PASS"
  else:
    assembly_status = "Warning: "

  if ~{coverage_trim} <= ~{coverage_threshold}:
      assembly_status += "Avg coverage < 95%; "
  if ~{meanbaseq_trim} <= ~{meanbaseq_threshold}:
      assembly_status += "Mean base quality < 30; "
  if ~{meanmapq_trim} <= ~{meanmapq_threshold}:
      assembly_status += "Mean map quality < 30"

  outstring="~{samplename},\
  ~{pangolin_lineage},~{pangolin_aLRT},\
  ~{nextclade_clade},~{nextclade_aa_subs},~{nextclade_aa_dels},\
  ~{fastqc_raw_pairs},~{seqy_pairs},~{seqy_percent},\
  ~{depth},~{depth_trim},~{coverage},~{coverage_trim},\
  ~{kraken_human},~{kraken_sc2},~{amp_fail},\
  ~{number_N},~{number_Degenerate},~{number_ATCG},~{number_Total},\
  ~{meanbaseq_trim},~{meanmapq_trim}," + assembly_status

  print(outstring)

  CODE
  >>>

  output {
    String  single_metrics = read_string(stdout())
  }

  runtime {
      docker:       "staphb/multiqc:1.7"
      memory:       "1 GB"
      cpu:          1
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}


task merge_metrics {

  input {
    Array[String]   all_metrics
  }

  command {
    echo "sample_id,\
    pangolin_lineage,pangolin_probability,\
    nextclade_lineage,nextclade_aaSubstitutions,nextclade_aaDeletions,\
    raw_pairs,pairs_after_cleaning,percent_kept_after_cleaning,\
    depth_before_trimming,depth_after_trimming,coverage_before_trimming,coverage_after_trimming,\
    %_human_reads,%_SARS-COV-2_reads,num_failed_amplicons,\
    num_N,num_degenerate,num_ACTG,num_total,\
    meanbaseq_trim,meanmapq_trim,assembly_status" >> run_results.csv

    echo "${sep="END" all_metrics}" >> run_results.csv
    sed -i "s/END/\n/g" run_results.csv
  }

  output {
    File    run_results = "run_results.csv"
  }

  runtime {
      docker:       "staphb/seqyclean:1.10.09"
      memory:       "1 GB"
      cpu:          1
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}
