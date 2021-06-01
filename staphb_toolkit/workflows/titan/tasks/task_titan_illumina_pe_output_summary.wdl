version 1.0

task sample_metrics {

  input {
    String    samplename
    String    pangolin_lineage
    String    pangolin_version
    String    nextclade_clade
    String    nextclade_version
    String    nextclade_aa_subs
    String    nextclade_aa_dels
    Int?       fastqc_raw_pairs
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
    String    docker = "staphb/multiqc:1.7"
    Int?      cpus = 1
    String?   memory = "1 GB"
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
  ~{fastqc_raw_pairs},\
  ~{depth},\
  ~{depth_trim},\
  ~{coverage},\
  ~{coverage_trim},\
  ~{kraken_human},\
  ~{kraken_sc2},\
  ~{amp_fail},\
  ~{number_N},\
  ~{pangolin_lineage},\
  ~{pangolin_version},\
  ~{nextclade_clade},\
  ~{nextclade_version},\
  ~{number_Degenerate},\
  ~{number_ATCG},\
  ~{number_Total},\
  ~{meanbaseq_trim},\
  ~{meanmapq_trim},\
  ~{nextclade_aa_subs},\
  ~{nextclade_aa_dels}" + assembly_status

  print(outstring)

  CODE
  >>>

  output {
    String  single_metrics = read_string(stdout())
  }

  runtime {
      docker:       "~{docker}"
      memory:       "~{memory}"
      cpu:          cpus
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}


task merge_metrics {

  input {
    Array[String]   all_metrics
    String          docker="staphb/trimmomatic:0.39"
    Int?            cpus = 1
    String?         memory = "1 GB"
  }

  command {
    echo "sample_id,\
    raw_read_pairs,\
    depth_before_trimming,\
    depth_after_trimming,\
    coverage_before_trimming,\
    coverage_after_trimming,\
    %_human_reads,\
    %_SARS-COV-2_reads,\
    num_failed_amplicons,\
    num_N,\
    pangolin_lineage,\
    pangolin_version,\
    nextclade_lineage,\
    nextclade_version,\
    num_degenerate,\
    num_ACTG,\
    num_total,\
    meanbaseq_trim,\
    meanmapq_trim,\
    assembly_status,\
    nextclade_aaSubstitutions,\
    nextclade_aaDeletions" >> run_results.csv

    echo "${sep="END" all_metrics}" >> run_results.csv
    sed -i "s/END/\n/g" run_results.csv
  }

  output {
    File    run_results = "run_results.csv"
  }

  runtime {
      docker:       "~{docker}"
      memory:       "~{memory}"
      cpu:          cpus
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}
