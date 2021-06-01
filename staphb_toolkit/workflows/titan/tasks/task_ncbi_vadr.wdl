version 1.0

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC. See https://github.com/ncbi/vadr/wiki/Coronavirus-annotation"
  }
  input {
    File      genome_fasta
    String    samplename
    String    vadr_opts="--glsearch -s -r --nomisc --mkey sarscov2 --alt_fail lowscore,fstukcnf,insertnn,deletinn --mdir /opt/vadr/vadr-models/"
    String    docker="staphb/vadr:1.2"
    Int?      cpus = 1
    String?   memory = "2 GB"
    Int       minlen=50
    Int       maxlen=30000
  }
  String out_base = basename(genome_fasta, '.fasta')
  command <<<
    set -e

    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      "~{genome_fasta}" \
      --minlen ~{minlen} \
      --maxlen ~{maxlen} \
      > "~{out_base}_trimmed.fasta"

    # run VADR
    v-annotate.pl \
      ~{vadr_opts} \
      "~{out_base}_trimmed.fasta" \
      "~{out_base}"


    # package everything for output
    tar -C "~{out_base}" -czvf "~{out_base}.vadr.tar.gz" .

    # prep alerts into a tsv file for parsing
    cat "~{out_base}/~{out_base}.vadr.alt.list" | cut -f 2 | tail -n +2 > "~{out_base}.vadr.alerts.tsv"
    cat "~{out_base}.vadr.alerts.tsv" | wc -l > NUM_ALERTS

    read -r num < NUM_ALERTS
    if [[ "$num" -lt 1 ]]; then
      echo true > vadr.result
    else
     echo false > vadr.result
    fi

  >>>
  output {
    File feature_tbl  = "~{out_base}/~{out_base}.vadr.pass.tbl"
    Int  num_alerts = read_int("NUM_ALERTS")
    File alerts_list = "~{out_base}/~{out_base}.vadr.alt.list"
    Array[Array[String]] alerts = read_tsv("~{out_base}.vadr.alerts.tsv")
    File outputs_tgz = "~{out_base}.vadr.tar.gz"
    Boolean vadr_result = read_boolean("vadr.result")
    String vadr_docker = docker
  }
  runtime {
    docker:           "~{docker}"
    memory:           "~{memory}"
    cpu:              cpus
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}
