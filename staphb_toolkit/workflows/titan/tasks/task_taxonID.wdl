version 1.0

task kraken2 {
  input {
  	File        read1
	  File? 		  read2
	  String      samplename
	  String?     kraken2_db = "/kraken2-db"
    Int?        cpus=4
    String?     memory = "8 GB"
    String      docker="staphb/kraken2:2.0.8-beta_hv"
  }

  command
  <<<
    # date and version control
    date | tee DATE

    kraken2 --version | head -n1 | tee VERSION

    if ! [ -z ~{read2} ]; then
      mode="--paired"
    fi

    kraken2 $mode \
      --classified-out cseqs#.fq \
      --threads ~{cpus} \
      --db ~{kraken2_db} \
      ~{read1} ~{read2} \
      --report ~{samplename}_kraken2_report.txt

    percentage_human=$(grep "Homo sapiens" ~{samplename}_kraken2_report.txt | cut -f 1)
    percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" ~{samplename}_kraken2_report.txt | cut -f1 )

    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_sc2" ] ; then percentage_sc2="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN
    echo $percentage_sc2 | tee PERCENT_SC2
  >>>

  output {
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
    File 	     kraken_report = "${samplename}_kraken2_report.txt"
    Float 	   percent_human = read_string("PERCENT_HUMAN")
    Float 	   percent_sc2   = read_string("PERCENT_SC2")
  }

  runtime {
    docker:       "~{docker}"
    memory:       "~{memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task pangolin {
  input {
    File        fasta
    String      samplename
    String      docker="staphb/pangolin:1.1.14"
    Int?        cpus = 4
    String?     memory = "8 GB"
  }

  command
  <<<
    # date and version control
    date | tee DATE
    pangolin --version | head -n1 | tee VERSION

    pangolin --outdir ~{samplename} ~{fasta}
    pangolin_lineage=$(tail -n 1 ~{samplename}/lineage_report.csv | cut -f 2 -d "," | grep -v "lineage")

    pangolin_aLRT=$(tail -n 1 ~{samplename}/lineage_report.csv | cut -f 3 -d "," )
    pangolin_stats=$(tail -n 1 ~{samplename}/lineage_report.csv | cut -f 4 -d "," )
    mv ~{samplename}/lineage_report.csv ~{samplename}_pango_lineage.csv

    echo $pangolin_lineage | tee PANGOLIN_LINEAGE
    echo $pangolin_aLRT | tee PANGOLIN_aLRT
    echo $pangolin_stats | tee PANGOLIN_STATS
  >>>

  output {
    String     date                 = read_string("DATE")
    String     version              = read_string("VERSION")
    String     pangolin_lineage     = read_string("PANGOLIN_LINEAGE")
    Float      pangolin_aLRT        = read_string("PANGOLIN_aLRT")
    Float      pangolin_stats       = read_string("PANGOLIN_STATS")
    File       pango_lineage_report = "${samplename}_pango_lineage.csv"
  }

  runtime {
    docker:       "~{docker}"
    memory:       "~{memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task pangolin2 {
  input {
    File        fasta
    String      samplename
    String      docker = "staphb/pangolin:2.3.2-pangolearn-2021-02-21"
    Int?        cpus = 4
    String?     memory = "8 GB"
  }

  command
  <<<
    # date and version control
    date | tee DATE
    echo "$(pangolin -v); $(pangolin -pv)" | tee VERSION
    set -e

    pangolin "~{fasta}" \
       --outfile "~{samplename}.pangolin_report.csv" \
       --verbose

    pangolin_lineage=$(tail -n 1 ~{samplename}.pangolin_report.csv | cut -f 2 -d "," | grep -v "lineage")

    pangolin_probability=$(tail -n 1 ~{samplename}.pangolin_report.csv | cut -f 3 -d "," )
    mv ~{samplename}.pangolin_report.csv ~{samplename}_pango2_lineage.csv

    echo $pangolin_lineage | tee PANGOLIN_LINEAGE
    echo $pangolin_probability | tee PANGOLIN_PROBABILITY
  >>>

  output {
    String     date                 = read_string("DATE")
    String     version              = read_string("VERSION")
    String     pangolin_lineage     = read_string("PANGOLIN_LINEAGE")
    String     pangolin_aLRT        = read_string("PANGOLIN_PROBABILITY")
    String     pangolin_docker      = docker
    File       pango_lineage_report = "${samplename}_pango2_lineage.csv"
  }

  runtime {
    docker:       "~{docker}"
    memory:       "~{memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task nextclade_one_sample {
    meta {
        description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        File    genome_fasta
        File?   root_sequence
        File?   auspice_reference_tree_json
        File?   qc_config_json
        File?   gene_annotations_json
        File?   pcr_primers_csv
        String  docker="neherlab/nextclade:0.14.2"
        Int?    cpus = 2
        String? memory = "3 GB"
    }
    String basename = basename(genome_fasta, ".fasta")
    command
    <<<
        set -e
        nextclade.js --version > VERSION
        nextclade.js \
            --input-fasta "~{genome_fasta}" \
            ~{"--input-root-seq " + root_sequence} \
            ~{"--input-tree " + auspice_reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_json} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json
        cp "~{basename}".nextclade.tsv input.tsv
        python3 <<CODE
        # transpose table
        with open('input.tsv', 'r', encoding='utf-8') as inf:
            with open('transposed.tsv', 'w', encoding='utf-8') as outf:
                for c in zip(*(l.rstrip().split('\t') for l in inf)):
                    outf.write('\t'.join(c)+'\n')
        CODE
        grep ^clade transposed.tsv | cut -f 2 | grep -v clade > NEXTCLADE_CLADE
        grep ^aaSubstitutions transposed.tsv | cut -f 2 | grep -v aaSubstitutions | sed 's/,/|/g' > NEXTCLADE_AASUBS
        grep ^aaDeletions transposed.tsv | cut -f 2 | grep -v aaDeletions | sed 's/,/|/g' > NEXTCLADE_AADELS
    >>>
    runtime {
        docker: "~{docker}"
        memory: "~{memory}"
        cpu:    cpus
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        String nextclade_version  = read_string("VERSION")
        File   nextclade_json     = "~{basename}.nextclade.json"
        File   auspice_json       = "~{basename}.nextclade.auspice.json"
        File   nextclade_tsv      = "~{basename}.nextclade.tsv"
        String nextclade_clade    = read_string("NEXTCLADE_CLADE")
        String nextclade_aa_subs  = read_string("NEXTCLADE_AASUBS")
        String nextclade_aa_dels  = read_string("NEXTCLADE_AADELS")
    }
}
