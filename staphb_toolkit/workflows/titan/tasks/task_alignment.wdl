version 1.0

task bwa {

  input {
    File        read1
    File        read2
    String      samplename
    String?     reference_genome="/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"
    Int?        cpus = 2
    String?     memory = "8 GB"
    String      docker="staphb/ivar:1.2.2_artic20200528"
  }

  command {
    # date and version control
    date | tee DATE
    echo "BWA $(bwa 2>&1 | grep Version )" | tee BWA_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # Map with BWA MEM
    bwa mem \
    -t ${cpus} \
    ${reference_genome} \
    ${read1} ${read2} |\
      samtools sort | samtools view -F 4 -o ${samplename}.sorted.bam

    # index BAMs
    samtools index ${samplename}.sorted.bam
  }

  output {
    String     bwa_version = read_string("BWA_VERSION")
    String     sam_version = read_string("SAMTOOLS_VERSION")
    File       sorted_bam = "${samplename}.sorted.bam"
    File       sorted_bai = "${samplename}.sorted.bam.bai"
  }

  runtime {
    docker:       "~{docker}"
    memory:       "~{memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task mafft {

  input {
    Array[File]   genomes
    Int?          cpus = 16
    String?       memory = "32 GB"
    String        docker="staphb/mafft:7.450"
  }

  command{
    # date and version control
    date | tee DATE
    mafft_vers=$(mafft --version)
    echo Mafft $(mafft_vers) | tee VERSION

    cat ${sep=" " genomes} | sed 's/Consensus_//;s/.consensus_threshold.*//' > assemblies.fasta
    mafft --thread -${cpus} assemblies.fasta > msa.fasta
  }

  output {
    String        date = read_string("DATE")
    String        version = read_string("VERSION")
    File          msa = "msa.fasta"
  }

  runtime {
    docker:       "~{docker}"
    memory:       "~{memory}"
    cpu:          cpus
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}
