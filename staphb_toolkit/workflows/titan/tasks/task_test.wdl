version 1.0

task pangolin2 {

input {
  File        fasta
  String      samplename
  Int?        cpus=40
  String      docker="staphb/pangolin:2.1.11-pangolearn-2021-02-05"

}

command{
  # date and version control
  date | tee DATE
  pangolin --version | head -n1 | tee VERSION
  set -e

  pangolin "~{fasta}" \
     --outfile "~{samplename}.pangolin_report.csv" \
     --verbose

  pangolin_lineage=$(tail -n 1 ${samplename}.pangolin_report.csv | cut -f 2 -d "," | grep -v "lineage")

  pangolin_probability=$(tail -n 1 ${samplename}.pangolin_report.csv | cut -f 3 -d "," )
  mv ${samplename}.pangolin_report.csv ${samplename}_pango2_lineage.csv

  echo $pangolin_lineage | tee PANGOLIN_LINEAGE
  echo $pangolin_probability | tee PANGOLIN_PROBABILITY
}

output {
  String     date                 = read_string("DATE")
  String     version              = read_string("VERSION")
  String     pangolin_lineage     = read_string("PANGOLIN_LINEAGE")
  String     pangolin_aLRT        = read_string("PANGOLIN_PROBABILITY")
  File       pango_lineage_report = "${samplename}_pango2_lineage.csv"
}

runtime {
  docker:     "~{docker}"
  memory:       "8 GB"
  cpu:          40
  disks:        "local-disk 100 SSD"
  preemptible:  0
}
}
