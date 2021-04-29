version 1.0

import "wf_titan_illumina_pe_cli.wdl" as titan_illumina_pe

workflow cli_wrapper {
  input {
    Array[Pair[Array[String], Pair[File,File]]] inputSamples
    Array[Array[String]] inputConfig
  }

  scatter (sample in inputSamples){
    call titan_illumina_pe.titan_illumina_pe{
      input:
        samplename
        seq_method = "Illumina paired-end"
        read1_raw = 
        read2_raw =
        primer_bed =
        pangolin_docker_image = "staphb/pangolin:2.3.2-pangolearn-2021-02-21"
    }
  }

}
