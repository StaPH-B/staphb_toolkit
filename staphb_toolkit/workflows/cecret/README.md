# CECRET

A workflow for getting consensus sequences out from paired-end fastq.gz or fastq reads from amplicon prepared libraries.

Required parameters :

1) bed file for primer sequences
2) reference file
3) gff file

* --reference_genome [MN908947.3.fasta](https://raw.githubusercontent.com/UPHL-BioNGS/Cecret/master/config/MN908947.3.fasta)
* --gff_file [MN908947.3.gff](https://raw.githubusercontent.com/UPHL-BioNGS/Cecret/master/config/MN908947.3.gff)
* --primer_bed [artic_V3_nCoV-2019.bed](https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.bed)

Default parameters worth noting are where the original paired-end fastq files are originally looked for and where the results will appear.
* params.reads = workflow.launchDir + '/Sequencing_reads/Raw'
* params.outdir = workflow.launchDir + "/cecret"


A full list of parameters and documentation can be found cecret's standalone repository : [https://github.com/UPHL-BioNGS/Cecret](https://github.com/UPHL-BioNGS/Cecret)
