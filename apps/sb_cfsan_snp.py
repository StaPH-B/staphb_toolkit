#!/usr/bin/env python3

#author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import os.path
import sys
import json
import pathlib
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

from core import fileparser
from core import sb_programs


class CFSAN_SNP():
    def __init__(self, runfiles=None, path=None, output_dir = None, configuration=None, reference=None):

        # set configuration file variables
        self.configuration = configuration
        if self.configuration:
            config_file_path = os.path.abspath(configuration)
        else:
            config_file_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4] + "/core/docker_config.json"

        with open(config_file_path, 'r') as config_file:
            self.config = json.load(config_file)

        # create output dir if it doesn't already exist
        if output_dir:
            self.output_dir = os.path.abspath(output_dir)
        else:
            self.output_dir = os.getcwd()

        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        # set reference variables from user input
        if not reference:
            raise ValueError("Must specify reference genome.")
        self.reference = reference
        self.reference_dir = pathlib.Path(self.reference).parent.absolute()
        self.reference_file = os.path.basename(self.reference)

        # set runfiles variable using fileparser if runfiles not provided
        if runfiles:
            self.runfiles = runfiles
        else:
            self.path = path
            self.runfiles = fileparser.ProcessFastqs(self.path, output_dir=output_dir)

        # set the read_file_path equal to the output_dir since reads have been copied/hard linked there
        if os.path.isdir(os.path.join(self.path, "AppResults")):
            self.path = output_dir
        else:

            self.runfiles.link_reads(output_dir=output_dir)
            self.path = output_dir

    def run(self):
        if not os.path.isfile(os.path.join(*[self.output_dir, "cfsan_output", "snpma.fasta"])):

            # create cfsansnp_output directory
            # cfsan_snp_out_dir = os.path.join(self.output_dir, "cfsan_snp_output")
            # if not os.path.isdir(cfsan_snp_out_dir):
            #     os.makedirs(cfsan_snp_out_dir)
            #     print("Directory for cfsansnp output made: ", cfsan_snp_out_dir)

            # mount info for datain and dataout dirs
            cfsan_snp_mounting = {self.output_dir: '/dataout', self.reference_dir: '/reference'}

            # command for creating the cfsan-snp
            cfsan_snp_configuration = self.config["parameters"]["cfsan-snp-pipeline"]
            cfsan_snp_parameters = cfsan_snp_configuration["params"]
            cfsan_snp_command = f"bash -c 'run_snp_pipeline.sh -m {cfsan_snp_parameters['mirrored_input']} -o /dataout/cfsan_snp_output --samples_dir /dataout/input_reads/ /reference/{self.reference_file}'"

            # create cfsan-snp sketch object
            cfsan_snp = sb_programs.Run(command=cfsan_snp_command, path=cfsan_snp_mounting, image=cfsan_snp_configuration["image"], tag=cfsan_snp_configuration["tag"])

            # run cfsan-snp
            cfsan_snp.run()
