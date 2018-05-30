#!/usr/bin/env python3

import os
import glob
import subprocess


class BSProject(object):

    def __init__(self, name):
        self.name = name  # name of BS project
        self.path = os.path.expanduser("~") + "/BaseSpace/Projects/" + name  # path to BS project
        self.reads = glob.glob(self.path + "/Samples/*/Files/*.fastq.gz")  # Read files for all samples in BS Project

        if not os.path.isdir(self.path):
            raise ValueError("BaseSpace project " + self.name + " not found. \nPlease check project id "
                                                                "and make sure that"
                                                                "BaseSpace is properly mounted")
        return

    def link_reads(self, output_dir=None, samples=None):
        if output_dir is None:
            raw_reads_dir = os.getcwd() + "/" + self.name + "/raw_reads/"
        else:
            raw_reads_dir = os.getcwd() + "/" + output_dir + "/raw_reads/"

        if not os.path.isdir(raw_reads_dir):
            os.makedirs(raw_reads_dir)
            print("Directory for raw reads made: ", raw_reads_dir)

        if samples is None:
            reads = self.reads
        else:
            samples = samples.split(", ")
            reads = list()
            for sample in samples:
                reads.extend(glob.glob(self.path + "/Samples/" + sample + "/Files/*.fastq.gz"))

        for read in reads:
            dest = raw_reads_dir + os.path.basename(read).replace("_001", "")
            if not os.path.isfile(dest):
                os.symlink(read, dest)
                print("Sym link for", read, "made at", dest)
            else:
                print("Sym link for", read, "already exists at", dest)
        return

    def post_to_basespace(self, report, report_dir):
        app_results = self.path + "/AppResults/"
        subprocess.call("cp -R" + " " + report_dir + " " + app_results, shell=True)

        os.chdir(app_results + report)
        subprocess.call("basemount-cmd mark-as-complete", shell=True)

        print(report_dir + " " + "uploaded to BaseSpace.")
