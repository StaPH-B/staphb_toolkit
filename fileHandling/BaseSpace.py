#!/usr/bin/env python3

import os
import glob
import subprocess
import sys
import argparse
import datetime
from time import gmtime, strftime


class BSProject(object):
    """

    """

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
        #TODO rename sample names using
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
        # TODO change date to date and time
        app_results_dir = self.path + "/AppResults/" + report + strftime("_%Y_%m_%d", gmtime())  + "/"
        subprocess.call("mkdir" + " " + app_results_dir + " " "&& cp -R" + " " + report_dir + "*" + " "
                        + app_results_dir, shell=True)

        os.chdir(app_results_dir)
        subprocess.call("basemount-cmd mark-as-complete", shell=True)

        print(report_dir + " " + "uploaded to BaseSpace.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="BaseSpace.py <input> [options]")
    parser.add_argument("input", type=str, help="Name of BaseSpace project")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir. Default will store output within "
                                                         "an output_dir with the same name as the <input> "
                                                         "BaseSpace Project")
    parser.add_argument("-link", action='store_true', help="Will link reads from <input> BaseSpace Project to "
                                                           "<output_dir>/raw_reads/")
    parser.add_argument("-postToBS", choices=["Belvidere", "Tredegar"], help="Will post the indicated report to the "
                                                                             "<input> BaseSpace Project")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.o
    link_reads = args.link
    reportToBS = args.postToBS

    if len(output_dir) > 0:
        pass
    else:
        output_dir = input_dir

    project = BSProject(input_dir)

    print("Selected BaseSpace Project: ", input_dir)

    if link_reads:
        project.link_reads(output_dir)

    if reportToBS:
        report_dir = output_dir + "/reports/" + reportToBS + "/"

        if not os.path.isdir(report_dir):
            raise ValueError("No " + reportToBS + " report was not found.")
        else:
            print(report_dir)
            print(reportToBS)
            print("Posting", " ", reportToBS, " ", "report to BaseSpace. . .")
            project.post_to_basespace(str(reportToBS), report_dir)


