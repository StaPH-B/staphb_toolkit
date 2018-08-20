#!/usr/bin/env python3

import os
import glob
import subprocess
import sys
import argparse
import re
from time import gmtime, strftime


class BSProject(object):
    """

    """

    def __init__(self, Project):
        self.name = Project  # name of BS project
        self.path = os.path.expanduser("~") + "/BaseSpace/Projects/" + self.name  # path to BS project
        self.reads = glob.glob(self.path + "/Samples/*/Files/*.fastq.gz")  # Read files for all samples in BS Project

        if not os.path.isdir(self.path):
            raise ValueError("BaseSpace project " + self.name + " not found. \nPlease check project id "
                                                                "and make sure that"
                                                                "BaseSpace is properly mounted")
        return

    def link_raw_reads(self, output_dir=None, samples=None):
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
            dest = raw_reads_dir + re.sub('S\d+_L\d+_R', "", os.path.basename(read))
            dest = dest.replace("_001","")
            print(dest)
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
    parser = argparse.ArgumentParser(usage="sb_basespace.py <input> [options]")
    parser.add_argument("input", type=str, help="Name of BaseSpace project")
    parser.add_argument("-o", default="", type=str, help="Name of output_dir. Default will store output within "
                                                         "an output_dir with the same name as the <input> "
                                                         "BaseSpace Project")
    parser.add_argument("-link_raw_reads", action='store_true', help="Will link reads from <input> BaseSpace Project "
                                                                     "to <output_dir>/raw_reads/")
    parser.add_argument("-postToBS", choices=["Belvidere", "Tredegar"], help="Will post the indicated report to the "
                                                                             "<input> BaseSpace Project")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.o
    link_raw_reads = args.link_raw_reads
    postToBS = args.postToBS

    if len(output_dir) > 0:
        pass
    else:
        output_dir = input_dir

    project = BSProject(input_dir)

    print("Selected BaseSpace Project: ", project.name)

    if link_raw_reads:
        print("linking raw reads. . .")
        project.link_raw_reads(output_dir)

    if postToBS:
        report_dir = output_dir + "/reports/" + postToBS + "/"

        if not os.path.isdir(report_dir):
            raise ValueError("No " + posttToBS + " report was not found.")
        else:
            print("Posting", " ", postToBS, " ", "report to BaseSpace. . .")
            project.post_to_basespace(str(postToBS), report_dir)


