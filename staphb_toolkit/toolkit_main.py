#!/usr/bin/env python3

#authors:
# Kelsey Florek (kelsey.florek@slh.wisc.edu)

import sys,os,re
import argparse
from shutil import copy
import json
import staphb_toolkit.core.container_handler as container
from staphb_toolkit.core.autopath import path_replacer
import staphb_toolkit.core.update_app as autoupdate
from datetime import date

def main():
    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            self.print_help()
            sys.stderr.write('\nerror: %s\n' % message)
            sys.exit(1)

    with open(os.path.abspath(os.path.dirname(__file__) + '/apps.json'),'r') as app_data_file:
        app_data = json.load(app_data_file)

    with open(os.path.abspath(os.path.dirname(__file__) + '/workflows.json'),'r') as workflows_data_path:
        workflows_data = json.load(workflows_data_path)

    parser = MyParser(description=f"StaPH-B ToolKit v{autoupdate.version}",usage="staphb-tk [optional arguments] <application/workflow> [application/workflow arguments]",add_help=True)
    parser.add_argument("command", nargs="*", help="application or workflow name")
    parser.add_argument("-ch","--command_help",default=False,action="store_true", help="get usage for the tool or workflow to run.")
    parser.add_argument("-v","--run_version",default="latest", metavar="<version>", help="version of tool or workflow to run. default: latest")
    parser.add_argument("-l","--list_tools",default=False,action="store_true", help="List all tools in the toolkit.")
    parser.add_argument("-w","--list_workflows",default=False,action="store_true", help="List all workflows in the toolkit.")
    parser.add_argument("-nv","--nextflow_version",nargs='?',const="get",metavar="<version>",help="Get or set the version of nextflow.")
    parser.add_argument("--update",default=False,action="store_true",help="Check for and install a ToolKit update.")
    parser.add_argument("--auto_update",default=False,action="store_true",help="Toggle automatic ToolKit updates. Default is off.")

    def print_tool_list():
        print("Available programs:")
        header = ["Command","Description","-------","-----------"]
        print(f"{header[0]:<25}{header[1]:^10}")
        print(f"{header[2]:<25}{header[3]:^10}")
        for key in app_data['apps']:
            print(f"{key:<25}{app_data['apps'][key]['description']:^10}")
        return

    def print_workflow_list():
        print("Available workflows:")
        header = ["Command","Description","-------","-----------"]
        print(f"{header[0]:<25}{header[1]:^10}")
        print(f"{header[2]:<25}{header[3]:^10}")
        for key in workflow_data['workflows']:
            print(f"{key:<25}{workflow_data['workflows'][key]['description']:^10}")
        return

    #handle the arguments and perform automatic path replacement
    parser_args = parser.parse_known_args()
    parser_args[0].command = parser_args[0].command + parser_args[1]
    print(parser_args[0])

    #check for updates
    if parser_args[0].update:
        autoupdate.check_for_updates()
        sys.exit(0)

    if parser_args[0].auto_update:
        #get current status
        update_status = autoupdate.check_update_status()
        if update_status:
            autoupdate.toggle_updater(False)
        else:
            autoupdate.toggle_updater(True)

    if autoupdate.check_update_status():
        autoupdate.check_for_updates()

    #display list of tools or workflows if needed
    if parser_args[0].list_tools:
        print_tool_list()
        sys.exit(0)

    if parser_args[0].list_workflows:
        print_workflow_list()
        sys.exit(0)

    if not parser_args[0].command:
        parser.print_help()
        sys.exit(0)

    #Check if command is a tool or workflow
    sb_tool = False
    sb_workflow = False
    toolName = parser_args[0].command[0]
    if toolName in app_data['apps'].keys():
        sb_tool = True
    elif toolName in workflows_data['workflows'].keys():
        sb_workflow = True
    else:
        print(f"The workflow or tool \"{toolName}\" is not in the StaPH-B Toolkit.")
        sys.exit(1)

    #Check if user needs command usage
    if parser_args[0].command_help:
        if sb_tool:
            print(app_data['apps'][toolName]['help'])
        if sb_workflow:
            print(workflows_data['workflows'][toolName]['help'])
        sys.exit(0)

    #Run autopathing
    arg_string,path_map = path_replacer(parser_args[0].command,os.getcwd())
    print(arg_string)
    sys.exit(0)
    #Run the program
    #-----------------------------------------
    program_object = container.Run(command=command, path=path_map, image=program_configuration["image"], tag=program_configuration["tag"])
    program_object.run()
