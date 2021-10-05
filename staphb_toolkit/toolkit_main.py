#!/usr/bin/env python3

#authors:
# Kelsey Florek (kelsey.florek@slh.wisc.edu)

import sys,os,re
import argparse
from shutil import copy
import json
import staphb_toolkit.core.container_handler as container
from staphb_toolkit.core.autopath import path_replacer
import staphb_toolkit.core.updates as updates
from datetime import date

def main():
    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            self.print_help()
            sys.stderr.write('\nerror: %s\n' % message)
            sys.exit(1)

    #get app metadata
    with open(os.path.abspath(os.path.dirname(__file__) + '/apps.json'),'r') as app_data_file:
        app_data = json.load(app_data_file)

    #get workflow metadata
    with open(os.path.abspath(os.path.dirname(__file__) + '/workflows.json'),'r') as workflows_data_path:
        workflow_data = json.load(workflows_data_path)

    #construct top level help menu
    parser = MyParser(description=f"StaPH-B ToolKit v{updates.tk_version}",usage="staphb-tk [optional arguments] <application/workflow> [application/workflow arguments]",add_help=True)
    subparser = parser.add_subparsers(title='application or workflow name',metavar='<application/workflow>',dest="app_name",parser_class=MyParser)
    parser.add_argument("-ch","--command_help",default=False,action="store_true", help="get usage for the tool or workflow to run.")
    parser.add_argument("-v","--run_version",default="latest", metavar="<version>", help="version of tool or workflow to run. default: latest")
    parser.add_argument("-l","--list_tools",default=False,action="store_true", help="List all tools in the toolkit.")
    parser.add_argument("-w","--list_workflows",default=False,action="store_true", help="List all workflows in the toolkit.")
    parser.add_argument("-nv","--nextflow_version",nargs='?',const="get",metavar="<version>",help="Get or set the version of nextflow.")
    parser.add_argument("--update",default=False,action="store_true",help="Check for and install a ToolKit update.")
    parser.add_argument("--auto_update",default=False,action="store_true",help="Toggle automatic ToolKit updates. Default is off.")

    #construct subparsers
    parser_data = {}
    for app in app_data['apps']:
        parser_data[app] = subparser.add_parser(app,add_help=False)
    for workflow in workflow_data['workflows']:
        parser_data[workflow] = subparser.add_parser(workflow,add_help=False)

    def print_tool_list():
        print("Available programs:")
        header = ["Command","Description","-------","-----------"]
        print(f"{header[0]:<25}{header[1]:^10}")
        print(f"{header[2]:<25}{header[3]:^10}")
        for app in app_data['apps']:
            print(f"{app:<25}{app_data['apps'][app]['description']:^10}")
        return

    def print_workflow_list():
        print("Available workflows:")
        header = ["Command","Description","-------","-----------"]
        print(f"{header[0]:<25}{header[1]:^10}")
        print(f"{header[2]:<25}{header[3]:^10}")
        for workflow in workflow_data['workflows']:
            print(f"{workflow:<25}{workflow_data['workflows'][workflow]['description']:^10}")
        return

    #handle the arguments and perform automatic path replacement
    parser_args = parser.parse_known_args()
    application = parser_args[0].app_name
    args = parser_args[1]

    #check for updates
    if parser_args[0].update:
        updates.check_for_updates()
        sys.exit(0)

    if parser_args[0].auto_update:
        #get current status
        update_status = updates.check_update_status()
        if update_status:
            updates.toggle_updater(False)
        else:
            updates.toggle_updater(True)

    if updates.check_update_status():
        updates.check_for_updates()

    #display list of tools or workflows if needed
    if parser_args[0].list_tools:
        print_tool_list()
        sys.exit(0)

    if parser_args[0].list_workflows:
        print_workflow_list()
        sys.exit(0)

    if application == None:
        parser.print_help()
        sys.exit(0)

    #Check if command is an app or workflow
    sb_tool = False
    sb_workflow = False
    toolName = application
    if toolName in app_data['apps'].keys():
        sb_tool = True
    elif toolName in workflows_data['workflows'].keys():
        sb_workflow = True
    else:
        print(f"The workflow or tool \"{toolName}\" is not in the StaPH-B Toolkit.")
        sys.exit(1)

    #Run autopathing if we are running an app
    if sb_tool:
        arg_string,path_map = path_replacer(args,os.getcwd())

    #Run the program
    #-----------------------------------------
    #app
    if sb_tool:

        try:
            e = app_data['apps'][application]['exec']
            if args == []:
                args.append(app_data['apps'][application]['help'])
        except KeyError:
            e = ''
            if args == []:
                print(app_data['apps'][application]['help'])
                sys.exit(0)
        command = e + " " + " ".join(args)
        image = app_data['apps'][application]['image']
        tag = parser_args[0].run_version
        program_object = container.Run(command=command, path=path_map, image=image, tag=tag)
        program_object.run()
    #workflow
    if sb_workflow:
        command = workflow_data['workflow'][application]['repo']
        sys.exit(0)
