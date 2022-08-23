#!/usr/bin/env python3

#authors:
# Kelsey Florek (kelsey.florek@slh.wisc.edu)

import sys,os,re
import argparse
from shutil import copy
from urllib.request import urlopen
from urllib import error as urlerror
from shutil import which
import json
import staphb_toolkit.lib.container_handler as container
from staphb_toolkit.lib.autopath import path_replacer
import staphb_toolkit.lib.updates as updates
import staphb_toolkit.lib.callnextflow as callnxf
from datetime import date
from pyfiglet import Figlet
from rich.console import Console
from rich.table import Table
from rich import print

#global app dir
appDir = os.path.abspath(os.path.dirname(__file__))

def main():
    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            self.print_help()
            sys.stderr.write('\nerror: %s\n' % message)
            sys.exit(1)
        def print_help(self, file=None):
            if file is None:
                pass
            f = Figlet(font='starwars')
            print("[bold cyan]"+f.renderText('StaPH-B ToolKit')+"[/bold cyan]")
            print('StaPH-B ToolKit')
            print(f"Version: {updates.tk_version}")
            print(f"Documentation: https://staphb.org/staphb_toolkit/")
            self._print_message(self.format_help(), file)

    #get app list and metadata
    app_list_url = "https://raw.githubusercontent.com/StaPH-B/docker-builds/master/staphb_toolkit_apps.json"
    try:
        response = urlopen(app_list_url)
    except (urlerror.HTTPError, urlerror.URLError):
        print("[bold red]Error: Cannot connect to GitHub to get app inventory.[/bold red]")
        sys.exit(1)
    app_data = json.loads(response.read())

    #get workflow list and metadata
    workflow_list_url = "https://raw.githubusercontent.com/StaPH-B/staphb_toolkit/main/workflows.json"
    try:
        response = urlopen(workflow_list_url)
        workflow_data = json.loads(response.read())
    except (urlerror.HTTPError, urlerror.URLError):
        print("[bold yellow]WARNING: Cannot connect to GitHub to get workflow inventory. Using local list instead.[/bold yellow]")
        with open(os.path.abspath( os.path.join(os.path.dirname(__file__), '..', 'workflows.json') ),'r') as workflow_data_path:
            workflow_data = json.load(workflow_data_path)

    #construct top level help menu
    parser = MyParser(description=f"",usage="staphb-tk [optional arguments] <application/workflow> [application/workflow arguments]",add_help=True)
    subparser = parser.add_subparsers(title='application or workflow name', metavar='<application/workflow>', dest="app_name", parser_class=MyParser)
    parser.add_argument("-l","--list_tools", default=False, action="store_true", help="List all tools in the toolkit.")
    parser.add_argument("-w","--list_workflows", default=False, action="store_true", help="List all workflows in the toolkit.")
    parser.add_argument("-wv","--workflow_version", default="latest", metavar="<version>", help="Version of tool or workflow to run. Default: latest")
    parser.add_argument("-c","--configuration", metavar="<config_file>", help="Specify a custom workflow configuration file.")
    parser.add_argument("-gc","--get_configuration", default=False, action="store_true", help="Get the configuration file for the specified workflow. Note: You may need to specify a version for the workflow using -wv to get the correct configuration file.")
    parser.add_argument("-nv","--nextflow_version", nargs='?', const="get", metavar="<version>", help="Get or set the version of nextflow.")
    parser.add_argument("--update", default=False, action="store_true", help="Check for and install a ToolKit update.")
    parser.add_argument("--auto_update", default=False, action="store_true", help="Toggle automatic ToolKit updates. Default is off.")

    #construct subparsers
    parser_data = {}
    for app in app_data['apps']:
        parser_data[app] = subparser.add_parser(app,add_help=False)
    for workflow in workflow_data['workflows']:
        parser_data[workflow] = subparser.add_parser(workflow,add_help=False)

    def print_tool_list():
        table = Table(title="Available programs:",title_justify="left",title_style="b",padding=(1,0,0,0))

        table.add_column("Command", justify="left",no_wrap=True)
        table.add_column("Description", justify="left")

        for app in app_data['apps']:
            table.add_row(app,app_data['apps'][app]['description'])

        console = Console()
        console.print(table)
        return

    def print_workflow_list():
        table = Table(title="Available workflows:",title_justify="left",title_style="b",padding=(1,0,0,0))

        table.add_column("Command", justify="left",no_wrap=True)
        table.add_column("Description", justify="left")

        for workflow in workflow_data['workflows']:
            table.add_row(workflow,workflow_data['workflows'][workflow]['description'])

        console = Console()
        console.print(table)
        return

    #handle the arguments and perform automatic path replacement for apps
    parser_args = parser.parse_known_args()
    application = parser_args[0].app_name
    args = parser_args[1]

    #check if nextflow is installed, if not install and set to latest
    if not which('nextflow') and not os.path.exists(os.path.join(appDir,'bin','nextflow')):
        updates.install_nextflow()
        updates.set_nf_version(version="latest")

    #set nextflow version
    if parser_args[0].nextflow_version == "get":
        updates.get_nf_version()
        sys.exit(0)
    elif parser_args[0].nextflow_version == "latest":
        updates.set_nf_version()
        sys.exit(0)
    elif parser_args[0].nextflow_version:
        updates.set_nf_version(parser_args[0].nextflow_version)
        sys.exit(0)
    else:
        pass

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
    elif toolName in workflow_data['workflows'].keys():
        sb_workflow = True
    else:
        print(f"[bold red]Error: The workflow or tool \"{toolName}\" is not in the StaPH-B Toolkit.[/bold red]")
        sys.exit(1)

    #Run the program
    #-----------------------------------------

    #---------------------
    #app
    if sb_tool:
        #run autopathing to mount correct file locations
        arg_string,path_map = path_replacer(args,os.getcwd())

        #see if we have an executable, also check if we need to display help
        try:
            e = app_data['apps'][application]['exec']
            if not arg_string:
                arg_string = arg_string + " " + app_data['apps'][application]['help']
        except KeyError:
            e = ''
            if not arg_string:
                print(app_data['apps'][application]['help'])
                sys.exit(0)

        #build command for running application
        command = e + " " + arg_string.lstrip()
        #get proper image and tag
        image = app_data['apps'][application]['image']
        tag = parser_args[0].workflow_version

        program_object = container.Run(command=command, path=path_map, image=image, tag=tag)
        program_object.run()

    #---------------------
    #workflow
    if sb_workflow:
        #get workflow information
        wf_version = parser_args[0].workflow_version
        if wf_version == 'latest':
            wf_version = workflow_data['workflows'][application]['default_branch']
        configFileName = workflow_data['workflows'][application]['configuration']
        repo = workflow_data['workflows'][application]['repo']

        #get workflow help if asked
        if (any(x in ['-h','--help','-help'] for x in parser_args[1]) or not parser_args[1]) and not parser_args[0].get_configuration:
            callnxf.get_workflow_help(repo,application,wf_version,workflow_data['workflows'][application]['schema'])
            sys.exit()

        #set configuration file
        config_file = None
        if parser_args[0].configuration:
            config_file = parser_args[0].configuration

        #get configuration file from repo if asked
        if parser_args[0].get_configuration:
            callnxf.get_configuration_file(repo,application,wf_version,configFileName)
            sys.exit(0)

        workflow_args = " ".join(parser_args[1])

        #run nextflow workflow
        callnxf.run(wf_version,repo,workflow_args,config_file)
