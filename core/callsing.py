#!/usr/bin/env python3

from spython.main import Client
import os
import argparse
import sys
import json
import shlex
#TODO add inidcator that container is being downloaded or updated

def shutdown():
    print('\nShutting down the running singularity containers and exiting...')
    #we arn't actually going to do anything

def call(container,command,cwd='',paths={},remove=True):
    ###load container
    container = 'docker://'+container
    Client.load(container)

    ###setup mount point paths
    #{"/path/outside":"/path/incontainer"}
    volumes = []
    if paths:
        for key in paths.keys():
            volumes.append(key+':'+paths[key])

    ###setup command as a list
    command_list = command.split()

    ###run the container
    output = Client.execute(command_list,bind=volumes,writable=True,options=['--pwd',cwd])
    #once container is finished return output as a string
    return output

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="callsing.py <docker_container> <command> [options]")
    parser.add_argument("docker_container", type=str, help="Name of the docker container to run using singularity")
    parser.add_argument("command", type=str, help="Command to execute in the singularity container")
    parser.add_argument("-d", default="", type=str, help="Path to working directory in container")
    parser.add_argument("-p", default="", type=str, help='Dictonary of paths to mount in the container e.x. \'{"/path/outside":"/path/incontainer"}\' note: the quoting')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    container_name = args.docker_container
    cmd = args.command
    cwd = args.d
    paths = ''

    if not cwd:
        cwd = os.getcwd()

    #check for path argument if empty set a default
    if args.p:
        try:
            paths = json.loads(args.p)
        except:
            print("Path incorrectly formatted.")
            sys.exit()
    else:
        paths = {cwd:"/data"}

    #run the container and display the output
    print(call(container_name,cmd,cwd=cwd,paths=paths))
