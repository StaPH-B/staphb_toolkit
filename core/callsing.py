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
    #we don't need to do anything since the containers arn't detached

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
