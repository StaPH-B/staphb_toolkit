#!/usr/bin/env python3

import docker
import os,sys
import argparse
import json
import shlex
import ast

def shutdown():
    print('\nShutting down the running docker containers and exiting...')
    client = docker.from_env()
    container_list = client.containers.list(filters={"label":"prog=sb_toolkit"})
    for container in container_list:
        print("shutting down container: ",container.name)
        container.kill()

def call(container,command,cwd='',paths={},remove=True):
    ###strip docker registery if it exists
    container = container.replace('docker://','')

    ###access docker environment
    client = docker.from_env()

    #update or pull image if necessary from docker hub
    low_client = docker.APIClient(base_url='unix://var/run/docker.sock')
    pull_stdout = low_client.pull(container,stream=True)
    out_lines = []
    for line in pull_stdout:
        out = ast.literal_eval(line.decode('utf-8'))
        out_lines.append(out['status'])
        if len(out_lines) == 4:
            print(f"Downloading container {container}, this could take some time.")
        elif len(out_lines) > 4:
            if "Downloaded newer image" in out['status']:
                print(out['status'])
        else:
            pass

    ###get the effectie user and group id's
    user = str(os.geteuid())+':'+str(os.getegid())

    ###setup mount point paths
    #{"/path/outside":"/path/incontainer"}
    volumes = {}
    if paths:
        for key in paths.keys():
            #TODO add more options than only read/write
            volumes[key] = {'bind':paths[key],'mode':'rw'}

    ###run the container
    #try block to run the container
    try:
        container_obj = client.containers.run(container,command,user=user,volumes=volumes,working_dir=cwd,remove=remove,detach=True,labels={"prog":"sb_toolkit"})
    except:
        #loop through output as it is streamed
        for line in container_obj.logs(stream=True):
            yield line.decode('utf-8')
    else:
        for line in container_obj.logs(stream=True):
            yield line.decode('utf-8').strip()
