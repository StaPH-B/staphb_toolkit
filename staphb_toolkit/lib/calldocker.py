#!/usr/bin/env python3

import docker
import os,sys
import argparse
import json
import shlex
import ast
from rich.progress import Progress
from rich import print

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
    low_client = docker.APIClient()
    pull_stdout = low_client.pull(container,stream=True)

    with Progress() as progress:
        task_total = 0
        task = progress.add_task(f"Pulling {container}",total = task_total,visible=False)
        for line in pull_stdout:
            out = ast.literal_eval(line.decode('utf-8'))
            status = out['status']
            if status == "Pulling fs layer":
                if task_total == 0:
                    task_total = 2
                    progress.update(task,total=task_total,visible=True)
                else:
                    task_total += 2
                    progress.update(task,total=task_total)

            elif status == "Download complete":
                progress.update(task,advance=1)

            elif status == "Pull complete":
                progress.update(task,advance=1)

            elif "Downloaded newer image" in status:
                    print(status)
                    break
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
        container_obj = client.containers.run(
            container,
            command,
            user=user,
            volumes=volumes,
            working_dir=cwd,
            remove=False,
            detach=True,
            labels={"prog":"sb_toolkit"}
        )
    except Exception as e:
        print(f"[bold red]There was an error trying to run the container {container}[/bold red]")
        print(e)

    for line in container_obj.logs(stream=True):
        yield line.decode('utf-8').strip()

    if remove == True:
        container_obj.remove()