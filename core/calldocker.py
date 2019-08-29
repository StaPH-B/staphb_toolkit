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
    #create empty variable for holding byte object output for the container logs
    output = b''
    #try block to run the container
    try:
        container_obj = client.containers.run(container,command,user=user,volumes=volumes,working_dir=cwd,remove=remove,detach=True,labels={"prog":"sb_toolkit"})
    except:
        #loop through output as it is streamed
        for line in container_obj.logs(stream=True):
            output += line
    else:
        for line in container_obj.logs(stream=True):
            output += line
    #once container is finished return output as a string
    return output.decode('utf-8')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="calldocker.py <docker_container> <command> [options]")
    parser.add_argument("docker_container", type=str, help="Name of the docker container to run")
    parser.add_argument("command", type=str, help="Command to execute in the docker container")
    parser.add_argument("-d", default="", type=str, help="Path to working directory in container")
    parser.add_argument("-p", default="", type=str, help='Dictonary of paths to mount in the container e.x. \'{"/path/outside":"/path/incontainer"}\' note: the quoting')
    parser.add_argument('-k', dest='k', default=True, action='store_false',help="Keep the container after the container has finished")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    container_name = args.docker_container
    cmd = args.command
    cwd = args.d
    remove = args.k
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
    print(call(container_name,cmd,cwd=cwd,paths=paths,remove=remove))
