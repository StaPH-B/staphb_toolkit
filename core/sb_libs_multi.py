#!/usr/bin/env python3

import os
import sys
import json
import multiprocessing as mp
from shutil import which
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

##Test to see if singularity or docker is installed
if which('docker'):
    from staphb_toolkit.core import calldocker as container_engine
elif which('singularity'):
    from staphb_toolkit.core import callsing as container_engine
else:
    print('Singularity or Docker is not installed or not in found in PATH')
    sys.exit(1)

##define and set signal handler, to handle detached containers
def handler(sig,frame):
    #shutdown containers
    container_engine.shutdown()
    sys.exit()

##attach signal handler
signal.signal(signal.SIGINT, handler)

class SB_lib_multi:
    def __init__(self, command_list=None, path=None, docker_image=None):
        self.path = path
        self.docker_image = docker_image
        self.command_list = command_list
        # TODO: Find a better way to grab json file
        docker_config = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))[:-4] + "/core/docker_config.json"
        docker_tag = ""

        with open(docker_config) as config_file:
            config_file = json.load(config_file)
            for image in config_file['images']:
                if image['name'] == docker_image:
                    docker_tag = image['tag']

        self.docker_tag = docker_tag
        self.parameters = parameters

    def run_lib(self,jobs):
        #create multiprocessing pool
        pool = mp.Pool(processes=jobs)

        #set signal handling to default during spawning of subthreads
        signal.signal(signal.SIGINT,signal.SIG_DFL)

        results = pool.starmap_async(cd.call,[[f"staphb/{self.docker_image}:{self.docker_tag}",cmd,'/data',self.path] for cmd in self.command_list])

        #reset signal handling to handle shutting down containers
        signal.signal(signal.SIGINT, handler)

        stdouts = results.get()
        for stdout in stdouts:
            print(stdout)
