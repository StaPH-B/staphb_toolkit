#!/usr/bin/env python3

import os
import sys
import json
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

class SB_lib:
    def __init__(self, parameters=None, path=None, docker_image=None, executable=None):
        self.path=path
        self.docker_image = docker_image
        self.executable = executable
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

    def run_lib(self):
        command = f"{self.executable} {self.parameters}"
        print(container_engine.call(f"staphb/{self.docker_image}:{self.docker_tag}", command, '/data', self.path))
