#!/usr/bin/env python3

# author: Kevin Libuit
#email: kevin.libuit@dgs.virginia.gov

import os
import sys
import json
import argparse
import re
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))
from staphB_ToolKit.core import fileparser
from staphB_ToolKit.core import calldocker


class SB_lib:
    def __init__(self, parameters=None, path=None, docker_image=None, executable=None):
        self.path=path
        self.docker_image = docker_image
        self.executable = executable
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
        # create paths for data
        mounting = {self.path: '/data'}
        out_dir = '/data'
        in_dir = '/data'

        command = f"{self.executable} {self.parameters}"
        print(calldocker.call(f"staphb/{self.docker_image}:{self.docker_tag}", command, '/dataout', mounting))
