#!/usr/bin/env python3

import os
import sys
from shutil import which
import signal,psutil
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

##Test to see if singularity or docker is installed
if which('docker'):
    from staphb_toolkit.core import calldocker as container_engine
elif which('singularity'):
    from staphb_toolkit.core import callsing as container_engine
else:
    print('Singularity or Docker is not installed or not in found in PATH')
    sys.exit(1)


class Run:
    def __init__(self, command, path, image, tag):
        self.path=path
        #self.command = "bash -c '" + command + "'"
        self.command = command
        self.image = image
        self.tag = tag

    def run(self):
        try:
            for line in container_engine.call(f"{self.image}:{self.tag}", self.command, '/data', self.path):
                print(line)
        except KeyboardInterrupt:
            container_engine.shutdown()
            sys.exit()
