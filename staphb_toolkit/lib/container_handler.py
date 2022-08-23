#!/usr/bin/env python3

import os
import sys
from shutil import which
import signal,psutil
from rich import print
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

##Test to see if singularity or docker is installed
if which('docker'):
    from staphb_toolkit.lib import calldocker as container_engine
elif which('singularity'):
    from staphb_toolkit.lib import callsing as container_engine
else:
    print('[bold red]Error: Singularity or Docker is not installed or not in found in PATH.[/bold red]')
    sys.exit(1)


class Run:
    def __init__(self, command, path, image, tag):
        self.path=path
        #look for bash pipes, if they exist use bash to run command otherwise run via docker
        if ">" in command or "<" in command or "|" in command or ";" in command or "&" in command:
            self.command = "bash -c '" + command + "'"
        else:
            self.command = command
        self.image = image
        self.tag = tag

    def run(self):
        #check if we got a local singularity image
        if os.path.isfile(self.image):
            container = self.image
        else:
            container = f"{self.image}:{self.tag}"

        try:
            for line in container_engine.call(container, self.command, '/data', self.path):
                print(line)
        except KeyboardInterrupt:
            container_engine.shutdown()
            sys.exit()
