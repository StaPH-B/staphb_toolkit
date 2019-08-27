#!/usr/bin/env python3

import os
import sys
import json
import multiprocessing as mp
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

    def run_lib(self,jobs):
        #initalize all workers to ignore signal int since we are handeling the keyboard interrupt ourself
        parent_id = os.getpid()
        def init_worker():
            #only set signal for docker since containers are detached
            if which('docker'):
                signal.signal(signal.SIGINT, signal.SIG_IGN)
            else:
                def sig_int(signal_num,frame):
                    parent = psutil.Process(parent_id)
                    for child in parent.children():
                        if child.pid != os.getpid():
                            child.kill()
                    psutil.Process(os.getpid()).kill()
                signal.signal(signal.SIGINT,sig_int)

        #create multiprocessing pool
        pool = mp.Pool(processes=jobs,initializer=init_worker)

        try:
            results = pool.starmap_async(container_engine.call,[[f"staphb/{self.docker_image}:{self.docker_tag}",cmd,'/data',self.path] for cmd in self.command_list])
            stdouts = results.get()

        except KeyboardInterrupt:
            pool.terminate()
            pool.join()
            #shutdown containers
            container_engine.shutdown()
            sys.exit()
        else:
            pool.close()
            pool.join()

        for stdout in stdouts:
            print(stdout)
