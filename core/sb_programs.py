#!/usr/bin/env python3

import os
import sys
from shutil import which
import multiprocessing as mp
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

class Run_multi:
    def __init__(self, command_list, path, image, tag):
        self.path = path
        self.command_list = command_list
        self.image = image
        self.tag = tag


    def run(self,jobs):
        #initalize all workers to ignore signal int since we are handeling the keyboard interrupt ourself
        parent_id = os.getpid()
        def init_worker():
            #set signal for docker since containers are detached and we will kill them separately
            if which('docker'):
                signal.signal(signal.SIGINT, signal.SIG_IGN)
            #for singularity we will kill the child processes when the main process gets a signal
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
            results = pool.starmap_async(container_engine.call,[[f"{self.image}:{self.tag}",cmd,'/data',self.path] for cmd in self.command_list])
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
