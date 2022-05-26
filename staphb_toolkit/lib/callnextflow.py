import os, sys
import staphb_toolkit.core.updates as updates
import pexpect
from shutil import which
import urllib.request
from datetime import date, datetime

DockerConfig = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config','docker.config'))
SingularityConfig = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config', 'singularity.config'))

def run(version,repo,args,config_path=None):
    nextflow_path = updates.nextflow_path
    nxf_version = updates.nxf_version
    env = os.environ.copy()
    env['NXF_VER'] = nxf_version
    env['NXF_HOME'] = os.path.expanduser("~")

    ##Test to see if singularity or docker is installed
    if config_path:
        config = config_path
    elif which('docker'):
        config = DockerConfig
    elif which('singularity'):
        config = SingularityConfig
    else:
        print('Singularity or Docker is not installed or not in found in PATH')
        sys.exit(1)

    #build command
    command = f'nextflow run -c {config} -r {version} {repo} {args}'
    cmd = nextflow_path + '/' + command

    #run command using nextflow in a subprocess
    child = pexpect.spawn(nextflow_path+'/nextflow pull '+repo,env=env)
    child.setecho(False)
    child.wait()
    child = pexpect.spawn(cmd,env=env)
    child.interact()

def get_configuration_file(url,workflow):
    configFileName = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+f"_{workflow}.config")
    print(f'Downloading {workflow} configuration file to: {configFileName}')
    urllib.request.urlretrieve(url, configFileName)
