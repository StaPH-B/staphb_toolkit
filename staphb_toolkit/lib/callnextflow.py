import os, sys
import staphb_toolkit.lib.updates as updates
import pexpect
from shutil import which
from urllib.request import urlopen
from urllib.request import urlretrieve
from urllib import error as urlerror
import json
from datetime import date, datetime
from rich.console import Console
from rich.table import Table
from rich import print

DockerConfig = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config','docker.config'))
SingularityConfig = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config', 'singularity.config'))

def recursive_lookup(key, d):
    def _lookup(key, d):
        if key in d:
            return d[key]

        for value in d.values():
            if isinstance(value, dict):
                result = _lookup(key, value)
                if result:
                    accumulator.update(result)

    accumulator = {}
    result = _lookup(key, d)
    if result:
        accumulator.update(result)
    return accumulator

def recursive_search(key, d):
    if key in d:
        return(d[key])
    for value in d.values():
        if isinstance(value, dict):
            return recursive_search(key, value)
    return None

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
        print('[bold red]Error: Singularity or Docker is not installed or not in found in PATH. [/bold red]')
        sys.exit(1)

    #build command
    command = f'nextflow run -c {config} -r {version} {repo} {args}'
    cmd = nextflow_path + '/' + command

    #run command using nextflow in a subprocess
    child = pexpect.spawn(nextflow_path+'/nextflow pull '+repo,env=env)
    child.interact()
    child.close()
    child = pexpect.spawn(cmd,env=env)
    child.interact()
    child.close()
    return

def get_configuration_file(repo,workflow,wf_version,remote_config_fname):
    config_url = f'https://raw.githubusercontent.com/{repo}/{wf_version}/{remote_config_fname}'
    configFileName = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+f"_{workflow}.config")
    print(f'Downloading {workflow} configuration file to: {configFileName}')
    try:
        urlretrieve(config_url, configFileName)
    except (urlerror.HTTPError, urlerror.URLError):
        print(schema_url)
        sys.exit(1)
    return

def get_workflow_help(repo,workflow,wf_version,schema_fname):
    if not schema_fname:
        print("[bold yellow]Sorry, no help is available for this workflow.[/bold yellow]")
        return
    schema_url = f"https://raw.githubusercontent.com/{repo}/{wf_version}/{schema_fname}"
    try:
        response = urlopen(schema_url)
    except (urlerror.HTTPError, urlerror.URLError):
        print(f"[bold red]Error: Cannot connect to GitHub repository {repo} to fetch parameters for {workflow}.[/bold red]")
        sys.exit(1)
    schema = json.loads(response.read())

    params = recursive_lookup("properties",schema)

    table = Table(title=f"{workflow} Parameters:",title_justify="left",title_style="b",padding=(1,0,0,0))

    table.add_column("Parameter", justify="left",no_wrap=True)
    table.add_column("Description", justify="left")

    for parameter in params.keys():
        table.add_row(parameter,params[parameter]['description'])

    console = Console()
    console.print(table)
    return
