import os, sys
from shutil import which
import re
import subprocess as sub
import shlex
import pexpect
from rich import print

#global app dir
appDir = os.path.abspath(os.path.dirname(__file__)+'/..')

# version number
versionPath = os.path.join(appDir,'lib','VERSION')
with open(versionPath,'r') as versionFile:
    tk_version = versionFile.readline().strip()

#nextflow install path
if which('nextflow'):
    nextflow_path = os.path.abspath(os.path.dirname(which('nextflow')))
else:
    nextflow_path = os.path.join(appDir,'bin')

#nextflow version tracker
nxfVersionPath = os.path.join(appDir,'bin','NXF_VERSION')
try:
    with open(nxfVersionPath,'r') as NXFversionFile:
        nxf_version = NXFversionFile.readline().strip()
except FileNotFoundError:
    nxf_version = ''

# selfupdate check file
selfupdate_status = os.path.join(appDir,'selfupdate')

# auto selfupdate toggle
def toggle_updater(update):
    if update:
        if not os.path.exists(selfupdate_status):
            cmd = "touch " + selfupdate_status
            sub.Popen(shlex.split(cmd)).wait()
            print('Auto-updates on')
    else:
        if os.path.exists(selfupdate_status):
            cmd = "rm " + selfupdate_status
            sub.Popen(shlex.split(cmd)).wait()
            print('Auto-updates off')

# check update status for toolkit
def check_update_status():
    return os.path.exists(selfupdate_status)

# check for toolkit updates
def check_for_updates():
    #regular expressing for checking python version
    re_pattern = r"python\s3.[7-9]|python\s3.[0-9]{2,}"
    cmd = ''
    #Determine pip command version
    if which('pip'):
        #get python version pip is using
        output = sub.check_output(['pip','-V']).decode('utf-8').strip()
        re_obj = re.search(re_pattern,fr"{output}")
        if re_obj:
            cmd = 'pip'
    if which('pip3') and not cmd:
        #get python version pip3 is using
        output = sub.check_output(['pip3','-V']).decode('utf-8')
        re_obj = re.search(re_pattern,fr"{output}")
        if re_obj:
            cmd = 'pip3'

    #check if we have a pip command
    if not cmd:
        print("[bold red]Error: Cannot find pip for python 3.7 or greater.[/bold red]")
        sys.exit(1)

    #run pip update
    cmd = cmd + " install staphb-toolkit --upgrade"
    print("Checking for updates...")
    p = sub.Popen(shlex.split(cmd),stdout=sub.PIPE,stderr=sub.PIPE)
    out, err = p.communicate()
    if err:
        print(err.decode("utf-8"))
    else:
        sout = out.decode("utf-8")
        if "Successfully installed staphb-toolkit" in sout:
            v = sout[-6:].strip()
            print(f"Updated to version {v}")
            print("Done.")
            sys.exit(0)
        else:
            print("No new updates.")

#get nextflow
def install_nextflow():
    #check if we have curl or wget
    if which('curl'):
        cmd = "curl -s https://get.nextflow.io | bash"
    elif which('wget'):
        cmd = "wget -qO- https://get.nextflow.io | bash"
    else:
        print('[bold red]Error: Installing NextFlow requires wget or curl to be installed and in the path.[/bold red]')
        sys.exit(1)

    #pexpect strip $
    def outputFilter(data):
        data = data.decode('utf-8')
        data = data.replace("$","")
        return data.encode('utf-8')

    #install nextflow
    print("Getting the most recent version of NextFlow...")
    try:
        os.makedirs(nextflow_path)
    except FileExistsError:
        pass
    child = pexpect.spawn('/bin/sh',cwd=nextflow_path)
    child.setecho(False)
    child.sendline(cmd)
    child.sendline('exit')
    child.interact(output_filter=outputFilter)

#get nextflow version
def get_nf_version():
    cmd = nextflow_path + "/nextflow -v"
    child = pexpect.spawn(cmd, cwd=os.path.expanduser("~"), env={'NXF_VER':nxf_version,'NXF_HOME':os.path.expanduser("~")})
    child.setecho(False)
    child.interact()

#set nextflow version
def set_nf_version(version='latest'):
    print("Switching NextFlow to version:",version)

    #update nextflow
    print("Performing Nextflow selfupdate this may take a moment...")
    cmd = nextflow_path + "/nextflow self-update"
    child = pexpect.spawn(cmd,cwd=os.path.expanduser("~"), env={'NXF_VER':version,'NXF_HOME':os.path.expanduser("~")})
    child.setecho(False)
    child.wait()

    #an empty version indicates latest
    if version == 'latest':
        version = ''

    #set version
    cmd = nextflow_path + "/nextflow -v"
    child = pexpect.spawn(cmd, cwd=os.path.expanduser("~"), env={'NXF_VER':version,'NXF_HOME':os.path.expanduser("~")})
    child.setecho(False)
    child.interact()
    
    #record new nextflow version
    with open(nxfVersionPath,'w') as NXFversionFile:
        NXFversionFile.write(version)