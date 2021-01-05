import os, sys
from shutil import which
import re
import subprocess as sub
import shlex

#version number
versionPath = os.path.abspath(os.path.dirname(__file__) + '/' + 'VERSION')
with open(versionPath,'r') as versionFile:
    version = versionFile.readline().strip()

#selfupdate check file
selfupdate_status = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))),'selfupdate')

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

def check_update_status():
    return os.path.exists(selfupdate_status)

def check_for_updates():
    #regular expressing for checking python version
    re_pattern = r"python\s3.[6-9]|python\s3.[0-9]{2,}"
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
        print("Cannot find pip for python 3.6 or greater...")
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
        else:
            print("No new updates.")
