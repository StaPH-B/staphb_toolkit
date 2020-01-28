#!/usr/bin/env python3
import subprocess, os, sys

print("Checking Dependencies...")

toolkit_dir = os.path.dirname(os.path.abspath(__file__))

python_version = ""
pip_cmd = ""
container = False
#check python version
try:
    stdout = subprocess.run(["python3","-V"],stdout=subprocess.PIPE)
    if "Python 3.7" in stdout.stdout.decode('utf-8'):
        python_version = "3.7"
    else:
        os.system("which python3")
        print("Check your python installation, python3 must be >= python version 3.7")
        sys.exit(1)
except FileNotFoundError:
    print("Check your python installation, python3 must be >= python version 3.7")
    sys.exit(1)

#check pip version
try:
    stdout = subprocess.run(["pip3","-V"],stdout=subprocess.PIPE)
    if "python 3.7" in stdout.stdout.decode('utf-8'):
        pip_cmd = "pip3"
except:
    pass

try:
    stdout = subprocess.run(["pip","-V"],stdout=subprocess.PIPE)
    if "python 3.7" in stdout.stdout.decode('utf-8'):
        pip_cmd = "pip"
except:
    pass

if not pip_cmd:
    print("Please verify the python 3.7 pip installation.")
    sys.exit(1)

if python_version and pip_cmd:
    print("Installing dependencies")
    requirements_path = os.path.join(toolkit_dir,"requirements.txt")
    stdout = subprocess.run([pip_cmd,"install","-r",requirements_path],stdout=subprocess.PIPE)
    print(stdout.stdout.decode('utf-8'))

#add to path for bash
#TODO make this more useful for other shells
print("Setting up path...")
stdout = subprocess.run(["which","staphb_toolkit"],stdout=subprocess.PIPE)
if not stdout.stdout.decode('utf-8'):
    homedir = os.environ['HOME']
    bashrc = os.path.join(homedir,".bashrc")
    with open(bashrc,'a') as outfile:
        outfile.write("#path for staphb_toolkit\n")
        outfile.write(f"export PATH=\"{toolkit_dir}:$PATH\"")
print("Done")
