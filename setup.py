#!/usr/bin/env python3
import subprocess, os, sys

print("Checking Dependencies...")

python_version = ""
pip_cmd = ""
container = False
#check python version
stdout = subprocess.run(["python3","-V"],stdout=subprocess.PIPE)
if "Python 3.7" in stdout.stdout.decode('utf-8'):
    python_cmd = "3.7"
else:
    os.system("which python3")
    print("Check your python installation, python3 must be >= python version 3.7")
    sys.exit(1)

#check pip version
stdout = subprocess.run(["pip3","-V"],stdout=subprocess.PIPE)
if "python 3.7" in stdout.stdout.decode('utf-8'):
    pip_cmd = "pip3"

stdout = subprocess.run(["pip","-V"],stdout=subprocess.PIPE)
if "python 3.7" in stdout.stdout.decode('utf-8'):
    pip_cmd = "pip"

if not pip_cmd:
    print("Please verify the python 3.7 pip installation.")
    sys.exit(1)

if python_version and pip_cmd:
    print("Installing dependencies")
    requirements_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"requirements.txt")
    stdout = subprocess.run([pip_cmd,"-r",requirements_path],stdout=subprocess.PIPE)
    print(stdout.stdout.decode('utf-8'))
