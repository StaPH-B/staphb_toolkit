#!/usr/bin/env python3
import setuptools, os

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("./staphb_toolkit/lib/VERSION",'r') as versionFile:
    version = versionFile.readline().strip()

setuptools.setup(
    name="staphb_toolkit",
    version=version,
    author="Kelsey Florek",
    author_email="kelsey.florek@slh.wisc.edu",
    description="A ToolKit of commonly used Public Health Bioinformatics Tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/StaPH-B/staphb_toolkit",
    packages=setuptools.find_packages(),
    package_data={
        'staphb_toolkit' : ['lib/VERSION'],
        '' : ['../requirements.txt','../workflows.json','../staphb_toolkit/config/docker.config','../staphb_toolkit/config/singularity.config','../staphb-tk']
    },
    entry_points={
        "console_scripts": [
            'staphb-tk = staphb_toolkit.toolkit_main:main'
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "spython>=0.0.73",
        "psutil>=5.6.3",
        "docker>=4.1.0",
        "pexpect>=4.8",
        "pyfiglet>=0.8.post1",
        "rich>=12.4.4"],
    python_requires='>=3.7',
)
