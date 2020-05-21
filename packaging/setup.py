#!/usr/bin/env python3
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="staphb_toolkit",
    version="1.0.5",
    author="Kelsey Florek, Kevin Libuit",
    author_email="kelsey.florek@slh.wisc.edu, kevin.libuit@dgs.virginia.gov",
    description="A ToolKit of commonly used Public Health Bioinformatics Tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/StaPH-B/staphb_toolkit",
    packages=setuptools.find_packages(),
    include_package_data=True,
    entry_points={
            "console_scripts":[
            'staphb-tk = staphb_toolkit.toolkit_apps:main',
            'staphb-wf = staphb_toolkit.toolkit_workflows:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "spython>=0.0.73",
        "psutil>=5.6.3",
        "docker>=4.1.0"],
    python_requires='>=3.6',
)
