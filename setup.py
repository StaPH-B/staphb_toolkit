#!/usr/bin/env python3
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="staphb_toolkit", # Replace with your own username
    version="0.1.0",
    author="Kelsey Florek, Kevin Libuit",
    author_email="kelsey.florek@slh.wisc.edu, kevin.libuit@dgs.virginia.gov",
    description="A ToolKit of commonly used Public Health Bioinformatics Tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/StaPH-B/staphb_toolkit",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
