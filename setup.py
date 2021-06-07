#!/usr/bin/env python

import setuptools
import subprocess
import sys
import re
from pathlib import Path


VERSIONFILE="sinto/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'sinto',
    version = verstr,
    description = "sinto: tools for single-cell data processing",
    long_description = long_description,
    long_description_content_type = "text/x-rst",
    author = 'Tim Stuart',
    install_requires=[
        'pysam>=0.14',
        'scipy',
        'numpy'
    ],
    entry_points = {'console_scripts': ['sinto = sinto.arguments:main']},
    author_email = 'tstuart@nygenome.org',
    url = 'https://github.com/timoast/sinto',
    packages = setuptools.find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research"
    ]
)
