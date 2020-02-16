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

def samtools():
    if (sys.version_info > (3, 0)):
        v = subprocess.check_output(['samtools', '--version']).decode().split()[1].split('.')
    else:
        v = subprocess.check_output(['samtools', '--version']).split()[1].split('.')
    major = int(v[0])
    if major >= 1:
        return True
    return False


if __name__ == "__main__":
    if not samtools():
        raise Exception("sinto requires samtools >= v1")

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
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    scripts = ["scripts/sinto"],
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
