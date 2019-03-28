#!/usr/bin/env python

from setuptools import setup
import versioneer
import subprocess
import sys


def samtools():
    if (sys.version_info > (3, 0)):
        v = subprocess.check_output(['samtools', '--version']).decode().split()[1].split('.')
    else:
        v = subprocess.check_output(['samtools', '--version']).split()[1].split('.')
    major = int(v[0])
    minor = int(v[1])
    if major >= 1:
        return True
    return False


if __name__ == "__main__":
    if not samtools():
        raise Exception("sinto requires samtools >= v1")


setup(
    name = 'sinto',
    version = versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="sinto: tools for single-cell data processing",
    author = 'Tim Stuart',
    install_requires = [
        'pysam>0.8',
    ],
    scripts = ["scripts/sinto"],
    author_email = 'tstuart@nygenome.org',
    url = 'https://github.com/timoast/sinto',
    packages = ['sinto']
)
