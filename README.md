# sinto
Tools for single-cell data processing

## Installation

Install from PyPI:

```bash
pip install sinto
```

Install from source:

```bash
git clone https://github.com/timoast/sinto.git
cd sinto
python setup.py install
```

## Features

### Filter barcodes

If you want to extract reads from a subset of cell barcodes, you can
do so using the `sinto filterbarcodes` command.

```
$ sinto filterbarcodes -h
usage: sinto filterbarcodes [-h] -b BAM -c CELLS -o OUTPUT [-t] [-s]
                            [-p NPROC] [-m MODE]

Filter reads based on input list of cell barcodes

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     Input bam file (must be indexed)
  -c CELLS, --cells CELLS
                        File or comma-separated list of cell barcodes. Can be
                        gzip compressed
  -o OUTPUT, --output OUTPUT
                        Name for output text file
  -t, --trim_suffix     Remove trail 2 characters from cell barcode in BAM
                        file
  -s, --sam             Output sam format (default bam output)
  -p NPROC, --nproc NPROC
                        Number of processors (default = 1)
  -m MODE, --mode MODE  Either tag (default) or readname. Some BAM file store
                        the cell barcode in the readname rather than under a
                        read tag
```