Sinto: single-cell analysis tools
=================================

.. image:: https://github.com/timoast/sinto/workflows/pytest/badge.svg
   :target: https://github.com/timoast/sinto/actions

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: http://bioconda.github.io/recipes/sinto/README.html

.. image:: https://badge.fury.io/py/sinto.svg
    :target: https://badge.fury.io/py/sinto

.. image:: https://pepy.tech/badge/sinto
    :target: https://pepy.tech/project/sinto

Sinto is a toolkit for processing aligned single-cell data. Sinto includes functions to:

- Subset reads from a BAM file by cell barcode
- Create a scATAC-seq fragments file from a BAM file
- Add read tags to a BAM file according to cell barcode information
- Add read groups based on read tags
- Copy or move read tags to another read tag
- Copy cell barcodes to/from read names or read tags
- Add cell barcodes to FASTQ read names

Read the documentation at https://timoast.github.io/sinto

