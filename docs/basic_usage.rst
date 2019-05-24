Basic usage
===========

Filter cell barcodes from BAM file
----------------------------------

Reads for a subset of cells can be extracted from a BAM file using the ``filterbarcodes`` command.
This requires a position-sorted, indexed BAM file, and a file containing a list of cell barcodes to retain.

.. code-block:: none

    sinto filterbarcodes [-h] -b BAM -c CELLS -o OUTPUT [-t] [-s]
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

Add read tags to BAM file
-------------------------

Read tags can be added to a BAM file according to which cell the read belongs to using the ``addtags`` command.
This requires a position-sorted and indexed BAM file, and a file specifying the tags to be added to each cell, for example:

.. code-block:: none

    TGGCAATGTTGAAGCG-1	CI	A
    GACCAATCACCATTCC-1	CI	A
    CAGGATTCAGAACTTC-1	CI	B
    GAACCTAAGAGAGGTA-1	CI	B
    ACATGGTGTAGACGCA-1	CI	C
    CCCTGATTCGGATAGG-1	CI	C

.. code-block:: none

    sinto addtags [-h] -b BAM -f TAGFILE -o OUTPUT [-t] [-s] [-p NPROC]
                        [-m MODE]

    Add read tags to reads from individual cells

    optional arguments:
    -h, --help            show this help message and exit
    -b BAM, --bam BAM     Input bam file (must be indexed)
    -f TAGFILE, --tagfile TAGFILE
                            Tab-delimited file containing cell barcode, tag to be
                            added, and tag identity. Can be gzip compressed
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


This will add a ``CI`` tag, with the tag set to A, B, or C depending on the cell barcode sequence.
