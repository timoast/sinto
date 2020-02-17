Basic usage
===========

Create scATAC-seq fragments file
--------------------------------

An ATAC-seq fragment file can be created from a BAM file using the ``fragments`` command.
The fragment file contains the position of each Tn5 integration site, the cell barcode 
associated with the fragment, and the number of times the fragment was sequenced. 
PCR duplicates are collapsed. See `10x Genomics <https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments>`_
for a further description of the fragment file format.

.. code-block:: none
    
    sinto fragments [-h] -b BAM -f FRAGMENTS [-m MIN_MAPQ] [-p NPROC]
                       [-t BARCODETAG] [-c CELLS]
                       [--barcode_regex BARCODE_REGEX] [--use_chrom USE_CHROM]
                       [--max_distance MAX_DISTANCE]

    Create ATAC-seq fragment file from BAM file

    optional arguments:
    -h, --help            show this help message and exit
    -b BAM, --bam BAM     Input bam file (must be indexed)
    -f FRAGMENTS, --fragments FRAGMENTS
                            Name and path for output fragments file. Note that the
                            output is not sorted or compressed. To sort the output
                            file use sort -k 1,1 -k2,2n
    -m MIN_MAPQ, --min_mapq MIN_MAPQ
                            Minimum MAPQ required to retain fragment (default =
                            30)
    -p NPROC, --nproc NPROC
                            Number of processors (default = 1)
    -t BARCODETAG, --barcodetag BARCODETAG
                            Read tag storing cell barcode information (default =
                            "CB")
    -c CELLS, --cells CELLS
                            Path to file containing cell barcodes to retain, or a
                            comma-separated list of cell barcodes. If None
                            (default), use all cell barocodes present in the BAM
                            file.
    --barcode_regex BARCODE_REGEX
                            Regular expression used to extract cell barcode from
                            read name. If None (default), extract cell barcode
                            from read tag. Use "[^:]*" to match all characters up
                            to the first colon.
    --use_chrom USE_CHROM
                            Regular expression used to match chromosomes to be
                            included in output. Default is "(?i)^chr" to match all
                            chromosomes starting with "chr", case insensitive
    --max_distance MAX_DISTANCE
                            Maximum distance between integration sites for the
                            fragment to be retained. Allows filtering of
                            implausible fragments that likely result from
                            incorrect mapping positions. Default is 5000 bp.

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
                            Name for output BAM file
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
                            Name for output BAM file
    -t, --trim_suffix     Remove trail 2 characters from cell barcode in BAM
                            file
    -s, --sam             Output sam format (default bam output)
    -p NPROC, --nproc NPROC
                            Number of processors (default = 1)
    -m MODE, --mode MODE  Either tag (default) or readname. Some BAM file store
                            the cell barcode in the readname rather than under a
                            read tag


This will add a ``CI`` tag, with the tag set to A, B, or C depending on the cell barcode sequence.
