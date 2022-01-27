Major version changes indicate new functionality
or breaking changes to existing functionality.

Minor version changes indicate bug fixes or
performance improvements to existing functionality
without breaking compatibility with previous versions.

Version 0.8
-----------

0.8.0
~~~~~

- Added ``tagtoname`` and ``nametotag`` commands to copy cell barcodes to/from read tags or read names

Version 0.7
-----------

0.7.6
~~~~~

- Added ``--collapse_within`` parameter to ``fragments`` function to enable only collapsing PCR duplicates if the cell barcode is the same (https://github.com/timoast/sinto/issues/36)
- Added tests for ``fragments`` function and associated small test dataset

0.7.5
~~~~~

- Added ``--shift_plus`` and ``--shift_minus`` parameters to configure Tn5 shift applied in ``fragments`` function (https://github.com/timoast/sinto/issues/33)

0.7.4
~~~~~

- Fixed bug causing some fragments at the end of contigs to be dropped (https://github.com/timoast/sinto/issues/31)

0.7.3
~~~~~

- Fixed error in soft clipping for ``fragments`` function (https://github.com/timoast/sinto/issues/29)

0.7.2
~~~~~

- Added ``min_distance`` parameter to ``fragments``
- Fixed bug in soft clipping for ``fragments`` function

0.7.1
~~~~~

- Code style update for ``tagtotag`` and ``tagtorg``
- Fix bug in ``filterbarcodes`` and ``addtags`` that caused lines in BAM header to be duplicated (https://github.com/timoast/sinto/issues/15)

0.7.0
~~~~~

- New ``tagtotag`` function to copy or move read tags

Version 0.6
-----------

0.6.1
~~~~~

- Bug fixes for ``barcode`` function
- Allow running ``barcode`` on unzipped FASTQ files

0.6.0
~~~~~

- New ``barcode`` function to add cell barcodes to read names in FASTQ file


Version 0.5
-----------

0.5.1
~~~~~

- Fix bug in ``utils.read_cell_barcode_file`` when same cell appears on multiple lines

0.5.0
~~~~~

- Add ``tagtorg`` command to add read groups to BAM according to cell barcode.

Version 0.4
-----------

0.4.3
~~~~~

- Throw error if file is not present for ``filterbarcodes`` and ``addtags``

0.4.2
~~~~~

- Fix bug when cell barcode is None for ``fragments`` function.

0.4.1
~~~~~

- Increase recursion limit to prevent error when running on genomes
  with >1000 scaffolds.

0.4.0
~~~~~

- Removed ``sam`` parameter from ``filterbarcodes``
- Allow multiple groups of cells to be specified in ``filterbarcodes``. 
  This will create a separate BAM file for each unique group of cells.

Version 0.3
-----------

0.3.4
~~~~~

- Memory improvements for ``fragments`` function

0.3.3
~~~~~

- Bug fix for ``fragments`` function when using chromosome containing zero fragments

0.3.2
~~~~~

- Added ``--barcodetag`` and ``--barcode_regex`` arguments to ``filterbarcodes``

0.3.1
~~~~~

- Better handling of BAM file opening/closing
- Add ``max_distance`` parameter to ``fragments`` to remove fragments over a certain length

0.3.0
~~~~~

- added ``fragments`` function to create scATAC fragment file from BAM file
- removed use of versioneer for version tracking


Version 0.2
-----------

- added ``addtags`` function to add read tags to BAM file for different groups of cells

Version 0.1
-----------

First release. Functionality:

- ``filterbarcodes``
