Major version changes indicate new functionality
or breaking changes to existing functionality.

Minor version changes indicate bug fixes or
performance improvements to existing functionality
without breaking compatibility with previous versions.

Version 0.4
-----------

0.4.1
~~~~~

- Increase recursion limit to prevent error when running on genomes
  with a >1000 scaffolds.

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
