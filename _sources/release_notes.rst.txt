Version 0.3
-----------

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
