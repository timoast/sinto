Processing scATAC-seq data
==========================

In this tutorial we'll demonstrate a simple pipeline for processing scATAC-seq data using Sinto,
bwa_, tabix_, and Genrich_.


Download data
-------------

For this example we'll download scATAC-seq data from
`Chen, Lake, and Zhang (2019) <https://www.nature.com/articles/s41587-019-0290-0>`_
This dataset actually contains both gene expression and DNA accessibility measurements for each
cell, but we can process the DNA accessibility data exactly as though it was a typical scATAC-seq
experiment. In this example we'll just process one replicate, but there are 12 replicates in
total in the original data.

.. code-block:: bash

    wget https://sra-pub-src-1.s3.amazonaws.com/SRR9672090/SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R1_001.fastq.gz.1
    wget https://sra-pub-src-1.s3.amazonaws.com/SRR9672090/SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R2_001.fastq.gz.1
    wget https://sra-pub-src-1.s3.amazonaws.com/SRR9672090/SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R3_001.fastq.gz.1

Attach cell barcodes
--------------------

In this experiment the cell barcodes were sequenced using the first 12 bp of read 1,
and the paired-end genomic DNA reads sequenced using read 2 and read 3. We can use the
``sinto barcode`` function to transfer the cell barcodes from read 1 into the read names
in read 2 and read 3. Storing the cell barcode in the read name is an easy way to track
which reads came from which cells.

This step is much faster when using unzipped fastq files, so we'll first unzip the fastq
files we downloaded earlier.

.. code-block:: bash

    gzip -S .gz.1 -d SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R1_001.fastq.gz.1
    gzip -S .gz.1 -d SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R2_001.fastq.gz.1
    gzip -S .gz.1 -d SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R3_001.fastq.gz.1

.. code-block:: bash

    sinto barcode -b 12 --barcode_fastq SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R1_001.fastq \
        --read1 SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R2_001.fastq \
        --read2 SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R3_001.fastq

.. code-block:: none

    Function run_barcode called with the following arguments:

    barcode_fastq	SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R1_001.fastq
    read1	SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R2_001.fastq
    read2	SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R3_001.fastq
    bases	12
    prefix	
    suffix	
    func	<function run_barcode at 0x7f28c93714c0>

    Function completed in  8.0 m 13.35 s

Map
---

We can map the paired-end genomic reads using any aligner we like. In this example we'll use
bwa mem. First we need to download the fasta file for the mm10 genome and create a genome index for bwa:


.. code-block:: bash

    wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    gzip -d mm10.fa.gz
    bwa index mm10.fa

Next we can map the reads with bwa mem:

.. code-block:: bash

    bwa mem -t 8 mm10.fa SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R2_001.barcoded.fastq \
        SC_20190509A_AdMsBrain_SNARE_chromatin_S16_R3_001.barcoded.fastq \
        | samtools view -b - > aln.bam

    # sort and index bam file
    samtools sort -@ 8 aln.bam -o aln.sort.bam
    samtools index -@ 8 aln.sort.bam

.. code-block:: bash

    # look at some of the aligned reads
    samtools view aln.sort.bam | head -n 3

.. code-block:: none

    AACCTCCAACTG:D00611:698:CDN4HANXX:5:1106:1580:94969	69	chr1	3000052	0	*	=	3000052	0	TCACATTCTCAGTGCACAATAGAACCCCTTACCTCCAATCCAGAGTAAACAAAGAGTCTACCACAAACACACATC	/</<///</BB/////<//////<///<//<<////<B<BBB///////<<//////<<//<////</<<FF<</	MC:Z:32M43S	AS:i:0	XS:i:0
    AACCTCCAACTG:D00611:698:CDN4HANXX:5:1106:1580:94969	137	chr1	3000052	0	32M43S	=	3000052	0	TTGAAGGTCTGGTAGAACTCTGCATTAAACCCCGACTCCAGTTGGAGGTTGTACTCTGCGTTGATACCACTTTTT	BBB/BFFFBFBFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFB<///////	NM:i:0	MD:Z:32	AS:i:32	XS:i:32
    TCCGCCGGAAAC:D00611:697:CD0V6ANXX:7:2116:7385:7732	117	chr1	3000098	0	*	=	3000098	0	AAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACATATACACTCAGCTTTAT	FFFFFFFFFFFFFFFFFFFFFFFFFFFF<FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBB	MC:Z:33M42S	AS:i:0	XS:i:0

Create a fragment file
----------------------

Finally we create a fragment file for the dataset using ``sinto fragments``

.. code-block:: bash

    sinto fragments -b aln.sort.bam -p 8 -f fragments.bed --barcode_regex "[^:]*"

.. code-block:: none

    Function run_fragments called with the following arguments:

    bam	aln.sort.bam
    fragments	fragments.bed
    min_mapq	30
    nproc	8
    barcodetag	CB
    cells	None
    barcode_regex	[^:]*
    use_chrom	(?i)^chr
    max_distance	5000
    chunksize	500000
    func	<function run_fragments at 0x7ff6976b5dc0>

    Function completed in  7.0 m 30.24 s

.. code-block:: bash

    # sort, compress, and index
    sort -k1,1 -k2,2n fragments.bed > fragments.sort.bed
    bgzip -@ 8 fragments.sort.bed
    tabix -p bed fragments.sort.bed.gz

    # clean up
    rm fragments.bed

    # take a look at the output
    gzip -dc fragments.sort.bed.gz | head

.. code-block:: none

    chr1	3003930	3003960	ATGGTGACTCAT	1
    chr1	3010068	3010136	CTCGGGTTTAGC	5
    chr1	3010538	3010589	GAAACGCAAAGT	8
    chr1	3011663	3011770	AACTGGGCTATC	2
    chr1	3011684	3011728	CTATTGTGCATA	12
    chr1	3012702	3012739	ATCAGAGACGGC	6
    chr1	3012724	3012775	TGCGGACGGTGG	3
    chr1	3012743	3012804	TGACATCCCCGA	2
    chr1	3019505	3019588	CACGCTCGTCTT	1
    chr1	3021571	3021620	GTATTCGTCGGG	9


Call peaks
----------

We can call peaks using all cells combined using the mapped bam file and Genrich_.
Genrich requires reads to be sorted by queryname, so we first re-sort the bam file and
then run Genrich using ATAC-seq mode (``-j``) to find peaks.

.. code-block:: bash

    samtools sort -n -@ 8 aln.sort.bam -o aln.qname.bam
    Genrich -j -t aln.qname.bam -o peaks.bed -v

.. code-block:: bash

    # look at some of the peaks
    head peaks.bed

.. code-block:: none

    chr1	3094859	3095377	peak_0	1000	.	913.693115	5.720295	-1	374
    chr1	3119581	3120170	peak_1	1000	.	1186.086792	6.389130	-1	420
    chr1	3120272	3120746	peak_2	1000	.	965.524841	6.166295	-1	285
    chr1	3121354	3121654	peak_3	1000	.	580.268433	4.948659	-1	129
    chr1	3292563	3292985	peak_4	1000	.	897.739197	5.319649	-1	257
    chr1	3299692	3299934	peak_5	1000	.	489.207764	5.154057	-1	98
    chr1	3309998	3310407	peak_6	1000	.	829.686340	6.153750	-1	210
    chr1	3322437	3322775	peak_7	1000	.	397.093872	3.806610	-1	128
    chr1	3360973	3361167	peak_8	1000	.	262.217590	3.911148	-1	120
    chr1	3369527	3369829	peak_9	1000	.	373.499207	4.094081	-1	149

Downstream analysis
-------------------

Downstream analysis steps, including quantifying counts in each peak for each cell,
can be performed using Signac_.

.. _bwa: https://github.com/lh3/bwa
.. _tabix: https://github.com/samtools/htslib
.. _Genrich: https://github.com/jsh58/Genrich
.. _Signac: https://satijalab.org/signac