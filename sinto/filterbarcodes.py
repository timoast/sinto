import pysam
from multiprocessing import Pool
import functools
import random
import string
from sinto import utils
from subprocess import call
import re
import os
from itertools import chain


def _iterate_reads(intervals, bam, cb, classes, trim_suffix, cellbarcode, readname_barcode):
    inputBam = pysam.AlignmentFile(bam, "rb")
    ident = "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(6)
    )
    filelist = [x + "_" + ident for x in classes]
    bamlist = [pysam.AlignmentFile(x, "wb", template=inputBam) for x in filelist]
    for i in intervals:
        for r in inputBam.fetch(i[0], i[1], i[2]):
            if readname_barcode is not None:
                re_match = readname_barcode.match(r.qname)
                cell_barcode = re_match.group()
            else:
                cell_barcode, _ = utils.scan_tags(r.tags, cb=cellbarcode)
            if cell_barcode is not None:
                if trim_suffix:
                    cell_barcode = cell_barcode[:-2]
                if cell_barcode in cb.keys():
                    cell_classes = cb[cell_barcode]
                    for i in cell_classes:
                        fileindex = filelist.index(i + "_" + ident)
                        bamlist[fileindex].write(r)
    for i in bamlist:
        i.close()
    inputBam.close()
    return ident


def mergeAll(idents, classes, nproc,  remove=True):
    """Merge all temp files for each class
    If remove is True, remove temp files after
    successful merging"""
    for i in classes:
        allfiles = [i + "_" + x for x in idents]
        output = i + '.bam'
        mergestring = (
            "samtools merge -@ " + str(nproc) + " " + output + " " + " ".join(allfiles)
        )
        call(mergestring, shell=True)
        if os.path.exists(output):
            [os.remove(i) for i in allfiles]
        else:
            raise Exception("samtools merge failed, temp files not deleted")


def filterbarcodes(
    cells, bam, readname_barcode, cellbarcode, sam=False, trim_suffix=True, nproc=1
):
    """Filter reads based on input list of cell barcodes

    Copy BAM entries matching a list of cell barcodes to a new BAM file.
    Output BAM files will be named according to the group name in the 
    file provided.

    Parameters
    ----------
    cells : str
        Path to file containing cell barcodes and the group associated with each barcode.
        File can be gzip compressed. A separate BAM file will be created for each 
        group of cells.
    bam : str
        Path to BAM file.
    trim_suffix: bool, optional
        Remove trailing 2 characters from cell barcode in bam file (sometimes needed to match 10x barcodes).
    nproc : int, optional
        Number of processors to use. Default is 1.
    cellbarcode : str
       Tag used for cell barcode. Default is CB (used by cellranger)
    readname_barcode : regex
        A regular expression for matching cell barcode in read name. If None (default),
        use the read tags.

    Raises
    ------
    Exception
        If samtools merge of temporary BAM files fails
    """
    nproc = int(nproc)
    cb = utils.read_cell_barcode_file(cells)
    unique_classes = list(set(chain.from_iterable(cb.values())))
    inputBam = pysam.AlignmentFile(bam, "rb")
    intervals = utils.chunk_bam(inputBam, nproc)
    inputBam.close()
    if readname_barcode is not None:
        readname_barcode = re.compile(readname_barcode)
    p = Pool(nproc)
    idents = p.map_async(
        functools.partial(
            _iterate_reads,
            bam=bam,
            cb=cb,
            classes=unique_classes,
            trim_suffix=trim_suffix,
            cellbarcode=cellbarcode,
            readname_barcode=readname_barcode
        ),
        intervals.values(),
    ).get(9999999)
    mergeAll(idents=idents, classes=unique_classes, nproc=nproc, remove=True)
