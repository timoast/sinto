import pysam
import gzip
from multiprocessing import Pool
import functools
import random
import string
import os
from sinto import utils
from subprocess import call


def _add_read_tags(intervals, bam, sam, output, cb, trim_suffix, mode):
    inputBam = pysam.AlignmentFile(bam, "rb")
    ident = "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(6)
    )
    if sam:
        outputBam = pysam.AlignmentFile(output + ident, "w", template=inputBam)
    else:
        outputBam = pysam.AlignmentFile(output + ident, "wb", template=inputBam)
    for i in intervals:
        for r in inputBam.fetch(i[0], i[1], i[2]):
            if mode == "tag":
                cell_barcode, _ = utils.scan_tags(r.tags)
            elif mode == "readname":
                cell_barcode = r.qname.split(":")[0]
            else:
                raise Exception("Unknown mode. Use either tag or readname")
            if cell_barcode is not None:
                if trim_suffix:
                    cell_barcode = cell_barcode[:-2]
                if cell_barcode in cb.keys():
                    r.tags += cb[cell_barcode]
            outputBam.write(r)
    outputBam.close()
    inputBam.close()
    return output + ident


def addtags(bam, tagfile, output, sam=False, trim_suffix=True, mode="tag", nproc=1):
    """Add tags to reads from individual cells

    Copies BAM entries to a new file, adding a read tag to cells matching an input table

    Parameters
    ----------
    bam : str
        Path to BAM file.
    tagfile : str
        Tab-delimited file containing cell barcode, read tag to be added, tag information
    output : str
        Name for output BAM file.
    sam : bool, optional
        Output SAM format. Default is BAM format.
    trim_suffix: bool, optional
        Remove trailing 2 characters from cell barcode in bam file (sometimes needed to match 10x barcodes).
    nproc : int, optional
        Number of processors to use. Default is 1.
    mode : str
        Either tag (default) or readname. Some BAM file store the cell barcode in the readname rather than under
        a read tag.

    Raises
    ------
    Exception
        If samtools merge of temporary BAM files fails
    """
    nproc = int(nproc)
    tags = utils.read_cell_barcode_tag_file(tagfile)
    inputBam = pysam.AlignmentFile(bam, "rb")
    intervals = utils.chunk_bam(inputBam, nproc)
    inputBam.close()
    p = Pool(nproc)
    tempfiles = p.map_async(
        functools.partial(
            _add_read_tags,
            bam=bam,
            sam=sam,
            output=output,
            cb=tags,
            trim_suffix=trim_suffix,
            mode=mode,
        ),
        intervals.values(),
    ).get(9999999)
    mergestring = (
        "samtools merge -@ " + str(nproc) + " " + output + " " + " ".join(tempfiles)
    )
    call(mergestring, shell=True)
    if os.path.exists(output):
        [os.remove(i) for i in tempfiles]
    else:
        raise Exception("samtools merge failed, temp files not deleted")
