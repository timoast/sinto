import pysam
from multiprocessing import Pool
import functools
import random
import string
from sinto import utils
from subprocess import call


def _iterate_reads(intervals, bam, sam, output, cb, trim_suffix, mode):
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
                    if cell_barcode[:-2] in cb:
                        outputBam.write(r)
                else:
                    if cell_barcode in cb:
                        outputBam.write(r)
    outputBam.close()
    inputBam.close()
    return output + ident


def filterbarcodes(
    cells, bam, output, sam=False, trim_suffix=True, nproc=1, mode="tag"
):
    """Filter reads based on input list of cell barcodes

    Copy BAM entries matching a list of cell barcodes to a new BAM file.

    Parameters
    ----------
    cells : str
        Path to file containing cell barcodes, or comma-separated list of cell barcodes. File can be gzip compressed.
    bam : str
        Path to BAM file.
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
    cb = utils.read_cells(cells)
    inputBam = pysam.AlignmentFile(bam, "rb")
    intervals = utils.chunk_bam(inputBam, nproc)
    inputBam.close()
    p = Pool(nproc)
    tempfiles = p.map_async(
        functools.partial(
            _iterate_reads,
            bam=bam,
            sam=sam,
            output=output,
            cb=cb,
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
