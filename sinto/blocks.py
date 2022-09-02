import pysam
from sinto import utils
import re
import os
from multiprocessing import Pool
import functools
import tempfile


def getBlocks(
    intervals,
    bam,
    min_mapq=30,
    cellbarcode="CB",
    umibarcode="UB",
    readname_barcode=None,
    cells=None,
    include_strand=False
):
    inputBam = pysam.AlignmentFile(bam, "rb")
    outfile = tempfile.NamedTemporaryFile(delete=False)
    outname = outfile.name
    outf = open(outname, "w+")
    for i in intervals:
        for r in inputBam.fetch(i[0], i[1], i[2]):
            if r.mapping_quality < min_mapq:
                continue
            if readname_barcode is not None:
                re_match = readname_barcode.match(r.qname)
                cell_barcode = re_match.group()
                _, umi = utils.scan_tags(r.tags, ub=umibarcode)
            else:
                cell_barcode, umi = utils.scan_tags(r.tags, cb=cellbarcode, ub=umibarcode)
            if cell_barcode is None or umi is None:
                continue
            if (cells is not None) and (cell_barcode not in cells):
                continue
            blocks = r.get_blocks()
            chrom = r.reference_name
            strand = "-" if r.is_reverse else "+"
            for i in blocks:
                outstr = "\t".join(map(str, [chrom, i[0], i[1], cell_barcode, umi]))
                if include_strand:
                    outstr += "\t" + strand
                outf.write(outstr + "\n")
    outf.close()
    inputBam.close()
    return(outname)


def blocks(
    bam,
    block_path,
    min_mapq=30,
    nproc=1,
    cellbarcode="CB",
    umibarcode="UB",
    readname_barcode=None,
    cells=None,
    include_strand=False,
):
    """Create scRNA-seq block file from BAM file

    Iterate over all aligned reads in a BAM file and extract the 
    coordinates of each aligned block. Each aligned block is written
    as a separate entry in a BED file, along with the cell barcode
    and UMI sequence for the read.

    Parameters
    ----------
    bam : str
        Path to BAM file
    block_path : str
        Path for output block file
    min_mapq : int
        Minimum MAPQ to retain fragment
    nproc : int, optional
        Number of processors to use. Default is 1.
    cellbarcode : str
        Tag used for cell barcode. Default is CB (used by cellranger)
    umibarcode : str
        Tag used for the UMI sequence. Defauly is UB (used by cellranger)
    readname_barcode : str, optional
        Regular expression used to match cell barocde stored in read name. 
        If None (default), use read tags instead. Use "[^:]*" to match all characters 
        before the first colon (":").
    cells : str
        File containing list of cell barcodes to retain. If None (default), use all cell barcodes
        found in the BAM file.
    include_strand : bool
        Include strand information
    """
    nproc = int(nproc)
    inputBam = pysam.AlignmentFile(bam, "rb")
    intervals = utils.chunk_bam(inputBam, nproc, unmapped=False)
    inputBam.close()
    if readname_barcode is not None:
        readname_barcode = re.compile(readname_barcode)
    cells = utils.read_cells(cells)
    p = Pool(nproc)
    tmpfiles = [
        p.map_async(
            functools.partial(
                getBlocks,
                bam=bam,
                min_mapq=int(min_mapq),
                cellbarcode=cellbarcode,
                umibarcode=umibarcode,
                readname_barcode=readname_barcode,
                cells=cells,
                include_strand=include_strand,
            ),
            intervals.values(),
        )
    ]
    filenames = [res.get() for res in tmpfiles]
    # cat files and write to output
    with open(block_path, "w") as outfile:
        for i in filenames:
            for j in i:
                with open(j, "r") as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(j)