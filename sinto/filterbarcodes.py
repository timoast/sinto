import pysam
import gzip
from multiprocessing import Pool
import functools
import random
import string
from glob import glob
import os
from subprocess import call


def chunk_bam(bamfile, nproc):
    """
    chunk file into n chunks for multicore
    """
    chrom_lengths = bamfile.lengths
    chunksize = sum(chrom_lengths) / int(nproc)
    intervals = []
    for x in range(1, nproc+1):
        position = chunksize*x
        intervals.append(find_chromosome_break(position, chrom_lengths, 0))
    return add_start_coords(intervals, chrom_lengths, bamfile)


def add_start_coords(intervals, chrom_lengths, bamfile):
    """
    given the intervals that will be handled by each core,
    break into genomic regions (removing chromosome-spanning intervals)
    """
    intervals = [[1,0]] + intervals
    ranges = [intervals[x-1] + intervals[x] for x in range(1, len(intervals))]
    d = {}
    x = 0
    for i in ranges:
        x += 1
        if i[0] == i[2]:
            d[x] = [(bamfile.get_reference_name(i[0]-1), i[1], i[3])]
        else:
            d[x] = [(bamfile.get_reference_name(i[0]-1), i[1], chrom_lengths[i[0] - 1])]
            nchrom = i[2] - i[0]
            for y in range(nchrom - 1):
                d[x].append((bamfile.get_reference_name(i[0] + y), 0, chrom_lengths[i[0] + y]))
            d[x].append((bamfile.get_reference_name(i[2] -1), 0, i[3]))
    return(d)


def find_chromosome_break(position, chromosomes, current_chrom):
    assert position <= sum(chromosomes), "position past end of genome"
    if position <= chromosomes[current_chrom]:
        return [current_chrom + 1, position]
    else:
        position = position - chromosomes[current_chrom]
        return find_chromosome_break(position, chromosomes, current_chrom + 1)


def iterate_reads(intervals, bam, sam, output, cb, trim_suffix, mode):
    inputBam = pysam.AlignmentFile(bam, 'rb')
    ident = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
    if sam:
        outputBam = pysam.AlignmentFile(output + ident, "w", template=inputBam)
    else:
        outputBam = pysam.AlignmentFile(output + ident, 'wb', template=inputBam)
    for i in intervals:
        for r in inputBam.fetch(i[0], i[1], i[2]):
            if mode == 'tag':
                cell_barcode, _ = scan_tags(r.tags)
            elif mode == 'readname':
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


def scan_tags(tags, cb="CB", ub="UB"):
    """
    Input bam tags
    Return UMI and cell barcode sequences

    Parameters
    ----------
    tags
        List of read tags
    cb : str
        Tag for cell barcode. Default is CB, as used by 10x
    ub : str
        Tag for UMI barcode. Default is UB, as used by 10x
    """
    cell_barcode = None
    umi = None
    for tag in tags:
        if tag[0] == cb:
            cell_barcode = tag[1]
        elif tag[0] == ub:
            umi = tag[1]
        else:
            pass
    return cell_barcode, umi


def chunk(seq, num):
    """
    cut list into n chunks
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out


def merge_thread_output(data):
    """
    merge multiple dictionaries of the same format into one
    """
    out = {}
    for d in data:
        for cell, counts in d.items():
            try:
                out[cell]
            except KeyError:
                out[cell] = counts
            else:
                out[cell] = [sum(x) for x in zip(counts, out[cell])]
    return out


def filterbarcodes(cells, bam, output, sam=False, trim_suffix=True, nproc=1, mode='tag'):
    """Filter reads based on input list of cell barcodes

    Copy BAM entries matching a list of cell barcodes to a new BAM file.

    Parameters
    ----------
    cells : str
        Path to file containing cell barcodes, or comma-separated list of cell barcodes. File can be gzip compressed.
    bam : str
        Path to BAM file.
    output : str
        Name for output file.
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
    if os.path.isfile(cells):
        if cells.endswith(".gz"):
            cb = [line.strip("\n") for line in gzip.open(cells, "b")]
        else:
            cb = [line.strip("\n") for line in open(cells, "r")]
    else:
        cb = cells.split(",")
    inputBam = pysam.AlignmentFile(bam, 'rb')
    intervals = chunk_bam(inputBam, nproc)
    inputBam.close()
    p = Pool(nproc)
    tempfiles = p.map_async(functools.partial(iterate_reads,
                            bam=bam, sam=sam, output=output,
                            cb=cb, trim_suffix=trim_suffix, mode=mode), intervals.values()).get(9999999)
    mergestring = 'samtools merge -@ ' + str(nproc) + ' ' + output + ' ' + ' '.join(tempfiles)
    call(mergestring, shell=True)
    if os.path.exists(output):
        [os.remove(i) for i in tempfiles]
    else:
        raise Exception("samtools merge failed, temp files not deleted")
