import functools
import time


def log_info(func):
    """Decorator that prints function arguments and runtime
    """
    @functools.wraps(func)
    def wrapper(args):
        print("Function {} called with the following arguments:\n".format(func.__name__))
        for arg in vars(args):
            print(str(arg) + '\t' + str(getattr(args, arg)))
        t1 = time.time()
        func(args)
        t2 = time.time()
        elapsed = [round(x, 2) for x in divmod(t2-t1, 60)]
        print("\nFunction completed in  {} m {} s\n".format(elapsed[0], elapsed[1]))
    return wrapper


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