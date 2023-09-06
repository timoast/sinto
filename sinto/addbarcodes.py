import gzip
import os

def a2b(s, inverse=False):
    if inverse:
        return s.decode('ascii')
    return bytes(s, encoding='ascii')

def correct_barcodes(barcodes, whitelist):
    from umi_tools import UMIClusterer

    counts = dict()

    for bc in set(barcodes):
        # add every barcode
        # convert to bytes as UMIClusterer expects
        counts[a2b(bc)] = 1

    # add whitelist to the counter, making it the most abundant so that
    # UMIClusterer makes them the "top"
    for bc in whitelist:
        if not 'N' in bc:
            counts[a2b(bc)] = 10000000000
    clusterer = UMIClusterer(cluster_method='directional')
    
    corrected = dict.fromkeys(barcodes)
    
    for entry in clusterer(counts, threshold=1):
        for bc in entry:
            # assign every sequence the first in the list
            corrected[a2b(bc, inverse=True)] = a2b(entry[0], inverse=True)

    # return the list of corrected
    return [corrected[bc] for bc in barcodes]

def addbarcodes(cb_position, fq1, fq2, fq3=None, prefix="", suffix="", wl=None):
    """Add cell barcode to read names

    Parameters
    ----------
    cb_position : int
        Base positions containing cell barcode information
    fq1 : str
        Fastq file containing cell barcode sequences
    fq2 : str
        Fastq file to add cell barcodes to
    fq3 : str
        Fastq file to add cell barcodes to (optional)
    prefix : str
        Prefix to append to cell barcodes
    suffix : str
        Suffix to append to cell barcodes
    """
    barcodes, whitelist = get_barcodes(f=fq1, bases=cb_position, prefix=prefix, suffix=suffix, wl=wl)
    if len(whitelist) > 0:
        barcodes = correct_barcodes(barcodes, whitelist)
    add_barcodes(f=fq2, cb=barcodes)
    if fq3 is not None:
        add_barcodes(f=fq3, cb=barcodes)


def get_barcodes(f, bases=12, prefix="", suffix="", wl=None):
    f_open = open_fastq(f)
    whitelist = []
    if f.endswith(".gz"):
        gz = True
    else:
        gz = False
    x = 0
    cb = []
    for i in f_open:
        if (x % 4 == 1):
            if gz:
                cb.append(prefix + i.decode("utf-8")[:bases] + suffix)
            else:
                cb.append(prefix + i[:bases] + suffix)
        x += 1
    f_open.close()
    if wl is not None:
        for line in open(wl):
            fields = line.split()
            whitelist.append(fields[0])
    return(cb, whitelist)

def add_barcodes(f, cb):
    f_open = open_fastq(f)
    if f.endswith(".gz"):
        gz = True
        o = f.replace(".fastq.gz", "").replace(".fq.gz", "") + ".barcoded.fastq.gz"
        outfile = gzip.GzipFile(o, mode = "wb")
    else:
        gz = False
        o = f.replace(".fastq", "").replace(".fq", "") + ".barcoded.fastq"
        outfile = open(o, "w+")
    x = 0
    y = 0
    for i in f_open:
        if (x % 4 ==0):
            if gz:
                rdname = i.decode("utf-8")
            else:
                rdname = i
            rdname = "@" + cb[y] + ":" + rdname[1:]
            if gz:
                outfile.write(rdname.encode("utf-8"))
            else:
                outfile.write(rdname)
            y += 1
        else:
            outfile.write(i)
        x += 1
    f_open.close()
    outfile.close()


def open_fastq(f):
    if os.path.isfile(f):
        if f.endswith(".gz"):
            f_open = gzip.open(f, "rb")
        else:
            f_open = open(f, "r")
        return(f_open)
    else:
        raise Exception("File not found")
