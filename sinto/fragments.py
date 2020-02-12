import pysam
from sinto import utils
from scipy import sparse
import numpy as np
from collections import Counter, defaultdict


def writeFragments(fragments, filepath):
    """Write fragments to file

    Parameters
    ----------
    fragments : list
        List of ATAC fragments
    filepath : str
        Path for output file
    """
    with open(filepath, "w") as outf:
        for i in fragments:
            outstr = "\t".join(map(str, i))
            outf.write(outstr + "\n")


def collapseFragments(fragments, retain_cross_bc_duplicates=False):
    """Collapse duplicate fragments

    retain_cross_bc_duplicates is not implemented (TODO)
    """
    fraglist = [list(x.values()) for x in list(fragments.values())]
    fragcoords_with_bc = [".".join(map(str, x)) for x in fraglist]
    fragcoords = [".".join(map(str, x[:3])) for x in fraglist]
    cellbarcodes = [x[3] for x in fraglist]
    counts = Counter(fragcoords_with_bc)
    # enumerate fragments and barcodes
    frag_id_lookup = id_lookup(l=fragcoords)
    bc_id_lookup = id_lookup(l=cellbarcodes)

    # get list of barcode index and fragment index from counts
    row = []
    col = []
    data = []
    for i in list(counts.items()):
        data.append(i[1])
        rowstr = i[0].split(".")[:3]
        bcstr = i[0].split(".")[3]
        row.append(frag_id_lookup[".".join(rowstr)])
        col.append(bc_id_lookup[bcstr])

    # create sparse matrix of fragment counts from fraglist (column, row, value)
    mat = sparse.coo_matrix(
        (data, (np.array(row), np.array(col))),
        shape=(max(frag_id_lookup.values()) + 1, max(bc_id_lookup.values()) + 1),
    )
    mat = mat.tocsr()

    # find which barcode contains the most counts for each fragment
    rowsums = mat.sum(axis=1)
    colsums = mat.sum(axis=0)
    rowmax = np.argmax(mat, axis=1)
    rowsum = mat.sum(axis=1).tolist()
    # maxval = (
    #     mat.max(axis=1).toarray().tolist()
    # )

    # # construct a new matrix with 1s in the max position
    # collapsed_mat = sparse.coo_matrix(
    #     (np.ones(mat.shape[0]), (np.arange(start=0, stop=mat.shape[0]), np.squeeze(np.array(rowmax))))
    # )
    # collapsed_mat = collapsed_mat.tocsr()

    # collapse back into a list of fragment coords and barcodes
    frag_inverse = dict(zip(frag_id_lookup.values(), frag_id_lookup.keys()))
    bc_inverse = dict(zip(bc_id_lookup.values(), bc_id_lookup.keys()))
    collapsed_frags = [frag_inverse[x] for x in np.arange(mat.shape[0])]
    collapsed_barcodes = [bc_inverse[x[0]] for x in rowmax.tolist()]
    collapsed = []
    for i in range(len(collapsed_barcodes)):
        frag = collapsed_frags[i].split(".")
        frag.append(collapsed_barcodes[i])
        frag.append(rowsum[i][0])
        collapsed.append(frag)
    return collapsed


def id_lookup(l):
    """Create dictionary where each unique item is key, value is the item numerical ID"""
    temp = defaultdict(lambda: len(temp))
    idx = [temp[x] for x in l]
    lookup = dict(zip(set(l), set(idx)))
    return lookup


def getFragments(bam, min_mapq=30, nproc=1, cellbarcode="CB"):
    """Extract ATAC fragments from BAM file

    Iterate over paired reads in a BAM file and extract the ATAC fragment coordinates

    Parameters
    ----------
    bam : str
        Path to BAM file
    min_mapq : int
        Minimum MAPQ to retain fragment
    nproc : int, optional
        Number of processors to use. Default is 1.
    cellbarcode : str
       Tag used for cell barcode. Default is CB (used by cellranger)
    """
    nproc = int(nproc)
    fragment_dict = dict()
    inputBam = pysam.AlignmentFile(bam, "rb")
    for i in inputBam.fetch():
        fragment_dict = updateFragmentDict(fragment_dict, i, min_mapq, cellbarcode)
    fragment_dict = filterFragmentDict(fragment_dict)
    return fragment_dict


def updateFragmentDict(fragments, segment, min_mapq, cellbarcode):
    """Update dictionary of ATAC fragments
    Takes a new aligned segment and adds information to the dictionary,
    returns a modified version of the dictionary

    Positions are 0-based
    Reads aligned to the + strand are shifted +4 bp
    Reads aligned to the - strand are shifted -5 bp

    Parameters
    ----------
    fragments : dict
        A dictionary containing ATAC fragment information
    segment : pysam.AlignedSegment
        An aligned segment
    min_mapq : int
        Minimum MAPQ to retain fragment
    cellbarcode : str
       Tag used for cell barcode. Default is CB (used by cellranger)
    """
    # because the cell barcode is not stored with each read pair (only one of the pair)
    # we need to look for each read separately rather than using the mate cigar / mate postion information
    cell_barcode, _ = utils.scan_tags(segment.tags, cb=cellbarcode)
    mapq = segment.mapping_quality
    chromosome = segment.reference_name
    qname = segment.query_name
    rstart = segment.reference_start
    rend = segment.reference_end
    qstart = segment.query_alignment_start
    is_reverse = segment.is_reverse
    if rend is None:
        return fragments
    # correct for soft clipping
    rstart = rstart + qstart
    # correct for 9 bp Tn5 shift
    if is_reverse:
        rend = rend - 5
    else:
        rstart = rstart + 4
    if mapq < min_mapq:
        return fragments
    elif qname in fragments.keys():
        if is_reverse:
            fragments[qname]["end"] = rend
        else:
            fragments[qname]["start"] = rstart
    else:
        fragments[qname] = {
            "chrom": chromosome,
            "start": None,
            "end": None,
            "cell": cell_barcode,
        }
        if is_reverse:
            fragments[qname]["end"] = rend
        else:
            fragments[qname]["start"] = rstart
    return fragments


def filterFragmentDict(fragments):
    """Remove invalid entries"""
    keepval = []
    for key, value in fragments.items():
        if all(fragments[key].values()):
            keepval.append(key)
    return {key: value for key, value in fragments.items() if key in keepval}


def fragments(
    bam,
    fragment_path,
    min_mapq=30,
    nproc=1,
    cellbarcode="CB",
    retain_cross_bc_duplicates=False,
):
    """Create ATAC fragment file from BAM file

    Iterate over reads in BAM file, extract fragment coordinates and cell barcodes.
    Collapse sequencing duplicates.

    Parameters
    ----------
    bam : str
        Path to BAM file
    fragment_path : str
        Path for output fragment file
    min_mapq : int
        Minimum MAPQ to retain fragment
    nproc : int, optional
        Number of processors to use. Default is 1.
    cellbarcode : str
       Tag used for cell barcode. Default is CB (used by cellranger)
    retain_cross_bc_duplicates : bool
        Retain fragments that have the same start and end position but originate from 
        different cell barcodes. Default is False
    """
    frags = getFragments(
        bam=bam, min_mapq=int(min_mapq), nproc=int(nproc), cellbarcode=cellbarcode
    )
    frags = collapseFragments(fragments=frags)
    writeFragments(fragments=frags, filepath=fragment_path)
