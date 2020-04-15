import pysam
from sinto import utils
from scipy import sparse
import numpy as np
from collections import Counter, defaultdict
from multiprocessing import Pool
import functools
import re
import gc
import tempfile
import os


def writeFragments(fragments, filepath):
    """Write fragments to file

    Parameters
    ----------
    fragments : list
        3-deep list of ATAC fragments
    filepath : str
        Path for output file
    """
    with open(filepath, "a") as outf:
        for i in fragments:
            outstr = "\t".join(map(str, i))
            outf.write(outstr + "\n")


def createPositionLookup(frags, pos=1):
    """Create dictionary where key is subsetted fragment coords
    (no start or no stop), value is the list of full fragments 
    that share the coordinates in the key.
    Only entries where >1 full fragments share the same coordinate
    are retained"""
    fragsplit = [x.split("|") for x in frags]
    posfrags = ["|".join([x[0], x[pos], x[3]]) for x in fragsplit]
    counts = Counter(posfrags)
    starts = dict()
    for i in range(len(frags)):
        if counts[posfrags[i]] > 1:
            try:
                starts[posfrags[i]]
            except KeyError:
                starts[posfrags[i]] = [frags[i]]
            else:
                starts[posfrags[i]].append(frags[i])
        else:
            pass
    return starts


def collapseOverlapFragments(counts, pos=1):
    """Collapse fragment counts that share a common start
    or end coordinate and are from the same cell"""
    keys = list(counts.keys())
    startfrags = createPositionLookup(frags=keys, pos=pos)
    for i in startfrags.keys():
        # extract counts for all fragments sharing start coord
        fullfrags = startfrags[i]
        countvec = [counts[x] for x in fullfrags]
        # find which fragment has the max
        # if ties, choose first in tie
        windex = countvec.index(max(countvec))
        winner = fullfrags[windex]
        # update count for the winner
        counts[winner] = sum(countvec)
        # remove losers
        del fullfrags[windex]
        for j in fullfrags:
            del counts[j]
    return counts


def collapseFragments(fragments):
    """Collapse duplicate fragments
    """
    fraglist = list(fragments.values())
    fragcoords_with_bc = ["|".join(map(str, x)) for x in fraglist]
    counts = Counter(fragcoords_with_bc)

    if len(fraglist) == 0:
        return list()

    # enumerate fragments and barcodes
    frag_id_lookup = id_lookup(l=["|".join(map(str, x[:3])) for x in fraglist])
    bc_id_lookup = id_lookup(l=[x[3] for x in fraglist])

    # collapse counts from the same cell barcode with partial overlap
    counts = collapseOverlapFragments(counts, pos=1)
    counts = collapseOverlapFragments(counts, pos=2)

    # get list of barcode index and fragment index from counts
    row = []
    col = []
    data = []
    for i in list(counts.items()):
        data.append(i[1])
        rowstr = i[0].split("|")[:3]
        bcstr = i[0].split("|")[3]
        row.append(frag_id_lookup["|".join(rowstr)])
        col.append(bc_id_lookup[bcstr])

    # free memory
    del fraglist
    del fragments
    del counts
    gc.collect()

    # create sparse matrix of fragment counts from fraglist (column, row, value)
    mat = sparse.coo_matrix(
        (data, (np.array(row), np.array(col))),
        shape=(max(frag_id_lookup.values()) + 1, max(bc_id_lookup.values()) + 1),
    )
    mat = mat.tocsr()

    # find which barcode contains the most counts for each fragment
    rowmax = np.argmax(mat, axis=1)
    rowsum = mat.sum(axis=1).tolist()

    # collapse back into a list of fragment coords and barcodes
    frag_inverse = dict(zip(frag_id_lookup.values(), frag_id_lookup.keys()))
    bc_inverse = dict(zip(bc_id_lookup.values(), bc_id_lookup.keys()))
    collapsed_frags = [frag_inverse[x] for x in np.arange(mat.shape[0])]
    collapsed_barcodes = [bc_inverse[x[0]] for x in rowmax.tolist()]
    collapsed = []
    for i in range(len(collapsed_barcodes)):
        if rowsum[i][0] < 1:
            continue
        else:
            frag = collapsed_frags[i].split("|")
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


def getFragments(
    interval,
    bam,
    min_mapq=30,
    cellbarcode="CB",
    readname_barcode=None,
    cells=None,
    max_distance=5000,
    chunksize=500000,
):
    """Extract ATAC fragments from BAM file

    Iterate over paired reads in a BAM file and extract the ATAC fragment coordinates

    Parameters
    ----------
    bam : str
        Path to BAM file
    min_mapq : int
        Minimum MAPQ to retain fragment
    cellbarcode : str
        Tag used for cell barcode. Default is CB (used by cellranger)
    readname_barcode : str, optional
        Regex to extract cell barcode from readname. If None,
        use the read tag instead.
    cells : list, optional
        List of cell barocodes to retain
    max_distance : int, optional
        Maximum distance between integration sites for the fragment to be retained.
        Allows filtering of implausible fragments that likely result from incorrect 
        mapping positions. Default is 5000 bp.
    chunksize : int
        Number of BAM entries to read through before collapsing and writing
        fragments to disk. Higher chunksize will use more memory but will be 
        faster.
    """
    fragment_dict = dict()
    inputBam = pysam.AlignmentFile(bam, "rb")
    outfile = tempfile.NamedTemporaryFile(delete=False)
    outname = outfile.name
    x = 0
    if readname_barcode is not None:
        readname_barcode = re.compile(readname_barcode)
    for i in inputBam.fetch(interval[0], 0, interval[1]):
        fragment_dict = updateFragmentDict(
            fragments=fragment_dict,
            segment=i,
            min_mapq=min_mapq,
            cellbarcode=cellbarcode,
            readname_barcode=readname_barcode,
            cells=cells,
            max_dist=max_distance,
        )
        x += 1
        if x > chunksize:
            current_position = i.reference_start
            complete = findCompleteFragments(
                fragments=fragment_dict,
                max_dist=max_distance,
                current_position=current_position,
                max_collapse_dist=20,
            )
            collapsed = collapseFragments(fragments=complete)
            writeFragments(fragments=collapsed, filepath=outname)
            x = 0
            gc.collect()
    # collapse and write the remaining fragments
    complete = findCompleteFragments(
        fragments=fragment_dict,
        max_dist=max_distance,
        current_position=i.reference_start,
        max_collapse_dist=-max_distance,  # make sure we get them all
    )
    collapsed = collapseFragments(fragments=complete)
    writeFragments(fragments=collapsed, filepath=outname)
    return outname


def findCompleteFragments(fragments, max_dist, current_position, max_collapse_dist=20):
    """Find complete fragments that are >max_dist bp away from
    the current BAM file position
    
    Parameters
    ----------
    fragments : dict
        A dictionary containing ATAC fragment information
    max_dist : int
        The maximum allowed distance between fragment start and 
        end positions
    current_position : int
        The current position being looked at in the position-sorted
        BAM file
    max_collapse_dist : int
        Maximum allowed distance for fragments from the same cell 
        barcode that share one Tn5 integration site (fragment 
        start or end coordinate) to be collapsed into a single 
        fragment.
    
    Moves completed fragments to a new dictionary
    Completed fragments will be deleted from the original dictionary
    """
    allkey = list(fragments.keys())
    completed = dict()
    d = max_dist + max_collapse_dist
    for key in allkey:
        if fragments[key][4]:  # complete fragment
            if (fragments[key][2] + d) < current_position:
                completed[key] = fragments[key][:-1]  # removes "completed" T/F information
                del fragments[key]
        else:
            # remove incomplete fragments that are
            # too far away to ever be complete
            if fragments[key][1] is None:
                if (fragments[key][2] + d) < current_position:
                    del fragments[key]
            elif fragments[key][2] is None:
                if (fragments[key][1] + d) < current_position:
                    del fragments[key]
            else:
                # start and end coordinates present without a cell barcode
                del fragments[key]
    return completed


def updateFragmentDict(
    fragments, segment, min_mapq, cellbarcode, readname_barcode, cells, max_dist
):
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
    readname_barcode : regex
        A compiled regex for matching cell barcode in read name. If None,
        use the read tags.
    cells : list
        List of cells to retain. If None, retain all cells found.
    max_dist : int
        Maximum allowed distance between fragment start and end sites
    """
    # because the cell barcode is not stored with each read pair (only one of the pair)
    # we need to look for each read separately rather than using the mate cigar / mate postion information
    if readname_barcode is not None:
        re_match = readname_barcode.match(segment.qname)
        cell_barcode = re_match.group()
    else:
        cell_barcode, _ = utils.scan_tags(segment.tags, cb=cellbarcode)
    if cells is not None and cell_barcode is not None:
        if cell_barcode not in cells:
            return fragments
    mapq = segment.mapping_quality
    if mapq < min_mapq:
        return fragments
    chromosome = segment.reference_name
    qname = segment.query_name
    rstart = segment.reference_start
    rend = segment.reference_end
    qstart = segment.query_alignment_start
    is_reverse = segment.is_reverse
    if (rend is None) or (rstart is None):
        return fragments
    # correct for soft clipping
    rstart = rstart + qstart
    # correct for 9 bp Tn5 shift
    if is_reverse:
        rend = rend - 5
    else:
        rstart = rstart + 4
    fragments = addToFragments(
        fragments, qname, chromosome, rstart, rend, cell_barcode, is_reverse, max_dist
    )
    return fragments


def addToFragments(
    fragments, qname, chromosome, rstart, rend, cell_barcode, is_reverse, max_dist
):
    """Add new fragment information to dictionary
    
    Parameters
    ----------

    fragments : dict
        A dictionary containing all the fragment information
    qname : str
        Read name
    chromosome : str
        Chromosome name
    rstart : int
        Alignment start position
    rend : int
        Alignment end position
    cell_barcode : str
        Cell barcode sequence
    is_reverse : bool
        Read is aligned to reverse strand
    max_dist : int
        Maximum allowed fragment size
    """
    if qname in fragments.keys():
        if is_reverse:
            current_coord = fragments[qname][1]
            if current_coord is None:
                # read aligned to the wrong strand, don't include
                del fragments[qname]
            elif ((rend - current_coord) > max_dist) or ((rend - current_coord) < 0):
                # too far away, don't include
                del fragments[qname]
            else:
                if cell_barcode is None and fragments[qname][3] is None:
                    # both fragment ends present but no cell barcode
                    del fragments[qname]
                else:
                    if fragments[qname][3] is None:
                        fragments[qname][3] = cell_barcode
                    fragments[qname][2] = rend
                    fragments[qname][4] = True
        else:
            current_coord = fragments[qname][2]
            if current_coord is None:
                del fragments[qname]
            elif ((current_coord - rstart) > max_dist) or (
                (current_coord - rstart) < 0
            ):
                del fragments[qname]
            else:
                if cell_barcode is None and fragments[qname][3] is None:
                    del fragments[qname]
                else:
                    if fragments[qname][3] is None:
                        fragments[qname][3] = cell_barcode
                    fragments[qname][1] = rstart
                    fragments[qname][4] = True
    else:
        # new read pair, add to dictionary
        fragments[qname] = [
            chromosome,     # chromosome    0
            None,           # start         1
            None,           # end           2
            cell_barcode,   # cell          3
            False,          # complete      4
        ]
        if is_reverse:
            fragments[qname][2] = rend
        else:
            fragments[qname][1] = rstart
    return fragments


def fragments(
    bam,
    fragment_path,
    min_mapq=30,
    nproc=1,
    cellbarcode="CB",
    chromosomes="(?i)^chr",
    readname_barcode=None,
    cells=None,
    max_distance=5000,
    chunksize=500000,
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
    chromosomes : str, optional
        Regular expression used to match chromosome names to include in the
        output file. Default is "(?i)^chr" (starts with "chr", case-insensitive).
        If None, use all chromosomes in the BAM file.
    readname_barcode : str, optional
        Regular expression used to match cell barocde stored in read name. 
        If None (default), use read tags instead. Use "[^:]*" to match all characters 
        before the first colon (":").
    cells : str
        File containing list of cell barcodes to retain. If None (default), use all cell barcodes
        found in the BAM file.
    max_distance : int, optional
        Maximum distance between integration sites for the fragment to be retained.
        Allows filtering of implausible fragments that likely result from incorrect 
        mapping positions. Default is 5000 bp.
    chunksize : int
        Number of BAM entries to read through before collapsing and writing
        fragments to disk. Higher chunksize will use more memory but will be 
        faster.
    """
    nproc = int(nproc)
    chrom = utils.get_chromosomes(bam, keep_contigs=chromosomes)
    cells = utils.read_cells(cells)
    p = Pool(nproc)
    frag_lists = [
        p.map_async(
            functools.partial(
                getFragments,
                bam=bam,
                min_mapq=int(min_mapq),
                cellbarcode=cellbarcode,
                readname_barcode=readname_barcode,
                cells=cells,
                max_distance=max_distance,
                chunksize=chunksize,
            ),
            list(chrom.items()),
        )
    ]
    filenames = [res.get() for res in frag_lists]
    # cat files and write to output
    with open(fragment_path, "w") as outfile:
        for i in filenames:
            for j in i:
                with open(j, "r") as infile:
                    for line in infile:
                        outfile.write(line)
                os.remove(j)
