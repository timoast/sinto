import pytest
import os
from sinto import filterbarcodes


@pytest.mark.parametrize("nproc", [1, 2])
@pytest.mark.parametrize("barcodes", ["data/barcodes.tsv", "data/barcodes.tsv.gz"])
def test_filterbarcodes(tmpdir, nproc, barcodes):
    basepath=os.path.dirname(os.path.realpath(__file__))
    bam=os.path.join(basepath, "data/test.bam")
    barcodes=os.path.join(basepath, barcodes)

    filterbarcodes.filterbarcodes(
        cells=barcodes,
        bam=bam,
        trim_suffix=False,
        nproc=nproc,
        readname_barcode=None,
        cellbarcode="CB",
        outdir=tmpdir,
        sam=True
    )

    exists = [os.path.exists(os.path.join(tmpdir, x)) for x in ['A.sam', 'B.sam', 'C.sam']]
    assert all(exists)
    
    # check contents
    with open(os.path.join(tmpdir, 'A.sam'), 'r') as inp:
        for l in inp:
            if l.startswith("@"):
                next
            else:
                break
    with open(os.path.join(basepath, "data/a_line1.sam"), 'r') as inp:
        for l2 in inp:
            next

    assert l == l2