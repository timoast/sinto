import pytest
import os
from sinto import addtags


@pytest.mark.parametrize("nproc", [1, 2])
def test_addtags(tmpdir, nproc):
    basepath=os.path.dirname(os.path.realpath(__file__))
    bam=os.path.join(basepath, "data/test.bam")
    tags=os.path.join(basepath, "data/tags.tsv")
    outf=os.path.join(tmpdir, "tags.sam")
    addtags.addtags(
        bam=bam,
        tagfile=tags,
        trim_suffix=False,
        output=outf,
        sam=True,
        nproc=nproc,
        mode="tag"
    )
    
    assert os.path.exists(outf)

    with open(outf, 'r') as inp:
        for l in inp:
            if l.startswith("@"):
                next
            else:
                break
    with open(os.path.join(basepath, "data/tags_line1.sam"), 'r') as inp:
        for l2 in inp:
            next
    assert l == l2