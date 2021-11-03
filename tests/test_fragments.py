import pytest
import os
from sinto import fragments


@pytest.mark.parametrize("collapse_within", [True, False])
def test_fragments(tmpdir, collapse_within):
    basepath=os.path.dirname(os.path.realpath(__file__))
    bam=os.path.join(basepath, "data/test.bam")
    outf=tmpdir / "frags.bed"
    fragments.fragments(
        bam=bam,
        fragment_path=outf,
        collapse_within=collapse_within,
    )
    output = [line.strip("\n") for line in open(outf, "r")]    
    if collapse_within:
        expected_path = os.path.join(basepath, "data/frags_within.bed")
    else:
        expected_path = os.path.join(basepath, "data/frags_across.bed")
    expected = [line.strip("\n") for line in open(expected_path, "r")]
    output.sort()
    expected.sort()
    assert output == expected
