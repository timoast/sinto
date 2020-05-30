import subprocess

import pytest


@pytest.mark.parametrize("delete", [True, False])
def test_copies_CB_tag_to_cb(tmpdir, delete):

    # given
    # modified from an example in the SAM specification v1
    sam = (
        "@HD VN:1.5 SO:coordinate\n"
        "@SQ SN:20 LN:63025520\n"
        "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:AAAA-1\n"
    ).replace(" ", "\t")
    in_sam_file = tmpdir / "in.sam"
    in_tag_file = tmpdir / "tags.txt"
    output = tmpdir / "output.sam"
    with open(in_sam_file, "w") as fh:
        fh.write(sam)
    with open(in_tag_file, "w") as fh:
        fh.write("AAAA-1")

    # when
    cmd = f"sinto tagtotag --from CB --to xx --bam - -o -"
    if delete:
        cmd += " --delete"
    cmd += f" < {in_sam_file} > {output}"
    subprocess.run(cmd, shell=True)

    # then
    out_sam = open(output).read()
    if delete:
        expected = (
            "@HD VN:1.5 SO:coordinate\n"
            "@SQ SN:20 LN:63025520\n"
            "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * xx:Z:AAAA-1\n"
        )

    else:
        expected = (
            "@HD VN:1.5 SO:coordinate\n"
            "@SQ SN:20 LN:63025520\n"
            "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:AAAA-1 xx:Z:AAAA-1\n"
        )
    expected = expected.replace(" ", "\t")
    assert expected == out_sam

