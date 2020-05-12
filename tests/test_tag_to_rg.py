import subprocess


def test_adds_CB_tag_to_rg(tmpdir):

    # given
    # modified from an example in the SAM specification v1
    sam = (
        "@HD VN:1.5 SO:coordinate\n"
        "@SQ SN:20 LN:63025520\n"
        "@RG ID:rg1 SM:sample_1 LB:1 PU:1 PL:ILLUMINA\n"
        "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:AAAA-1 RG:Z:rg1\n"
    ).replace(" ", "\t")
    in_sam_file = tmpdir / "in.sam"
    in_tag_file = tmpdir / "tags.txt"
    output = tmpdir/'output.sam'
    with open(in_sam_file, "w") as fh:
        fh.write(sam)
    with open(in_tag_file, "w") as fh:
        fh.write("AAAA-1")

    # when
    subprocess.run(
        f"sinto tagtorg --tagfile {in_tag_file} --tag CB -b - -o - < {in_sam_file} > {output}",
        shell=True,
    )

    # then
    out_sam = open(output).read()
    expected = (
        "@HD VN:1.5 SO:coordinate\n"
        "@SQ SN:20 LN:63025520\n"
        "@RG ID:rg1 SM:sample_1 LB:1 PU:1 PL:ILLUMINA\n"
        "@RG ID:rg1:AAAA-1 SM:sample_1:AAAA-1 LB:1 PU:1 PL:ILLUMINA\n"
        "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:AAAA-1 RG:Z:rg1:AAAA-1\n"
    ).replace(" ", "\t")
    assert expected == out_sam

def test_ignores_unknown_cb_tag(tmpdir):

    # given
    # modified from an example in the SAM specification v1
    sam = (
        "@HD VN:1.5 SO:coordinate\n"
        "@SQ SN:20 LN:63025520\n"
        "@RG ID:rg1 SM:sample_1 LB:1 PU:1 PL:ILLUMINA\n"
        "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:AAAA-1 RG:Z:rg1\n"
        "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:CCCC-1 RG:Z:rg1\n"
    ).replace(" ", "\t")
    in_sam_file = tmpdir / "in.sam"
    in_tag_file = tmpdir / "tags.txt"
    output = tmpdir/'output.sam'
    with open(in_sam_file, "w") as fh:
        fh.write(sam)
    with open(in_tag_file, "w") as fh:
        fh.write("AAAA-1")

    # when
    subprocess.run(
        f"sinto tagtorg --tagfile {in_tag_file} --tag CB -b - -o - < {in_sam_file} > {output}",
        shell=True,
    )

    # then
    out_sam = open(output).read()
    expected = (
        "@HD VN:1.5 SO:coordinate\n"
        "@SQ SN:20 LN:63025520\n"
        "@RG ID:rg1 SM:sample_1 LB:1 PU:1 PL:ILLUMINA\n"
        "@RG ID:rg1:AAAA-1 SM:sample_1:AAAA-1 LB:1 PU:1 PL:ILLUMINA\n"
        "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:AAAA-1 RG:Z:rg1:AAAA-1\n"
        "r002 0 20 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA * CB:Z:CCCC-1 RG:Z:rg1\n"
    ).replace(" ", "\t")
    assert expected == out_sam
