import sys
from contextlib import closing

import pysam


def header_line_to_str(line):
    return "\t".join(f"{k}:{v}" for k, v in line.items())


def build_header(header, tag_vals):
    new_rg_lines = []
    for val in tag_vals:
        for rg in header["RG"]:
            rg = rg.copy()
            rg["SM"] += f":{val}"
            rg["ID"] += f":{val}"
            new_rg_lines.append("@RG\t" + header_line_to_str(rg))
    return str(header) + "\n".join(new_rg_lines) + "\n"


OUT_FORMAT_CONVERSION = {
    "t": "",
    "b": "b",
    "u": "bu",
}


def tagtorg(bam, tag, output, tag_value_file, out_format="t"):
    """Add tags to reads from individual cells

    Copies BAM entries to a new file, adding a read tag to cells matching an input table

    Parameters
    ----------
    bam : str
        Path to SAM/BAM file, or "-" to read from stdin.
    tag_value_file : str
        Text file with one line per expected tag value.
    output : str
        Name for output SAM/BAM file, or "-" to write to stdout.
    out_format : str
        One of "t" (SAM), "b" (BAM), or "u" (uncompressed BAM) ("t" default)
    """

    infile = pysam.AlignmentFile(bam, "r")
    with open(tag_value_file) as fh:
        tag_vals = {val.rstrip() for val in fh.readlines()}
    new_header = build_header(infile.header, tag_vals)

    outfile = pysam.AlignmentFile(
        output, "w" + OUT_FORMAT_CONVERSION[out_format], text=new_header
    )
    with closing(outfile) as outfile:
        for rec in infile:
            try:
                tag_val = rec.get_tag(tag)
            except KeyError:
                pass
            else:
                if tag_val in tag_vals:
                    rec.set_tag(
                        "RG",
                        f"{rec.get_tag('RG')}:{tag_val}",
                        value_type="Z",
                        replace=True,
                    )
            finally:
                outfile.write(rec)
