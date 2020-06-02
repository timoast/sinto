from contextlib import closing

import pysam

from sinto.constants import OUT_FORMAT_CONVERSION


def tagtotag(
    bam: str,
    output: str,
    from_tag: str,
    to_tag: str,
    delete: bool,
    out_format: str = "t",
) -> None:
    """Copy/move tags

    Copies BAM entries to a new file while copying a read tag to another read tag
    and optionally deleting the originating tag.

    Parameters
    ----------
    bam : str
        Path to SAM/BAM file, or "-" to read from stdin.
    output : str
        Name for output SAM/BAM file, or "-" to write to stdout.
    from_tag: str
        Read tag to copy from.
    to_tag: str
        Read tag to copy to
    delete: bool
        Delete from_tag from read if True.
    out_format : str
        One of "t" (SAM), "b" (BAM), or "u" (uncompressed BAM) ("t" default)
    """
    assert len(from_tag) == 2
    assert len(to_tag) == 2

    infile = pysam.AlignmentFile(bam)
    outfile = pysam.AlignmentFile(
        output, "w" + OUT_FORMAT_CONVERSION[out_format], template=infile
    )

    for rec in infile:
        try:
            tag_val, tag_type = rec.get_tag(from_tag, with_value_type=True)
        except KeyError:
            pass
        else:
            rec.set_tag(to_tag, tag_val, value_type=tag_type, replace=True)
            if delete:
                rec.set_tag(from_tag, None)
        finally:
            outfile.write(rec)
