import pysam
import re
from sinto import utils
from sinto.constants import OUT_FORMAT_CONVERSION
from contextlib import closing


def move(
    bam: str,
    output: str,
    cb_tag: str = "CB",
    cb_regex: str = "[^:]*",
    from_tag: bool = True,
    out_format: str = "t",
) -> None:
    """Copy cell barcode from tag to read name or vice versa
    
    Copies BAM entries to a new file. Note that cell barcodes are copied, the original tag/readname is not removed.

    Parameters
    ----------
    bam : str
        Path to BAM file.
    output : str
        Name of output BAM or SAM file
    cb_tag : str
        Name of tag containing cell barcode, or name to store cell barcode under if copying from read name
    cb_regex : str
        Regular expression used to extract cell barcode from read name
    from_tag : bool
        Copy cell barcodes from read tag to read name. If False, copy from read name to tag
    out_format : str
        One of "t" (SAM), "b" (BAM), or "u" (uncompressed BAM) ("t" default)
    """
    assert len(cb_tag) == 2
    
    inputBam = pysam.AlignmentFile(bam)
    outputBam = pysam.AlignmentFile(
        output, "w" + OUT_FORMAT_CONVERSION[out_format], template=inputBam
    )
    readname_barcode = re.compile(cb_regex)
    with closing(inputBam) as inputBam:
        with closing(outputBam) as outputBam:
            for r in inputBam.fetch(until_eof=True):
                tags = r.get_tags()
                if from_tag:
                    cell_barcode, _ = utils.scan_tags(tags, cb=cb_tag)
                    r.query_name = cell_barcode + ":" + r.query_name
                else:
                    re_match = readname_barcode.match(r.query_name)
                    cell_barcode = re_match.group()
                    r.set_tag(cb_tag, cell_barcode, value_type="Z", replace=True)
                outputBam.write(r)
