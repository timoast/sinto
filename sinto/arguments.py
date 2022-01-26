import sys
from argparse import ArgumentParser
import pkg_resources
from sinto import cli


version = pkg_resources.require("sinto")[0].version
parser = ArgumentParser(description="Tools for single-cell data processing")
parser.add_argument(
    "-v", "--version", action="version", version="%(prog)s " + str(version)
)
subparsers = parser.add_subparsers(title="Subcommands")

# filterbarcodes
parser_filterbarcodes = subparsers.add_parser(
    "filterbarcodes", description="Filter reads based on input list of cell barcodes"
)
parser_filterbarcodes.add_argument(
    "-b", "--bam", help="Input bam file (must be indexed)", required=True, type=str
)
parser_filterbarcodes.add_argument(
    "-c",
    "--cells",
    help="File or comma-separated list of cell barcodes. Can be gzip compressed",
    required=True,
    type=str,
)
parser_filterbarcodes.add_argument(
    "-t",
    "--trim_suffix",
    help="Remove trail 2 characters from cell barcode in BAM file",
    action="store_true",
    default=False,
)
parser_filterbarcodes.add_argument(
    "-p",
    "--nproc",
    help="Number of processors (default = 1)",
    required=False,
    default=1,
    type=int,
)
parser_filterbarcodes.add_argument(
    "--barcode_regex",
    help="""
    Regular expression used to extract cell barcode from read name. If None (default), extract cell barcode from read tag.
    Use "[^:]*" to match all characters up to the first colon.
    """,
    required=False,
    type=str,
    default=None,
)
parser_filterbarcodes.add_argument(
    "--barcodetag",
    help='Read tag storing cell barcode information (default = "CB")',
    required=False,
    type=str,
    default="CB",
)
parser_filterbarcodes.set_defaults(func=cli.run_filterbarcodes)

# addtags
parser_addtags = subparsers.add_parser(
    "addtags", description="Add read tags to reads from individual cells"
)
parser_addtags.add_argument(
    "-b", "--bam", help="Input bam file (must be indexed)", required=True, type=str
)
parser_addtags.add_argument(
    "-f",
    "--tagfile",
    help="Tab-delimited file containing cell barcode, tag to be added, and tag identity. Can be gzip compressed",
    required=True,
    type=str,
)
parser_addtags.add_argument(
    "-o", "--output", help="Name for output text file", required=True, type=str
)
parser_addtags.add_argument(
    "-t",
    "--trim_suffix",
    help="Remove trail 2 characters from cell barcode in BAM file",
    action="store_true",
    default=False,
)
parser_addtags.add_argument(
    "-s",
    "--sam",
    help="Output sam format (default bam output)",
    required=False,
    action="store_true",
    default=False,
)
parser_addtags.add_argument(
    "-p",
    "--nproc",
    help="Number of processors (default = 1)",
    required=False,
    default=1,
    type=int,
)
parser_addtags.add_argument(
    "-m",
    "--mode",
    help="Either tag (default) or readname. Some BAM file store the cell barcode in the readname rather than under a read tag",
    required=False,
    default="tag",
    type=str,
)
parser_addtags.set_defaults(func=cli.run_addtags)

# tagtorg
parser_tagtorg = subparsers.add_parser(
    "tagtorg",
    description="Append a read tag to the read group ID of each read. "
    "Also appends the read tag to the SM field of the read group.",
)
parser_tagtorg.add_argument(
    "-b",
    "--bam",
    help="Input SAM/BAM file, '-' reads from stdin",
    required=True,
    type=str,
)
parser_tagtorg.add_argument(
    "--tag",
    help="Read tag to extract the value from that is appended to the read group. "
    "Default is 'CB', the tag that is used in 10x sequencing to identify cells.",
    default="CB",
    type=str,
)
parser_tagtorg.add_argument(
    "-f",
    "--tagfile",
    help="List of expected tag values. Reads with tag values that are not in this list are not altered.",
    required=True,
    type=str,
)
parser_tagtorg.add_argument(
    "-o",
    "--output",
    help="Output SAM/BAM file, '-' outputs to stdout (default '-')",
    default="-",
    type=str,
)
parser_tagtorg.add_argument(
    "-O",
    "--outputformat",
    help="Output format. One of 't' (SAM), 'b' (BAM),"
    " or 'u' (uncompressed BAM) ('t' default)",
    default="t",
)


parser_tagtorg.set_defaults(func=cli.run_tagtorg)

# tagtotag
parser_tagtotag = subparsers.add_parser(
    "tagtotag",
    description="Copies BAM entries to a new file while copying a read tag"
                " to another read tag"
                " and optionally deleting the originating tag.",
)
parser_tagtotag.add_argument(
    "-b",
    "--bam",
    help="Input SAM/BAM file, '-' reads from stdin",
    required=True,
    type=str,
)
parser_tagtotag.add_argument(
    "--from",
    help="Read tag to copy from.",
    required=True,
    dest='from_',
    type=str,
)
parser_tagtotag.add_argument(
    "--to",
    help="Read tag to copy to.",
    required=True,
    type=str,
)
parser_tagtotag.add_argument(
    "--delete",
    help="Delete originating tag after copy (i.e. move).",
    action='store_true',
)
parser_tagtotag.add_argument(
    "-o",
    "--output",
    help="Output SAM/BAM file, '-' outputs to stdout (default '-')",
    default="-",
    type=str,
)
parser_tagtotag.add_argument(
    "-O",
    "--outputformat",
    help="Output format. One of 't' (SAM), 'b' (BAM),"
    " or 'u' (uncompressed BAM) ('t' default)",
    default="t",
)


parser_tagtotag.set_defaults(func=cli.run_tagtotag)

# fragments
parser_fragments = subparsers.add_parser(
    "fragments", description="Create ATAC-seq fragment file from BAM file"
)
parser_fragments.add_argument(
    "-b", "--bam", help="Input bam file (must be indexed)", required=True, type=str
)
parser_fragments.add_argument(
    "-f",
    "--fragments",
    help="""
    Name and path for output fragments file. Note that the output is not sorted or compressed.
    To sort the output file use sort -k 1,1 -k2,2n
    """,
    required=True,
    type=str,
)
parser_fragments.add_argument(
    "-m",
    "--min_mapq",
    help="Minimum MAPQ required to retain fragment (default = 30)",
    required=False,
    default=30,
    type=int,
)
parser_fragments.add_argument(
    "-p",
    "--nproc",
    help="Number of processors (default = 1)",
    required=False,
    type=int,
    default=1,
)
parser_fragments.add_argument(
    "-t",
    "--barcodetag",
    help='Read tag storing cell barcode information (default = "CB")',
    required=False,
    type=str,
    default="CB",
)
parser_fragments.add_argument(
    "-c",
    "--cells",
    help="""
    Path to file containing cell barcodes to retain, or a
    comma-separated list of cell barcodes.
    If None (default), use all cell barocodes present in the BAM file.
    """,
    required=False,
    type=str,
    default=None,
)
parser_fragments.add_argument(
    "--barcode_regex",
    help="""
    Regular expression used to extract cell barcode from read name. If None (default), extract cell barcode from read tag.
    Use "[^:]*" to match all characters up to the first colon.
    """,
    required=False,
    type=str,
    default=None,
)
parser_fragments.add_argument(
    "--use_chrom",
    help="""
    Regular expression used to match chromosomes to be included in output.
    Default is "(?i)^chr" to match all chromosomes starting with "chr",
    case insensitive
    """,
    required=False,
    default="(?i)^chr",
    type=str,
)
parser_fragments.add_argument(
    "--max_distance",
    help="""
    Maximum distance between integration sites for the fragment to be retained.
    Allows filtering of implausible fragments that likely result from incorrect
    mapping positions. Default is 5000 bp.
    """,
    required=False,
    default=5000,
    type=int,
)
parser_fragments.add_argument(
    "--min_distance",
    help="""
    Minimum distance between integration sites for the fragment to be retained.
    Allows filtering of implausible fragments that likely result from incorrect
    mapping positions. Default is 10 bp.
    """,
    required=False,
    default=10,
    type=int,
)
parser_fragments.add_argument(
    "--chunksize",
    help="""
    Number of BAM file entries to iterate over before collapsing the fragments and writing
    to disk. Higher chunksize will use more memory but will be faster.
    """,
    required=False,
    default=500000,
    type=int,
)
parser_fragments.add_argument(
    "--shift_plus",
    help="""
    Number of bases to shift Tn5 insertion position by on the forward strand.
    """,
    required=False,
    default=4,
    type=int,
)
parser_fragments.add_argument(
    "--shift_minus",
    help="""
    Number of bases to shift Tn5 insertion position by on the reverse strand.
    """,
    required=False,
    default=-5,
    type=int,
)
parser_fragments.add_argument(
    "--collapse_within",
    help="""
    Take cell barcode into account when collapsing duplicate fragments. Setting this
    flag means that fragments with the same coordinates can be identified provided 
    they originate from a different cell barcode.
    """,
    action="store_true",
    default=False,
)
parser_fragments.set_defaults(func=cli.run_fragments)

# barcodes
parser_barcode = subparsers.add_parser(
    "barcode", description="Add cell barcode sequences to read names in FASTQ file."
)
parser_barcode.add_argument(
    "--barcode_fastq", help="FASTQ file containing cell barcode sequences", required=True, type=str
)
parser_barcode.add_argument(
    "--read1", help="FASTQ file containing read 1", required=True, type=str
)
parser_barcode.add_argument(
    "--read2", help="FASTQ file containing read 2", required=False, type=str, default=None
)
parser_barcode.add_argument(
    "-b", "--bases", help="Number of bases to extract from barcode-containing FASTQ", required=True, type=int
)
parser_barcode.add_argument(
    "--prefix", help="Prefix to add to cell barcodes", required=False, type=str, default=""
)
parser_barcode.add_argument(
    "--suffix", help="Suffix to add to cell barcodes", required=False, type=str, default=""
)
parser_barcode.set_defaults(func=cli.run_barcode)

# tagtoname
parser_tagtoname = subparsers.add_parser(
    "tagtoname",
    description="""
    Copy cell barcode sequences from tag to read names. Cell barcodes will be added
    as a readname prefix, followed by ":"
    """
)
parser_tagtoname.add_argument(
    "-b",
    "--bam",
    help="Input SAM/BAM file, '-' reads from stdin",
    required=True,
    type=str,
)
parser_tagtoname.add_argument(
    "-o",
    "--output",
    help="Output SAM/BAM file, '-' outputs to stdout (default '-')",
    default="-",
    type=str,
)
parser_tagtoname.add_argument(
    "-O",
    "--outputformat",
    help="Output format. One of 't' (SAM), 'b' (BAM),"
    " or 'u' (uncompressed BAM) ('t' default)",
    default="t",
)
parser_tagtoname.add_argument(
    "--tag",
    help="Read tag to copy from.",
    required=False,
    default="CB",
    type=str,
)
parser_tagtoname.set_defaults(func=cli.run_tagtoname)

# nametotag
parser_nametotag = subparsers.add_parser(
    "nametotag", description="Copy cell barcode sequences from read name to read tag"
)
parser_nametotag.add_argument(
    "-b",
    "--bam",
    help="Input SAM/BAM file, '-' reads from stdin",
    required=True,
    type=str,
)
parser_nametotag.add_argument(
    "-o",
    "--output",
    help="Output SAM/BAM file, '-' outputs to stdout (default '-')",
    default="-",
    type=str,
)
parser_nametotag.add_argument(
    "-O",
    "--outputformat",
    help="Output format. One of 't' (SAM), 'b' (BAM),"
    " or 'u' (uncompressed BAM) ('t' default)",
    default="t",
)
parser_nametotag.add_argument(
    "--barcode_regex",
    help="""
    Regular expression used to extract cell barcode from read name.
    Default ("[^:]*") matches all characters up to the first colon.
    """,
    required=False,
    type=str,
    default="[^:]*",
)
parser_nametotag.add_argument(
    "--tag",
    help="Read tag to copy to.",
    required=False,
    default="CB",
    type=str,
)
parser_nametotag.set_defaults(func=cli.run_nametotag)


def main():
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    else:
        options = parser.parse_args()
        sys.setrecursionlimit(200000)
        options.func(options)

if __name__ == "__main__":
    main()
