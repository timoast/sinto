from sinto import utils
from sinto import filterbarcodes
from sinto import addtags


@utils.log_info
def run_filterbarcodes(options):
    """Wraps the sctools.filterbarcodes function for use on the command line
    """
    filterbarcodes.filterbarcodes(
        cells=options.cells,
        bam=options.bam,
        trim_suffix=options.trim_suffix,
        output=options.output,
        sam=options.sam,
        nproc=options.nproc,
        mode=options.mode,
    )


@utils.log_info
def run_addtags(options):
    """Wraps the sctools.addtags function for use on the command line
    """
    addtags.addtags(
        bam=options.bam,
        tagfile=options.tagfile,
        trim_suffix=options.trim_suffix,
        output=options.output,
        sam=options.sam,
        nproc=options.nproc,
        mode=options.mode,
    )
