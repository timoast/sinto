from sinto import utils
from sinto import filterbarcodes


@utils.log_info
def run_filterbarcodes(options):
    """Wraps the sctools.filterbarcodes function for use on the command line
    """
    filterbarcodes.filterbarcodes(
        cells=options.cells, bam=options.bam, trim_suffix=options.trim_suffix,
        output=options.output, sam=options.sam, nproc=options.nproc, mode=options.mode
        )
