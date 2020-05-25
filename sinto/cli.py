from sinto import utils, filterbarcodes, addtags, fragments, tagtorg, addbarcodes


@utils.log_info
def run_filterbarcodes(options):
    """Wraps the sctools.filterbarcodes function for use on the command line
    """
    filterbarcodes.filterbarcodes(
        cells=options.cells,
        bam=options.bam,
        trim_suffix=options.trim_suffix,
        nproc=options.nproc,
        readname_barcode=options.barcode_regex,
        cellbarcode=options.barcodetag,
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


@utils.log_info
def run_fragments(options):
    """Wraps the sctools.fragments function for use on the command line
    """
    fragments.fragments(
        bam=options.bam,
        fragment_path=options.fragments,
        min_mapq=options.min_mapq,
        nproc=options.nproc,
        cellbarcode=options.barcodetag,
        readname_barcode=options.barcode_regex,
        chromosomes=options.use_chrom,
        cells=options.cells,
        max_distance=options.max_distance,
        chunksize=options.chunksize,
    )

@utils.log_info
def run_tagtorg(options):
    tagtorg.tagtorg(
        bam=options.bam,
        tag_value_file=options.tagfile,
        tag=options.tag,
        output=options.output,
        out_format=options.outputformat,
    )

@utils.log_info
def run_barcode(options):
    addbarcodes.addbarcodes(
        cb_position=options.bases,
        fq1=options.barcode_fastq,
        fq2=options.read1,
        fq3=options.read2,
        prefix=options.prefix,
        suffix=options.suffix
    )
