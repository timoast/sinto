from sinto import filterbarcodes
import functools
import time


def log_info(func):
    """Decorator that prints function arguments and runtime
    """
    @functools.wraps(func)
    def wrapper(args):
        print("Function {} called with the following arguments:\n".format(func.__name__))
        for arg in vars(args):
            print(str(arg) + '\t' + str(getattr(args, arg)))
        t1 = time.time()
        func(args)
        t2 = time.time()
        elapsed = [round(x, 2) for x in divmod(t2-t1, 60)]
        print("\nFunction completed in  {} m {} s\n".format(elapsed[0], elapsed[1]))
    return wrapper


@log_info
def run_filterbarcodes(options):
    """Wraps the sctools.filterbarcodes function for use on the command line
    """
    filterbarcodes(cells=options.cells, bam=options.bam, trim_suffix=options.trim_suffix,
                   output=options.output, sam=options.sam, nproc=options.nproc)