import argparse
import sys
from multiprocessing.spawn import freeze_support
import logging
from Mikado.version import __version__


def main(call_args=None):

    """
    Main launcher function for the pipeline.
    :param call_args: optional argument string to be passed to execute the commands.
    Otherwise, the string will be derived from sys.argv[1:]
    """

    if call_args is None:
        call_args = sys.argv[1:]
    from Mikado.subprograms import configure, prepare, serialise, pick, compare, util

    parser = argparse.ArgumentParser(prog="Mikado",
                                     description="""Mikado is a program to analyse RNA-Seq data
and determine the best transcript for each locus in accordance to user-specified criteria.""")

    parser.add_argument("--version", default=False, action="store_true",
                        help="Print Mikado current version and exit.")

    subparsers = parser.add_subparsers(
        title="Components",
        help="""These are the various components of Mikado:

""")
    subparsers.add_parser("configure",
                          help="This utility guides the user through \
the process of creating a configuration file for Mikado.")
    subparsers.choices["configure"] = configure.configure_parser()
    subparsers.choices["configure"].prog = "Mikado configure"

    subparsers.add_parser("prepare", help ="Mikado prepare analyses an input GTF file and \
prepares it for the picking analysis by sorting its transcripts \
and performing some simple consistency checks.")
    subparsers.choices["prepare"] = prepare.prepare_parser()
    subparsers.choices["prepare"].prog = "Mikado prepare"

    subparsers.add_parser("serialise", help="Mikado serialise creates the database used \
by the pick program. It handles Junction and ORF BED12 files as well as BLAST XML results.")
    subparsers.choices["serialise"] = serialise.serialise_parser()
    subparsers.choices["serialise"].prog = "Mikado serialise"

    subparsers.add_parser("pick", help="Mikado pick analyses a sorted GTF/GFF files in order \
to identify its loci and choose the best transcripts according to user-specified criteria. \
It is dependent on files produced by the \"prepare\" and \"serialise\" components.")
    subparsers.choices["pick"] = pick.pick_parser()
    subparsers.choices["pick"].prog = "Mikado pick"

    subparsers.add_parser("compare", help="Mikado compare produces a detailed comparison of \
reference and prediction files. It has been directly inspired \
by Cufflinks's cuffcompare and ParsEval.")
    subparsers.choices["compare"] = compare.compare_parser()
    subparsers.choices["compare"].prog = "Mikado compare"

    subparsers.add_parser("util", help="Miscellaneous utilities")
    subparsers.choices["util"] = util.util_parser()
    subparsers.choices["util"].prog = "Mikado util"

    try:
        args = parser.parse_args(call_args)
        if hasattr(args, "func"):
            args.func(args)
        elif len(call_args) > 0 and call_args[0] == "util":
            util.util_parser().print_help()
        elif args.version is True:
            print("Mikado v{}".format(__version__))
            sys.exit(0)
        else:
            parser.print_help()
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except BrokenPipeError:
        pass
    except Exception as exc:
        logger = logging.getLogger("main")
        logger.error("Mikado crashed, cause:")
        logger.exception(exc)
        import multiprocessing as mp
        for child in mp.active_children():
            child.terminate()

        sys.exit(1)


if __name__ == '__main__':
    # __spec__ = "Mikado"
    freeze_support()
    sys.exit(main())
