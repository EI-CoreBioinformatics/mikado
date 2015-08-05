#!/usr/bin/env python3
# coding: utf-8

"""Launcher of the Mikado pipeline."""

import argparse
import sys
import mikado_lib.loci_objects


if "line_profiler" not in dir():
    def profile(function):
        """
        Mock wrapper to imitate the profile decorator
        :param function: the function to be wrapped
        :return:
        """
        def inner(*args, **kwargs):
            """
            Returns the wrapped function
            :param args: arguments to be passed
            :param kwargs: keyword arguments to be passed
            :return:
            """
            return function(*args, **kwargs)
        return inner


@profile
def main():
    """
    Main script function.
    """
    parser = argparse.ArgumentParser("Launcher of the Mikado pipeline.")
    parser.add_argument("-p", "--procs", type=int, default=None,
                        help="Number of processors to use. Default: look in the configuration file (1 if undefined)")
    parser.add_argument("--json_conf", type=argparse.FileType(), required=True,
                        help="JSON/YAML configuration file for scoring transcripts.")
    parser.add_argument("--subloci_out", type=str, default=None)
    parser.add_argument("--monoloci_out", type=str, default=None)
    parser.add_argument("--loci_out", type=str, default=None,
                        help="""This output file is mandatory.
                        If it is not specified in the configuration file, it must be provided here.""")
    parser.add_argument("--no_cds", action="store_true", default=None,
                        help="Flag. If set, not CDS information will be printed out in the GFF output files.")
    parser.add_argument('--source', type=str, default=None,
                        help='Source field to use for the output files.')
    parser.add_argument('--purge', action='store_true', default=False,
                        help='''Flag. If set, the pipeline will suppress any loci whose transcripts
                        do not pass the requirements set in the JSON file.'''
                        )
    parser.add_argument('--cache', action='store_true', default=False,
                        help='''Flag. If set, the Mikado DB will be pre-loaded into memory for faster access.
                        WARNING: this option will increase memory usage and the preloading might be quite slow.'''
                        )
    parser.add_argument("--single", action="store_true", default=False,
                        help="""Flag. If set, Creator will be launched with a single process.
                        Useful for debugging purposes only.""")
    log_options = parser.add_argument_group("Log options")
    log_options.add_argument("-l", "--log", default=None,
                             help="File to write the log to. Default: decided by the configuration file.")
    verbosity = log_options.add_mutually_exclusive_group()
    verbosity.add_argument("--verbose", action="store_true",
                           default=False, help="Flag. If set, the debug mode will be activated.")
    verbosity.add_argument("--noverbose", action="store_true",
                           default=False, help="Flag. If set, the debug mode will be activated.")
    log_options.add_argument("-lv", "--log-level", dest="log_level",
                             choices=["DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"], default=None,
                             help="Logging level. Default: retrieved by the configuration file.")
    parser.add_argument("gff", type=argparse.FileType(), nargs="?", default=None)

    args = parser.parse_args()

    args.json_conf.close()
    args.json_conf = mikado_lib.json_utils.to_json(args.json_conf.name)

    if args.procs is not None:
        args.json_conf["run_options"]["threads"] = args.procs

    if args.log == "stderr":
        args.json_conf["log_settings"]['log'] = None
    elif args.log is not None:
        args.json_conf["log_settings"]['log'] = args.log

    if args.cache is True:
        args.json_conf["run_options"]["preload"] = True

    args.json_conf["single_thread"] = args.single

    if args.log_level is not None:
        args.json_conf["log_settings"]['log_level'] = args.log_level
    elif args.verbose is True:
        args.json_conf["log_settings"]['log_level'] = "DEBUG"
    elif args.noverbose is True:
        args.json_conf["log_settings"]['log_level'] = "ERROR"

    if args.monoloci_out is not None:
        args.json_conf["monoloci_out"] = args.monoloci_out
    if args.subloci_out is not None:
        args.json_conf["subloci_out"] = args.subloci_out
    if args.loci_out is not None:
        args.json_conf["loci_out"] = args.loci_out
    if args.no_cds is not None:
        args.json_conf["run_options"]["exclude_cds"] = True
    if args.source is not None:
        args.json_conf["source"] = args.source
    if args.purge is not None:
        args.json_conf["run_options"]["purge"] = True

    if args.gff is not None:
        args.gff.close()
        args.gff = args.gff.name
        args.json_conf["input"] = args.gff

    creator = mikado_lib.loci_objects.Creator.Creator(args.json_conf, commandline=" ".join(sys.argv))
    creator()  # Run


if __name__ == "__main__":
    main()
