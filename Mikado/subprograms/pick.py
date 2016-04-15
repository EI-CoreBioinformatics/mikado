#!/usr/bin/env python3
# coding: utf-8

"""Launcher of the Mikado pick step."""

import argparse
import sys
from ..picking import Picker
from ..configuration.configurator import to_json
# from . import to_gff
from ..exceptions import UnsortedInput


def check_log_settings(args):

    """
    Quick method to check the consistency of log settings
    from the namespace.
    :param args: a Namespace
    :return: args
    """

    if args.log == "stderr":
        args.json_conf["log_settings"]['log'] = None
    elif args.log is not None:
        args.json_conf["log_settings"]['log'] = args.log

    if args.log_level is not None:
        args.json_conf["log_settings"]['log_level'] = args.log_level
    elif args.verbose is True:
        args.json_conf["log_settings"]['log_level'] = "DEBUG"
    elif args.noverbose is True:
        args.json_conf["log_settings"]['log_level'] = "ERROR"

    return args


def check_run_options(args):
    """
    Quick method to check the consistency of run option settings
    from the namespace.
    :param args: a Namespace
    :return: args
    """

    if args.start_method is not None:
        args.json_conf["multiprocessing_method"] = args.start_method

    if args.procs is not None:
        args.json_conf["pick"]["run_options"]["procs"] = args.procs

    if args.shm_db is not None or args.shm is True:
        args.shm = True
        args.json_conf["pick"]["run_options"]["shm"] = True
        # I will deal with it being None or not in Creator
        args.json_conf["pick"]["run_options"]["shm_db"] = args.shm_db

    if args.preload is True:
        args.json_conf["pick"]["run_options"]["preload"] = True

    args.json_conf["pick"]["run_options"]["single_thread"] = args.single

    if args.no_cds is not None:
        args.json_conf["pick"]["run_options"]["exclude_cds"] = True
    if args.purge is not None:
        args.json_conf["pick"]["run_options"]["purge"] = True

    if args.output_dir is not None:
        args.json_conf["pick"]["files"]["output_dir"] = args.output_dir

    if args.source is not None:
        args.json_conf["pick"]["output_format"]["source"] = args.source
    if args.prefix is not None:
        args.json_conf["pick"]["output_format"]["id_prefix"] = args.prefix

    if args.sqlite_db is not None:
        args.json_conf["db_settings"]["db"] = args.sqlite_db
        args.json_conf["db_settings"]["dbtype"] = "sqlite"

    for key in ["loci_out", "gff", "monoloci_out", "subloci_out", "log"]:
        if getattr(args, key):
            if key == "gff":
                args.json_conf["pick"]["files"]["input"] = getattr(
                    args,
                    key,
                    args.json_conf["pick"]["files"]["input"])
            else:
                val = getattr(args, key, args.json_conf["pick"]["files"][key])
                if key in ("loci_out", "monoloci_out", "subloci_out"):
                    if val.split(".")[-1] not in ("gff3", "gff"):
                        val = "{0}.gff3".format(val)

                args.json_conf["pick"]["files"][key] = val

    return args


def pick(args):

    """
    This function launches the pick step, using the data derived from the Namespace.
    :param args: argparse Namespace with the configuration for the run.

    """

    args = check_log_settings(args)
    args = check_run_options(args)

    creator = Picker(args.json_conf, commandline=" ".join(sys.argv))
    try:
        creator()  # Run
    except UnsortedInput as err:
        print(err, file=sys.stderr)
        sys.exit(1)


def pick_parser():
    """
    Parser for the picking step.
    """
    parser = argparse.ArgumentParser("Launcher of the Mikado pipeline.")
    parser.add_argument("--start-method", dest="start_method",
                        choices=["fork", "spawn", "forkserver"],
                        default=None, help="Multiprocessing start method.")
    parser.add_argument("-p", "--procs", type=int, default=None,
                        help="""Number of processors to use.
                        Default: look in the configuration file (1 if undefined)""")
    parser.add_argument("--json-conf", dest="json_conf",
                        type=to_json, required=True,
                        help="JSON/YAML configuration file for scoring transcripts.")
    parser.add_argument("--subloci_out", type=str, default=None)
    parser.add_argument("--monoloci_out", type=str, default=None)
    parser.add_argument("--loci_out", type=str, default=None,
                        help="""This output file is mandatory.
                        If it is not specified in the configuration file,
                        it must be provided here.""")
    parser.add_argument("--prefix", type=str, default=None,
                        help="Prefix for the genes. Default: Mikado")
    parser.add_argument("--no_cds", action="store_true", default=None,
                        help="""Flag. If set, not CDS information
                        will be printed out in the GFF output files.""")
    parser.add_argument('--source', type=str, default=None,
                        help='Source field to use for the output files.')
    parser.add_argument('--purge', action='store_true', default=False,
                        help='''Flag. If set, the pipeline will suppress any loci
                        whose transcripts do not pass the requirements set in the JSON file.''')
    parser.add_argument("-shm", "--shared-memory", dest="shm", default=False, action="store_true",
                        help="Flag. If set, the DB will be copied into memory.")
    parser.add_argument("-shmdb", "--shared-memory-db", dest="shm_db", default=None, type=str,
                        help="""Name of the shared memory DB.
                        WARNING: if set, the DB copy will be persistently copied
                        into memory, so that multiple pickers can share.""")
    parser.add_argument('--preload', action='store_true', default=False,
                        help='''Flag. If set, the Mikado DB will be pre-loaded
                        into memory for faster access. WARNING: this option will
                        increase memory usage and the preloading might be quite slow.''')
    parser.add_argument("-db", "--sqlite-db", dest="sqlite_db",
                        default=None, type=str,
                        help="Location of an SQLite database to overwrite what is specified \
                             in the configuration file.")
    parser.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=None,
                        help="Output directory. Default: current working directory")
    parser.add_argument("--single", action="store_true", default=False,
                        help="""Flag. If set, Creator will be launched with a single process.
                        Useful for debugging purposes only.""")
    log_options = parser.add_argument_group("Log options")
    log_options.add_argument("-l", "--log", default=None,
                             help="""File to write the log to.
                             Default: decided by the configuration file.""")
    verbosity = log_options.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="store_true",
                           default=False, help="Flag. If set, the debug mode will be activated.")
    verbosity.add_argument("-nv", "--noverbose", action="store_true",
                           default=False, help="Flag. If set, the debug mode will be activated.")
    log_options.add_argument("-lv", "--log-level", dest="log_level",
                             choices=["DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"], default=None,
                             help="Logging level. Default: retrieved by the configuration file.")
    parser.add_argument("gff", nargs="?", default=None)
    parser.set_defaults(func=pick)
    return parser
