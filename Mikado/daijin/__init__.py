#!/usr/bin/env python3

"""
The Daijin module is built to manage multiple alignments and assemblies of
RNA-Seq data, and subsequently merge them with Mikado.
"""


import sys
import os
import argparse
import datetime
import time
import json
import yaml
import snakemake
from snakemake.utils import min_version
from ..utilities.log_utils import create_default_logger
from ..configuration.daijin_configurator import create_daijin_config, check_config
import shutil
import pkg_resources

# import logging
# import logging.handlers

min_version("3.5")

TIME_START = time.time()
NOW = datetime.datetime.fromtimestamp(TIME_START).strftime('%Y-%m-%d_%H:%M:%S')

DAIJIN_DIR = pkg_resources.resource_filename("Mikado", "daijin")
assert pkg_resources.resource_exists("Mikado", "daijin")


# noinspection PyPep8Naming
def get_sub_commands(SCHEDULER):
    res_cmd = ""
    sub_cmd = ""

    if SCHEDULER == "LSF":
        sub_cmd = "bsub"
        res_cmd = " ".join([" -R rusage[mem={cluster.memory}]span[ptile={threads}] -n {threads}",
                            "-q {cluster.queue} -oo /dev/null",
                            "-J {rule} -oo daijin_logs/{rule}_%j.out"])
    elif SCHEDULER == "PBS":
        sub_cmd = "qsub"
        res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q {cluster.queue}"
    elif SCHEDULER == "SLURM":
        sub_cmd = "sbatch"
        res_cmd = " ".join([" -N 1 -n 1 -c {threads} -p {cluster.queue} --mem={cluster.memory}",
                            "-J {rule} -o daijin_logs/{rule}_%j.out -e daijin_logs/{rule}_%j.err"])
    return res_cmd, sub_cmd


def create_parser():

    """
    Function to create the command-line parser for Daijin.
    :return:
    """

    parser = argparse.ArgumentParser("""Execute pipeline""")
    # parser.add_argument("config",
    #                     help="Configuration file to use for running daijin.")
    parser.add_argument("-c", "--hpc_conf", default=pkg_resources.resource_filename(
        "Mikado", os.path.join("daijin", "hpc.yaml")),
                        help="""Configuration file that allows the user to override
                        resource requests for each rule when running under a scheduler
                        in a HPC environment.""")
    parser.add_argument("-d", "--dryrun", action="store_true", default=False,
                        help="Do a dry run for testing.")
    parser.add_argument("--jobs", "-J", action="store", metavar="N", type=int, default="10",
                        help="Maximum number of cluster jobs to execute concurrently.")
    parser.add_argument("--cores", "-C", action="store", nargs="?", metavar="N", type=int, default="1000",
                        help="Use at most N cores in parallel (default: 1000).")
    parser.add_argument("--threads", "-t", action="store", metavar="N", type=int, default=None,
                        help="""Maximum number of threads per job.
                        Default: None (set in the configuration file)""")
    parser.add_argument("--no_drmaa", "-nd", action='store_true', default=False,
                        help="Use this flag if you wish to run without DRMAA, for example, \
if running on a HPC and DRMAA is not available, or if running locally on your own machine or server.")
    parser.add_argument("--rerun-incomplete", "--ri", action='store_true', default=False,
                        dest="rerun_incomplete",
                        help="Re-run all jobs the output of which is recognized as incomplete.")
    parser.add_argument("--forcerun", "-R", nargs="+", metavar="TARGET",
                        help="Force the re-execution or creation of the given rules or files. \
                        Use this option if you changed a rule and want to have all its output in your \
                        workflow updated.")
    parser.add_argument("--detailed-summary", "-D", action='store_true', default=False,
                        dest="detailed_summary",
                        help="Print detailed summary of all input and output files")
    parser.add_argument("--list", "-l", action='store_true', default=False,
                        help="List resources used in the workflow")
    parser.add_argument("--dag", action='store_true', default=False,
                        help="Do not execute anything and print the redirected acylic graph \
                        of jobs in the dot language.")
    return parser


def create_config_parser():

    """
    Function to create the configuration file for Daijin.
    :return:
    """

    parser = argparse.ArgumentParser("""Configure the pipeline""")
    runtime = parser.add_argument_group("Options related to how to run Daijin - threads, cluster configuration, etc.")
    runtime.add_argument("-c", "--cluster_config",
                        type=str, default=None,
                        help="Cluster configuration file to write to.")
    runtime.add_argument("--threads", "-t", action="store", metavar="N", type=int, default=4,
                        help="""Maximum number of threads per job. Default: %(default)s""")
    runtime.add_argument("-od", "--out-dir", dest="out_dir", default=None, required=False,
                        help="Output directory. Default if unspecified: chosen name.")
    runtime.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"),
                    help="Output file. If the file name ends in \"json\", the file will be in JSON format; \
                        otherwise, Daijin will print out a YAML file. Default: STDOUT.")
    runtime.add_argument("--scheduler", default="", choices=["", "SLURM", "LSF", "PBS"],
                        help="Scheduler to use. Default: None - ie, either execute everything on the local machine or use DRMAA to submit and control jobs (recommended).")
    reference = parser.add_argument_group("Arguments related to the reference species.")
    reference.add_argument("--name", default="Daijin", help="Name of the species under analysis.")
    reference.add_argument("--genome", "-g", required=True,
                        help="Reference genome for the analysis, in FASTA format. Required.")
    reference.add_argument("--transcriptome", help="Reference annotation, in GFF3 or GTF format.",
                        default="")
    paired_reads = parser.add_argument_group("Arguments related to the input paired reads.")
    paired_reads.add_argument("-r1", "--left_reads", dest="r1",
                        nargs="+",
                        default=[], required=False,
                        help="Left reads for the analysis. Required.")
    paired_reads.add_argument("-r2", "--right_reads", dest="r2",
                        nargs="+",
                        default=[], required=False,
                        help="Right reads for the analysis. Required.")
    paired_reads.add_argument("-s", "--samples",
                        nargs="+",
                        default=[], required=False,
                        help="Sample names for the analysis. Required.")
    paired_reads.add_argument(
        "-st", "--strandedness", nargs="+",
        default=[], required=False, choices=["fr-unstranded", "fr-secondstrand", "fr-firststrand"],
        help="Strandedness of the reads. Specify it 0, 1, or number of samples times. Choices: %(choices)s.")
    parser.add_argument("-al", "--aligners", choices=["gsnap", "star", "hisat", "tophat2"], required=True,
                        default=[], nargs="*", help="Aligner(s) to use for the analysis. Choices: %(choices)s")
    parser.add_argument("-as", "--assemblers", dest="asm_methods", required=True,
                        choices=["class", "cufflinks", "stringtie", "trinity"],
                        default=[], nargs="*", help="Assembler(s) to use for the analysis. Choices: %(choices)s")
    mikado = parser.add_argument_group("Options related to the Mikado phase of the pipeline.")
    # scoring = parser.add_argument_group("Options related to the scoring system")
    mikado.add_argument("--scoring", type=str, default=None,
                         choices=pkg_resources.resource_listdir(
                             "Mikado", os.path.join("configuration", "scoring_files")),
                         help="Available scoring files.")
    mikado.add_argument("--copy-scoring", default=False,
                         type=str, dest="copy_scoring",
                         help="File into which to copy the selected scoring file, for modification.")
    mikado.add_argument("-m", "--modes", default=["permissive"], nargs="+",
                        choices=["nosplit", "split", "permissive", "stringent", "lenient"],
                        required=False,
                        help="Mikado pick modes to run. Choices: %(choices)s")
    mikado.add_argument("--flank", default=None, type=int,
                        required=False,
                        help="Amount of flanking for grouping transcripts in superloci during the pick phase of Mikado.")
    mikado.add_argument("--prot-db", dest="prot_db", default=[], nargs="+",
                        help="Protein database to compare against, for Mikado.")
    parser.set_defaults(func=create_daijin_config)
    return parser


# pylint: disable=too-many-locals
def assemble_transcripts_pipeline(args):

    """
    This section of Daijin is focused on creating the necessary configuration for
    driving the pipeline.
    :param args:
    :return:
    """

    if args.config.endswith("json"):
        loader = json.load
    else:
        loader = yaml.load

    with open(args.config, 'r') as _:
        doc = loader(_)

    # Check the configuration
    check_config(doc)

    # pylint: disable=invalid-name
    LABELS = doc["short_reads"]["samples"]
    R1 = doc["short_reads"]["r1"]
    R2 = doc["short_reads"]["r2"]
    READS_DIR = doc["out_dir"] + "/1-reads"
    SCHEDULER = doc["scheduler"] if doc["scheduler"] else ""
    CWD = os.path.abspath(".")
    # pylint: enable=invalid-name

    res_cmd, sub_cmd = get_sub_commands(SCHEDULER)

    # Create log folder
    if not os.path.exists("daijin_logs"):
        os.makedirs("daijin_logs")
    elif not os.path.isdir("daijin_logs"):
        raise OSError("{} is not a directory!".format("daijin_logs"))

    if (len(R1) != len(R2)) and (len(R1) != len(LABELS)):
        print("R1, R2 and LABELS lists are not the same length.  Please check and try again")
        exit(1)

    if not os.path.exists(READS_DIR):
        os.makedirs(READS_DIR)

    for read1, read2, label in zip(R1, R2, LABELS):
        suffix = read1.split(".")[-1]
        if suffix not in ("gz", "bz2"):
            suffix = ""
        else:
            suffix = ".{}".format(suffix)

        r1out = READS_DIR + "/" + label + ".R1.fq{}".format(suffix)
        r2out = READS_DIR + "/" + label + ".R2.fq{}".format(suffix)
        if not os.path.islink(r1out):
            os.symlink(os.path.abspath(read1), r1out)

        if not os.path.islink(r2out):
            os.symlink(os.path.abspath(read2), r2out)

    # Launch using SnakeMake
    assert pkg_resources.resource_exists("Mikado",
                                         os.path.join("daijin", "tr.snakefile"))

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    snakemake.snakemake(
        pkg_resources.resource_filename("Mikado",
                                        os.path.join("daijin", "tr.snakefile")),
        dryrun=args.dryrun,
        cores=args.cores,
        nodes=args.jobs,
        configfile=args.config,
        config=additional_config,
        workdir=CWD,
        cluster_config=args.hpc_conf,
        cluster=sub_cmd + res_cmd if args.no_drmaa else None,
        drmaa=res_cmd if (not args.no_drmaa and res_cmd is not None) else None,
        printshellcmds=True,
        snakemakepath=shutil.which("snakemake"),
        stats="daijin_tr_" + NOW + ".stats",
        force_incomplete=args.rerun_incomplete,
        detailed_summary=args.detailed_summary,
        list_resources=args.list,
        latency_wait=60 if SCHEDULER else 1,
        printdag=args.dag,
        forceall=args.dag,
        forcerun=args.forcerun)
# pylint: enable=too-many-locals


def mikado_pipeline(args):

    """
    This function launches the sub-section dedicated to the Mikado pipeline.
    :param args:
    :return:
    """

    if args.config.endswith("json"):
        loader = json.load
    else:
        loader = yaml.load

    with open(args.config, 'r') as _:
        doc = loader(_)

    check_config(doc)

    # pylint: disable=invalid-name
    SCHEDULER = doc["scheduler"] if ("scheduler" in doc and doc["scheduler"]) else ""
    CWD = os.path.abspath(".")
    # pylint: enable=invalid-name

    res_cmd, sub_cmd = get_sub_commands(SCHEDULER)

    if not os.path.exists("daijin_logs"):
        os.makedirs("daijin_logs")
    elif not os.path.isdir("daijin_logs"):
        raise OSError("{} is not a directory!".format("daijin_logs"))

    # Launch using SnakeMake
    assert pkg_resources.resource_exists("Mikado",
                                         os.path.join("daijin", "mikado.snakefile"))

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    snakemake.snakemake(
        pkg_resources.resource_filename("Mikado",
                                        os.path.join("daijin", "mikado.snakefile")),
        cores=args.cores,
        dryrun=args.dryrun,
        nodes=args.jobs,
        configfile=args.config,
        config=additional_config,
        workdir=CWD,
        cluster_config=args.hpc_conf,
        cluster=sub_cmd if args.no_drmaa else None,
        drmaa=res_cmd if (not args.no_drmaa and res_cmd is not None) else None,
        printshellcmds=True,
        snakemakepath=shutil.which("snakemake"),
        stats="daijin_tr_" + NOW + ".stats",
        force_incomplete=args.rerun_incomplete,
        detailed_summary=args.detailed_summary,
        list_resources=args.list,
        latency_wait=60 if not SCHEDULER == "" else 1,
        printdag=args.dag,
        forceall=args.dag,
        forcerun=args.forcerun)


def main(call_args=None):

    """
    Main call function.
    :param call_args: Arguments to use to launch the pipeline. If unspecified, the default behaviour
    (using CL arguments) will be adopted.
    :return:
    """

    if call_args is None:
        call_args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        """A Directed Acyclic pipeline for gene model reconstruction from RNA seq data.
        Basically, a pipeline for driving Mikado. It will first align RNAseq reads against
        a genome using multiple tools, then creates transcript assemblies using multiple tools,
        and find junctions in the alignments using Portcullis.
        This input is then passed into Mikado.""")

    subparsers = parser.add_subparsers(
        title="Pipelines",
        help="""These are the pipelines that can be executed via daijin.""")

    subparsers.add_parser("configure",
                          help="Creates the configuration files for Daijin execution.")
    subparsers.choices["configure"] = create_config_parser()
    subparsers.choices["configure"].prog = "daijin configure"
    subparsers.choices["configure"].set_defaults(func=create_daijin_config)

    subparsers.add_parser("assemble",
                          description="Creates transcript assemblies from RNAseq data.",
                          help="""A pipeline that generates a variety of transcript assemblies
                          using various aligners and assemblers, as well a producing
                          a configuration file suitable for driving Mikado.""")
    subparsers.choices["assemble"] = create_parser()
    subparsers.choices["assemble"].add_argument(
        "config",
        help="Configuration file to use for running the transcript assembly pipeline.")
    subparsers.choices["assemble"].prog = "daijin assemble"
    subparsers.choices["assemble"].set_defaults(func=assemble_transcripts_pipeline)

    subparsers.add_parser("mikado",
                          description="Run full mikado pipeline",
                          help="""Using a supplied configuration file that describes
                          all input assemblies to use, it runs the Mikado pipeline,
                          including prepare, BLAST, transdecoder, serialise and pick.""")
    subparsers.choices["mikado"] = create_parser()
    subparsers.choices["mikado"].add_argument(
        "config",
        help="Configuration file to use for running the Mikado step of the pipeline.")
    subparsers.choices["mikado"].prog = "daijin mikado"
    subparsers.choices["mikado"].set_defaults(func=mikado_pipeline)

    try:
        args = parser.parse_args(call_args)
        if hasattr(args, "func"):
            args.func(args)
        else:
            parser.print_help()
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except BrokenPipeError:
        pass
    except Exception as exc:
        logger = create_default_logger("main")
        logger.error("Daijin crashed, cause:")
        logger.exception(exc)
        sys.exit(1)

    # args = parser.parse_args(call_args)
    # pylint: disable=broad-except
    # try:
    #
    #     args = parser.parse_args(call_args)
    #     if len(call_args) == 0:
    #         parser.print_help()
    #         sys.exit(1)
    #
    #     if call_args[0] == "assemble":
    #         assemble_transcripts_pipeline(args)
    #     elif call_args[0] == "mikado":
    #         mikado_pipeline(args)
    #     elif call_args[0] == "configure":
    #         daijin_config(args)
    #     else:
    #         raise ValueError("Invalid subprogram specified!")
    #
    # except KeyboardInterrupt:
    #     raise KeyboardInterrupt
    # except BrokenPipeError:
    #     pass
    # except Exception as exc:
    #     logger = create_default_logger("main")
    #     logger.error("daijin crashed, cause:")
    #     logger.exception(exc)
    #     sys.exit(2)
    # pylint: enable=broad-except

if __name__ == '__main__':
    # pylint: disable=redefined-builtin
    # noinspection PyShadowingBuiltins
    __spec__ = "Mikado"
    # pylint: enable=redefined-builtin
    main()
