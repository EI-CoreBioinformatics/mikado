#!/usr/bin/env python3

"""
The DAGGER module is built to manage multiple alignments and assemblies of
RNA-Seq data, and subsequently merge them with Mikado.
"""


import sys
import os
import argparse
import datetime
import time
import yaml
import snakemake
from snakemake.utils import min_version
from ..utilities.log_utils import create_default_logger
import shutil
import pkg_resources

# import logging
# import logging.handlers

min_version("3.5")

TIME_START = time.time()
NOW = datetime.datetime.fromtimestamp(TIME_START).strftime('%Y-%m-%d_%H:%M:%S')

# DAGGER_DIR = os.path.dirname(os.path.realpath(__file__))
DAGGER_DIR = pkg_resources.resource_filename("Mikado", "dagger")
assert pkg_resources.resource_exists("Mikado", "dagger")


# noinspection PyPep8Naming

def create_parser():

    """
    Function to create the command-line parser for DAGGER.
    :return:
    """

    parser = argparse.ArgumentParser("""Execute pipeline""")
    # parser.add_argument("config",
    #                     help="Configuration file to use for running dagger.")
    parser.add_argument("-c", "--hpc_conf", default=DAGGER_DIR + "/hpc.json",
                        help="""Configuration file that allows the user to override
                        resource requests for each rule when running under a scheduler
                        in a HPC environment.""")
    parser.add_argument("-N", "--max_nodes", type=int, default="10",
                        help="Maximum number of nodes to use concurrently")
    parser.add_argument("-n", "--max_cores", type=int, default="1000",
                        help="Maximum number of cores to use concurrently")
    #TODO: Implement this.
    parser.add_argument("-t", "--threads", type=int, default=None,
                        help="""Maximum number of threads per job.
                        Default: None (set in the configuration file)""")
    parser.add_argument("-d", "--no_drmaa", action='store_true', default=False,
                        help="Use this flag if running on a HPC and DRMAA is not available")
    parser.add_argument("--force_incomplete", action='store_true', default=False,
                        help="Force snakemake to rerun incomplete steps")
    parser.add_argument("--detailed_summary", action='store_true', default=False,
                        help="Print detailed summary of all input and output files")
    parser.add_argument("--list_resources", action='store_true', default=False,
                        help="List resources used in the workflow")
    parser.add_argument("--make_dag", action='store_true', default=False,
                        help="Produce a DAG rather than execute the workflow")
    return parser


# pylint: disable=too-many-locals
def assemble_transcripts_pipeline(args):

    """
    This section of DAGGER is focused on creating the necessary configuration for
    :param args:
    :return:
    """

    with open(args.config, 'r') as _:
        doc = yaml.load(_)

    # pylint: disable=invalid-name
    LABELS = doc["samples"]
    R1 = doc["r1"]
    R2 = doc["r2"]
    READS_DIR = doc["out_dir"] + "/1-reads"
    SCHEDULER = doc["scheduler"] if doc["scheduler"] else ""
    CWD = os.path.abspath(".")
    # pylint: enable=invalid-name

    res_cmd = ""
    sub_cmd = ""

    if SCHEDULER == "LSF":
        sub_cmd = "bsub"
        res_cmd = " ".join(["-R rusage[mem={cluster.memory}]span[ptile={threads}]",
                            "-n {threads} -q {cluster.queue} -oo /dev/null"])
    elif SCHEDULER == "PBS":
        sub_cmd = "qsub"
        res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q {cluster.queue}"
    elif SCHEDULER == "SLURM":
        sub_cmd = "sbatch"
        res_cmd = " -N 1 -n 1 -c {threads} -p {cluster.queue} --mem={cluster.memory}"

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
                                         os.path.join("dagger", "tr.snakefile"))

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    snakemake.snakemake(
        pkg_resources.resource_filename("Mikado",
                                        os.path.join("dagger", "tr.snakefile")),
        cores=args.max_cores,
        nodes=args.max_nodes,
        configfile=args.config,
        config=additional_config,
        workdir=CWD,
        cluster_config=args.hpc_conf,
        cluster=sub_cmd + res_cmd if args.no_drmaa else None,
        drmaa=res_cmd if not args.no_drmaa else None,
        printshellcmds=True,
        snakemakepath=shutil.which("snakemake"),
        stats="dagger_tr_" + NOW + ".stats",
        force_incomplete=args.force_incomplete,
        detailed_summary=args.detailed_summary,
        list_resources=args.list_resources,
        latency_wait=60 if SCHEDULER else 1,
        printdag=args.make_dag,
        forceall=args.make_dag)
# pylint: enable=too-many-locals


def mikado_pipeline(args):

    """
    This function launches the sub-section dedicated to the Mikado pipeline.
    :param args:
    :return:
    """

    with open(args.config, 'r') as _:
        doc = yaml.load(_)

    # pylint: disable=invalid-name
    SCHEDULER = doc["scheduler"] if ("scheduler" in doc and doc["scheduler"]) else ""
    CWD = os.path.abspath(".")
    # pylint: enable=invalid-name

    res_cmd = ""
    sub_cmd = ""

    if SCHEDULER == "LSF":
        sub_cmd = "bsub"
        res_cmd = " ".join([" -R rusage[mem={cluster.memory}]span[ptile={threads}] -n {threads}",
                            "-q {cluster.queue} -oo /dev/null"])
    elif SCHEDULER == "PBS":
        sub_cmd = "qsub"
        res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q {cluster.queue}"
    elif SCHEDULER == "SLURM":
        sub_cmd = "sbatch"
        res_cmd = " -N 1 -n 1 -c {threads} -p {cluster.queue} --mem={cluster.memory}"

    # Launch using SnakeMake
    assert pkg_resources.resource_exists("Mikado",
                                         os.path.join("dagger", "mikado.snakefile"))

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    snakemake.snakemake(
        pkg_resources.resource_filename("Mikado",
                                        os.path.join("dagger", "mikado.snakefile")),
        cores=args.max_cores,
        nodes=args.max_nodes,
        configfile=args.config,
        config=additional_config,
        workdir=CWD,
        cluster_config=args.hpc_conf,
        cluster=sub_cmd + res_cmd if args.no_drmaa else None,
        drmaa=res_cmd if not args.no_drmaa else None,
        printshellcmds=True,
        snakemakepath=shutil.which("snakemake"),
        stats="dagger_tr_" + NOW + ".stats",
        force_incomplete=args.force_incomplete,
        detailed_summary=args.detailed_summary,
        list_resources=args.list_resources,
        latency_wait=60 if not SCHEDULER == "" else 1,
        printdag=args.make_dag,
        forceall=args.make_dag)


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
        """A Directed Acylic Graph pipeline for GEne model Reconstruction from RNA seq data.
        Basically, a pipeline for driving Mikado. It will first align RNAseq reads against
        a genome using multiple tools, then creates transcript assemblies using multiple tools,
        and find junctions in the alignments using portcullis.
        This input is then passed into Mikado.""")

    subparsers = parser.add_subparsers(
        title="Pipelines",
        help="""These are the pipelines that can be executed via dagger.""")

    subparsers.add_parser("assemble",
                          description="Creates transcript assemblies from RNAseq data.",
                          help="""A pipeline that generates a variety of transcript assemblies
                          using various aligners and assemblers, as well a producing
                          a configuration file suitable for driving Mikado.""")
    subparsers.choices["assemble"] = create_parser()
    subparsers.choices["assemble"].add_argument(
        "config",
        help="Configuration file to use for running the transcript assembly pipeline.")
    subparsers.choices["assemble"].prog = "dagger assemble"

    subparsers.add_parser("mikado",
                          description="Run full mikado pipeline",
                          help="""Using a supplied configuration file that describes
                          all input assemblies to use, it runs the Mikado pipeline,
                          including prepare, BLAST, transdecoder, serialise and pick.""")
    subparsers.choices["mikado"] = create_parser()
    subparsers.choices["mikado"].add_argument(
        "config",
        help="Configuration file to use for running the Mikado step of the pipeline.")
    subparsers.choices["mikado"].prog = "dagger mikado"

    # pylint: disable=broad-except
    try:

        args = parser.parse_args(call_args)
        if len(call_args) == 0:
            parser.print_help()
            sys.exit(1)

        if call_args[0] == "assemble":
            assemble_transcripts_pipeline(args)
        elif call_args[0] == "mikado":
            mikado_pipeline(args)

    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except BrokenPipeError:
        pass
    except Exception as exc:
        logger = create_default_logger("main")
        logger.error("dagger crashed, cause:")
        logger.exception(exc)
        sys.exit(2)
    # pylint: enable=broad-except

if __name__ == '__main__':
    # pylint: disable=redefined-builtin
    # noinspection PyShadowingBuiltins
    __spec__ = "Mikado"
    # pylint: enable=redefined-builtin
    main()
