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
import functools
import inspect
import toml
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader
# import logging
# import logging.handlers

system_hpc_yaml = pkg_resources.resource_filename("Mikado", os.path.join("daijin", "hpc.yaml"))

min_version("3.5")

TIME_START = time.time()
NOW = datetime.datetime.fromtimestamp(TIME_START).strftime('%Y-%m-%d_%H:%M:%S')

DAIJIN_DIR = pkg_resources.resource_filename("Mikado", "daijin")
assert pkg_resources.resource_exists("Mikado", "daijin")


# noinspection PyPep8Naming
def get_sub_commands(SCHEDULER, prefix, additional):
    res_cmd = ""
    sub_cmd = ""

    if SCHEDULER == "LSF":
        sub_cmd = "bsub"
        res_cmd = " ".join([" -R rusage[mem={{cluster.memory}}]span[ptile={{threads}}] -n {{threads}}",
                            "-q {{cluster.queue}} -oo /dev/null",
                            "-J {prefix}_{rule} -oo daijin_logs/{prefix}_{{rule}}_%j.out"]).format(
            prefix=prefix)
    elif SCHEDULER == "PBS":
        sub_cmd = "qsub"
        res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q {cluster.queue}"
    elif SCHEDULER == "SLURM":
        sub_cmd = "sbatch"
        res_cmd = " ".join([" -N 1 -n 1 -c {{threads}} -p {{cluster.queue}} --mem={{cluster.memory}}",
                            "-J {prefix}_{{rule}} -o daijin_logs/{prefix}_{{rule}}_%j.out -e daijin_logs/{prefix}_{{rule}}_%j.err"]).format(prefix=prefix)

    res_cmd = "{} {}".format(res_cmd, additional)
    return res_cmd, sub_cmd


def create_parser():

    """
    Function to create the command-line parser for Daijin.
    :return:
    """

    parser = argparse.ArgumentParser("""Execute pipeline""")
    # parser.add_argument("config",
    #                     help="Configuration file to use for running daijin.")
    parser.add_argument("-c", "--hpc_conf",
                        default="daijin_hpc.yaml",
                        help="""Configuration file that allows the user to override resource requests for each rule \
when running under a scheduler in a HPC environment.""")
    parser.add_argument("--latency-wait", default=None, type=int, dest="latency_wait",
                        help="Latency wait for Daijin. Default: 1s if local, 60s for scheduler jobs.")
    parser.add_argument("-d", "--dryrun", action="store_true", default=False,
                        help="Do a dry run for testing.")
    parser.add_argument("--jobs", "-J", action="store", metavar="N", type=int, default="10",
                        help="Maximum number of cluster jobs to execute concurrently.")
    parser.add_argument("--cores", "-C", action="store", nargs="?", metavar="N", type=int, default="1000",
                        help="Use at most N cores in parallel (default: 1000).")
    parser.add_argument("--threads", "-t", action="store", metavar="N", type=int, default=None,
                        help="""Maximum number of threads per job. \
Default: None (set in the configuration file)""")
    parser.add_argument("--exe", default="daijin_exe.yaml",
                        help="""Configuration file containing the information on the software versions to be used. \
Default: None, Daijin presumes that all needed programs are already present in the environment.""")
    parser.add_argument("--no_drmaa", "-nd", action='store_true', default=False,
                        help="Use this flag if you wish to run without DRMAA, for example, \
if running on a HPC and DRMAA is not available, or if running locally on your own machine or server.")
    parser.add_argument("-ad", "--additional-drmaa", default="", type=str,
                        dest="additional_drmaa", help="Additional parameters to be added to the DRMAA command.")
    parser.add_argument("--rerun-incomplete", "--ri", action='store_true', default=False,
                        dest="rerun_incomplete",
                        help="Re-run all jobs the output of which is recognized as incomplete.")
    parser.add_argument("--use-conda", action="store_true", default=False, dest="use_conda")
    parser.add_argument("--forcerun", "-R", nargs="+", metavar="TARGET",
                        help="Force the re-execution or creation of the given rules or files. \
                        Use this option if you changed a rule and want to have all its output in your \
                        workflow updated.")
    parser.add_argument("--cleanup-metadata", dest="cleanup_metadata",
                        nargs="+",
                        default=[],
                        help="List of files that are complete even if Daijin has lost track of them because eg. they were generated by external processes.")
    parser.add_argument("--prefix", default="daijin",
                        help="Optional prefix to prepend to job names while using DRMAA.")
    parser.add_argument("--detailed-summary", "-D", action='store_true', default=False,
                        dest="detailed_summary",
                        help="Print detailed summary of all input and output files")
    parser.add_argument("--list", "-l", action='store_true', default=False,
                        help="List resources used in the workflow")
    parser.add_argument("--dag", action='store_true', default=False,
                        help="Do not execute anything and print the redirected acylic graph \
                        of jobs in the dot language.")
    parser.add_argument("--nolock", action="store_true", default=False,
                        help="Do not lock the working directory. Use with caution!")
    return parser


def create_config_parser():

    """
    Function to create the configuration file for Daijin.
    :return:
    """

    parser = argparse.ArgumentParser("""Configure the pipeline""")
    runtime = parser.add_argument_group(""""Options related to how to run Daijin - threads, cluster configuration, etc.
Please note that the development of the assembly pipeline has now been DISCONTINUED. We will continue to provide support
for the "mikado" part of daijin.""")
    runtime.add_argument("-c", "--cluster_config",
                        type=str, default=None,
                        help="Cluster configuration file to write to.")
    parser.add_argument("--full", action="store_true", default=False)
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed number.")
    output_format = parser.add_mutually_exclusive_group()
    output_format.add_argument("-j", "--json", action="store_true", default=False,
                               help="Output will be in JSON (default: inferred by filename, with TOML as fallback).")
    output_format.add_argument("-y", "--yaml", action="store_true", default=False,
                               help="Output will be in YAML (default: inferred by filename, with TOML as fallback).")
    output_format.add_argument("--toml", action="store_true", default=False,
                               help="Output will be in TOML (default: inferred by filename, with TOML as fallback).")
    runtime.add_argument("--threads", "-t", action="store", metavar="N", type=int, default=4,
                        help="""Maximum number of threads per job. Default: %(default)s""")
    runtime.add_argument("-od", "--out-dir", dest="out_dir", default=None, required=False,
                        help="Output directory. Default if unspecified: chosen name.")
    runtime.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"),
                        help="Output file. If not specified, it will be printed to STDOUT.\
Daijin will try to infer the type of configuration (TOML, YAML, JSON) from the output file name, with TOML as the\
default. If one of --json, --yaml, --toml flags is specified, it will override the filename inference.")
    runtime.add_argument("--scheduler", default="", choices=["local", "SLURM", "LSF", "PBS"],
                        help="Scheduler to use. Default: None - ie, either execute everything on the local machine or use DRMAA to submit and control jobs (recommended).")
    runtime.add_argument("--exe", default="daijin_exe.yaml",
                         help="Configuration file for the executables.")
    runtime.add_argument("-q", "--queue",
                         default=None,
                         help="Name of queue to be used in the HPC. Required if a scheduler has been selected.")
    reference = parser.add_argument_group("Arguments related to the reference species.")
    reference.add_argument("--name", default="Daijin", help="Name of the species under analysis.")
    reference.add_argument("--genome", "-g", required=True,
                           help="Reference genome for the analysis, in FASTA format. Required.")
    reference.add_argument("--transcriptome", help="Reference annotation, in GFF3 or GTF format.",
                           default="")
    paired_reads = parser.add_argument_group("Arguments related to the input paired reads.")
    paired_input = paired_reads.add_mutually_exclusive_group()
    paired_input.add_argument("--sample-sheet", dest="sample_sheet", default=None,
                              help="""Sample sheet containing the details of the input reads, separated by tabs, in the form:
                              - Read1
                              - Read2 (optional)
                              - Sample name
                              - Strandedness
                              - Boolean flag, True if it is a long read sample, False otherwise. Optional, default if "False".""")
    paired_input.add_argument("-r1", "--left_reads", dest="r1",
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
        default=[], required=False, choices=["fr-unstranded", "fr-secondstrand", "fr-firststrand", "f", "r"],
        help="Strandedness of the reads. Specify it 0, 1, or number of samples times. Choices: %(choices)s.")
    parser.add_argument("-al", "--aligners", choices=["gsnap", "star", "hisat", "tophat"], required=False,
                        default=[], nargs="*", help="Aligner(s) to use for the analysis. Choices: %(choices)s")
    parser.add_argument("-lal", "--long-read-aligners", dest="long_aln_methods", choices=["star", "gmap"],
                        required=False, default=[], nargs="*",
                        help="Aligner(s) to use for long reads. Choices: %(choices)s")
    parser.add_argument("-as", "--assemblers", dest="asm_methods", required=False,
                        choices=["class", "cufflinks", "stringtie", "trinity", "scallop"],
                        default=[], nargs="*", help="Assembler(s) to use for the analysis. Choices: %(choices)s")
    mikado = parser.add_argument_group("Options related to the Mikado phase of the pipeline.")
    # scoring = parser.add_argument_group("Options related to the scoring system")
    scoring_file = mikado.add_mutually_exclusive_group(required=True)
    scoring_file.add_argument("--scoring", type=str, default=None,
                         help="""Available scoring files. Either provide your own of choose from
                         one of the pre-packaged scoring files:
                             {}""".format("\n".join(["- {}".format(_) for _ in pkg_resources.resource_listdir(
                             "Mikado", os.path.join("configuration", "scoring_files"))]
                                                )))
    scoring_file.add_argument("--custom-scoring", type=str, default=None,
                              help="""Pre-created scoring file. If the file does not exist,
                              the skeleton of a scoring file will be written out at the provided location.""")
    mikado.add_argument("--copy-scoring", default=False,
                         type=str, dest="copy_scoring",
                         help="File into which to copy the selected scoring file, for modification.")
    mikado.add_argument("-m", "--modes", default=["stringent"], nargs="+",
                        choices=["nosplit", "split", "permissive", "stringent", "lenient"],
                        required=False,
                        help="Mikado pick modes to run. Choices: %(choices)s")
    mikado.add_argument("-i", "--intron-range",
                        dest="intron_range", type=int, nargs=2,
                        default=None,
                        help="""Range into which intron lengths should fall, as a couple of integers.
                                Transcripts with intron lengths outside of this range will be penalised.
                                Default: (60, 900)""")
    mikado.add_argument("--flank", default=None, type=int,
                        required=False,
                        help="Amount of flanking for grouping transcripts in superloci during the pick phase of Mikado.")
    mikado.add_argument("--prot-db", dest="prot_db", default=[], nargs="+",
                        help="Protein database to compare against, for Mikado.")
    mikado.add_argument("--use-blast", dest="use_blast", action="store_true",
                        default=False, help="Flag. If set, Daijin will use BLAST instead of DIAMOND.")
    mikado.add_argument("--use-transdecoder", dest="use_transdecoder", action="store_true",
                        default=False, help="Flag. If set, Daijin will use TransDecoder instead of Prodigal.")
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
    elif args.config.endswith("yaml"):
        loader = functools.partial(yaml.load, Loader=yLoader)
    else:
        loader = functools.partial(toml.load)

    with open(args.config, 'r') as _:
        doc = loader(_)

    if args.exe and os.path.exists(args.exe):
        if args.exe.endswith("json"):
            loader = json.load
        else:
            loader = functools.partial(yaml.load, Loader=yLoader)
        with open(args.exe) as _:
            doc["load"] = loader(_)

    # print(doc["load"])
    # Check the configuration
    check_config(doc)

    # pylint: disable=invalid-name

    if not "short_reads" in doc and not "long_reads" in doc:
        print("No short reads section or long reads sections was present in the configuration.  Please include your samples and try again")
        exit(1)

    LABELS = []
    R1 = []
    R2 = []
    LR_LABELS = []
    LR_FILES = []
    
    if "short_reads" in doc:
        LABELS = doc["short_reads"]["samples"]
        R1 = doc["short_reads"]["r1"]
        R2 = doc["short_reads"]["r2"]
    
    if "long_reads" in doc:
        LR_LABELS = doc["long_reads"]["samples"]
        LR_FILES = doc["long_reads"]["files"]
    READS_DIR = doc["out_dir"] + "/1-reads"
    SCHEDULER = doc["scheduler"] if doc["scheduler"] else ""
    CWD = os.path.abspath(".")
    # pylint: enable=invalid-name

    res_cmd, sub_cmd = get_sub_commands(SCHEDULER, args.prefix, args.additional_drmaa)

    # Create log folder
    if not os.path.exists("daijin_logs"):
        os.makedirs("daijin_logs")
    elif not os.path.isdir("daijin_logs"):
        raise OSError("{} is not a directory!".format("daijin_logs"))

    if (len(R1) != len(R2)) and (len(R1) != len(LABELS)):
        print("R1, R2 and LABELS lists are not the same length.  Please check and try again")
        exit(1)

    if len(LR_LABELS) != len(LR_FILES):
        print("long read samples and file arrays in the configuration file are not the same length.  Please check and try again")
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
    
    for lr_file, label in zip(LR_FILES, LR_LABELS):
        suffix = lr_file.split(".")[-1]
        compress = ""
        if suffix in ("gz", "bz2"):
            compress = "." + suffix[:]
            suffix = lr_file.split(".")[-2]

        if suffix in ("fa", "fna", "fasta"):
            suffix = ".fa" + compress
        elif suffix in ("fq", "fastq"):
            suffix = ".fq" + compress
        else:
            suffix = ".{}".format(suffix)

        out = READS_DIR + "/" + label + ".long{}".format(suffix)
        if not os.path.islink(out):
            os.symlink(os.path.abspath(lr_file), out)


    # Launch using SnakeMake
    assert pkg_resources.resource_exists("Mikado",
                                         os.path.join("daijin", "assemble.smk"))

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    cluster_var = None
    if args.no_drmaa is True and sub_cmd:
        cluster_var = sub_cmd + res_cmd

    drmaa_var = None
    if args.no_drmaa is False and res_cmd:
        try:
            import drmaa
            _ = drmaa.Session()
        except (RuntimeError,ImportError,AttributeError):
            print("WARNING: DRMAA not installed or not configured properly. Switching to local/cluster mode. Please use the \"-nd\" flag to run Daijin if you do not plan to use DRMAA.", file=sys.stderr)
            drmaa_var = None
            args.no_drmaa = True
        else:
            drmaa_var = res_cmd

    if SCHEDULER == "local":
        hpc_conf = None
        drmaa_var = None
        cluster_var = None
    elif drmaa_var or cluster_var:
        if os.path.exists(args.hpc_conf):
            hpc_conf = args.hpc_conf
        else:
            hpc_conf = system_hpc_yaml
    else:
        hpc_conf = None

    yaml_file = open("daijin.{}.yaml".format(NOW), "wt")
    yaml.dump(doc, yaml_file)
    yaml_file.flush()
    shutil.copystat(args.config, yaml_file.name)

    if args.latency_wait is not None:
        latency = abs(args.latency_wait)
    elif SCHEDULER not in ('', "local"):
        latency = 60
    else:
        latency = 1

    kwds = {
        "dryrun": args.dryrun,
        "cores": args.cores,
        "nodes": args.jobs,
        "config": additional_config,
        "workdir": CWD,
        "cluster_config": hpc_conf,
        "cluster": cluster_var,
        "drmaa": drmaa_var,
        "printshellcmds": True,
        "snakemakepath": shutil.which("snakemake"),
        "use_conda": args.use_conda,
        "stats": "daijin_tr_" + NOW + ".stats",
        "force_incomplete": args.rerun_incomplete,
        "detailed_summary": args.detailed_summary,
        "list_resources": args.list,
        "latency_wait": latency,
        "printdag": args.dag,
        "forceall": args.dag,
        "forcerun": args.forcerun,
        "lock": (not args.nolock),
        "printreason": True
    }

    if "configfile" in inspect.getfullargspec(snakemake.snakemake).args:
        kwds["configfile"] = yaml_file.name
    elif "configfiles" in inspect.getfullargspec(snakemake.snakemake).args:
        kwds["configfiles"] = [yaml_file.name]
    else:
        raise KeyError("No configfile key found")

    snakemake.snakemake(
        pkg_resources.resource_filename("Mikado",
                                        os.path.join("daijin", "assemble.smk")),
        **kwds
        )
# pylint: enable=too-many-locals


def mikado_pipeline(args):

    """
    This function launches the sub-section dedicated to the Mikado pipeline.
    :param args:
    :return:
    """

    if args.config.endswith("json"):
        loader = json.load
    elif args.config.endswith("yaml"):
        loader = functools.partial(yaml.load, Loader=yLoader)
    else:
        loader = functools.partial(toml.load)
    with open(args.config, 'r') as _:
        doc = loader(_)

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    if args.exe and os.path.exists(args.exe):
        if args.exe.endswith("json"):
            loader = json.load
        else:
            loader = functools.partial(yaml.load, Loader=yLoader)
        with open(args.exe) as _:
            doc["load"] = loader(_)

    check_config(doc)

    # pylint: disable=invalid-name
    SCHEDULER = doc["scheduler"] if ("scheduler" in doc and doc["scheduler"]) else ""
    CWD = os.path.abspath(".")
    # pylint: enable=invalid-name

    res_cmd, sub_cmd = get_sub_commands(SCHEDULER, args.prefix, args.additional_drmaa)

    if not os.path.exists("daijin_logs"):
        os.makedirs("daijin_logs")
    elif not os.path.isdir("daijin_logs"):
        raise OSError("{} is not a directory!".format("daijin_logs"))

    # Launch using SnakeMake
    assert pkg_resources.resource_exists("Mikado",
                                         os.path.join("daijin", "mikado.smk"))

    cluster_var = None
    if args.no_drmaa is True and sub_cmd:
        cluster_var = sub_cmd + res_cmd

    drmaa_var = None
    if args.no_drmaa is False and res_cmd:
        try:
            import drmaa
            _ = drmaa.Session()
        except (RuntimeError,ImportError,AttributeError):
            print("WARNING: DRMAA not installed or not configured properly. Switching to local/cluster mode. Please use the \"-nd\" flag to run Daijin if you do not plan to use DRMAA.", file=sys.stderr)
            drmaa_var = None
            args.no_drmaa = True
        else:
            drmaa_var = res_cmd        

    if drmaa_var or cluster_var:
        if os.path.exists(args.hpc_conf):
            hpc_conf = args.hpc_conf
        else:
            hpc_conf = system_hpc_yaml
    else:
        hpc_conf = None

    yaml_file = open("daijin.{}.yaml".format(NOW), "wt")
    yaml.dump(doc, yaml_file)
    yaml_file.flush()
    shutil.copystat(args.config, yaml_file.name)

    if SCHEDULER == "local":
        hpc_conf = None
        drmaa_var = None
        cluster_var = None
    elif drmaa_var or cluster_var:
        if os.path.exists(args.hpc_conf):
            hpc_conf = args.hpc_conf
        else:
            hpc_conf = system_hpc_yaml
    else:
        hpc_conf = None

    if args.latency_wait:
        latency = abs(args.latency_wait)
    elif SCHEDULER not in ('', "local"):
        latency = 60
    else:
        latency = 1

    BLASTX_CHUNKS = max(int(doc["blastx"]["chunks"]), doc["threads"])
    if BLASTX_CHUNKS > doc["blastx"]["chunks"]:
        print("INFO: Increasing the number of chunks for DIAMOND/BLASTX to match the requested threads, \
as Mikado serialise relies on having a number of chunks equal or greater than the number of requested threads.")

    kwds = {
        "ignore_ambiguity": False,
        "cores": args.cores,
        "dryrun": args.dryrun,
        "nodes": args.jobs,
        "config": additional_config,
        "workdir": CWD,
        "cluster_config": hpc_conf,
        "cluster": cluster_var,
        "drmaa": drmaa_var,
        "latency_wait": latency,
        "printshellcmds": True,
        "use_conda": args.use_conda,
        "snakemakepath": shutil.which("snakemake"),
        "stats": "daijin_tr_" + NOW + ".stats",
        "force_incomplete": args.rerun_incomplete,
        "detailed_summary": args.detailed_summary,
        "list_resources": args.list,
        "printdag": args.dag,
        "forceall": args.dag,
        "forcerun": args.forcerun,
        "lock": (not args.nolock),
        "printreason": True
    }

    if "configfile" in inspect.getfullargspec(snakemake.snakemake).args:
        kwds["configfile"] = yaml_file.name
    elif "configfiles" in inspect.getfullargspec(snakemake.snakemake).args:
        kwds["configfiles"] = [yaml_file.name]
    else:
        raise KeyError("No configfile key found")

    snakemake.snakemake(
        pkg_resources.resource_filename("Mikado",
                                        os.path.join("daijin", "mikado.smk")),
        **kwds
    )

    
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
        This input is then passed into Mikado.
        
        WARNING: the "assemble" part of this pipeline will be soon DEPRECATED. 
        """)

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
                          a configuration file suitable for driving Mikado.
                          WARNING: this part of the Daijin pipeline will be DEPRECATED in future releases
                          as a new, more complete annotation pipeline is currently in development.""")
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


if __name__ == '__main__':
    # pylint: disable=redefined-builtin
    # noinspection PyShadowingBuiltins
    __spec__ = "Mikado"
    # pylint: enable=redefined-builtin
    main()
