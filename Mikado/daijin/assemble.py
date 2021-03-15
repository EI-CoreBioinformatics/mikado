import rapidjson as json
import yaml
import functools
import toml
import os
from ..configuration import DaijinConfiguration
from . import get_sub_commands, system_hpc_yaml, NOW
import pkg_resources
import shutil
import sys
import inspect
try:
    import drmaa
    _ = drmaa.Session()
    drmaa_available = True
except (RuntimeError, ImportError, AttributeError):
    drmaa_available = False
import snakemake
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader


def assemble_transcripts_pipeline(args):

    """
    This section of Daijin is focused on creating the necessary configuration for
    driving the pipeline.
    :param args: CLI arguments from argparse
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
    _ = DaijinConfiguration.Schema().load(doc)

    # pylint: disable=invalid-name
    if not ("short_reads" in doc or "long_reads" in doc):
        print("No short reads section or long reads sections was present in the configuration. Please include your "
              "samples and try again", file=sys.stderr)
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
        print("long read samples and file arrays in the configuration file are not the same length. Please check and "
              "try again")
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
    assert pkg_resources.resource_exists("Mikado", os.path.join("daijin", "assemble.smk"))

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    cluster_var = None
    if args.no_drmaa is True and sub_cmd:
        cluster_var = sub_cmd + res_cmd

    drmaa_var = None
    if args.no_drmaa is False and res_cmd:
        if drmaa_available is False:
            print("WARNING: DRMAA not installed or not configured properly. Switching to local/cluster mode. Please "
                  "use the \"-nd\" flag to run Daijin if you do not plan to use DRMAA.", file=sys.stderr)
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

    yaml_file = open(os.path.join(doc["out_dir"], "daijin.{}.yaml".format(NOW)), "wt")
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
