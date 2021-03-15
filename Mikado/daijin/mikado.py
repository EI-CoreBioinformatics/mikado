import rapidjson as json
import yaml
import functools
import toml
import os
from ..configuration import DaijinConfiguration
from . import get_sub_commands, system_hpc_yaml, NOW
import pkg_resources
from dataclasses import asdict
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


def mikado_pipeline(args):

    """
    This function launches the sub-section dedicated to the Mikado pipeline.
    :param args: argparse Namespace
    :return:
    """

    if args.config.endswith("json"):
        loader = json.load
    elif args.config.endswith("yaml"):
        loader = functools.partial(yaml.load, Loader=yLoader)
    else:
        loader = functools.partial(toml.load)
    with open(args.config, 'r') as _:
        daijin_config = loader(_)

    additional_config = {}
    if args.threads is not None:
        additional_config["threads"] = args.threads

    if args.exe and os.path.exists(args.exe):
        if args.exe.endswith("json"):
            loader = json.load
        else:
            loader = functools.partial(yaml.load, Loader=yLoader)
        with open(args.exe) as _:
            daijin_config["load"] = loader(_)

    daijin_config = DaijinConfiguration.Schema().load(daijin_config)

    # pylint: disable=invalid-name
    SCHEDULER = daijin_config.scheduler if daijin_config.scheduler else ""
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
        if drmaa_available is False:
            msg = "WARNING: DRMAA not installed or not configured properly. Switching to local/cluster mode. Please " \
                  "use the \"-nd\" flag to run Daijin if you do not plan to use DRMAA."
            print(msg, file=sys.stderr)
            drmaa_var = None
            args.no_drmaa = True
        else:
            drmaa_var = res_cmd

    yaml_file = open(os.path.join(daijin_config.out_dir, "daijin.{}.yaml".format(NOW)), "wt")
    yaml.dump(asdict(daijin_config), yaml_file)
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

    BLASTX_CHUNKS = max(int(daijin_config.blastx.chunks), daijin_config.threads)
    if BLASTX_CHUNKS > daijin_config.blastx.chunks:
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
