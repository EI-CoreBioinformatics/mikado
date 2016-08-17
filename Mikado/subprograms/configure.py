#!/usr/bin/env python3

"""Stub of pre-configurer for Mikado"""


import yaml
import itertools
import re
import os
from pkg_resources import resource_listdir, resource_stream
import argparse
import sys
from ..configuration import configurator
from ..exceptions import InvalidJson
from ..utilities import comma_split  # , merge_dictionaries
import json
from collections import Counter
import tempfile

__author__ = 'Luca Venturini'


def check_has_requirements(dictionary, schema, key=None, first_level=True):

    """
    Method to find all keys that
    :param dictionary:
    :param schema:
    :param key:
    :return:
    """

    required = []

    for new_key, value in dictionary.items():
        if isinstance(value, dict):
            assert "properties" in schema[new_key], new_key
            if "SimpleComment" in schema[new_key]:
                required.append((key, new_key, "SimpleComment"))
            if "required" in schema[new_key]:
                for req in schema[new_key]["required"]:
                    required.append((key, new_key, req))

            for k in check_has_requirements(value, schema[new_key]["properties"],
                                            key=new_key,
                                            first_level=False):
                if k is None:
                    continue
                nkey = [key]
                nkey.extend(k)
                nkey = tuple(nkey)
                required.append(nkey)
        elif first_level is True:
            if new_key in ("Comment", "SimpleComment"):
                continue
            elif new_key in schema:
                # if "SimpleComment" in schema[new_key]:
                #     required.append((key, new_key, "SimpleComment"))

                if "required" in schema[new_key] and schema[new_key]["required"] is True:
                    required.append([new_key])
        else:
            continue

    return required


def get_key(new_dict, key, default):

    """
    Recursive method to get a nested key from inside the "default" dict
    and transfer it, keeping the tree structure, inside the
    new_dict
    :param new_dict: dictionary to transfer the key to
    :param key: composite key
    :param default: dictionary to extract the key from
    :return: new_dict (with updated structure)
    """

    if isinstance(default[key[0]], dict):
        assert len(key) > 1
        new_dict.setdefault(key[0], new_dict.get(key[0], dict()))
        new_dict = get_key(new_dict[key[0]], key[1:], default[key[0]])
    else:
        assert len(key) == 1
        new_dict[key[0]] = default[key[0]]
    return new_dict


def create_simple_config():

    """
    Method to create a stripped down configuration dictionary
    containing only SimpleComments and required fields.
    :return:
    """

    default = configurator.to_json("", simple=True)
    validator = configurator.create_validator(simple=True)

    del default["scoring"]
    del default["requirements"]
    # del default["soft_requirements"]

    new_dict = dict()
    composite_keys = [(ckey[1:]) for ckey in
                      check_has_requirements(default,
                                             validator.schema["properties"])]

    # Sort the composite keys by depth
    for ckey in sorted(composite_keys, key=len, reverse=True):
        defa = default
        # Get to the latest position
        for key in ckey:
            try:
                defa = defa[key]
            except KeyError:
                raise KeyError(key, defa)
        val = defa
        for k in reversed(ckey):
            val = {k: val}

        new_dict = configurator.merge_dictionaries(new_dict, val)

    return new_dict


def print_config(output, out):

    """
    Function to print out the prepared configuration.
    :param output: prepared output, a huge string.
    :type output: str

    :param out: output handle.
    """

    comment = []
    comment_level = -1

    for line in output.split("\n"):
        # comment found
        if line.lstrip().startswith(("Comment", "SimpleComment")) or comment:
            level = sum(1 for _ in itertools.takewhile(str.isspace, line))
            line = re.sub("Comment:", "", re.sub("SimpleComment:", "", line))
            if comment:
                if level > comment_level or line.lstrip().startswith("-"):
                    comment.append(line.strip())
                else:
                    for comment_line in iter(_ for _ in comment if _ != ''):
                        print("{spaces}#  {comment}".format(spaces=" "*comment_level,
                                                            comment=re.sub(
                                                                "'", "", re.sub("^- ", "",
                                                                                comment_line))),
                              file=out)
                    if level < comment_level:
                        print("{0}{{}}".format(" " * comment_level), file=out)
                    comment = []
                    comment_level = -1

                    print(line.rstrip(), file=out)
            else:
                comment = [re.sub("(Comment|SimpleComment):", "", line.strip())]
                comment_level = level
        else:
            print(line.rstrip(), file=out)

    if comment:
        for comment_line in comment:
            print("{spaces}#{comment}".format(spaces=" "*comment_level, comment=comment_line),
                  file=out)


def create_config(args):
    """
    Utility to create a default configuration file.
    :param args:
    :return:
    """

    if args.full is True:
        default = configurator.to_json("")
        del default["scoring"]
        del default["requirements"]
        config = default
    else:
        config = create_simple_config()

    if args.external is not None:
        if args.external.endswith("json"):
            loader = json.load
        else:
            loader = yaml.load
        with open(args.external) as external:
            external_conf = loader(external)
        # Overwrite values specific to Mikado
        if "mikado" in external_conf:
            mikado_conf = dict((key, val) for key, val in external_conf["mikado"].items() if key in config)
            config = configurator.merge_dictionaries(config, mikado_conf)
        # Leave all other values, including those in a "mikado" section that are not also present in the default config
        for key in config:
            if key in external_conf:
                del external_conf[key]
            if "mikado" in external_conf and isinstance(external_conf["mikado"], dict):
                __mikado_keys = [key for key in external_conf["mikado"] if key in config]
                for key in __mikado_keys:
                    del external_conf["mikado"][key]
        config = configurator.merge_dictionaries(config, external_conf)

    if args.reference is not None:
        config["reference"]["genome"] = args.reference

    if args.junctions is not None:
        config["serialise"]["files"]["junctions"] = args.junctions

    if args.blast_targets is not None:
        config["serialise"]["files"]["blast_targets"] = args.blast_targets

    if args.gff:
        args.gff = args.gff.split(",")
        __gff_counter = Counter()
        __gff_counter.update(args.gff)

        if __gff_counter.most_common()[0][1] > 1:
            raise InvalidJson(
                "Repeated elements among the input GFFs! Duplicated files: {}".format(
                    ", ".join(_[0] for _ in __gff_counter.most_common() if _[1] > 1)
                ))

        config["prepare"]["files"]["gff"] = args.gff

        if args.labels != '':
            args.labels = args.labels.split(",")
            if not len(args.labels) == len(args.gff):
                raise ValueError("""Length mismatch between input files and labels!
                GFFs: {0} (length {1})
                Labels: {2} (length {3})""".format(
                    args.gff, len(args.gff),
                    args.labels, len(args.labels)))
            config["prepare"]["files"]["labels"] = args.labels

        if args.strand_specific_assemblies != "":
            args.strand_specific_assemblies = args.strand_specific_assemblies.split(",")
            if (len(args.strand_specific_assemblies) > len(args.gff) or
                    any([(_ not in args.gff) for _ in args.strand_specific_assemblies])):
                raise InvalidJson("Invalid strand-specific assemblies specified")
            config["prepare"]["files"]["strand_specific_assemblies"] = args.strand_specific_assemblies

    elif args.list:
        config["pick"]["source_score"] = dict()
        with open(args.list) as list_file:
            files, labels, strandedness, scores = [], [], [], []
            files_counter = Counter()
            for line in list_file:
                try:
                    _fields = line.rstrip().split("\t")
                    filename, label, stranded =  _fields[:3]

                    if not os.path.exists(filename):
                        raise ValueError("Invalid file name: {}".format(filename))
                    files.append(filename)
                    if label in labels:
                        raise ValueError("Non-unique label specified: {}".format(label))
                    labels.append(label)
                    strandedness.append(stranded)
                    if len(_fields) > 3:
                        score = float(_fields[3])
                    else:
                        score = 0
                    scores.append(score)
                except ValueError as exc:
                    raise ValueError("Malformed inputs file. Error:\n{}".format(exc))
            files_counter.update(files)
            if files_counter.most_common()[0][1] > 1:
                raise InvalidJson(
                    "Repeated elements among the input GFFs! Duplicated files: {}".format(
                        ", ".join(_[0] for _ in files_counter.most_common() if _[1] > 1)))
            if any([_ not in ("True", "False") for _ in strandedness]):
                raise InvalidJson("Invalid values for strandedness in the list file.")
            config["prepare"]["files"]["labels"] = list(labels)
            config["prepare"]["files"]["gff"] = list(files)
            config["prepare"]["files"]["strand_specific_assemblies"] = [files[_[0]] for _ in enumerate(strandedness)
                                                                        if _[1] == "True"]
            for source, score in zip(labels, scores):
                config["pick"]["source_score"][source] = score

    elif args.no_files is True:
        for stage in ["pick", "prepare", "serialise"]:
            if "files" in config[stage]:
                del config[stage]["files"]
            # except KeyError:
            #     raise KeyError(stage)
        del config["reference"]
        del config["db_settings"]

    if args.scoring is not None:
        if args.copy_scoring is not False:
            with open(args.copy_scoring, "wt") as out:
                with resource_stream("Mikado", os.path.join("configuration",
                                                            "scoring_files",
                                                            args.scoring)) as original:
                    for line in original:
                        print(line.decode().rstrip(), file=out)
            args.scoring = args.copy_scoring

        config["pick"]["scoring_file"] = args.scoring

    if args.mode is not None:
        if args.mode == "nosplit":
            config["pick"]["chimera_split"]["execute"] = False
        else:
            config["pick"]["chimera_split"]["execute"] = True
            if args.mode == "split":
                config["pick"]["chimera_split"]["blast_check"] = False
            else:
                config["pick"]["chimera_split"]["blast_check"] = True
                config["pick"]["chimera_split"]["blast_params"]["leniency"] = args.mode.upper()


    # Check that the configuration file is correct
    tempcheck = tempfile.NamedTemporaryFile("wt", suffix=".yaml")
    output = yaml.dump(config, default_flow_style=False)
    print_config(output, tempcheck)
    tempcheck.flush()
    try:
        configurator.to_json(tempcheck.name)
    except InvalidJson as exc:
        raise InvalidJson("Created an invalid configuration file! Error:\n{}".format(exc))

    if args.json is True or args.out.name.endswith("json"):
        json.dump(config, args.out, sort_keys=True, indent=4)
    else:
        output = yaml.dump(config, default_flow_style=False)
        print_config(output, args.out)


def configure_parser():
    """
    Parser for the configuration utility.
    :return: the parser.
    :rtype: argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser(description="Configuration utility for Mikado",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--full", action="store_true", default=False)
    scoring = parser.add_argument_group("Options related to the scoring system")
    scoring.add_argument("--scoring", type=str, default=None,
                         help="Scoring file to use. Mikado provides the following: {}".format(
                             ",".join(resource_listdir("Mikado", os.path.join("configuration", "scoring_files"))
                         )))
    scoring.add_argument("--copy-scoring", default=False,
                         type=str, dest="copy_scoring",
                         help="File into which to copy the selected scoring file, for modification.")
    parser.add_argument("--strand-specific", default=False,
                        action="store_true",
                        help="""Boolean flag indicating whether all the assemblies are strand-specific.""")
    files = parser.add_mutually_exclusive_group()
    files.add_argument("--no-files", dest="no_files",
                       help="""Remove all files-specific options from the printed configuration file.
                       Invoking the "--gff" option will disable this flag.""",
                       default=False, action="store_true")
    files.add_argument("--gff", help="Input GFF/GTF file(s), separated by comma", type=str)
    files.add_argument("--list", help="""List of the inputs, one by line, in the form:
<file1>  <label>  <strandedness (true/false)>""")
    parser.add_argument("--reference", help="Fasta genomic reference.", default=None)
    serialisers = parser.add_argument_group(
        "Options related to the serialisation step")
    serialisers.add_argument("--junctions", type=comma_split, default=[])
    serialisers.add_argument("-bt", "--blast_targets", type=comma_split, default=None)
    parser.add_argument("--strand-specific-assemblies", type=str, default="",
                        dest="strand_specific_assemblies",
                        help="""List of strand-specific assemblies among the inputs.""")
    parser.add_argument("--labels", type=str, default="",
                        help="""Labels to attach to the IDs of the transcripts of the input files,
        separated by comma.""")
    parser.add_argument("--external", help="""External configuration file to overwrite/add values from.
    Parameters specified on the command line will take precedence over those present in the configuration file.""")
    parser.add_argument("--mode", default=None,
                        choices=["nosplit", "stringent", "lenient", "permissive", "split"],
                        help="""Mode in which Mikado will treat transcripts with multiple ORFs.
- nosplit: keep the transcripts whole.
- stringent: split multi-orf transcripts if two consecutive ORFs have both BLAST hits
             and none of those hits is against the same target.
- lenient: split multi-orf transcripts as in stringent, and additionally, also when
           either of the ORFs lacks a BLAST hit (but not both).
- permissive: like lenient, but also split when both ORFs lack BLAST hits
- split: split multi-orf transcripts regardless of what BLAST data is available.""")
    parser.add_argument("-j", "--json", action="store_true", default=False,
                        help="Output will be in JSON instead of YAML format.")
    parser.add_argument("out", nargs='?', default=sys.stdout, type=argparse.FileType('w'))
    parser.set_defaults(func=create_config)
    return parser
