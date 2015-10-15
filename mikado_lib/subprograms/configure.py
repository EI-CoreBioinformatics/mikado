#!/usr/bin/env python3

__author__ = 'Luca Venturini'

import yaml
import itertools
import re
import argparse
import sys
import mikado_lib
import jsonschema

"""Stub of pre-configurer for mikado_lib"""


def check_has_requirements(dictionary, schema, key=None):

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
            assert "properties" in schema[new_key]
            if "SimpleComment" in schema[new_key]:
                required.append((key, new_key, "SimpleComment"))
            if "required" in schema[new_key]:
                for req in schema[new_key]["required"]:
                    required.append((key, new_key, req))

            for k in check_has_requirements(value, schema[new_key]["properties"], key=new_key):
                if k is None:
                    continue
                nkey = [key]
                nkey.extend(k)
                nkey = tuple(nkey)
                required.append(nkey)
        else:
            continue

    return required


def get_key(new_dict, key, default):

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

    default = mikado_lib.configuration.configurator.to_json("", simple=True)
    validator = mikado_lib.configuration.configurator.create_validator(simple=True)

    del default["scoring"]
    del default["requirements"]
    del default["soft_requirements"]

    new_dict = dict()
    composite_keys = [(ckey[1:]) for ckey in
                      check_has_requirements(default,
                                             validator.schema["properties"])]

    # Sort the composite keys by depth
    for ckey in sorted(composite_keys, key=len, reverse=True):
        defa = default
        # Get to the latest position
        for pos, k in enumerate(ckey):
            try:
                defa = defa[k]
            except KeyError:
                raise KeyError(k, defa)
        val = defa
        for k in reversed(ckey):
            val = {k: val}

        new_dict = mikado_lib.configuration.configurator.merge_dictionaries(new_dict, val)

    return new_dict


def create_config(args):
    """
    Utility to create a default configuration file.
    :param args:
    :return:
    """

    if args.full is True:
        default = mikado_lib.configuration.configurator.to_json("")
        del default["scoring"]
        del default["requirements"]
        config = default
    else:
        config = create_simple_config()

    if args.gff:
        args.gff = args.gff.split(",")
        config["prepare"]["gff"] = args.gff

        if args.labels != '':
            args.labels = args.labels.split(",")
            if not len(args.labels) == len(args.gff):
                raise ValueError("""Length mismatch between input files and labels!
                GFFs: {0} (length {1})
                Labels: {2} (length {3})""".format(
                    args.gff, len(args.gff),
                    args.labels, len(args.labels)))
            config["prepare"]["labels"] = args.labels

    comment = []
    comment_level = -1
    output = yaml.dump(config, default_flow_style=False)

    for line in output.split("\n"):
        # comment found
        if (line.lstrip().startswith("Comment") or
                line.lstrip().startswith("SimpleComment") or comment):
            level = sum(1 for _ in itertools.takewhile(str.isspace, line))
            if comment:
                if level > comment_level or line.lstrip().startswith("-"):
                    comment.append(line.strip())
                else:
                    for l in iter(_ for _ in comment if _ != ''):
                        print("{spaces}#  {comment}".format(spaces=" "*comment_level,
                                                            comment=re.sub("'", "", re.sub("^- ", "", l))),
                              file=args.out)
                    if level < comment_level:
                        print("{0}{{}}".format(" " * comment_level), file=args.out)
                    comment = []
                    comment_level = -1

                    print(line.rstrip(), file=args.out)
            else:
                if args.full is True:
                    comment = [re.sub("Comment:", "", line.strip())]
                else:
                    comment = [re.sub("SimpleComment:", "", line.strip())]
                comment_level = level
        else:
            print(line.rstrip(), file=args.out)

    if comment:
        for l in comment:
            print("{spaces}#{comment}".format(spaces=" "*comment_level, comment=l),
                  file=args.out)


def configure_parser():
    parser = argparse.ArgumentParser("Configuration utility")
    parser.add_argument("--full", action="store_true", default=False)
    parser.add_argument("--labels", type=str, default="",
                        help="""Labels to attach to the IDs of the transcripts of the input files,
                        separated by comma.""")
    parser.add_argument("--gff", help="Input GFF/GTF file(s), separated by comma", type=str)
    parser.add_argument("out", nargs='?', default=sys.stdout, type=argparse.FileType('w'))
    parser.set_defaults(func=create_config)
    return parser
