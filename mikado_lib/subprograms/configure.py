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

    if args.simple is False:
        default = mikado_lib.configuration.configurator.to_json("")
        del default["scoring"]
        del default["requirements"]
        output = yaml.dump(default, default_flow_style=False)
    else:
        output = yaml.dump(create_simple_config(), default_flow_style=False)

    comment = []
    comment_level = -1

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
                    comment = []
                    comment_level = -1
                    print(line.rstrip(), file=args.out)
            else:
                if args.simple is True:
                    comment = [re.sub("SimpleComment:", "", line.strip())]
                else:
                    comment = [re.sub("Comment:", "", line.strip())]
                comment_level = level
        else:
            print(line.rstrip(), file=args.out)

    if comment:
        for l in comment:
            print("{spaces}#{comment}".format(spaces=" "*comment_level, comment=l),
                  file=args.out)


def configure_parser():
    parser = argparse.ArgumentParser("Configuration utility")
    parser.add_argument("--simple", action="store_true", default=False)
    parser.add_argument("out", nargs='?', default=sys.stdout, type=argparse.FileType('w'))
    parser.set_defaults(func=create_config)
    return parser
