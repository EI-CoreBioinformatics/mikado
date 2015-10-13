#!/usr/bin/env python3

__author__ = 'Luca Venturini'

import yaml
import itertools
import re
import argparse
import sys
import mikado_lib

"""Stub of pre-configurer for mikado_lib"""


def create_config(args):
    """
    Utility to create a default configuration file.
    :param args:
    :return:
    """
    default = mikado_lib.configuration.configurator.to_json("")

    del default["scoring"]
    del default["requirements"]

    output = yaml.dump(default, default_flow_style=False)

    comment = []
    comment_level = -1

    for line in output.split("\n"):
        if line.lstrip().startswith("Comment") or comment: # comment found
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
    parser.add_argument("out", nargs='?', default=sys.stdout, type=argparse.FileType('w'))
    parser.set_defaults(func=create_config)
    return parser
