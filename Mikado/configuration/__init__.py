"""
This module defines the functions needed to check the sanity of the configuration file,
plus the JSON schemas for the configuration and scoring files.
"""


from . import configurator
import itertools
import re


__author__ = 'Luca Venturini'


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
