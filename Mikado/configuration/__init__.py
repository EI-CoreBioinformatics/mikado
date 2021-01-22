"""
This module defines the functions needed to check the sanity of the configuration file,
plus the JSON schemas for the configuration and scoring files.
"""


from .configuration import MikadoConfiguration
from .daijin_configuration import DaijinConfiguration
from . import configurator
import itertools
import re
import textwrap
import pkg_resources
import json


__author__ = 'Luca Venturini'


def print_toml_config(output, out):

    import sys

    schema = json.load(pkg_resources.resource_stream("Mikado.configuration", "configuration_blueprint.json"))
    daijin_schema = json.load(pkg_resources.resource_stream("Mikado.configuration", "daijin_schema.json"))
    daijin_found = False

    lines = []

    level = schema
    for line in output.split("\n"):
        if line.startswith("["):
            keys = line.rstrip().replace("[", "").replace("]", "").split(".")
            level = schema
            for key in keys:
                if key not in level["properties"]:
                    level = daijin_schema
                    daijin_found = True
                if key not in level["properties"]:
                    raise KeyError("Unknown key found: {}".format(key))
                level = level["properties"][key]
            # print("We are at", key, "with level: ", level, file=sys.stderr)
            comment = []
            # title = level.get("title", None)
            # if title:
            #     comment += ["# " + _ for _ in textwrap.wrap(title)]
            description = level.get("description", None)
            if description:
                comment += ["# " + _ for _ in textwrap.wrap(description)]
            lines.append(line)
            lines.extend(comment)
        else:
            comment = []
            if "=" in line:
                key = line.split("=")[0].strip()
                if "Comment" in key:
                    raise KeyError((line, key, level))
                description = level["properties"].get(key, dict()).get("description", None)
                if description:
                    comment += ["# " + key + ": " + _ for _ in textwrap.wrap(description)]
                lines.extend(comment)
            lines.append(line.rstrip())

    if daijin_found is True:
        print(*["# " + _ for _ in textwrap.wrap(daijin_schema["title"])], sep="\n", file=out)
    else:
        print(*["# " + _ for _ in textwrap.wrap(schema["title"])], sep="\n", file=out)

    print(*lines, sep="\n", file=out)


def print_config(output, out, format="yaml"):

    """
    Function to print out the prepared configuration.
    :param output: prepared output, a huge string.
    :type output: str

    :param out: output handle.
    """

    comment = []
    comment_level = -1

    if format == "json":
        import json
        print(json.dumps(output, indent=4, sort_keys=True), file=out)

    # elif format == "yaml":

    for line in output.split("\n"):
        # comment found
        if line.lstrip().startswith(("Comment", "SimpleComment")) or comment:
            level = sum(1 for _ in itertools.takewhile(str.isspace, line))
            line = re.sub("Comment:", "", re.sub("SimpleComment:", "", line))
            line = re.sub("^- ", "", line)
            if comment:
                if level > comment_level or line.lstrip().startswith("-"):
                    comment.append(line.strip())
                else:
                    comment_line = " ".join([_ for _ in comment if _ != ''])
                    comment_line = re.sub("'", "", re.sub(" - ", "\n- ", comment_line))
                    for part in comment_line.split("\n"):
                        for part_line in textwrap.wrap(part.rstrip(), 80, replace_whitespace=True,
                                                      initial_indent=" "*comment_level + "# ",
                                                      subsequent_indent=" "*comment_level + "# "):
                            part_line = re.sub("# - ", "# ", part_line)
                            print(part_line, file=out)
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
            if new_key in schema:
                if "required" in schema[new_key] and schema[new_key]["required"] is True:
                    required.append([new_key])
        else:
            continue

    return required
