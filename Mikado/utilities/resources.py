"""Helpers for accessing package data files."""

import os
from importlib import resources


def _resource_parts(*parts):
    for part in parts:
        for chunk in str(part).replace(os.sep, "/").split("/"):
            if chunk not in ("", "."):
                yield chunk


def resource(package, *path_parts):
    return resources.files(package).joinpath(*_resource_parts(*path_parts))


def resource_file(package, *parts):
    return str(resource(package, *parts))


def resource_binary_stream(package, *parts):
    return resource(package, *parts).open("rb")


def resource_exists(package, *parts):
    return resource(package, *parts).exists()


def resource_listdir(package, *parts):
    return [entry.name for entry in resource(package, *parts).iterdir()]
