#!/usr/bin/env python3
# coding: utf-8

"""Script to serialise junctions into the database."""

import argparse

from Mikado.serializers import junction
from Mikado.configuration import json_utils


def main():
    """
    Main script function.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--fai", default=None)
    parser.add_argument("--max-objects", dest="max_objects", default=10 ** 5, type=int)
    parser.add_argument("--json-conf", default=None, dest="json_conf", type=json_utils.to_json)
    parser.add_argument("bed12")
    parser.add_argument("db", nargs="?", default=":memory:")
    args = parser.parse_args()

    serializer = junction.JunctionSerializer(args.bed12, args.db, fai=args.fai, json_conf=args.json_conf,
                                             maxobjects=args.max_objects)
    serializer()


if __name__ == "__main__":
    main()
