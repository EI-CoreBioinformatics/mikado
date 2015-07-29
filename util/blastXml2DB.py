#!/usr/bin/env python3
# coding: utf-8

"""
Script to conver the BlastXML output into a tabular format.
"""

import os
import argparse
from mikado_lib.serializers.blast_utils import XmlSerializer
from mikado_lib import json_utils
from Bio import SeqIO


def to_seqio(string):
    """
    Convert a string to a SeqIO index.

    :param string
    :type string: str
    """
    assert os.path.exists(string) and os.path.isfile(string) and os.stat(string).st_size > 0
    return SeqIO.index(string, "fasta")


def main():
    """
    Main script function.
    """

    parser = argparse.ArgumentParser("Script to conver the BlastXML output into a tabular format.")
    parser.add_argument("--max_target_seqs", type=int, default=float("Inf"),
                        help="Maximum number of target sequences.")
    parser.add_argument("--maxobjects", type=int, default=10**5, help="Maximum number of objects to cache in memory.")
    parser.add_argument("--definition", action="store_true", default=False,
                        help="Use query def instead of ID for the output.")
    parser.add_argument("--query_seqs", default=None, type=to_seqio, help="Query sequences")
    parser.add_argument("--target_seqs", default=None, type=to_seqio, help="Target sequences")
    parser.add_argument("--json-conf", dest="json_conf", default=None, type=json_utils.to_json)
    parser.add_argument("xml", type=str, help="XML file to parse.")
    parser.add_argument("dbout", type=str, default=":memory:",
                        nargs='?', help="Optional output file. Default: :memory:")

    args = parser.parse_args()

    XmlSerializer(args.dbout, args.xml,
                  keep_definition=args.definition,
                  max_target_seqs=args.max_target_seqs,
                  maxobjects=args.maxobjects,
                  target_seqs=args.target_seqs,
                  query_seqs=args.query_seqs,
                  json_conf=args.json_conf
                  )()


if __name__ == '__main__':
    main()
