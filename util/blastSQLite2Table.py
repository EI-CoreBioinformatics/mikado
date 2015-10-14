#!/usr/bin/env python3
# coding: utf-8

"""
Script to convert the BlastXML output into a tabular format.
"""

import sys
import argparse
from mikado_lib.serializers.blast_serializer import Hit
from mikado_lib.utilities.dbutils import connect
from mikado_lib.configuration.configurator import to_json
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker


def main():
    """
    Main script function. Just a thin wrapper over XmlSerializer.
    """

    parser = argparse.ArgumentParser(
        "Script to convert the BlastXML output into a tabular format.")
    db_settings = parser.add_mutually_exclusive_group("DB connection settings")
    db_settings.add_argument("--db", type=str, help="SQLite DB file to parse.")
    db_settings.add_argument("--json", type=to_json, default=None,
                             type=argparse.FileType,
                             help="Configuration file to use for connection settings.")
    parser.add_argument("out", type=argparse.FileType("w"), default=sys.stdout,
                        nargs='?', help="Optional output file. Default: %(default)s")

    args = parser.parse_args()

    if args.db is not None:
        engine = create_engine("sqlite:///{0}".format(args.db))
    else:
        engine = connect(args.json)
    session = sessionmaker()
    session.configure(bind=engine)
    current_session = session()

    for hit in current_session.query(Hit):
        print(hit)
        for hsp in hit.hsps:
            print("\t", hsp)

if __name__ == '__main__':
    main()
