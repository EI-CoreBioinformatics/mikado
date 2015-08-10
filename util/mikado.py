#!/usr/bin/env python3

import argparse
import mikado_lib.subprograms


def all():
    """
    Mock. It will launch the whole pipeline.
    """
    parser = argparse.ArgumentParser()
    return parser


def main():

    """
    Main launcher function for the pipeline.
    """

    parser = argparse.ArgumentParser(prog="Mikado")
    subparsers = parser.add_subparsers(help="Available programs in the suite.")
    subparsers.choices["compare"] = mikado_lib.subprograms.compare.compare_parser()
    subparsers.choices["pick"] = mikado_lib.subprograms.pick.pick_parser()
    subparsers.choices["prepare"] = mikado_lib.subprograms.prepare.prepare_parser()
    subparsers.choices["serialise"] = mikado_lib.subprograms.serialise.serialise_parser()
    subparsers.choices["all"] = all()

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_usage()

if __name__ == '__main__':
    main()
