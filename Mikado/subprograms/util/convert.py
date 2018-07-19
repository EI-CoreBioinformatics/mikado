#!/usr/bin/env python3
import argparse
import sys
from ...parsers import to_gff
from ...loci import Transcript, Gene


def launch(args):

    args.gf.close()
    args.gf = args.gf.name

    parser = to_gff(args.gf)
    current = None
    if parser.__annot_type__ == "gtf":
        out_format = "gff3"
    elif parser.__annot_type__ in ("bed12", "gff3"):
        out_format = "gtf"
    else:
        raise TypeError("Invalid annotation type: {}".format(parser.__annot_type__))

    if (args.out_format is not None and
                args.out_format != out_format):
        if args.out_format != parser.__annot_type__:
            out_format = args.out_format
        else:
            print("Input and output format are the same, aborting.", file=sys.stderr)
            sys.exit(1)

    if out_format == "gff3":
        print("##gff-version\t3", file=args.out)

    for line in parser:
        if line.header is True:
            continue
        if current is None or ((line.is_gene or line.is_transcript) and line.gene != current.id):
            if current:
                print(current.format(out_format), file=args.out)

            if hasattr(line, "feature") and "superlocus" in line.feature:  # Hack for Mikado files
                continue
            elif line.is_gene is True:
                current = Gene(line)
            else:
                if line.parent is None:
                    line.parent = "{}.gene".format(line.id)  # Hack for BED12 files
                current = Gene(Transcript(line))
        elif parser.__annot_type__ == "gtf" and line.is_exon is True and (
                        current is None or current.id != line.gene):
            if current:
                print(current.format(out_format), file=args.out)
            current = Gene(Transcript(line))
        elif (parser.__annot_type__ == "gtf" and
                      line.is_exon is True and line.transcript not in current.transcripts):
            current.add(Transcript(line))
        elif line.is_transcript is True:
            current.add(Transcript(line))
        elif line.is_exon is True:
            current.add_exon(line)
    if current is not None:
        print(current.format(out_format), file=args.out)


def convert_parser():
    """
    Parser for the command line
    :return:
    """

    parser = argparse.ArgumentParser(
        "Utility to covert across GTF, GFF3 and BED12.")
    parser.add_argument("-of", "--out-format", dest="out_format",
                        choices=["bed12", "gtf", "gff3"], default=None)
    parser.add_argument("gf", type=argparse.FileType())
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?", default=sys.stdout)
    parser.set_defaults(func=launch)
    return parser
