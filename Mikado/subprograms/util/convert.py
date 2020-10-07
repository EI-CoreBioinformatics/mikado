#!/usr/bin/env python3
import argparse
import sys


def launch(args):

    from ...parsers import to_gff
    from ...parsers.bed12 import BED12
    from ...loci import Transcript, Gene

    if args.gf == "-":
        if args.in_format is None:
            raise ValueError("I need a format if it cannot be inferred from the string")
        parser = to_gff(sys.stdin, input_format=args.in_format)
    else:
        parser = to_gff(args.gf, input_format=args.in_format)

    current = None
    if parser.__annot_type__ == "gtf":
        out_format = "gff3"
    elif parser.__annot_type__ in ("bed12", "gff3", "bam"):
        out_format = "gtf"
    else:
        raise TypeError("Invalid annotation type: {}".format(parser.__annot_type__))

    if (args.out_format is not None and
                args.out_format != out_format):
        if args.out_format != parser.__annot_type__ or args.transcriptomic is True:
            out_format = args.out_format
        else:
            print("Input and output format are the same, aborting.", file=sys.stderr)
            sys.exit(1)

    if out_format == "gff3":
        print("##gff-version\t3", file=args.out)

    mock_gene_counter = 0

    for line in parser:
        if line.header is True:
            continue
        if parser.__annot_type__ == "bam":
            # BAM file. We need to do things a little bit differently
            if line.is_unmapped is True:
                continue

            mock_gene_counter += 1
            gene = "gene_{mock_gene_counter}".format(**locals())
            transcript = Transcript(line)
            transcript.parent = gene
            print(Gene(transcript).format(out_format, transcriptomic=args.transcriptomic), file=args.out)
            continue

        elif isinstance(line, BED12) or (line.is_exon is False and line.gene is None):
            mock_gene_counter += 1
            gene = "gene_{mock_gene_counter}".format(**locals())
            line.parent = gene

        if current is None or ((line.is_gene or line.is_transcript) and line.gene != current.id):
            if current:
                print(current.format(out_format, transcriptomic=args.transcriptomic), file=args.out)
            if hasattr(line, "feature") and "superlocus" in line.feature:  # Hack for Mikado files
                continue
            elif line.is_gene is True:
                current = Gene(line)
            elif hasattr(line, "feature") and line.feature in ("chromosome", "region"):
                continue
            else:
                if line.parent is None:
                    line.parent = "{}.gene".format(line.id)  # Hack for BED12 files
                current = Gene(Transcript(line))
        elif parser.__annot_type__ == "gtf" and line.is_exon is True and (
                        current is None or current.id != line.gene):
            if current:
                print(current.format(out_format, transcriptomic=args.transcriptomic), file=args.out)
            current = Gene(Transcript(line))
        elif (parser.__annot_type__ == "gtf" and
                      line.is_exon is True and line.transcript not in current.transcripts):
            current.add(Transcript(line))
        elif line.is_transcript is True:
            current.add(Transcript(line))
        elif line.is_exon is True:
            current.add_exon(line)

    if current is not None:
        print(current.format(out_format, transcriptomic=args.transcriptomic), file=args.out)


def convert_parser():
    """
    Parser for the command line
    :return:
    """

    parser = argparse.ArgumentParser(
        "Utility to covert across GTF, GFF3 and BED12.")
    parser.add_argument("-of", "--out-format", dest="out_format",
                        choices=["bed12", "gtf", "gff3"], default=None)
    parser.add_argument("-if", "--in-format", dest="in_format",
                        choices=["bed12", "gtf", "gff3", "bam"], default=None)
    parser.add_argument("-t", "--transcriptomic", default=False, action="store_true",
                        help="Flag. If on, the file will be converted to a transcriptomic version.")
    parser.add_argument("gf")
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?", default=sys.stdout)
    parser.set_defaults(func=launch)
    return parser
