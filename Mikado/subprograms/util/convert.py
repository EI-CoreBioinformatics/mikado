#!/usr/bin/env python3
import argparse
import sys
from ...exceptions import InvalidTranscript, InvalidParsingFormat
from ...parsers import parser_factory
from ...parsers.GFF import GFF3
from ...parsers.GTF import GTF
from ...parsers.bed12 import BED12, Bed12Parser
from ...loci import Transcript, Gene


def _convert_bam(parser, args, out_format):
    current = None
    mock_gene_counter = 0
    _line_counter = 0
    for line in parser:
        if line.is_unmapped is True:
            continue

        mock_gene_counter += 1
        gene = "gene_{mock_gene_counter}".format(**locals())
        try:
            transcript = Transcript(line)
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            raise InvalidTranscript("Invalid BAM line at position {} for file {}. Error: {}".format(
                _line_counter, parser.name, exc
            ))
        transcript.parent = gene
        print(Gene(transcript).format(out_format, transcriptomic=args.transcriptomic), file=args.out)
    return


def _convert_gtf(parser: GTF, args, out_format):
    genes = dict()
    current = None
    found_ids = set()

    for line_counter, line in enumerate(parser, start=1):
        if line.header is True:
            continue
        elif hasattr(line, "feature") and "superlocus" in line.feature or line.feature in ("chromosome", "region"):
            continue
        if line.gene in genes:
            if line.is_transcript is True and line.transcript in genes[line.gene]:
                if args.assume_sorted is True:
                    raise InvalidParsingFormat(
                        f"The transcript at line no. {line_counter} has already appeared: the input is not sorted.")
                gene = genes[line.gene]
                assert isinstance(gene, Gene)
                assert line.strand == gene.strand
                assert line.chrom == gene.chrom
                gene.transcripts[line.transcript].start, gene.transcripts[line.transcript].end = line.start, line.end
                gene.transcripts[line.transcript]._set_expandable(False)
                gene.transcripts[line.transcript].attributes.update(line.attributes)
            elif line.is_transcript is True:
                if args.assume_sorted is True and line.transcript in found_ids:
                    raise InvalidParsingFormat(
                        f"The transcript at line no. {line_counter} has already appeared: the input is not sorted.")
                genes[line.gene].add(Transcript(line))
            elif line.is_gene:
                if line.gene in genes or line.gene in found_ids:
                    if args.assume_sorted is True:
                        raise InvalidParsingFormat(
                            f"The transcript at line no. {line_counter} has already appeared: the input is not sorted.")
                genes[line.gene].start, genes[line.gene].end = line.start, line.end
                genes[line.gene].attributes.update(line.attributes)
            elif line.is_exon:
                if line.transcript in genes[line.gene].transcripts:
                    genes[line.gene].add_exon(line)
                else:
                    genes[line.gene].add(Transcript(line))
        else:
            if line.is_gene:
                if args.assume_sorted is True and line.gene in genes:
                    raise InvalidParsingFormat(
                        f"The gene at line no. {line_counter} has already appeared: the input is not sorted.")
                genes[line.gene] = Gene(line)
            else:
                if args.assume_sorted is True and line.gene in found_ids or line.transcript in found_ids:
                    raise InvalidParsingFormat(
                        f"The gene at line no. {line_counter} has already appeared: the input is not sorted.")
                genes[line.gene] = Gene(Transcript(line))

        found_ids.add(line.transcript)
        if current and current != line.gene and args.assume_sorted is True:
            print(genes[current].format(out_format, transcriptomic=args.transcriptomic), file=args.out)
            del genes[current]

        current = line.gene

    if args.assume_sorted is False:
        for gid, gene in genes.items():
            print(gene.format(out_format, transcriptomic=args.transcriptomic), file=args.out)
    elif current is not None:
        print(genes[current].format(out_format, transcriptomic=args.transcriptomic), file=args.out)

    return


def _convert_bed12(parser: Bed12Parser, args, out_format):
    mock_gene_counter = 0
    for line in parser:
        if line.header is True:
            continue
        gid = "gene_{mock_gene_counter}".format(**locals())
        line.parent = gid
        gene = Gene(Transcript(line))
        print(gene.format(out_format, transcriptomic=args.transcriptomic), file=args.out)

    return


def _convert_gff(parser: GFF3, args, out_format):
    orphaned = dict()
    tid2gene = dict()
    genes = dict()

    current = None
    for line in parser:
        if line.header is True:
            continue
        elif line.feature in ("region", "chrom", "superlocus", "protein"):
            continue
        # We need to have better checks here regarding sorted/unsorted inputs here
        if line.is_gene is False and args.assume_sorted is True:
            if line.is_transcript is False:
                if current is None or line.parent[0] not in genes[current].transcripts:
                    raise InvalidParsingFormat("The provided GFF is not sorted.")
                genes[current].add_exon(line)
            else:
                assert genes[current].id in line.parent[0]
                genes[current].add(Transcript(line))
        elif line.is_gene is True:
            if args.assume_sorted is True and current is not None:
                print(genes[current].format(out_format, transcriptomic=args.transcriptomic), file=args.out)
                del genes[current]
            genes[line.gene] = Gene(line)
            current = line.gene
            if args.assume_sorted is False:
                orph = list(orphaned.keys())
                for tid in orph:
                    if orphaned[tid].parent and orphaned[tid].parent[0] == current:
                        genes[current].add(orphaned[tid])
                        del orphaned[tid]
        elif line.is_transcript is True and args.assume_sorted is False:
            transcript = Transcript(line)
            if transcript.parent:
                tid2gene[transcript.id] = transcript.parent[0]
            if transcript.id in orphaned:
                assert transcript.chrom == orphaned[transcript.id].chrom
                assert transcript.strand == orphaned[transcript.id].strand
                orphaned[transcript.id].start, orphaned[transcript.id].end = transcript.start, transcript.end
                orphaned[transcript.id].attributes.update(transcript.attributes)
                orphaned[transcript.id]._set_expandable(False)
                orphaned[transcript.id].parent = transcript.gene
                if orphaned[transcript.id].parent and orphaned[transcript.id].parent[0] in genes:
                    genes[orphaned[transcript.id].parent[0]].add(orphaned[transcript.id])
                    del orphaned[transcript.id]
            elif transcript.parent and transcript.parent[0] in genes:
                genes[transcript.parent[0]].add(transcript)
            else:
                orphaned[transcript.id] = transcript
        else:
            assert args.assume_sorted is False
            # Exons might have multiple parents ...
            if any(_ in tid2gene for _ in line.parent):
                gids = set([tid2gene[_] for _ in line.parent])
                assert len(gids) == 1, line
                gid = gids.pop()
                if gid in genes:
                    genes[gid].add_exon(line)
                    for parent in line.parent:
                        if parent not in genes[gid]:
                            genes[gid].add(Transcript(line))
                else:
                    for parent in line.parent:
                        if parent in orphaned:
                            orphaned[parent].add_exon(line)
                        else:
                            orphaned[parent] = Transcript(line)
            else:
                for parent in line.parent:
                    if parent in orphaned:
                        orphaned[parent].add_exon(line)
                    else:
                        orphaned[parent] = Transcript(line)

    if current and args.assume_sorted is True:
        print(genes[current].format(out_format, transcriptomic=args.transcriptomic), file=args.out)

    if args.assume_sorted is False:
        for tid in orphaned:
            if not orphaned[tid].parent:
                orphaned[tid].parent = "{}.gene".format(orphaned[tid].id)
                genes[orphaned[tid].parent[0]] = orphaned[tid]
            else:
                for parent in orphaned[tid].parent:
                    if parent in genes:
                        genes[parent].add(orphaned[tid])
                    else:
                        genes[parent] = Gene(orphaned[tid])

        for gid, gene in genes.items():
            print(gene.format(out_format, transcriptomic=args.transcriptomic), file=args.out)

    return


def launch(args):

    if args.gf == "-":
        if args.in_format is None:
            raise ValueError("I need a format if it cannot be inferred from the string")
        parser = parser_factory(sys.stdin, input_format=args.in_format)
    else:
        parser = parser_factory(args.gf, input_format=args.in_format)

    if parser.__annot_type__ == "gtf":
        out_format = "gff3"
    elif parser.__annot_type__ in ("bed12", "gff3", "bam"):
        out_format = "gtf"
    else:
        raise TypeError("Invalid annotation type: {}".format(parser.__annot_type__))

    if args.out_format is not None and args.out_format != out_format:
        if args.out_format != parser.__annot_type__ or args.transcriptomic is True:
            out_format = args.out_format
        else:
            print("Input and output format are the same, aborting.", file=sys.stderr)
            sys.exit(1)

    if out_format == "gff3":
        print("##gff-version\t3", file=args.out)

    if parser.__annot_type__ == "bam":
        _convert_bam(parser, args, out_format)
    elif parser.__annot_type__ == "gtf":
        _convert_gtf(parser, args, out_format)
    elif parser.__annot_type__ == "bed12":
        _convert_bed12(parser, args, out_format)
    else:
        _convert_gff(parser, args, out_format)

    if not args.out == sys.stdout:
        args.out.flush()
        args.out.close()


def convert_parser():
    """
    Parser for the command line
    :return:
    """

    parser = argparse.ArgumentParser(
        "Utility to covert across GTF, GFF3 and BED12.")
    parser.add_argument("-as", "--assume-sorted", default=False, action="store_true")
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
