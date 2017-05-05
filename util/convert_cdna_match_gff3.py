#!/usr/bin/env python3

import Mikado
import argparse
import sys
from collections import defaultdict
import operator

__doc__ = """Script to convert a cDNA_match GFF3 from exonerate into either a match/match_part or a GTF file."""

def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-of", "--out-format", choices=["gtf", "match"], default="gtf", dest="format")
    parser.add_argument("gff", nargs="?", default="-")
    parser.add_argument("out", nargs="?", default=sys.stdout)
    args = parser.parse_args()

    if args.gff == "-":
        args.gff = sys.stdin
        pass
    else:
        args.gff = open(args.gff)
        pass
    
    if args.out == "-":
        args.out = sys.stdout
    elif not args.out is sys.stdout:
        args.out = open(args.out, "wt")
    
    matches = defaultdict(list)
    positions = dict()
    
    
    for line in args.gff:
        line = Mikado.parsers.GFF.GffLine(line)
        if line.header is True:
            continue
        matches[line.id].append(line)
        if line.id not in positions:
            positions[line.id] = (line.chrom, line.start, line.end)
        else:
            if positions[line.id][0] != line.chrom:
                raise ValueError("Chromosome mismatch!")
            positions[line.id] = (line.chrom,
                                  min(line.start, positions[line.id][1]),
                                  max(line.end, positions[line.id][2]))
        continue

    rev_positions = defaultdict(list)
    for tid, val in positions.items():
        rev_positions[val].append(tid)
    
    for pos in sorted(rev_positions.keys(), key=operator.itemgetter(0, 1, 2)):
        for match in rev_positions[pos]:
            transcript = Mikado.transcripts.Transcript()
            transcript.chrom, transcript.strand, transcript.id = matches[match][0].chrom, matches[match][0].strand, match
            transcript.attributes["Target"] = matches[match][0].attributes["Target"]
            transcript.add_exons(matches[match])
            transcript.finalize()
            if args.format == "gtf":
                transcript.parent = transcript.id
                print(transcript.format("gtf", with_cds=False), file=args.out)
                pass
            else:
                for line in transcript.format("gff", with_cds=False).split("\n"):
                    line = line.split("\t")
                    if line[2] == "exon":
                        line[2] = "match_part"
                    else:
                        line[2] = "match"
                    print(*line, sep="\t", file=args.out)
                    continue
                pass
            continue
        continue

    args.gff.close()
    args.out.close()
    return

main()
