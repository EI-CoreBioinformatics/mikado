#!/usr/bin/env python3

from Mikado.parsers.GFF import GFF3
from Mikado.loci import Gene, Transcript
import sys
import argparse
from collections import defaultdict
from Mikado.utilities.log_utils import create_default_logger


__doc__ = """Script that tries to emulate GenomeTools sort"""


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("gff3", type=GFF3,
                        help="Input GFF3 file")
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?",
                        default=sys.stdout)
    args = parser.parse_args()

    genes = dict()
    indices = defaultdict(dict)
    tid2gid = dict()

    logger = create_default_logger("sorter", level="ERROR")
    previous = None

    print("##gff-version\t3", file=args.out)

    for row in args.gff3:
        if row.header:
            continue
        elif not (row.is_exon or row.is_gene or row.is_transcript):
            continue
        elif row.is_gene:
            gene = Gene(row, logger=logger)
            assert gene.id not in indices[row.chrom]
            genes[gene.id] = gene
            if previous is not None:
                genes[previous].finalize()
                key = (genes[previous].start, genes[previous].end, genes[previous].strand)
                if key not in indices[genes[previous].chrom]:
                    indices[genes[previous].chrom][key] = []
                indices[genes[previous].chrom][key].append(previous)
            previous = gene.id

        elif row.is_transcript:
            transcript = Transcript(row, logger=logger, trust_orf=True)
            if row.parent[0] not in genes:
                gene = Gene(transcript, logger=logger)
                genes[gene.id] = gene
                previous = gene.id
            tid2gid[transcript.id] = transcript.parent[0]
            genes[row.parent[0]].add(transcript)
        elif row.is_exon:
            found = False
            for parent in row.parent:
                if parent not in tid2gid:
                    continue
                genes[tid2gid[parent]][parent].add_exon(row)
                found = True
            if not found:
                raise KeyError("I have not found the parent for\n\t{}".format(row))

    if previous is not None:
        genes[previous].finalize()
        key = (genes[previous].start, genes[previous].end, genes[previous].strand)
        if key not in indices[genes[previous].chrom]:
            indices[genes[previous].chrom][key] = []
        indices[genes[previous].chrom][key].append(previous)

    for chrom in iter(sorted(indices.keys())):
        for pos in iter(sorted(indices[chrom].keys())):
            for gene in indices[chrom][pos]:
                print(genes[gene].format("gff3"), file=args.out)
                print("###", file=args.out)

    args.out.close()
    args.gff3.close()
    return

main()