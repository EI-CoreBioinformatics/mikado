#!/usr/bin/env python3

from Mikado.loci import Transcript, Gene
import argparse
import sys
from Mikado.parsers.GFF import GFF3
from Mikado.utilities import to_gff
from intervaltree import Interval
from Mikado.utilities.log_utils import create_default_logger

__doc__ = """This script truncates transcript with UTR exons separated by long introns."""


def remove_introns_from_transcr(transcript, args):
    transcript.finalize()
    if transcript.monoexonic is True:
        return transcript
    to_remove = set()

    cstart, cend = sorted([transcript.combined_cds_start, transcript.combined_cds_end])
    candidate_exons = [exon for exon in transcript.exons if
                       exon[1] < transcript.combined_cds_start or
                       exon[0] > transcript.combined_cds_end]
    
    # print(candidate_exons)
    if len(candidate_exons) == 0:
        return transcript
    
    for intron in transcript.introns - transcript.combined_cds_introns:
        if intron[1] - intron[0] + 1 >= args.max_intron:
            # print(transcript.id, "intron", intron, intron[1] - intron[0] + 1, file=sys.stderr)
            if intron[1] < cstart:
                to_remove.update([exon for exon in candidate_exons if exon[1] < intron[0]])
            elif intron[0] > cend:
                to_remove.update([exon for exon in candidate_exons if exon[0] > intron[1]])
            else:
                raise ValueError

    if len(to_remove) > 0:
        # cds = [_ for _ in transcript.segments if _[0] == "CDS"]
        # transcript.strip_cds()
        # transcript.combined_utr = []
        # transcript.segments = []
        transcript.finalized = False
        # transcript.internal_orfs[0].extend(cds)
        for exon in to_remove:
            transcript.remove_exon(Interval(exon[0], exon[1]))
        # print(cds)
        # for exon in cds:
        #     transcript.add_exon(exon, "CDS")
        # transcript.add_exons([(_[1][0], _[1][1]) for _ in cds], features="CDS")
        # transcript.internal_orfs[0].extend(cds)
        transcript.start = min(_[0] for _ in transcript.exons)
        transcript.end = max(_[1] for _ in transcript.exons)
        # print(transcript.internal_orfs)
        transcript.finalize()
        assert transcript.combined_cds_length > 0
    return transcript

def remove_introns(current, args):
    assert len(current.transcripts) > 0
    for tid in current.transcripts:
        transcript = current.transcripts[tid]
        if transcript.monoexonic:
            continue
        current.transcripts[tid] = remove_introns_from_transcr(transcript, args)
        # assert isinstance(current.transcripts[tid], Transcript)

    assert len(current.transcripts) > 0
    current.start = min(_.start for _ in current)
    current.end = max(_.end for _ in current)
    current.finalize()
    return current
    

def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-mi", "--max-intron", default=10000, dest="max_intron",
                        type=int, help="Maximum intron length for UTR introns.")
    parser.add_argument("gff", type=to_gff)
    parser.add_argument("out", default=sys.stdout, type=argparse.FileType("wt"), nargs="?")
    args = parser.parse_args()

    if args.max_intron < 0:
        raise ValueError("Max intron length <0 specified! {0}".format(args.max_intron))

    ref_gff = isinstance(args.gff, GFF3)
    if ref_gff:
        form = "gff3"
    else:
        form = "gtf"

    current = None
    current_transcript = None

    last_header = []
    for record in args.gff:
        if record.header is True:
            # print(record, file=sys.stderr)
            if current is not None:
                current = remove_introns(current, args)
                print(current.format(form), file=args.out)
                print(*last_header, sep="\n", file=args.out, end='')
                current = None
                current_transcript = None
            print(*last_header, sep="\n", end="")
            last_header = [record]
            continue
        if record.feature not in ("gene", "mRNA", "CDS", "exon"):
            continue
        
        if record.is_gene is True and ref_gff:
            print(record, file=sys.stderr)
            last_header = []
            if current is not None:
                # current = remove_introns(current, args)
                print(current.format(form), file=args.out)
                print(*last_header, sep="\n", file=args.out, end='')
                current = None
                current_transcript = None
        if record.is_transcript:
            if ref_gff is False:
                if current_transcript is not None:
                    # current_transcript = remove_introns_from_transcr(current_transcript,
                    #                                                 args)
                    assert current_transcript.combined_cds_length > 0
                    print(current_transcript, file=args.out)
                    print(*last_header, sep="\n", file=args.out)
                    last_header = []
            elif ref_gff is True:
                if current_transcript is not None:
                    if current is None:
                        current = Gene(current_transcript)
                        current.add(current_transcript)
                    else:
                        assert current_transcript.parent[0] != current.id
                        current.add(current_transcript)
                    # if current.id == current_transcript.parent[0]:
                        
                    # else:
                    #     current = remove_introns(current, args)
                    #     print(current.format(form), file=args.out)
                    #     print("###", file=args.out)
                    #     current = None
                
                # elif current_transcript is not None:
                #     current = Gene(current_transcript)

            current_transcript = Transcript(record)
        elif record.is_exon:
            if record.feature not in ("CDS", "exon"):
                continue
            current_transcript.add_exon(record)
        else:
            continue
        continue
    
    if ref_gff and current is not None:
        print(*last_header, sep="\n", file=args.out)
        last_header = []
        current = remove_introns(current, args)
        print(current.format(form), file=args.out)
    elif not ref_gff and current_transcript is not None:
        current_transcript = remove_introns_from_transcr(current_transcript,
                                                         args)
        print(current_transcript.format(form), file=args.out)
        print(*last_header, sep="\n", file=args.out, end='')
        
main()
