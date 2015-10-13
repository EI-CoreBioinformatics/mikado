#!/usr/bin/env python3
# coding: utf-8

"""Script to trim down the terminal exons of multiexonic transcripts"""

import sys
import argparse
import Mikado.parsers
import Mikado.loci_objects.transcript


def strip_terminal(transcript, args) -> Mikado.loci_objects.transcript.Transcript:
    """This function will take as input a transcript and then:
    - return immediately if the transcript is monoexonic
    - trim the terminal exons to a length of max(args.max_length (if longer), CDS boundary)

    :param transcript: the input transcript
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param args: argparse configuration namespace
    :type args: argparse.Namespace

    :rtype : Mikado.loci_objects.transcript.Transcript
    """

    transcript.finalize()
    # Return immediately if the transcript is monoexonic
    if transcript.monoexonic is True:
        return transcript

    transcript.finalized = False

    if transcript.selected_cds_length == 0:
        # Non-coding transcript: trim the terminal exons and the return
        first = list(transcript.exons[0])
        if (first[1] - first[0] + 1) > args.max_length:
            newfirst = list(first)
            newfirst[0] = first[1] - args.max_length
            transcript.start = newfirst[0]
            transcript.exons[0] = tuple(newfirst)
        last = transcript.exons[-1]
        if (last[1] - last[0] + 1) > args.max_length:
            newlast = list(last[:])
            newlast[1] = last[0] + args.max_length
            transcript.end = newlast[1]
            transcript.exons[-1] = tuple(newlast)
            if transcript.selected_cds_length > 0:
                last_utr = list(filter(lambda uu: uu[1] == last[1], transcript.combined_utr))[0]
                transcript.combined_utr.remove(last_utr)
                transcript.combined_utr.append(tuple(newlast))

    else:
        # Coding transcript
        # Order cds_start and end irrespectively of strand
        cds_start, cds_end = sorted([transcript.selected_cds_end, transcript.selected_cds_start])
        exons = transcript.exons
        first = exons[0]
        last = exons[-1]
        if first[1] < cds_start:  # first exon is only UTR
            if first[1] - first[0] + 1 > args.max_length:
                # First exon longer than max_length; trim it, remove from UTR and exons, substitute with trimmed exon
                newfirst = (first[1] - args.max_length, first[1])
                transcript.combined_utr.remove(first)
                transcript.exons.remove(first)
                transcript.combined_utr.append(newfirst)
                transcript.exons.append(newfirst)
            else:
                # Leave as it is
                newfirst = first
        elif first[0] < cds_start <= first[1]:
            # Partly UTR partly CDS
            if first[1] - first[0] + 1 > args.max_length:
                newfirst = (min(cds_start, first[1] - args.max_length), first[1])
                u = list(filter(lambda uu: uu[0] == first[0], transcript.combined_utr))[
                    0]  # Retrieve and remove UTR segment
                transcript.combined_utr.remove(u)
                transcript.exons.remove(first)
                transcript.exons.append(newfirst)
                if newfirst[0] < cds_start:
                    # Create new UTR segment
                    newu = (newfirst[0], cds_start - 1)
                    transcript.combined_utr.append(newu)
            else:
                newfirst = first
        else:  # Beginning of first exon == beginning of CDS
            newfirst = first
        transcript.start = newfirst[0]

        if last[0] > cds_end:
            if last[1] - last[0] + 1 > args.max_length:
                newlast = (last[0], last[0] + args.max_length)
                transcript.combined_utr.remove(last)
                transcript.exons.remove(last)
                transcript.combined_utr.append(newlast)
                transcript.exons.append(newlast)
            else:
                newlast = last
        elif last[0] <= cds_end < last[1]:
            if last[1] - last[0] + 1 > args.max_length:
                newlast = (last[0], max(cds_end, last[0] + args.max_length))
                u = list(filter(lambda uu: uu[1] == last[1], transcript.combined_utr))[0]
                transcript.combined_utr.remove(u)
                transcript.exons.remove(last)
                transcript.exons.append(newlast)
                if newlast[1] > cds_end:
                    newu = (cds_end + 1, newlast[1])
                    transcript.combined_utr.append(newu)
            else:
                newlast = last
        else:
            newlast = last

        transcript.end = newlast[1]

    transcript.finalize()
    return transcript


def main():
    """
    Main script function.
    """
    def to_ann(string):
        """
         Function to select the right parser for the input file.
        :param string: the filename
        :type string: str
        """
        if string.endswith("gtf"):
            return Mikado.parsers.GTF.GTF(string)
        elif string.endswith("gff") or string.endswith("gff3"):
            return Mikado.parsers.GFF.GFF3(string)
        else:
            raise ValueError("Unrecognized format")

    parser = argparse.ArgumentParser("Script to trim down the terminal exons of multiexonic transcripts")
    parser.add_argument("-ml", "--max_length", type=int, default=50, help="Maximal length of trimmed terminal exons")
    parser.add_argument("--as_gtf", default=False, action="store_true",
                        help="Flag. If set, the output will be in GTF rather than GFF3 format.")
    parser.add_argument("ann", type=to_ann, help="Reference GTF/GFF output file.")
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType('w'))
    args = parser.parse_args()

    if args.as_gtf is False:
        print("##gff-version 3", file=args.out)

    transcript = None

    for record in args.ann:
        if record.is_transcript is True:
            if transcript is not None:
                print(strip_terminal(transcript, args).__str__(to_gtf=args.as_gtf), file=args.out)
            transcript = Mikado.loci_objects.transcript.Transcript(record)
        elif record.is_exon is True:
            transcript.add_exon(record)

    print(strip_terminal(transcript, args).__str__(to_gtf=args.as_gtf), file=args.out)


if __name__ == "__main__":
    main()
