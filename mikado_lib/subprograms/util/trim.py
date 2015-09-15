#!/usr/bin/env python3
# coding: utf-8

"""Script to trim down the terminal exons of multiexonic transcripts"""

__author__ = 'Luca Venturini'

import sys
import argparse
import mikado_lib.parsers
import mikado_lib.loci_objects.transcript
from mikado_lib.subprograms import to_gff


def trim_noncoding(transcript, args):
    """
    Function to trim the terminal exons for non-coding transcripts,
     i.e. the simplest case.
    :param transcript:
    :return: transcript
    """

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
            last_utr = [segment for segment in transcript.combined_utr if
                        segment[1] == last[1]][0]
            transcript.combined_utr.remove(last_utr)
            transcript.combined_utr.append(tuple(newlast))
    return transcript


def trim_start(transcript, cds_start, max_length=0):

    """
    Method to trim the beginning of a transcript, given both
    the cds_start and the maximal allowed length after trimming.
    :param transcript:
    :param cds_start:
    :param max_length:
    :return: transcript
    """

    first = transcript.exons[0]
    if first[1] < cds_start:  # first exon is only UTR
        if first[1] - first[0] + 1 > max_length:
            # First exon longer than max_length; trim it,
            # remove from UTR and exons, substitute with trimmed exon
            newfirst = (first[1] - max_length, first[1])
            transcript.combined_utr.remove(first)
            transcript.exons.remove(first)
            transcript.combined_utr.append(newfirst)
            transcript.exons.append(newfirst)
        else:
            # Leave as it is
            newfirst = first
    elif first[0] < cds_start <= first[1]:
        # Partly UTR partly CDS
        if first[1] - first[0] + 1 > max_length:
            newfirst = (min(cds_start, first[1] - max_length), first[1])
            # Retrieve and remove offending UTR segment
            utr_segment = [segment for segment in transcript.combined_utr
                           if segment[0] == first[0]][0]
            transcript.combined_utr.remove(utr_segment)
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
    return transcript


def trim_end(transcript, cds_end, max_length=0):
    """
    Method to trim the end of a transcript, given both
    the cds_end and the maximal allowed length after trimming.
    :param transcript:
    :param cds_start:
    :param max_length:
    :return: transcript
    """

    last = transcript.exons[-1]
    if last[0] > cds_end:
        if last[1] - last[0] + 1 > max_length:
            newlast = (last[0], last[0] + max_length)
            transcript.combined_utr.remove(last)
            transcript.exons.remove(last)
            transcript.combined_utr.append(newlast)
            transcript.exons.append(newlast)
        else:
            newlast = last
    elif last[0] <= cds_end < last[1]:
        if last[1] - last[0] + 1 > max_length:
            newlast = (last[0], max(cds_end, last[0] + max_length))
            utr_segment = [segment for segment in transcript.combined_utr if
                           segment[1] == last[1]][0]
            transcript.combined_utr.remove(utr_segment)
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
    return transcript


def trim_coding(transcript, args):
    """
    Function to trim the terminal exons for coding transcripts,
     i.e. the more complex case as we have to account for CDS
     as well.
    :param transcript:
    :return: transcript
    """
    # Coding transcript
    # Order cds_start and end irrespectively of strand
    cds_start, cds_end = sorted([transcript.selected_cds_end, transcript.selected_cds_start])

    transcript = trim_start(transcript, cds_start,
                            max_length=args.max_length)
    transcript = trim_end(transcript, cds_end,
                          max_length=args.max_length)

    return transcript


def strip_terminal(transcript, args) -> mikado_lib.loci_objects.transcript.Transcript:
    """This function will take as input a transcript and then:
    - return immediately if the transcript is monoexonic
    - trim the terminal exons to a length of max(args.max_length (if longer), CDS boundary)

    :param transcript: the input transcript
    :type transcript: mikado_lib.loci_objects.transcript.Transcript

    :param args: argparse configuration namespace
    :type args: argparse.Namespace

    :rtype : mikado_lib.loci_objects.transcript.Transcript
    """

    transcript.finalize()
    # Return immediately if the transcript is monoexonic
    if transcript.monoexonic is True:
        return transcript

    transcript.finalized = False

    if transcript.selected_cds_length == 0:
        transcript = trim_noncoding(transcript, args)
    else:
        transcript = trim_coding(transcript, args)

    transcript.finalize()
    return transcript


def launch(args):

    """
    Launcher of the program.
    """

    if args.as_gtf is False:
        print("##gff-version 3", file=args.out)

    transcript = None

    for record in args.ann:
        if record.is_transcript is True:
            if transcript is not None:
                print(strip_terminal(transcript, args).__str__(to_gtf=args.as_gtf), file=args.out)
            transcript = mikado_lib.loci_objects.transcript.Transcript(record)
        elif record.is_exon is True:
            transcript.add_exon(record)

    print(strip_terminal(transcript, args).__str__(to_gtf=args.as_gtf), file=args.out)


def trim_parser():

    """
    Command line parser for the utility.
    """

    parser = argparse.ArgumentParser(
        "Script to trim down the terminal exons of multiexonic transcripts")
    parser.add_argument("-ml", "--max_length", type=int, default=50,
                        help="Maximal length of trimmed terminal exons")
    parser.add_argument("--as_gtf", default=False, action="store_true",
                        help="Flag. If set, the output will be in GTF rather than GFF3 format.")
    parser.add_argument("ann", type=to_gff, help="Reference GTF/GFF output file.")
    parser.add_argument("out", nargs="?", default=sys.stdout,
                        type=argparse.FileType('w'))
    parser.set_defaults(func=launch)
    return parser
