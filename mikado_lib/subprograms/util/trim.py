#!/usr/bin/env python3
# coding: utf-8

"""Script to trim down the terminal exons of multiexonic transcripts"""


import sys
import argparse
from intervaltree import Interval
from ...loci_objects import Transcript
from .. import to_gff
from ...exceptions import InvalidTranscript
from ...utilities.log_utils import create_default_logger

__author__ = 'Luca Venturini'


def trim_noncoding(transcript, max_length=0):
    """
    Function to trim the terminal exons for non-coding transcripts,
     i.e. the simplest case.
    :param transcript: the transcript to be trimmed.
    :type transcript: Transcript

    :param max_length: maximum length of terminal exons.
    :type max_length: int

    :return: transcript
    :rtype: Transcript
    """

    # Non-coding transcript: trim the terminal exons and the return

    first = transcript.exons[0]
    if (first[1] - first[0] + 1) > max_length:
        # newfirst = list(first)
        # newfirst[0] = first[1] - max_length
        newfirst = Interval(first.end - max_length, first.end)
        transcript.start = newfirst[0]
        transcript.exons[0] = newfirst
    last = transcript.exons[-1]
    if (last[1] - last[0] + 1) > max_length:
        newlast = list(last[:])
        newlast[1] = last[0] + max_length
        transcript.end = newlast[1]
        transcript.exons[-1] = Interval(*newlast)
        # if transcript.selected_cds_length > 0:
        #     last_utr = [segment for segment in transcript.combined_utr if
        #                 segment[1] == last[1]][0]
        #     transcript.combined_utr.remove(last_utr)
        #     transcript.combined_utr.append(tuple(newlast))
    return transcript


def trim_start(transcript, cds_start, max_length=0):

    """
    Method to trim the beginning of a transcript, given both
    the cds_start and the maximal allowed length after trimming.
    :param transcript: the transcript to be trimmed
    :type transcript: Transcript

    :param cds_start: the position at which the CDS starts for this transcript.
    :type cds_start: int

    :param max_length: maximum length of terminal exons.
    :type max_length: int

    :rtype: Transcript
    """

    first = transcript.exons[0]
    if first[1] < cds_start:  # first exon is only UTR
        if first[1] - first[0] + 1 > max_length:
            # First exon longer than max_length; trim it,
            # remove from UTR and exons, substitute with trimmed exon
            newfirst = Interval(first[1] - max_length, first[1])
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
            newfirst = Interval(min(cds_start, first[1] - max_length), first[1])
            # Retrieve and remove offending UTR segment
            utr_segment = [segment for segment in transcript.combined_utr
                           if segment[0] == first[0]][0]
            transcript.combined_utr.remove(utr_segment)
            transcript.exons.remove(first)
            transcript.exons.append(newfirst)
            if newfirst[0] < cds_start:
                # Create new UTR segment
                newu = Interval(newfirst[0], cds_start - 1)
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
    :param transcript: the transcript to analyse
    :type transcript: Transcript

    :param cds_end: position of CDS end
    :type cds_end: int

    :param max_length: maximum length of terminal exons.
    :type max_length: int

    :return: transcript
    :rtype: mikado_lib.loci_objects.transcript.Transcript
    """

    # We have to sort b/c we appended a new exon in trim_start
    last = sorted(transcript.exons)[-1]
    if last.begin > cds_end:
        if last[1] - last[0] + 1 > max_length:
            newlast = Interval(last[0], last[0] + max_length)
            transcript.combined_utr.remove(last)
            transcript.exons.remove(last)
            transcript.combined_utr.append(newlast)
            transcript.exons.append(newlast)
        else:
            newlast = last
    elif last[0] <= cds_end < last[1]:
        if last[1] - last[0] + 1 > max_length:
            newlast = Interval(last[0], max(cds_end, last[0] + max_length))
            utr_segment = [segment for segment in transcript.combined_utr if
                           segment[1] == last[1]][0]
            transcript.combined_utr.remove(utr_segment)
            transcript.exons.remove(last)
            transcript.exons.append(newlast)
            if newlast[1] > cds_end:
                newu = Interval(cds_end + 1, newlast[1])
                transcript.combined_utr.append(newu)
        else:
            newlast = last
    else:
        newlast = last

    transcript.end = newlast[1]
    assert all([(exon.end <= newlast.end for exon in transcript.exons)])
    return transcript


def trim_coding(transcript, logger, max_length=0):
    """
    Function to trim the terminal exons for coding transcripts,
     i.e. the more complex case as we have to account for CDS
     as well.
    :param transcript: the transcript to analyse.
    :type transcript: Transcript

    :param logger: a logger instance
    :type logger: logging.Logger

    :param max_length: maximum length of terminal exons.
    :type max_length: int

    :return: transcript
    :rtype: Transcript
    """
    # Coding transcript
    # Order cds_start and end irrespectively of strand
    cds_start, cds_end = sorted([transcript.selected_cds_end,
                                 transcript.selected_cds_start])

    if cds_start >= cds_end:
        logger.warning("{0} has coincident CDS start and end coordinates. Stripping CDS".format(
            transcript.id))
        transcript.strip_cds()
        transcript = trim_noncoding(transcript, max_length=max_length)

    # assert cds_start < cds_end, "\n".join([str(_) for _ in (cds_start, cds_end, transcript)])
    transcript = trim_start(transcript, cds_start,
                            max_length=max_length)
    transcript = trim_end(transcript, cds_end,
                          max_length=max_length)

    return transcript


def strip_terminal(transcript, args) -> Transcript:
    """This function will take as input a transcript and then:
    - return immediately if the transcript is monoexonic
    - trim the terminal exons to a length of max(args.max_length (if longer), CDS boundary)

    :param transcript: the input transcript
    :type transcript: Transcript

    :param args: argparse configuration namespace
    :type args: argparse.Namespace

    :rtype : Transcript
    """

    transcript.finalize()
    # Return immediately if the transcript is monoexonic
    if transcript.monoexonic is True:
        return transcript

    transcript.finalized = False
    transcript.segments = []

    # noinspection PyUnresolvedReferences
    max_l = args.max_length
    try:
        if transcript.selected_cds_length == 0:
            transcript = trim_noncoding(transcript, max_length=max_l)
        else:
            transcript = trim_coding(transcript, args.logger, max_length=max_l)
    except InvalidTranscript as exc:
        args.logger.exception(exc)
        return None

    # noinspection PyUnresolvedReferences
    transcript.finalize()
    assert transcript.start == transcript.exons[0][0], (transcript.start, transcript.exons)
    assert transcript.end == transcript.exons[-1][1], (transcript.end, transcript.exons)
    return transcript


def launch(args):

    """
    Launcher of the program.
    :param args: the argparse Namespace.
    """

    args.logger = create_default_logger("trimmer")
    args.logger.setLevel("WARN")

    if args.as_gtf is False:
        print("##gff-version 3",
              file=args.out)

    transcript = None

    if args.as_gtf is True:
        format_name = "gtf"
    else:
        format_name = "gff3"

    if format_name == "gff3":
        print("WARNING: proper GFF3 format not implemented yet!",
              file=sys.stderr)

    for record in args.ann:
        if record.is_transcript is True:
            if transcript is not None:
                trimmed = strip_terminal(transcript, args)
                if trimmed is not None:
                    print(trimmed.format(format_name), file=args.out)
            transcript = Transcript(record)
        elif record.is_exon is True:
            transcript.add_exon(record)

    if transcript is not None:
        trimmed = strip_terminal(transcript, args)
        print(trimmed.format(format_name), file=args.out)


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
