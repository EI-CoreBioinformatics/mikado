#!/usr/bin/env python3
# coding: utf-8

"""Script to trim down the terminal exons of multiexonic transcripts"""

import argparse
import sys
from ...utilities.log_utils import create_default_logger
from ...exceptions import InvalidTranscript


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

    transcript.unfinalize()
    exons = sorted(transcript.exons)
    first, last = exons[0], exons[-1]
    if (first[1] - first[0] + 1) > max_length:
        # newfirst = list(first)
        # newfirst[0] = first[1] - max_length
        newfirst = tuple([first[1] - max_length, first[1]])
        transcript.start = newfirst[0]
        transcript.remove_exon(first)
        transcript.add_exon(newfirst)
        # transcript.segments.remove(("exon", first))
        # transcript.segments.append(("exon", newfirst))

    # last = transcript.exons[-1]
    if (last[1] - last[0] + 1) > max_length:
        newlast = tuple([last[0], last[0] + max_length])
        transcript.end = last[0] + max_length
        transcript.remove_exon(last)
        transcript.add_exon(newlast)
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
            newfirst = tuple([first[1] - max_length, first[1]])
            transcript.remove_exon(first)
            transcript.add_exon(newfirst)
        else:
            # Leave as it is
            newfirst = first
    elif first[0] < cds_start <= first[1]:
        # Partly UTR partly CDS
        if first[1] - first[0] + 1 > max_length:
            newfirst = tuple([min(cds_start, first[1] - max_length), first[1]])
            # Retrieve and remove offending UTR segment
            transcript.remove_exon(first)
            transcript.add_exon(newfirst)
            transcript.add_exon((cds_start, newfirst[1]), feature="CDS")
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
    :rtype: Mikado.loci_objects.transcript.Transcript
    """

    # We have to sort b/c we appended a new exon in trim_start
    last = sorted(transcript.exons)[-1]
    if last[0] > cds_end:  # Only UTR
        if last[1] - last[0] + 1 > max_length:

            newlast = tuple([last[0], last[0] + max_length])
            if last != newlast:
                transcript.remove_exon(last)
                transcript.add_exon(newlast)
        else:
            newlast = last
    elif last[0] <= cds_end < last[1]:
        if last[1] - last[0] + 1 > max_length:
            newlast = tuple([last[0], max(cds_end, last[0] + max_length)])
            transcript.remove_exon(last)
            transcript.add_exon(newlast)
            transcript.add_exon((newlast[0], cds_end), feature="CDS")
        else:
            newlast = last
    else:
        newlast = last

    transcript.end = newlast[1]
    assert all([(exon[1] <= newlast[1] for exon in transcript.exons)])
    assert all([(segment[1][1] <= newlast[1] for segment in transcript.segments)])
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

    # Wrong CDS

    if transcript.selected_cds_start is None:
        logger.debug("Non-coding transcript %s submitted, calling the correct function.",
                     transcript.id)
        return trim_noncoding(transcript, max_length=max_length)

    cds_start, cds_end = sorted([transcript.selected_cds_end,
                                 transcript.selected_cds_start])

    transcript.unfinalize()
    assert transcript.internal_orfs == []
    if cds_start == cds_end:
        logger.warning("{0} has coincident CDS start and end coordinates. Stripping CDS".format(
            transcript.id))
        transcript.strip_cds()
        transcript.feature = "transcript"
        # Now we have to re-deinitialise the transcript
        transcript.finalized = False
        transcript = trim_noncoding(transcript,
                                    max_length=max_length)
        return transcript

    assert cds_start < cds_end, "\n".join([str(_) for _ in (cds_start, cds_end, transcript)])
    transcript = trim_start(transcript, cds_start,
                            max_length=max_length)
    transcript = trim_end(transcript, cds_end,
                          max_length=max_length)

    return transcript


def strip_terminal(transcript, args):
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

    if transcript.selected_cds_start is not None:
        coding = True
    else:
        coding = False

    # noinspection PyUnresolvedReferences
    max_l = args.max_length
    try:
        if coding:
            assert sorted([transcript.selected_cds_start,
                           transcript.selected_cds_end])
            transcript = trim_coding(transcript, args.logger, max_length=max_l)
        else:
            transcript = trim_noncoding(transcript, max_length=max_l)
    except InvalidTranscript as exc:
        args.logger.exception(exc)
        return None

    # noinspection PyUnresolvedReferences
    transcript.finalize()
    assert transcript.start == transcript.exons[0][0], (transcript.start, transcript.exons)
    assert transcript.end == transcript.exons[-1][1], (transcript.end, transcript.exons)
    return transcript

def trim_gene(current_gene, format_name, args):
    """"""

    if current_gene is not None:
        current_gene.finalize()
        new_start, new_end = float("inf"), float("-inf")
        for tid in list(current_gene.transcripts)[:]:
            trimmed = strip_terminal(current_gene[tid], args)
            current_gene.remove(tid)  # We have to remove the old transcript anyway
            if trimmed is not None:
                current_gene.add(trimmed)
                new_start, new_end = min(new_start,
                                         trimmed.start), max(new_end,
                                                             trimmed.end)
        # Reset them manually to avoid warnings
        current_gene.start, current_gene.end = new_start, new_end
        current_gene.finalize()
        if len(current_gene.transcripts) > 0:
            print(current_gene.format(format_name), file=args.out)
    return


def launch(args):

    """
    Launcher of the program.
    :param args: the argparse Namespace.
    """

    from ...parsers import to_gff
    from ...loci import Transcript
    from ...transcripts.reference_gene import Gene

    args.gff = to_gff(args.gff)
    args.logger = create_default_logger("trimmer")
    args.logger.setLevel("WARN")

    if args.as_gtf is False:
        print("##gff-version 3",
              file=args.out)

    if args.as_gtf is True:
        format_name = "gtf"
    else:
        format_name = "gff3"

    if args.ann.file_format == "gff3":
        current_gene = None
        for record in args.ann:
            if record.is_gene is True:
                trim_gene(current_gene, format_name, args)
                current_gene = Gene(record, logger=args.logger)
            elif record.is_transcript is True:
                assert current_gene is not None
                assert record.parent[0] == current_gene.id
                current_gene.add(Transcript(record))
            elif record.is_exon is True:
                current_gene.add_exon(record)
            elif record.is_derived is True:
                for derived in record.derived_from:
                    current_gene[derived].add_derived_child(record.id)

        trim_gene(current_gene, format_name, args)
    else:
        # When analysing a GTF, we can blissfully ignore everything
        current_gene = None
        for record in args.ann:
            if current_gene is not None and current_gene.id == record.gene:
                if record.is_exon is True:
                    if record.transcript in current_gene:
                        current_gene.add_exon(record)
                    else:
                        current_gene.add(Transcript(record))
                elif record.is_transcript is True:
                    assert record.transcript not in current_gene
                    current_gene.add(Transcript(record))
            else:
                trim_gene(current_gene, format_name, args)
                current_gene = Gene(Transcript(record), logger=args.logger)
            continue
        trim_gene(current_gene, format_name, args)


def trim_parser():

    """
    Command line parser for the utility.
    """

    parser = argparse.ArgumentParser(
        "Script to trim down the terminal exons of multiexonic transcripts")
    parser.add_argument("-ml", "--max_length", type=int, default=50,
                        help="Maximal length of trimmed terminal exons")
    parser.add_argument("--as-gtf", default=False, action="store_true",
                        dest="as_gtf",
                        help="Flag. If set, the output will be in GTF rather than GFF3 format.")
    parser.add_argument("ann", type=str, help="Reference GTF/GFF output file.")
    parser.add_argument("out", nargs="?", default=sys.stdout,
                        type=argparse.FileType('w'))
    parser.set_defaults(func=launch)
    return parser
