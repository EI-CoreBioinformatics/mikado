#!/usr/bin/env python3
# coding: utf-8

"""
Subprogram that constitutes the first step of the Mikado pipeline.
"""

import sys
import os
import argparse
import operator
import collections
import logging
import itertools
from mikado_lib import exceptions
import copy
from mikado_lib.loci_objects.transcriptchecker import TranscriptChecker
# from mikado_lib.parsers import GTF, GFF
from mikado_lib.subprograms import to_gff
from Bio import SeqIO
import functools
import multiprocessing
import multiprocessing.connection
import multiprocessing.sharedctypes

__author__ = 'Luca Venturini'


def grouper(iterable, num, fillvalue=(None, None)):
    """Collect data into fixed-length chunks or blocks.
    Source: itertools standard library documentation
    https://docs.python.org/3/library/itertools.html?highlight=itertools#itertools-recipes
    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * num
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def create_transcript(lines,
                      fasta_seq,
                      lenient=False,
                      strand_specific=False):
    """Function to create the checker.

    :param lines: all the exon lines for an object
    :type lines: list[GTF.GtfLine]

    :param fasta_seq: str

    :type lenient: bool
    :type strand_specific: bool

    """
    logger = logging.getLogger("main")
    logger.debug("Starting with %s", lines[0].transcript)

    transcript_line = copy.deepcopy(lines[0])
    transcript_line.feature = "transcript"
    transcript_line.start = min(r.start for r in lines)
    transcript_line.end = max(r.end for r in lines)
    transcript_object = TranscriptChecker(transcript_line,
                                          fasta_seq, lenient=lenient,
                                          strand_specific=strand_specific)
    transcript_object.logger = logger
    for line in lines:
        transcript_object.add_exon(line)
    try:
        transcript_object.finalize()
        transcript_object.check_strand()
    except exceptions.IncorrectStrandError:
        logger.info("Discarded %s because of incorrect fusions of splice junctions",
                    lines[0].transcript)
        transcript_object = None
    except exceptions.InvalidTranscript:
        logger.info("Discarded generically invalid transcript %s",
                    lines[0].transcript)
        transcript_object = None
    logger.debug("Finished with %s", lines[0].transcript)
    return transcript_object


def store_transcripts(exon_lines, fasta, logger, min_length=0):

    """
    Function that analyses the exon lines from the original file
    and organises the data into a proper dictionary.
    :param exon_lines: dictionary of exon lines, ordered by TID
    :type exon_lines: dict
    :return: transcripts: dictionary which will be the final output
    :rtype: transcripts
    """

    logger.info("Starting to organise %d transcripts", len(exon_lines))
    transcripts = collections.defaultdict(dict)
    for tid in exon_lines:
        tlines = exon_lines[tid]
        tlength = sum(exon.end + 1 - exon.start for exon in exon_lines[tid])
        # Discard transcript under a certain size
        if tlength < min_length:
            logger.warn("Discarding %s because its size (%d) is under the minimum of %d",
                        tid, tlength, min_length)
            continue
        start, end = min(x.start for x in tlines), max(x.end for x in tlines)
        chrom = tlines[0].chrom
        if (start, end) not in transcripts[chrom]:
            transcripts[chrom][(start, end)] = []
        transcripts[chrom][(start, end)].append(tid)

    logger.info("Starting to sort %d transcripts", len(exon_lines))
    keys = []
    for chrom in sorted(transcripts.keys()):
        for key in sorted(transcripts[chrom].keys(),
                          key=operator.itemgetter(0, 1)):
            seq = fasta[chrom][key[0]-1:key[1]]
            keys.extend([tid, seq] for tid in transcripts[chrom][key])
    logger.info("Finished to sort %d transcripts", len(exon_lines))

    return keys


def perform_check(keys, exon_lines, args, logger):

    """
    This is the most important method. After preparing the data structure,
    this function creates the real transcript instances and checks that
    they are correct when looking at the underlying genome sequence.
    This is also the point at which we start using multithreading, if
    so requested.
    :param keys: sorted list of [tid, sequence]
    :param exon_lines: dictionary hosting the exon lines, indexed by TID
    :param args: the namespace
    :param logger: logger
    :return:
    """

    counter = 0

    # pylint: disable=no-member
    pool = multiprocessing.Pool(args.threads)
    # pylint: enable=no-member

    # Use functools to pre-configure the function
    # with all necessary arguments aside for the lines
    partial_checker = functools.partial(create_transcript,
                                        lenient=args.lenient,
                                        strand_specific=args.strand_specific)

    logger.info("Starting to analyse the transcripts looking at the underlying sequence")
    if args.single is True:
        for tid, seq in keys:
            transcript_object = partial_checker(exon_lines[tid], seq)
            if transcript_object is None:
                continue
            counter += 1
            if counter >= 10**4 and counter % (10**4) == 0:
                logger.info("Retrieved %d transcript positions", counter)
            elif counter >= 10**3 and counter % (10**3) == 0:
                logger.debug("Retrieved %d transcript positions", counter)
            print(transcript_object.__str__(to_gtf=True), file=args.out)
            print(transcript_object.fasta, file=args.out_fasta)
    else:
        for group in grouper(keys, 100):
            if group is None:
                continue
            results = [
                pool.apply_async(partial_checker, args=(exon_lines[tid], seq))
                for (tid, seq) in iter(x for x in group if x != (None, None))
            ]
            for transcript_object in results:
                transcript_object = transcript_object.get()
                if transcript_object is None:
                    continue
                counter += 1
                if counter >= 10**4 and counter % (10**4) == 0:
                    logger.info("Retrieved %d transcript positions", counter)
                elif counter >= 10**3 and counter % (10**3) == 0:
                    logger.debug("Retrieved %d transcript positions", counter)
                print(transcript_object.__str__(to_gtf=True), file=args.out)
                print(transcript_object.fasta, file=args.out_fasta)
        pool.close()
        pool.join()

    logger.info("Finished to analyse %d transcripts (%d retained)",
                len(exon_lines), counter)
    return


def prepare(args):
    """Main script function."""

    if args.labels != '':
        args.labels = args.labels.split(",")
        # Checks labels are unique
        assert len(set(args.labels)) == len(args.labels)
        assert not any([True for _ in args.labels if _.strip() == ''])
        if len(args.labels) != len(args.gff):
            raise ValueError("Incorrect number of labels specified")
    else:
        args.labels = [""] * len(args.gff)

    if isinstance(args.out_fasta, str):
        args.out_fasta = open(args.out_fasta, 'w')
    if isinstance(args.out, str):
        args.out = open(args.out, 'w')

    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "{asctime}:{levelname} - {filename}:{lineno} - {funcName} - {message}",
        style="{")
    handler.setFormatter(formatter)
    logger = logging.getLogger("main")
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
    elif args.quiet is True:
        logger.setLevel(logging.WARN)

    to_seqio = functools.partial(to_seqio_complete,
                                 logger_instance=logger)
    args.fasta.close()
    args.fasta = to_seqio(args.fasta.name)

    exon_lines = collections.defaultdict(list)

    previous_file_ids = collections.defaultdict(set)

    for label, gff_handle in zip(args.labels, args.gff):
        found_ids = set.union(set(), *previous_file_ids.values())
        for row in gff_handle:
            if row.is_exon is False:
                continue
            # This convoluted block is due to the fact that transcript is a list for
            # GFF files (as it is a parent), while it is a string for GTF files

            if isinstance(row.transcript, list):
                if label != '':
                    row.transcript = ["{0}_{1}".format(label, tid) for tid in row.transcript]
                for tid in row.transcript:
                    if tid in found_ids:
                        assert label != ''
                        raise ValueError("""{0} has already been found in another file,
                        this will cause unsolvable collisions. Please rerun preparation using
                        labels to tag each file.""")
                    previous_file_ids[gff_handle.name].add(tid)
                    exon_lines[tid].append(row)
            else:
                if label != '':
                    row.transcript = "{0}_{1}".format(label, row.transcript)
                if row.transcript in found_ids:
                    assert label != ''
                    raise ValueError("""{0} has already been found in another file,
                        this will cause unsolvable collisions. Please rerun preparation using
                        labels to tag each file.""")
                exon_lines[row.transcript].append(row)

    logger.info("Finished loading exon lines")

    # Prepare the sorted data structure
    keys = store_transcripts(exon_lines,
                             args.fasta,
                             logger,
                             min_length=args.minimum_length)

    perform_check(keys, exon_lines, args, logger)

    logger.info("Finished")


def to_seqio_complete(string, logger_instance=None):
    """
    Function to index a FASTA file using SeqIO.
    :param string: a vaild file name
    :type string: str

    :param logger_instance: a logging.Logger instance

    # :param manager_instance: a multiprocessing.manager_instance instance
    # :type manager_instance: None | multiprocessing.Manager
    """

    logger_instance.info("Loading reference file")
    if (not os.path.exists(string) or
            not os.path.isfile(string) or
            not os.stat(string).st_size > 0):
        exc = ValueError("Invalid input file.")
        logger_instance.exception(exc)
        raise exc
    seqdict = SeqIO.to_dict(SeqIO.parse(open(string), 'fasta'))

    seqdict = dict((seq, str(seqdict[seq].seq)) for seq in seqdict)
    logger_instance.info("Finished loading reference file")
    return seqdict


def prepare_parser():
    """
    This function defines the parser for the command line interface
    of the program.
    :return: an argparse.Namespace object
    :rtype: argparse.Namespace
    """

    def to_cpu_count(string):
        """
        :param string: cpu requested
        :rtype: int
        """
        try:
            string = int(string)
        except:
            raise
        return max(1, string)

    def positive(string):
        return abs(int(string))

    parser = argparse.ArgumentParser("""Script to prepare a GTF for the pipeline;
    it will perform the following operations:
    1- add the "transcript" feature
    2- sort by coordinates
    3- check the strand""")
    parser.add_argument("--fasta", type=argparse.FileType(), required=True,
                        help="Genome FASTA file. Required.")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="store_true", default=False)
    verbosity.add_argument("-q", "--quiet", action="store_true", default=False)
    parser.add_argument("-s", "--strand-specific", dest="strand_specific",
                        action="store_true", default=False,
                        help="""Flag. If set, monoexonic transcripts
                        will be left on their strand rather than being
                        moved to the unknown strand.""")
    parser.add_argument("-l", "--lenient", action="store_true", default=False,
                        help="""Flag. If set, transcripts with mixed +/-
                        splices will not cause exceptions but rather
                        be annotated as problematic.""")
    parser.add_argument("-m", "--minimum_length", default=0, type=positive,
                        help="Minimum length for transcripts. Default: no filtering.")
    parser.add_argument("-t", "--threads",
                        help="Number of processors to use (default %(default)s)",
                        type=to_cpu_count, default=1)
    parser.add_argument("--labels", type=str, default="",
                        help="""Labels to attach to the IDs of the transcripts of the input files,
                        separated by comma.""")
    parser.add_argument("--single", action="store_true", default=False,
                        help="Disable multi-threading. Useful for debugging.")
    parser.add_argument("-o", "--out", default='mikado.gtf', type=argparse.FileType('w'),
                        help="Output file. Default: mikado.fasta.")
    parser.add_argument("-of", "--out_fasta", default="mikado.fasta",
                        type=argparse.FileType('w'),
                        help="Output file. Default: mikado.fasta.")
    parser.add_argument("gff", type=to_gff, help="Input GFF/GTF file(s).", nargs="+")
    parser.set_defaults(func=prepare)
    return parser
