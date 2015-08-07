#!/usr/bin/env python3
# coding: utf-8

"""Script to prepare a GTF for the pipeline; it will perform the following operations:
    1- add the "transcript" feature
    2- sort by coordinates
    3- check the strand
"""

import sys
import os
import argparse
import operator
import collections
import logging
# import time
from mikado_lib import exceptions
import copy
from mikado_lib.loci_objects.transcriptchecker import TranscriptChecker
from mikado_lib.parsers import GTF
from Bio import SeqIO
import functools
import multiprocessing
import multiprocessing.connection
import multiprocessing.sharedctypes


def create_transcript(lines,
                      fasta_seq,
                      lenient=False,
                      strand_specific=False,
                      ):
    """Function to create the checker.

    :param lines: all the exon lines for an object
    :type lines: list[GTF.GtfLine]

    :param fasta_seq: str

    :type lenient: bool
    :type strand_specific: bool

    :param sender: Pipe sender
    :type sender: multiprocessing.connection.Connection | None

    """
    logger = logging.getLogger("main")
    logger.debug("Starting with {0}".format(lines[0].transcript))

    transcript_line = copy.deepcopy(lines[0])
    transcript_line.feature = "transcript"
    transcript_line.start = min(r.start for r in lines)
    transcript_line.end = max(r.end for r in lines)
    transcript_object = TranscriptChecker(transcript_line, fasta_seq, lenient=lenient,
                                          strand_specific=strand_specific)
    transcript_object.logger = logger
    for line in lines:
        transcript_object.add_exon(line)
    try:
        transcript_object.finalize()
        transcript_object.check_strand()
    except exceptions.IncorrectStrandError:
        logger.info("Discarded {0} because of incorrect fusions of splice junctions".format(lines[0].transcript))
        transcript_object = None
    except exceptions.InvalidTranscript:
        logger.info("Discarded generically invalid transcript {0}".format(lines[0].transcript))
        transcript_object = None
    logger.debug("Finished with {0}".format(lines[0].transcript))
    return transcript_object


def main():
    """Main script function."""

    def to_gff(string):
        """
        Function to verify the input file.
        :param string: filename
        """
        if string[-4:] != ".gtf":
            raise TypeError("This script takes as input only GTF files.")

        gff = GTF.GTF(string)
        for record in gff:
            if record.header is False:
                gff.close()
                return GTF.GTF(string)
        raise ValueError("Empty GTF file provided.")

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
        if not os.path.exists(string) or not os.path.isfile(string) or not os.stat(string).st_size > 0:
            exc = ValueError("Invalid input file.")
            logger_instance.exception(exc)
            raise exc
        seqdict = SeqIO.to_dict(SeqIO.parse(open(string), 'fasta'))

        seqdict = dict((seq, str(seqdict[seq].seq)) for seq in seqdict)
        # if manager_instance is not None:
        #     logger_instance.debug("Loading index into shared memory object")
        #     # seqdict = dict((x, multiprocessing.Array("u", seqdict[x], lock=False)) for x in seqdict)
        #     seqdict = manager_instance.dict(seqdict, lock=False)
        logger_instance.info("Finished loading reference file")
        return seqdict

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

    handler = logging.StreamHandler()
    formatter = logging.Formatter("{asctime}:{levelname} - {filename}:{lineno} - {funcName} - {message}",
                                  style="{")
    handler.setFormatter(formatter)
    logger = logging.getLogger("main")
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

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
    parser.add_argument("-s", "--strand-specific", dest="strand_specific", action="store_true", default=False,
                        help="""Flag. If set, monoexonic transcripts will be left on their strand
                        rather than being moved to the unknown strand.""")
    parser.add_argument("-l", "--lenient", action="store_true", default=False,
                        help="""Flag. If set, transcripts with mixed +/- splices will not cause exceptions
                        but rather be annotated as problematic.""")
    parser.add_argument("-t", "--threads", help="Number of processors to use (default %(default)s)",
                        type=to_cpu_count, default=1)
    parser.add_argument("--single", action="store_true", default=False,
                        help="Disable multi-threading. Useful for debugging.")
    parser.add_argument("gff", type=to_gff, help="Input GTF file.")
    parser.add_argument("out", default=sys.stdout, nargs='?', type=argparse.FileType('w'),
                        help="Output file. Default: STDOUT.")
    args = parser.parse_args()

    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
    elif args.quiet is True:
        logger.setLevel(logging.WARN)

    if args.single is False:
        manager = multiprocessing.Manager()
    else:
        manager = None
    to_seqio = functools.partial(to_seqio_complete,
                                 logger_instance=logger)
    args.fasta.close()
    args.fasta = to_seqio(args.fasta.name)

    exon_lines = dict()

    for row in args.gff:
        if row.feature != "exon":
            continue
        if row.transcript not in exon_lines:
            exon_lines[row.transcript] = []
        exon_lines[row.transcript].append(row)

    logger.info("Finished loading exon lines")
    transcripts = collections.defaultdict(dict)

    pool = multiprocessing.Pool(args.threads)

    # Use functools to pre-configure the function with all necessary arguments
    # aside for the lines
    partial_checker = functools.partial(create_transcript,
                                        lenient=args.lenient,
                                        strand_specific=args.strand_specific)

    logger.info("Starting to sort {0} transcripts".format(len(exon_lines)))

    for tid in exon_lines:
        tlines = exon_lines[tid]
        start, end = min(x.start for x in tlines), max(x.end for x in tlines)
        chrom = tlines[0].chrom
        if (start, end) not in transcripts[chrom]:
            transcripts[chrom][(start, end)] = []
        transcripts[chrom][(start, end)].append(tid)

    counter = 0
    keys = []
    for chrom in sorted(transcripts.keys()):
        for key in sorted(transcripts[chrom].keys(), key=operator.itemgetter(0, 1)):
            seq = args.fasta[chrom][key[0]:key[1]]
            if args.single is False:
                for tid in transcripts[chrom][key]:
                    keys.append(pool.apply_async( partial_checker, args=(exon_lines[tid], seq)))
            else:
                keys.extend([tid, seq] for tid in transcripts[chrom][key])

    pool.close()

    if args.single is True:
        for tid,seq in keys:
            transcript_object = partial_checker(exon_lines[tid], seq)
            if transcript_object is None:
                continue
            counter += 1
            if counter >= 10**4 and counter % (10**4) == 0:
                logger.info("Retrieved {0} transcript positions".format(counter))
            elif counter >= 10**3 and counter % (10**3) == 0:
                logger.debug("Retrieved {0} transcript positions".format(counter))
            print(transcript_object.__str__(to_gtf=True), file=args.out)

    else:
        for key in keys:
            transcript_object = key.get()
            if transcript_object is None:
                continue
            counter += 1
            if counter >= 10**4 and counter % (10**4) == 0:
                logger.info("Retrieved {0} transcript positions".format(counter))
            elif counter >= 10**3 and counter % (10**3) == 0:
                logger.debug("Retrieved {0} transcript positions".format(counter))
            print(transcript_object.__str__(to_gtf=True), file=args.out)

    logger.info("Finished to analyse {0} transcripts ({1} retained)".format(
        len(exon_lines),
        counter
        ))

if __name__ == "__main__":
    main()
