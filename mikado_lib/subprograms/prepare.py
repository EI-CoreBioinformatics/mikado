#!/usr/bin/env python3
# coding: utf-8

"""
Subprogram that constitutes the first step of the mikado_lib pipeline.
"""

import sys
import io
import os
import argparse
import operator
import collections
import logging
import itertools
from .. import exceptions
import copy
from ..loci_objects.transcriptchecker import TranscriptChecker
from ..configuration.configurator import to_json
from . import to_gff
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

    :param iterable: the iterable to be considered for grouping.
    :param num: length of the chunks.
    :type num: int

    :param fillvalue: the default filler for missing positions while grouping.
    :type fillvalue: tuple
    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * num
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def create_transcript(lines,
                      fasta_seq,
                      lenient=False,
                      strand_specific=False,
                      canonical_splices=(("GT", "AG"),
                                         ("GC", "AG"),
                                         ("AT", "AC"))):
    """Function to create the checker.

    :param lines: all the exon lines for an object
    :type lines: list[GTF.GtfLine]

    :param fasta_seq: str

    :type lenient: bool
    :type strand_specific: bool

    :param canonical_splices: the splices considered as canonical for the species.
    :type canonical_splices: list[tuple]

    """
    logger = logging.getLogger("main")
    logger.debug("Starting with %s", lines[0].transcript)

    transcript_line = copy.deepcopy(lines[0])
    transcript_line.feature = "transcript"
    transcript_line.start = min(r.start for r in lines)
    transcript_line.end = max(r.end for r in lines)
    transcript_object = TranscriptChecker(transcript_line,
                                          fasta_seq, lenient=lenient,
                                          strand_specific=strand_specific,
                                          canonical_splices=canonical_splices)

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


def store_transcripts(exon_lines, fasta, logger, min_length=0, cache=False):

    """
    Function that analyses the exon lines from the original file
    and organises the data into a proper dictionary.
    :param exon_lines: dictionary of exon lines, ordered by TID
    :type exon_lines: dict

    :param fasta: the FASTA sequence to use for the transcript.

    :param logger: logger instance.
    :type logger: logging.Logger

    :param min_length: minimal length of the transcript.
    If it is not met, the transcript will be discarded.
    :type min_length: int

    :param cache: whether the FASTA file of the genome has been preloaded into RAM.
    :type cache: bool

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
            logger.warning("Discarding %s because its size (%d) is under the minimum of %d",
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
        if cache is False:
            logger.debug("Copying %s into memory", chrom)
            chrom_seq = str(fasta[chrom].seq)
        else:
            chrom_seq = fasta[chrom]

        logger.debug("Starting with %s", chrom)
        for key in sorted(transcripts[chrom].keys(),
                          key=operator.itemgetter(0, 1)):
            seq = chrom_seq[key[0]-1:key[1]]
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
    partial_checker = functools.partial(
        create_transcript,
        lenient=args.json_conf["prepare"]["lenient"],
        strand_specific=args.json_conf["prepare"]["strand_specific"])

    logger.info("Starting to analyse the transcripts looking at the underlying sequence")
    if args.json_conf["prepare"]["single"] is True:
        for tid, seq in keys:
            transcript_object = partial_checker(exon_lines[tid], seq)
            if transcript_object is None:
                continue
            counter += 1
            if counter >= 10**4 and counter % (10**4) == 0:
                logger.info("Retrieved %d transcript positions", counter)
            elif counter >= 10**3 and counter % (10**3) == 0:
                logger.debug("Retrieved %d transcript positions", counter)
            print(transcript_object.__str__(to_gtf=True),
                  file=args.json_conf["prepare"]["out"])
            print(transcript_object.fasta,
                  file=args.json_conf["prepare"]["out_fasta"])
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
                print(transcript_object.__str__(to_gtf=True),
                      file=args.json_conf["prepare"]["out"])
                print(transcript_object.fasta,
                      file=args.json_conf["prepare"]["out_fasta"])
        pool.close()
        pool.join()

    logger.info("Finished to analyse %d transcripts (%d retained)",
                len(exon_lines), counter)
    return


def setup(args):
    """Method to set up the analysis using the JSON configuration
    and the command line options.

    :param args: the ArgumentParser-derived namespace.
    """

    if args.json_conf["prepare"]["log"]:
        handler = logging.FileHandler(
            args.json_conf["prepare"]["log"],
            "w")
    else:
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

    if args.gff:
        args.json_conf["prepare"]["gff"] = args.gff
    else:
        assert args.json_conf["prepare"]["gff"]

    if args.labels != '':
        args.labels = args.labels.split(",")
        # Checks labels are unique
        assert len(set(args.labels)) == len(args.labels)
        assert not any([True for _ in args.labels if _.strip() == ''])
        if len(args.labels) != len(args.json_conf["prepare"]["gff"]):
            raise ValueError("Incorrect number of labels specified")
        args.json_conf["prepare"]["labels"] = args.labels
    else:
        if not args.json_conf["prepare"]["labels"]:
            args.labels = [""] * len(args.json_conf["prepare"]["gff"])
            args.json_conf["prepare"]["labels"] = args.labels

    for option in ["cache",
                   "out", "out_fasta", "fasta",
                   "minimum_length", "threads", "single"]:
        if ((getattr(args, option) or getattr(args, option) == 0) and
                getattr(args, option) is not False):
            args.json_conf["prepare"][option] = getattr(args, option)

    return args, logger


def load_from_gff(exon_lines, gff_handle, label, found_ids, strip_cds=False):
    """
    Method to load the exon lines from GFF3 files.
    :param exon_lines: the defaultdict which stores the exon lines.
    :param gff_handle: The handle for the GTF to be parsed.
    :param label: label to be attached to all transcripts.
    :type label: str
    :param found_ids: set of IDs already found in other files.
    :type found_ids: set
    :param strip_cds: boolean flag. If true, all CDS lines will be ignored.
    :type strip_cds: bool
    :return:
    """
    transcript2genes = dict()
    new_ids = set()
    for row in gff_handle:
        if row.is_transcript is True:
            if label != '':
                row.id = "{0}_{1}".format(label, row.id)
            transcript2genes[row.id] = row.parent[0]
            continue
        elif not row.is_exon:
            continue
        elif row.is_exon is True:
            if not row.is_cds or (row.is_cds is True and strip_cds is False):
                if label != '':
                    row.transcript = ["{0}_{1}".format(label, tid) for tid in row.transcript]
                parents = row.transcript[:]
                for tid in parents:
                    if tid in found_ids:
                        assert label == ''
                        raise ValueError("""{0} has already been found in another file,
                        this will cause unsolvable collisions. Please rerun preparation using
                        labels to tag each file.""")
                    new_ids.add(tid)
                    new_row = copy.deepcopy(row)
                    new_row.parent = new_row.id = tid
                    new_row.attributes["gene_id"] = transcript2genes[tid]
                    new_row.name = tid
                    exon_lines[tid].append(new_row)
            else:
                continue
    return exon_lines, new_ids


def load_from_gtf(exon_lines, gff_handle, label, found_ids, strip_cds=False):
    """
    Method to load the exon lines from GTF files.
    :param exon_lines: the defaultdict which stores the exon lines.
    :param gff_handle: The handle for the GTF to be parsed.
    :param label: label to be attached to all transcripts.
    :type label: str
    :param found_ids: set of IDs already found in other files.
    :type found_ids: set
    :param strip_cds: boolean flag. If true, all CDS lines will be ignored.
    :type strip_cds: bool
    :return:
    """

    new_ids = set()
    for row in gff_handle:
        if row.is_exon is False or (row.is_cds is True and strip_cds is True):
            continue
        if label != '':
            row.transcript = "{0}_{1}".format(label, row.transcript)
        if row.transcript in found_ids:
            assert label != ''
            raise ValueError("""{0} has already been found in another file,
                this will cause unsolvable collisions. Please rerun preparation using
                labels to tag each file.""")
        exon_lines[row.transcript].append(row)
        new_ids.add(row.transcript)
    return exon_lines, new_ids


def load_exon_lines(args, logger):

    """This function loads all exon lines from the GFF inputs into a
     defaultdict instance.
    :param args: the Namespace from the command line.
    :param logger: the logger instance.
    :type logger: logging.Logger
    :return: exon_lines
    :rtype: collections.defaultdict[list]
    """

    exon_lines = collections.defaultdict(list)
    previous_file_ids = collections.defaultdict(set)
    for label, gff_name in zip(args.json_conf["prepare"]["labels"],
                               args.json_conf["prepare"]["gff"]):
        logger.info("Starting with %s", gff_name)
        gff_handle = to_gff(gff_name)
        found_ids = set.union(set(), *previous_file_ids.values())
        strip_cds = args.json_conf["prepare"]["strip_cds"]

        if gff_handle.__annot_type__ == "gff3":
            exon_lines, new_ids = load_from_gff(exon_lines, gff_handle,
                                                label, found_ids, strip_cds=strip_cds)
        else:
            exon_lines, new_ids = load_from_gtf(exon_lines, gff_handle,
                                                label, found_ids, strip_cds=strip_cds)

        previous_file_ids[gff_handle.name] = new_ids

    return exon_lines


def prepare(args):
    """Main script function.

    :param args: the ArgumentParser-derived namespace.
    """

    args, logger = setup(args)

    if not isinstance(args.json_conf["prepare"]["fasta"], io.TextIOWrapper):
        if not (isinstance(args.json_conf["prepare"]["fasta"], str) and
                os.path.exists(args.json_conf["prepare"]["fasta"])):
            logger.critical("Invalid FASTA file: %s",
                            args.json_conf["prepare"]["fasta"])
            sys.exit(1)
        else:
            pass
    else:
        args.json_conf["prepare"]["fasta"].close()
        args.json_conf["prepare"]["fasta"] = args.json_conf["prepare"]["fasta"].name

    assert len(args.json_conf["prepare"]["gff"]) > 0
    assert len(args.json_conf["prepare"]["gff"]) == len(args.json_conf["prepare"]["labels"])

    args.json_conf["prepare"]["out_fasta"] = open(args.json_conf["prepare"]["out_fasta"], 'w')
    args.json_conf["prepare"]["out"] = open(args.json_conf["prepare"]["out"], 'w')

    to_seqio = functools.partial(to_seqio_complete,
                                 cache=args.json_conf["prepare"]["cache"],
                                 logger_instance=logger)

    args.json_conf["prepare"]["fasta"] = to_seqio(args.json_conf["prepare"]["fasta"])

    logger.info("Started loading exon lines")

    exon_lines = load_exon_lines(args, logger)

    logger.info("Finished loading exon lines")

    # Prepare the sorted data structure
    keys = store_transcripts(exon_lines,
                             args.json_conf["prepare"]["fasta"],
                             logger,
                             min_length=args.json_conf["prepare"]["minimum_length"],
                             cache=args.json_conf["prepare"]["cache"])

    perform_check(keys, exon_lines, args, logger)

    logger.info("Finished")


def to_seqio_complete(string, cache=False, logger_instance=None):
    """
    Function to index a FASTA file using SeqIO.
    :param string: a vaild file name
    :type string: str

    :param logger_instance: a logging.Logger instance

    :param cache: boolean flag. If set to True,the genome will be preloaded in memory.
    Fast but expensive!
    :type cache: bool
    """

    logger_instance.info("Loading reference file")
    if (not os.path.exists(string) or
            not os.path.isfile(string) or
            not os.stat(string).st_size > 0):
        exc = ValueError("Invalid input file.")
        logger_instance.exception(exc)
        raise exc
    if cache is True:
        seqdict = SeqIO.to_dict(SeqIO.parse(open(string), 'fasta'))
        seqdict = dict((seq, str(seqdict[seq].seq)) for seq in seqdict)
    else:
        seqdict = SeqIO.index(string, "fasta")

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
        """
        Simple function to return the absolute value of the integer of the input string.
        :param string:
        :return:
        """

        return abs(int(string))

    parser = argparse.ArgumentParser("""Script to prepare a GTF for the pipeline;
    it will perform the following operations:
    1- add the "transcript" feature
    2- sort by coordinates
    3- check the strand""")
    parser.add_argument("--fasta", type=argparse.FileType(),
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
    parser.add_argument("-m", "--minimum_length", default=200, type=positive,
                        help="Minimum length for transcripts. Default: 200 bps.")
    parser.add_argument("-t", "--threads",
                        help="Number of processors to use (default %(default)s)",
                        type=to_cpu_count, default=1)
    parser.add_argument("-scds", "--strip_cds", action="store_true", default=False,
                        help="Boolean flag. If set, ignores any CDS/UTR segment.")
    parser.add_argument("--labels", type=str, default="",
                        help="""Labels to attach to the IDs of the transcripts of the input files,
                        separated by comma.""")
    parser.add_argument("--single", action="store_true", default=False,
                        help="Disable multi-threading. Useful for debugging.")
    parser.add_argument("-o", "--out", default=None,
                        help="Output file. Default: mikado_prepared.fasta.")
    parser.add_argument("-of", "--out_fasta", default=None,
                        help="Output file. Default: mikado_prepared.fasta.")
    parser.add_argument("--cache", default=False, action="store_true",
                        help="Whether to load the whole genome in memory or not.")
    parser.add_argument("--json-conf", dest="json_conf",
                        type=to_json, default="",
                        help="Configuration file.")
    parser.add_argument("gff", help="Input GFF/GTF file(s).", nargs="*")
    parser.set_defaults(func=prepare)
    return parser
