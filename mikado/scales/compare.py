#!/usr/bin/env python3
# coding: utf-8

"""
Launcher for the mikado compare utility.
"""

import collections
import logging
import multiprocessing
import os
import re
import sys
from logging import handlers as log_handlers
from ..loci_objects.reference_gene import Gene
from .accountant import Accountant
from .assigner import Assigner
from ..loci_objects.transcript import Transcript
from ..parsers.GFF import GFF3
from ..utilities.log_utils import create_default_logger

__author__ = 'Luca Venturini'


def setup_logger(args, manager):

    """
    Function to setup the logger for the compare function.
    :param args:
    :param manager:
    :return:
    """

    logger = create_default_logger("main_compare")
    args.log_queue = manager.Queue()
    args.queue_handler = log_handlers.QueueHandler(args.log_queue)
    log_queue_listener = log_handlers.QueueListener(args.log_queue, logger)
    log_queue_listener.propagate = False
    log_queue_listener.start()

    if args.log is not None:
        if os.path.exists(args.log):
            os.remove(args.log)
        handler = logging.FileHandler(args.log)
        handler.setFormatter(logger.handlers[0].formatter)
        # Remove stream handler
        logger.removeHandler(logger.handlers[0])
        logger.addHandler(handler)
    else:
        handler = logger.handlers[0]

    if args.verbose is False:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)

    logger.propagate = False

    queue_logger = logging.getLogger("main_queue")
    for handler in queue_logger.handlers:
        queue_logger.removeHandler(handler)
    if args.verbose is False:
        queue_logger.setLevel(logging.INFO)
    else:
        queue_logger.setLevel(logging.DEBUG)
    main_queue_handler = log_handlers.QueueHandler(args.log_queue)
    queue_logger.propagate = False
    queue_logger.addHandler(main_queue_handler)

    return args, handler, logger, log_queue_listener, queue_logger


def finalize_reference(genes, positions, queue_logger, args):

    """
:param genes:
:param positions:
:param queue_logger:
:param args:
:return:
"""

    non_coding_to_remove = set()
    genes_to_remove = set()
    for gid in genes:
        genes[gid].set_logger(queue_logger)
        genes[gid].finalize(exclude_utr=args.exclude_utr)
        if len(genes[gid]) == 0:
            genes_to_remove.add(gid)
            continue
        if args.protein_coding is True:
            to_remove = []
            for tid in genes[gid].transcripts:
                if genes[gid].transcripts[tid].combined_cds_length == 0:
                    to_remove.append(tid)
                    queue_logger.debug("No CDS for %s", tid)
            if len(to_remove) == len(genes[gid].transcripts):
                non_coding_to_remove.add(gid)
                queue_logger.debug("Noncoding gene: %s", gid)
                continue
            elif len(to_remove) > 0:
                for tid in to_remove:
                    genes[gid].remove(tid)
        key = (genes[gid].start, genes[gid].end)
        if key not in positions[genes[gid].chrom]:
            positions[genes[gid].chrom][key] = []
        positions[genes[gid].chrom][key].append(genes[gid])

    for gid in genes_to_remove:
        queue_logger.warn("Removed from reference: %s; error: %s",
                          gid, genes[gid].exception_message)
        del genes[gid]
    for gid in non_coding_to_remove:
        del genes[gid]
    return genes, positions


def prepare_reference(args, queue_logger, ref_gff=False):

    """
    Method to prepare the data structures that hold the reference
    information for the parsing.
    :param args:
    :param queue_logger:
    :param ref_gff:
    :return:
    """

    genes = dict()
    positions = collections.defaultdict(dict)
    transcript2gene = dict()

    for row in args.reference:
        # Assume we are going to use GTF for the moment
        if row.is_transcript is True:
            queue_logger.debug("Transcript\n%s", str(row))
            transcript = Transcript(row)
            transcript2gene[row.id] = row.gene
            if row.gene not in genes:
                genes[row.gene] = Gene(transcript, gid=row.gene)
            genes[row.gene].add(transcript)
            assert transcript.id in genes[row.gene].transcripts
        elif row.is_exon is True:
            if ref_gff is True:
                found = False
                for transcript in row.transcript:
                    if transcript in transcript2gene:
                        # We have to perform the check because there are some GFFs
                        # e.g. TAIR
                        # where CDSs are defined within a spurious "Protein" feature
                        found = True
                        gid = transcript2gene[transcript]
                        genes[gid][transcript].add_exon(row)
                if found is False:
                    queue_logger.warn("This feature has no corresponding transcript! %s",
                                      str(row))
            else:
                if row.gene in genes and row.transcript in genes[row.gene].transcripts:
                    genes[row.gene][row.transcript].add_exon(row)
                else:
                    if row.gene not in genes:
                        exc = KeyError("Gene {0} not encountered yet, line: {1}".format(
                            row.gene, str(row)
                        ))
                    elif row.transcript not in genes[row.gene]:
                        exc = KeyError("TID {0} not found for {1}; found tids: {2}".format(
                            row.transcript, row.gene,
                            ",".join(genes[row.gene].transcripts.keys())
                        ))
                    else:
                        exc = KeyError("Weird error; both {0} and {1} are in the hash..".format(
                            row.gene, row.transcript
                        ))
                    queue_logger.exception(exc)
                    raise exc

    genes, positions = finalize_reference(genes, positions, queue_logger, args)

    if len(genes) == 0:
        raise KeyError("No genes remained for the reference!")
    return genes, positions


def parse_prediction(args, genes, positions, queue_logger):

    """
    This function performs the real comparison between the reference and the prediction.
     It needs the following inputs:
    :param args: the Namespace with the necessary parameters
    :param genes: Dictionary with the reference genes
    :param positions: Dictionary with the positions of the reference genes
    :param queue_logger: Logger
    :return:
    """

    # start the class which will manage the statistics
    accountant_instance = Accountant(genes, args)
    # genes: a dictionary with the reference annotation, indexed by GID
    # positions: a dictionary of the form dict[chrom][(start,end)] = [gene object]
    assigner_instance = Assigner(genes, positions, args, accountant_instance)

    transcript = None
    for row in args.prediction:
        if row.header is True:
            continue
        #         queue_logger.debug("Row:\n{0:>20}".format(str(row)))
        if row.is_transcript is True:
            queue_logger.debug("Transcript row:\n%s", str(row))
            if transcript is not None:
                if re.search(r"\.orf[0-9]+$", transcript.id) and \
                        (not transcript.id.endswith("orf1")):
                    pass
                else:
                    assigner_instance.get_best(transcript)
            transcript = Transcript(row)
        elif row.is_exon is True:
            queue_logger.debug("Adding exon to transcript %s: %s",
                               transcript.id, row)
            transcript.add_exon(row)
        else:
            continue

    if transcript is not None:
        if re.search(r"\.orf[0-9]+$", transcript.id) and not transcript.id.endswith("orf1"):
            pass
        else:
            assigner_instance.get_best(transcript)

    # Finish everything, including printing refmap and stats
    assigner_instance.finish()
    args.prediction.close()

def compare(args):
    """
    This function performs the comparison between two different files.

    :param args: the argparse Namespace
    """

    # Flags for the parsing
    ref_gff = isinstance(args.reference, GFF3)

    if os.path.dirname(args.out) != ''\
            and os.path.dirname(args.out) != os.path.dirname(os.path.abspath(".")):
        dirname = os.path.dirname(args.out)
        if os.path.exists(dirname):
            assert os.path.isdir(dirname)
        else:
            os.makedirs(dirname)
    # pylint: disable=no-member
    context = multiprocessing.get_context()
    manager = context.Manager()
    # pylint: enable=no-member

    args, handler, logger, log_queue_listener, queue_logger = setup_logger(
        args, manager)

#    queue_logger.propagate = False
    queue_logger.info("Start")
    args.commandline = " ".join(sys.argv)
    queue_logger.info("Command line: %s", args.commandline)

    queue_logger.info("Starting parsing the reference")

    logger.handlers[0].flush()

    genes, positions, = prepare_reference(args,
                                          queue_logger,
                                          ref_gff=ref_gff)

    # Needed for refmap
    queue_logger.info("Finished preparation; found %d reference gene%s",
                      len(genes), "s" if len(genes) > 1 else "")
    queue_logger.debug("Gene names (first 20): %s",
                       "\n\t".join(list(genes.keys())[:20]))

    try:
        parse_prediction(args, genes, positions, queue_logger)
    except Exception as err:
        queue_logger.exception(err)
        log_queue_listener.enqueue_sentinel()
        handler.close()
        log_queue_listener.stop()
        args.queue_handler.close()
        raise

    queue_logger.info("Finished")
    log_queue_listener.enqueue_sentinel()
    handler.close()
    log_queue_listener.stop()
    args.queue_handler.close()
    return
