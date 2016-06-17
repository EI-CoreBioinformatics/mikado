#!/usr/bin/env python3
# coding: utf-8

"""
Launcher for the Mikado compare utility.
Launcher for the Mikado compare utility.
"""

import collections
import logging
import multiprocessing
import os
import re
import sys
import csv
from logging import handlers as log_handlers
from ..loci.reference_gene import Gene
from .accountant import Accountant
from .assigner import Assigner
from .resultstorer import ResultStorer
from ..loci.transcript import Transcript
from ..parsers.GFF import GFF3
from ..utilities.log_utils import create_default_logger
from ..exceptions import CorruptIndex
from ..utilities import to_gff
import json
import gzip
import itertools

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


def finalize_reference(genes, positions, queue_logger, args) \
        -> (dict, collections.defaultdict(dict)):

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
        genes[gid].logger = queue_logger
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
        positions[genes[gid].chrom][key].append(gid)

    for gid in genes_to_remove:
        queue_logger.warn("Removed from reference: %s; error: %s",
                          gid, genes[gid].exception_message)
        del genes[gid]
    for gid in non_coding_to_remove:
        del genes[gid]
    return genes, positions


def prepare_reference(args, queue_logger, ref_gff=False) -> (dict, collections.defaultdict(dict)):

    """
    Method to prepare the data structures that hold the reference
    information for the parsing.
    :param args:
    :param queue_logger:
    :param ref_gff:
    :return: genes, positions
    """

    genes = dict()
    positions = collections.defaultdict(dict)
    transcript2gene = dict()

    for row in args.reference:
        # Assume we are going to use GTF for the moment
        if row.is_transcript is True:
            queue_logger.debug("Transcript\n%s", str(row))
            transcript = Transcript(row, logger=queue_logger)
            transcript2gene[row.id] = row.gene
            if row.gene not in genes:
                genes[row.gene] = Gene(transcript,
                                       gid=row.gene,
                                       logger=queue_logger)
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
                        genes[row.gene] = Gene(None, gid=row.gene, logger=queue_logger)
                    if row.transcript not in genes[row.gene]:
                        transcript = Transcript(row, logger=queue_logger)
                        transcript2gene[row.id] = row.gene
                        genes[row.gene].add(transcript)
                    genes[row.gene][row.transcript].add_exon(row)

    genes, positions = finalize_reference(genes, positions, queue_logger, args)

    if len(genes) == 0:
        raise KeyError("No genes remained for the reference!")
    return genes, positions


def parse_prediction(args, genes, positions, queue_logger):

    """
    This function performs the real comparison between the reference and the prediction.
     It needs the following inputs:
    :param args: the Namespace with the necessary parameters
    :param genes: Dictionary with the reference genes, of the form
    dict[chrom][(start,end)] = [gene object]
    :param positions: Dictionary with the positions of the reference genes, of the form
    dict[chrom][IntervalTree]
    :param queue_logger: Logger
    :return:
    """

    # start the class which will manage the statistics
    accountant_instance = Accountant(genes, args)
    assigner_instance = Assigner(genes, positions, args, accountant_instance)

    transcript = None
    if hasattr(args, "self") and args.self is True:
        args.prediction = to_gff(args.reference.name)
    ref_gff = isinstance(args.prediction, GFF3)
    __found_with_orf = set()

    for row in args.prediction:
        if row.header is True:
            continue
        #         queue_logger.debug("Row:\n{0:>20}".format(str(row)))
        if row.is_transcript is True:
            queue_logger.debug("Transcript row:\n%s", str(row))
            if transcript is not None:
                if re.search(r"\.orf[0-9]+$", transcript.id):
                    __name = re.sub(r"\.orf[0-9]+$", "", transcript.id)
                    if __name not in __found_with_orf:
                        __found_with_orf.add(__name)
                        assigner_instance.get_best(transcript)
                    else:
                        pass
                else:
                    assigner_instance.get_best(transcript)
            transcript = Transcript(row, logger=queue_logger)
        elif row.is_exon is True:
            if ref_gff is False:
                if transcript is None or (transcript is not None and transcript.id != row.transcript):
                    if transcript is not None:
                        if re.search(r"\.orf[0-9]+$", transcript.id) and \
                                (not transcript.id.endswith("orf1")):
                            pass
                        else:
                            assigner_instance.get_best(transcript)
                    queue_logger.debug("New transcript: %s", row.transcript)
                    transcript = Transcript(row, logger=queue_logger)

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


def parse_self(args, genes, queue_logger):

    """
    This function is called when we desire to compare a reference against itself.
    :param args:
    :param genes:
    :param queue_logger:
    :return:
    """

    queue_logger.info("Performing a self-comparison.")

    tmap_out = open("{0}.tmap".format(args.out), 'wt')
    tmap_rower = csv.DictWriter(tmap_out, ResultStorer.__slots__, delimiter="\t")
    tmap_rower.writeheader()

    for gid, gene in genes.items():
        assert isinstance(gene, Gene)
        if len(gene.transcripts) == 1:
            continue

        combinations = itertools.combinations(gene.transcripts.keys(), 2)
        for combination in combinations:
            result, _ = Assigner.compare(gene.transcripts[combination[0]],
                                         gene.transcripts[combination[1]])
            tmap_rower.writerow(result.as_dict())

    queue_logger.info("Finished.")


def load_index(args, queue_logger):

    """
    Function to load the genes and positions from the indexed GFF.
    :param args:
    :param queue_logger:
    :return: genes, positions
    :rtype: ((None|collections.defaultdict),(None|collections.defaultdict))
    """

    genes, positions = None, None
    with gzip.open("{0}.midx".format(args.reference.name), "rt") as index:
        positions = collections.defaultdict(dict)
        try:
            cp_genes = json.load(index)
            genes = dict()
            for gid, gobj in cp_genes.items():
                gene = Gene(None, logger=queue_logger)
                gene.load_dict(gobj,
                               exclude_utr=args.exclude_utr,
                               protein_coding=args.protein_coding)
                # Necessary for when we are excluding non-coding
                if len(gene.transcripts) > 0:
                    genes[gid] = gene
                    if (gene.start, gene.end) not in positions[gene.chrom]:
                        positions[gene.chrom][(gene.start, gene.end)] = []
                    positions[gene.chrom][(gene.start, gene.end)].append(gene.id)
                else:
                    queue_logger.warning("No transcripts for %s", gid)

            if not (isinstance(genes, dict) and
                    isinstance(positions, collections.defaultdict) and
                    positions.default_factory is dict):
                raise EOFError
        except (EOFError, json.decoder.JSONDecodeError) as exc:
            queue_logger.exception(exc)
            raise CorruptIndex("Invalid index file")
        except (TypeError, ValueError) as exc:
            genes, positions = None, None
            queue_logger.exception(exc)
            raise CorruptIndex("Corrupted index file; deleting and rebuilding.")

    return genes, positions


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

    logger.handlers[0].flush()

    genes, positions = None, None

    if args.index is True:
        queue_logger.info("Starting to create an index for %s", args.reference.name)
        if os.path.exists("{0}.midx".format(args.reference.name)):
            queue_logger.warning("Removing the old index")
            os.remove("{0}.midx".format(args.reference.name))
        args.protein_coding = False
        args.exclude_utr = False
        genes, positions = prepare_reference(args,
                                             queue_logger,
                                             ref_gff=ref_gff)
        with gzip.open("{0}.midx".format(args.reference.name), "wt") as index:
                cp_genes = dict()
                for gid in genes:
                    cp_genes[gid] = genes[gid].as_dict()
                json.dump(cp_genes, index)
        queue_logger.info("Finished to create an index for %s, with %d genes",
                          args.reference.name, len(genes))
    else:
        if (os.path.exists("{0}.midx".format(args.reference.name)) and
            os.stat(args.reference.name).st_mtime >= os.stat(
                    "{0}.midx".format(args.reference.name)).st_mtime):
            queue_logger.warning("Reference index obsolete, deleting and rebuilding.")
            os.remove("{0}.midx".format(args.reference.name))
        elif os.path.exists("{0}.midx".format(args.reference.name)):
            queue_logger.info("Starting loading the indexed reference")
            try:
                genes, positions = load_index(args, queue_logger)
            except CorruptIndex:
                queue_logger.warning("Reference index corrupt, deleting and rebuilding.")
                os.remove("{0}.midx".format(args.reference.name))
                genes, positions = None, None
        if genes is None:
            queue_logger.info("Starting parsing the reference")
            exclude_utr, protein_coding = args.exclude_utr, args.protein_coding
            if args.no_save_index is False:
                args.exclude_utr, args.protein_coding = False, False
            genes, positions = prepare_reference(args,
                                                 queue_logger,
                                                 ref_gff=ref_gff)
            if args.no_save_index is False:
                with gzip.open("{0}.midx".format(args.reference.name), "wt") as index:
                    cp_genes = dict()
                    for gid in genes:
                        cp_genes[gid] = genes[gid].as_dict()
                    json.dump(cp_genes, index)
                if exclude_utr is True or protein_coding is True:
                    args.exclude_utr, args.protein_coding = exclude_utr, protein_coding
                    positions = collections.defaultdict(dict)
                    finalize_reference(genes, positions, queue_logger, args)

        assert isinstance(genes, dict)

        # Needed for refmap
        queue_logger.info("Finished preparation; found %d reference gene%s",
                          len(genes), "s" if len(genes) > 1 else "")
        queue_logger.debug("Gene names (first 20): %s",
                           "\n\t".join(list(genes.keys())[:20]))

        try:
            if hasattr(args, "internal") and args.internal is True:
                parse_self(args, genes, queue_logger)
            else:
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
