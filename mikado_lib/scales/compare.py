#!/usr/bin/env python3
# coding: utf-8

import sys
import argparse
import os
import re
import collections
import multiprocessing
import logging
from logging import handlers as log_handlers
from mikado_lib.loci_objects.transcript import Transcript
from mikado_lib.scales.assigner import Assigner
from mikado_lib.scales.reference_gene import Gene
from mikado_lib.parsers.GFF import GFF3
from mikado_lib.scales.accountant import Accountant

__author__ = 'Luca Venturini'


def compare(args):
    """
    This function performs the comparison between two different files.

    :param args: the argparse Namespace
    :type args: argparse.Namespace

    """

    # Flags for the parsing
    if type(args.reference) is GFF3:
        ref_gff = True
    else:
        ref_gff = False

    if os.path.dirname(args.out) != '' and os.path.dirname(args.out) != os.path.dirname(os.path.abspath(".")):
        dirname = os.path.dirname(args.out)
        if os.path.exists(dirname):
            assert os.path.isdir(dirname)
        else:
            os.makedirs(dirname)

    context = multiprocessing.get_context()  # @UndefinedVariable
    manager = context.Manager()

    logger = logging.getLogger("main")
    formatter = logging.Formatter("{asctime} - {name} - {levelname} - {message}", style="{")
    args.log_queue = manager.Queue()
    args.queue_handler = log_handlers.QueueHandler(args.log_queue)
    log_queue_listener = log_handlers.QueueListener(args.log_queue, logger)
    log_queue_listener.propagate = False
    log_queue_listener.start()

    if args.log is None:
        handler = logging.StreamHandler()
    else:
        if os.path.exists(args.log):
            os.remove(args.log)

        handler = logging.FileHandler(args.log)
    handler.setFormatter(formatter)

    if args.verbose is False:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)

    queue_logger = logging.getLogger("main_queue")
    if args.verbose is False:
        queue_logger.setLevel(logging.INFO)
    else:
        queue_logger.setLevel(logging.DEBUG)
    main_queue_handler = log_handlers.QueueHandler(args.log_queue)
    queue_logger.addHandler(main_queue_handler)

    queue_logger.propagate = False
    queue_logger.info("Start")
    args.commandline = " ".join(sys.argv)
    queue_logger.info("Command line: {0}".format(args.commandline))

    queue_logger.info("Starting parsing the reference")

    genes = dict()
    positions = collections.defaultdict(dict)

    transcript2gene = dict()
    logger.handlers[0].flush()

    for row in args.reference:
        # Assume we are going to use GTF for the moment
        if row.header is True:
            continue
        if row.is_transcript is True:
            queue_logger.debug("Transcript\n{0}".format(str(row)))
            tr = Transcript(row)
            transcript2gene[row.id] = row.gene
            if row.gene not in genes:
                genes[row.gene] = Gene(tr, gid=row.gene)
            genes[row.gene].add(tr)
            assert tr.id in genes[row.gene].transcripts
        elif row.is_exon is True:
            if ref_gff is True:
                for tr in row.transcript:
                    gid = transcript2gene[tr]
                    genes[gid][tr].add_exon(row)
            else:
                try:
                    genes[row.gene][row.transcript].add_exon(row)
                except KeyError as exc:
                    assert row.gene in genes
                    queue_logger.exception(exc)
                    queue_logger.exception("Keys for {0}: {1}".format(row.gene, genes[row.gene].transcripts.keys()))
                    raise
        else:
            continue
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
                    logger.debug("No CDS for {0}".format(tid))
            if len(to_remove) == len(genes[gid].transcripts):
                non_coding_to_remove.add(gid)
                logger.debug("Noncoding gene: {0}".format(gid))
                continue
            elif len(to_remove) > 0:
                for tid in to_remove:
                    genes[gid].remove(tid)
        key = (genes[gid].start, genes[gid].end)
        if key not in positions[genes[gid].chrom]:
            positions[genes[gid].chrom][key] = []
        positions[genes[gid].chrom][key].append(genes[gid])

    for gid in genes_to_remove:
        queue_logger.warn("Removed from reference: {0}; error: {1}".format(gid, genes[gid].exception_message))
        del genes[gid]
    for gid in non_coding_to_remove:
        del genes[gid]

    if len(genes) == 0:
        raise KeyError("No genes remained for the reference!")

    # Needed for refmap
    queue_logger.info("Finished preparation; found {0} reference genes".format(len(genes)))
    queue_logger.debug("Gene names (first 20): {0}".format("\n\t".join(list(genes.keys())[:20])))

    accountant_instance = Accountant(genes, args)  # start the class which will manage the statistics
    # genes: a dictionary with the reference annotation, indexed by GID
    # positions: a dictionary of the form dict[chrom][(start,end)] = [gene object]
    assigner_instance = Assigner(genes, positions, args, accountant_instance)

    transcript = None
    for row in args.prediction:
        if row.header is True:
            continue
        #         queue_logger.debug("Row:\n{0:>20}".format(str(row)))
        if row.is_transcript is True:
            queue_logger.debug("Transcript row:\n{0}".format(str(row)))
            if transcript is not None:
                if re.search("\.orf[0-9]+$", transcript.id) and not transcript.id.endswith("orf1"):
                    pass
                else:
                    try:
                        assigner_instance.get_best(transcript)
                    except Exception as err:
                        queue_logger.exception(err)
                        log_queue_listener.enqueue_sentinel()
                        handler.close()
                        log_queue_listener.stop()
                        args.queue_handler.close()
                        raise
            transcript = Transcript(row)
        elif row.is_exon is True:
            try:
                transcript.add_exon(row)
            except Exception as err:
                queue_logger.exception(err)
                break
        else:
            continue

    if transcript is not None:
        if re.search("\.orf[0-9]+$", transcript.id) and not transcript.id.endswith("orf1"):
            pass
        else:
            try:
                assigner_instance.get_best(transcript)
            except Exception as err:
                queue_logger.exception(err)
                log_queue_listener.enqueue_sentinel()
                handler.close()
                log_queue_listener.stop()
                args.queue_handler.close()
                raise

    assigner_instance.finish()
    queue_logger.info("Finished")
    log_queue_listener.enqueue_sentinel()
    handler.close()
    log_queue_listener.stop()
    args.queue_handler.close()
    return
