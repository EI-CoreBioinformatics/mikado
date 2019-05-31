#!/usr/bin/env python3
# coding: utf-8

"""
Launcher for the Mikado compare utility.
Launcher for the Mikado compare utility.
"""

import collections
import csv
import logging
import multiprocessing
import os
import re
import sys
from logging import handlers as log_handlers
from ..transcripts.transcript import Transcript
from .assigner import Assigner
from .resultstorer import ResultStorer
from ..exceptions import CorruptIndex
from ..loci.reference_gene import Gene
from ..parsers.GFF import GFF3
from ..parsers.bed12 import Bed12Parser
from ..parsers.bam_parser import BamParser
from ..parsers import to_gff
from ..utilities.log_utils import create_default_logger, formatter
import magic
import sqlite3
import multiprocessing as mp
import tempfile
from .gene_dict import check_index
from ..exceptions import InvalidTranscript
import json
import gzip
import itertools
import functools
import msgpack
from .gene_dict import GeneDict


__author__ = 'Luca Venturini'


# Hack to give the ujson library this exception class
# This becomes necessary when we happen to have a corrupted index
if not hasattr(json, "decoder"):

    class Decoder:
        class JSONDecodeError(TypeError):
            pass
    json.decoder = Decoder


def setup_logger(args, manager):

    """
    Function to setup the logger for the compare function.
    :param args:
    :param manager:
    :return:
    """

    args.log_queue = manager.Queue()
    args.queue_handler = log_handlers.QueueHandler(args.log_queue)

    if args.log is not None:
        _log_folder = os.path.dirname(args.log)
        if _log_folder and not os.path.exists(_log_folder):
            os.makedirs(_log_folder)
        handler = logging.FileHandler(args.log, mode="w")
        logger = logging.getLogger("main_compare")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.propagate = False
    else:
        logger = create_default_logger("main_compare")
        handler = logger.handlers[0]

    if args.verbose is False:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)

    logger.propagate = False
    log_queue_listener = log_handlers.QueueListener(args.log_queue, logger)
    log_queue_listener.propagate = False
    log_queue_listener.start()

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


def finalize_reference(genes, positions, queue_logger, exclude_utr=False, protein_coding=False) \
        -> (dict, collections.defaultdict):

    """
:param genes:
:param positions:
:param queue_logger:
:param exclude_utr:
:param protein_coding:
:return:
"""

    non_coding_to_remove = set()
    genes_to_remove = set()
    for gid in genes:
        genes[gid].logger = queue_logger
        genes[gid].finalize(exclude_utr=exclude_utr)
        if len(genes[gid]) == 0:
            genes_to_remove.add(gid)
            continue
        if protein_coding is True:
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


def prepare_reference(reference, queue_logger, ref_gff=False,
                      exclude_utr=False, protein_coding=False) -> (dict, collections.defaultdict):

    """
    Method to prepare the data structures that hold the reference
    information for the parsing.
    :param reference:
    :param queue_logger:
    :param ref_gff:
    :return: genes, positions
    :param exclude_utr:
    :param protein_coding:
    """

    genes = dict()
    positions = collections.defaultdict(dict)
    transcript2gene = dict()

    for row in reference:
        # Assume we are going to use GTF for the moment
        if row.header is True:
            continue
        elif reference.__annot_type__ == Bed12Parser.__annot_type__:
            transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
            if transcript.parent:
                gid = transcript.parent[0]
            else:
                transcript.parent = gid = row.id
            transcript2gene[row.id] = gid
            if gid not in genes:
                genes[gid] = Gene(transcript, gid=gid, logger=queue_logger)
            genes[gid].add(transcript)

        elif row.is_transcript is True or (ref_gff is True and row.feature == "match"):
            queue_logger.debug("Transcript\n%s", str(row))
            transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
            if row.feature == "match":
                gid = row.id
            else:
                gid = row.gene

            transcript2gene[row.id] = gid
            if gid not in genes:
                genes[gid] = Gene(transcript, gid=gid, logger=queue_logger)
            genes[gid].add(transcript)
            assert transcript.id in genes[gid].transcripts
        elif row.is_exon is True:
            if ref_gff is True:
                if "cDNA_match" in row.feature:
                    row.parent = row.id
                    # row.gene = row.id
                    if row.id not in transcript2gene:
                        genes[row.id] = Gene(None, gid=row.id, logger=queue_logger)
                        transcript2gene[row.id] = row.id
                        transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
                        genes[row.id].add(transcript)
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
                        transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
                        transcript2gene[row.id] = row.gene
                        genes[row.gene].add(transcript)
                    genes[row.gene][row.transcript].add_exon(row)

    genes, positions = finalize_reference(genes, positions, queue_logger, exclude_utr=exclude_utr,
                                          protein_coding=protein_coding)

    if len(genes) == 0:
        raise KeyError("No genes remained for the reference!")
    return genes, positions


class Assigners(mp.Process):

    def __init__(self, index, args, queue, returnqueue, counter, dump_dbname):
        super().__init__()
        self.__connection = sqlite3.connect(dump_dbname, timeout=600)
        self.__cursor = self.__connection.cursor()

        # self.accountant_instance = Accountant(genes, args, counter=counter)
        if hasattr(args, "fuzzymatch"):
            self.__fuzzymatch = args.fuzzymatch
        else:
            self.__fuzzymatch = 0

        self.assigner_instance = Assigner(index, args,
                                          printout_tmap=False,
                                          counter=counter,
                                          fuzzymatch=self.__fuzzymatch)

        self.queue = queue
        self.returnqueue = returnqueue
        self.__args = args
        self.__counter = counter

    def run(self):
        while True:
            transcr = self.queue.get()
            #
            if transcr == "EXIT":
                self.queue.put_nowait("EXIT")
                # out = os.path.join(self.__args.out + str(self.__counter) + ".pickle")
                self.assigner_instance.dump()
                self.returnqueue.put(self.assigner_instance.db.name)
                break
            else:
                dumped = self.__cursor.execute("SELECT json FROM dump WHERE idx=?", (transcr,)).fetchone()
                dumped = json.loads(dumped[0])
                transcr = Transcript()
                transcr.load_dict(dumped, trust_orf=True)
                self.assigner_instance.get_best(transcr)

        self.__connection.close()


class FinalAssigner(mp.Process):

    def __init__(self, index: str, args, queue):

        super().__init__()
        self.index, self.args = index, args
        self.queue = queue

    def run(self):
        self.assigner = Assigner(self.index, self.args, printout_tmap=True)
        while True:
            dbname = self.queue.get()
            if dbname == "EXIT":
                break
            self.assigner.load_result(dbname)
        self.assigner.finish()


def transmit_transcript(index, transcript: Transcript, connection: sqlite3.Connection):
    transcript.finalize()
    connection.execute("INSERT INTO dump VALUES (?, ?)",
                       (index, json.dumps(transcript.as_dict(remove_attributes=True))))


def get_best_result(_, transcript, assigner_instance: Assigner):
    transcript.finalize()
    assigner_instance.get_best(transcript)


def parse_prediction(args, index, queue_logger):

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
    transcript = None
    if hasattr(args, "self") and args.self is True:
        args.prediction = to_gff(args.reference.name)
    ref_gff = isinstance(args.prediction, GFF3)
    __found_with_orf = set()
    is_bed12 = isinstance(args.prediction, Bed12Parser)
    queue = mp.JoinableQueue(-1)
    returnqueue = mp.JoinableQueue(-1)
    dump_dbhandle = tempfile.NamedTemporaryFile(delete=True, prefix=".compare_dump", suffix=".db", dir=".")
    dump_db = sqlite3.connect(dump_dbhandle.name, check_same_thread=False, timeout=60)
    dump_db.execute("CREATE TABLE dump (idx INTEGER, json BLOB)")
    dump_db.execute("CREATE UNIQUE INDEX dump_idx ON dump(idx)")
    dump_db.commit()

    queue_logger.info("Starting to parse the prediction")
    if args.processes > 1:
        procs = [Assigners(index, args, queue, returnqueue, counter, dump_dbhandle.name)
                            for counter in range(1, args.processes)]
        final_proc = FinalAssigner(index, args, returnqueue)
        [proc.start() for proc in procs]
        final_proc.start()
        transmitter = functools.partial(transmit_transcript, connection=dump_db)
        assigner_instance = None
    else:
        procs = []
        assigner_instance = Assigner(index, args, printout_tmap=True, )
        transmitter = functools.partial(get_best_result, assigner_instance=assigner_instance)

    constructor = functools.partial(Transcript,
                                    logger=queue_logger, trust_orf=True, accept_undefined_multi=True)

    done = 0
    invalids = set()
    name_counter = collections.Counter()  # This is needed for BAMs

    lastdone = 1

    for row in args.prediction:
        if args.prediction.__annot_type__ == BamParser.__annot_type__:
            if row.is_unmapped is True:
                continue
            if transcript is not None:
                done += 1
                if done and done % 10000 == 0:
                    queue_logger.info("Parsed %s transcripts", done)
                    if assigner_instance is None:
                        dump_db.commit()
                        [queue.put(_) for _ in range(lastdone, done)]
                    lastdone = done

                transmitter(done, transcript)
            transcript = Transcript(row, accept_undefined_multi=True)
            if name_counter.get(row.query_name):
                name = "{}_{}".format(row.query_name, name_counter.get(row.query_name))
            else:
                name = row.query_name
            transcript.id = transcript.name = transcript.alias = name
            transcript.parent = transcript.attributes["gene_id"] = "{0}.gene".format(name)
        elif row.header is True:
            continue
        elif (row.is_transcript is True or
              args.prediction.__annot_type__ == Bed12Parser.__annot_type__ or
              row.feature == "match"):
            queue_logger.debug("Transcript row:\n%s", str(row))
            if transcript is not None:
                if re.search(r"\.orf[0-9]+$", transcript.id):
                    __name = re.sub(r"\.orf[0-9]+$", "", transcript.id)
                    if __name not in __found_with_orf:
                        __found_with_orf.add(__name)
                        done += 1
                        if done and done % 10000 == 0:
                            queue_logger.info("Parsed %s transcripts", done)
                            if assigner_instance is None:
                                dump_db.commit()
                                [queue.put(_) for _ in range(lastdone, done)]
                            lastdone = done
                        transmitter(done, transcript)
                    else:
                        pass
                else:
                    done += 1
                    if done and done % 10000 == 0:
                        queue_logger.info("Parsed %s transcripts", done)
                        if assigner_instance is None:
                            dump_db.commit()
                            [queue.put(_) for _ in range(lastdone, done)]
                        lastdone = done
                    try:
                        transmitter(done, transcript)
                    except AssertionError:
                        raise AssertionError((transcript.id, str(row)))
            try:
                transcript = constructor(row)
            except (InvalidTranscript, AssertionError, TypeError, ValueError):
                queue_logger.warning("Row %s is invalid, skipping.", row)
                transcript = None
                invalids.add(row.id)
                continue

            if is_bed12:
                transcript.parent = transcript.id

        elif row.is_exon is True:
            # Case 1: we are talking about cDNA_match and GFF
            if any(_ in invalids for _ in row.parent):
                # Skip children of invalid things
                continue
            if ref_gff is True and "match" not in row.feature:
                if transcript is None:
                    raise TypeError("Transcript not defined inside the GFF; line:\n{}".format(row))
                else:
                    queue_logger.debug("Adding exon to transcript %s: %s", transcript.id, row)
                    transcript.add_exon(row)
            elif ref_gff is True and "match" in row.feature:
                if transcript is not None and row.id == transcript.id:
                    transcript.add_exon(row)
                elif transcript is not None and transcript.id in row.parent:
                    transcript.add_exon(row)
                elif transcript is None or (transcript is not None and row.id != transcript.id):
                    if transcript is not None:
                        if re.search(r"\.orf[0-9]+$", transcript.id) and \
                                (not transcript.id.endswith("orf1")):
                            pass
                        else:
                            done += 1
                            if done and done % 10000 == 0:
                                queue_logger.info("Done %s transcripts", done)
                                if assigner_instance is None:
                                    dump_db.commit()
                                    [queue.put(_) for _ in range(lastdone, done)]
                                lastdone = done

                            transmitter(done, transcript)
                    queue_logger.debug("New transcript: %s", row.transcript)
                    transcript = constructor(row)
                assert len(transcript.exons) > 0, (row,)
            elif ref_gff is False:
                if transcript is None or (transcript is not None and transcript.id != row.transcript):
                    if transcript is not None:
                        if re.search(r"\.orf[0-9]+$", transcript.id) and \
                                (not transcript.id.endswith("orf1")):
                            pass
                        else:
                            done += 1
                            if done and done % 10000 == 0:
                                queue_logger.info("Done %s transcripts", done)
                                if assigner_instance is None:
                                    dump_db.commit()
                                    [queue.put(_) for _ in range(lastdone, done)]
                                lastdone = done
                            transmitter(done, transcript)
                    queue_logger.debug("New transcript: %s", row.transcript)
                    transcript = constructor(row)
                transcript.add_exon(row)
            else:
                raise TypeError("Unmatched exon: {}".format(row))

        elif row.header:
            continue
        else:
            queue_logger.debug("Skipped row: {}".format(row))

    if transcript is not None:
        if re.search(r"\.orf[0-9]+$", transcript.id) and not transcript.id.endswith("orf1"):
            pass
        else:
            done += 1
            if done and done % 10000 == 0:
                queue_logger.info("Done %s transcripts", done)
                if assigner_instance is None:
                    [queue.put(_) for _ in range(lastdone, done)]
                lastdone = done
            transmitter(done, transcript)

    if assigner_instance is None:
        dump_db.commit()
        [queue.put(_) for _ in range(lastdone, done + 1)]
    queue_logger.info("Finished parsing, %s transcripts in total", done)
    dump_dbhandle.close()
    if assigner_instance is None:
        queue.put("EXIT")
        [proc.join() for proc in procs]
        returnqueue.put("EXIT")
        final_proc.join()
    else:
        returnqueue.close()
        assert not procs, procs
        assert isinstance(assigner_instance, Assigner)
        assigner_instance.finish()
    queue.close()


def parse_self(args, gdict: GeneDict, queue_logger):

    """
    This function is called when we desire to compare a reference against itself.
    :param args:
    :param gdict: gene MIDX database
    :param queue_logger:
    :return:
    """

    queue_logger.info("Performing a self-comparison.")

    if args.gzip is False:
        tmap_out = open("{0}.tmap".format(args.out), 'wt')
    else:
        tmap_out = gzip.open("{0}.tmap.gz".format(args.out), 'wt')
    tmap_rower = csv.DictWriter(tmap_out, ResultStorer.__slots__, delimiter="\t")
    tmap_rower.writeheader()

    for gid, gene in gdict.items():
        assert isinstance(gene, Gene)
        if len(gene.transcripts) == 1:
            continue

        combinations = itertools.combinations(gene.transcripts.keys(), 2)
        for combination in combinations:
            result, _ = Assigner.compare(gene.transcripts[combination[0]],
                                         gene.transcripts[combination[1]],
                                         fuzzy_match=args.fuzzymatch)
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

    # genes, positions = None, None

    # New: now we are going to use SQLite for a faster experience
    wizard = magic.Magic(mime=True)

    if wizard.from_file("{0}.midx".format(args.reference.name)) == b"application/gzip":
        queue_logger.warning("Old index format detected. Starting to generate a new one.")
        raise CorruptIndex("Invalid index file")

    try:
        conn = sqlite3.connect("{0}.midx".format(args.reference.name))
        cursor = conn.cursor()
        tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
        if sorted(tables) != sorted([("positions",), ("genes",)]):
            raise CorruptIndex("Invalid database file")
        # Integrity check
        res = cursor.execute("PRAGMA integrity_check;").fetchone()
        if res[0] != "ok":
            raise CorruptIndex("Corrupt database, integrity value: {}".format(res[0]))

    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid database file")

    positions = dict()
    try:
        for counter, obj in enumerate(cursor.execute("SELECT * from positions")):
            chrom, start, end, gid = obj
            if chrom not in positions:
                positions[chrom] = collections.defaultdict(list)
            positions[chrom][(start, end)].append(gid)
    except sqlite3.DatabaseError:
        raise CorruptIndex("Invalid index file. Rebuilding.")

    genes = dict()
    for gid, obj in cursor.execute("SELECT * from genes"):
        try:
            gene = Gene(None, logger=queue_logger)
            gene.load_dict(msgpack.loads(obj, raw=False),
                           exclude_utr=args.exclude_utr,
                           protein_coding=args.protein_coding,
                           trust_orf=True)
            if len(gene.transcripts) > 0:
                genes[gid] = gene
            else:
                queue_logger.warning("No transcripts for %s", gid)
        except (EOFError, json.decoder.JSONDecodeError) as exc:
            queue_logger.exception(exc)
            raise CorruptIndex("Invalid index file")
        except (TypeError, ValueError) as exc:
            queue_logger.exception(exc)
            raise CorruptIndex("Corrupted index file; deleting and rebuilding.")

    return genes, positions


def create_index(reference, queue_logger, index_name, ref_gff=False,
                 exclude_utr=False, protein_coding=False):

    """Method to create the simple indexed database for features."""

    import tempfile
    import shutil

    temp_db = tempfile.mktemp(suffix=".db")

    conn = sqlite3.connect(temp_db)
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE positions (chrom text, start integer, end integer, gid text)")
    cursor.execute("CREATE INDEX pos_idx ON positions (chrom, start, end)")
    genes, positions = prepare_reference(reference, queue_logger, ref_gff=ref_gff,
                                         exclude_utr=exclude_utr, protein_coding=protein_coding)

    for chrom in positions:
        for key in positions[chrom]:
            start, end = key
            for gid in positions[chrom][key]:
                cursor.execute("INSERT INTO positions VALUES (?, ?, ?, ?)", (chrom,
                                                                             start,
                                                                             end,
                                                                             gid))
    cursor.execute("CREATE TABLE genes (gid text, json blob)")
    cursor.execute("CREATE INDEX gid_idx on genes(gid)")
    for gid, gobj in genes.items():
        cursor.execute("INSERT INTO genes VALUES (?, ?)", (gid, msgpack.dumps(gobj.as_dict()))) #, cls=NumpyEncoder)))
    cursor.close()
    conn.commit()
    conn.close()

    shutil.move(temp_db, index_name)

    return


def compare(args):
    """
    This function performs the comparison between two different files.

    :param args: the argparse Namespace
    """

    # Flags for the parsing

    ref_gff = isinstance(args.reference, GFF3)

    # pylint: disable=no-member
    context = multiprocessing.get_context()
    manager = context.Manager()
    # pylint: enable=no-member

    if isinstance(args.out, str):
        _out_folder = os.path.dirname(args.out)
        if _out_folder and not os.path.exists(_out_folder):
            os.makedirs(_out_folder)

    args, handler, logger, log_queue_listener, queue_logger = setup_logger(
        args, manager)

#    queue_logger.propagate = False
    queue_logger.info("Start")
    args.commandline = " ".join(sys.argv)
    queue_logger.info("Command line: %s", args.commandline)

    logger.handlers[0].flush()

    index_name = os.path.abspath("{0}.midx".format(args.reference.name))
    if args.index is True:
        queue_logger.info("Starting to create an index for %s", args.reference.name)
        if os.path.exists("{0}.midx".format(args.reference.name)):
            queue_logger.warning("Removing the old index")
            try:
                os.remove("{0}.midx".format(args.reference.name))
            except (OSError, PermissionError) as exc:
                queue_logger.critical(exc)
                queue_logger.critical(
                    "I cannot delete the old index, due to permission errors. Please investigate and relaunch.")
                sys.exit(1)
        args.protein_coding = False
        args.exclude_utr = False
        create_index(args.reference, queue_logger, index_name, ref_gff=ref_gff,
                     protein_coding=args.protein_coding, exclude_utr=args.exclude_utr)
        queue_logger.info("Finished to create an index for %s", args.reference.name)
    else:
        if os.path.dirname(args.out) and os.path.dirname(args.out) != os.path.dirname(os.path.abspath(".")):
            dirname = os.path.dirname(args.out)
            if os.path.exists(dirname):
                assert os.path.isdir(dirname)
            else:
                os.makedirs(dirname)
        if (os.path.exists(index_name) and os.stat(args.reference.name).st_mtime >= os.stat(index_name).st_mtime):
            queue_logger.warning("Reference index obsolete, deleting and rebuilding.")
            try:
                os.remove("{0}.midx".format(args.reference.name))
            except (OSError, PermissionError) as exc:
                queue_logger.error(
                    "I cannot delete the old index due to permission errors. I will create a temporary one instead.")
                __index = tempfile.NamedTemporaryFile(suffix=".midx")
                index_name = __index.name
            create_index(args.reference, queue_logger, index_name, ref_gff=ref_gff,
                         protein_coding=args.protein_coding, exclude_utr=args.exclude_utr)
        elif os.path.exists(index_name):
            # queue_logger.info("Starting loading the indexed reference")
            queue_logger.info("Index found")
            try:
                check_index(args.reference.name, queue_logger)
                queue_logger.info("Index valid, proceeding.")
            except CorruptIndex as exc:
                queue_logger.warning(exc)
                queue_logger.warning("Reference index corrupt, deleting and rebuilding.")
                try:
                    os.remove("{0}.midx".format(args.reference.name))
                except (OSError, PermissionError) as exc:
                    queue_logger.error(
                        "I cannot delete the old index due to permission errors. I will create a temporary one instead.")
                    __index = tempfile.NamedTemporaryFile(suffix=".midx")
                    index_name = __index.name
                create_index(args.reference, queue_logger, index_name, ref_gff=ref_gff,
                             protein_coding=args.protein_coding, exclude_utr=args.exclude_utr)
        else:
            if args.no_save_index is True:
                __index = tempfile.NamedTemporaryFile(suffix=".midx")
                index_name = __index.name
            create_index(args.reference, queue_logger, index_name, ref_gff=ref_gff,
                         protein_coding=args.protein_coding, exclude_utr=args.exclude_utr)
        try:
            if hasattr(args, "internal") and args.internal is True:
                # raise NotImplementedError()
                parse_self(args, GeneDict(index_name), queue_logger)
            else:
                parse_prediction(args, index_name, queue_logger)
        except Exception as err:
            queue_logger.exception(err)
            log_queue_listener.enqueue_sentinel()
            handler.close()
            log_queue_listener.stop()
            args.queue_handler.close()
            raise

    queue_logger.info("Finished")
    log_queue_listener.enqueue_sentinel()
    log_queue_listener.stop()
    args.queue_handler.close()
    [_.close() for _ in logger.handlers]
    handler.flush()
    handler.close()
    args.reference.close()
    if hasattr(args.prediction, "close"):
        args.prediction.close()

    return
