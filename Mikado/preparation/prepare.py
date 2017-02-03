import os
import tempfile
import gc
from .checking import create_transcript, CheckingProcess
from .annotation_parser import AnnotationParser, load_from_gtf, load_from_gff
from ..parsers import to_gff
import operator
import collections
import io
from .. import exceptions
import logging.handlers
import random
import functools
import multiprocessing
import multiprocessing.connection
import multiprocessing.sharedctypes
import pyfaidx
import logging
from ..utilities import path_join, merge_partial
from collections import Counter
import sqlite3
try:
    import ujson as json
except ImportError:
    import json

__author__ = 'Luca Venturini'


def __cleanup(args, shelves):
    """Private function to close the opened handles."""

    if hasattr(args.json_conf["reference"]["genome"], "close"):
        args.json_conf["reference"]["genome"].close()
        args.json_conf["reference"]["genome"] = args.json_conf["reference"]["genome"].filename

    for frole in ("out", "out_fasta"):
        if hasattr(args.json_conf["prepare"]["files"][frole], "close"):
            args.json_conf["prepare"]["files"][frole].close()
            args.json_conf["prepare"]["files"][frole] = args.json_conf["prepare"]["files"][frole].name

    [os.remove(fname) for fname in shelves if os.path.exists(fname)]


def store_transcripts(shelf_stacks, logger, keep_redundant=False):

    """
    Function that analyses the exon lines from the original file
    and organises the data into a proper dictionary.
    :param shelf_stacks: dictionary containing the name and the handles of the shelf DBs
    :type shelf_stacks: dict

    :param logger: logger instance.
    :type logger: logging.Logger

    :param keep_redundant: boolean flag. If set to True, redundant transcripts will be kept in the output.
    :type keep_redundant: bool

    :return: transcripts: dictionary which will be the final output
    :rtype: transcripts
    """

    transcripts = collections.defaultdict(dict)

    for shelf_name in shelf_stacks:

        for values in shelf_stacks[shelf_name]["cursor"].execute("SELECT chrom, start, end, strand, tid FROM dump"):
            chrom, start, end, strand, tid = values
            if (start, end) not in transcripts[chrom]:
                transcripts[chrom][(start, end)] = []
            transcripts[chrom][(start, end)].append((tid, shelf_name))

    for chrom in sorted(transcripts.keys()):

        logger.debug("Starting with %s (%d positions)",
                     chrom,
                     len(transcripts[chrom]))

        for key in sorted(transcripts[chrom].keys(),
                          key=operator.itemgetter(0, 1)):
            tids = transcripts[chrom][key]
            if len(tids) > 1:
                exons = collections.defaultdict(list)
                for tid, shelf in tids:
                    strand, features = next(shelf_stacks[shelf]["cursor"].execute(
                        "select strand, features from dump where tid = ?", (tid,)))
                    features = json.loads(features)
                    exon_set = tuple(sorted([(exon[0], exon[1], strand) for exon in
                                            features["features"]["exon"]],
                                            key=operator.itemgetter(0, 1)))

                    # exon_set = tuple(sorted(
                    #     [(exon[0], exon[1], shelf_stacks[shelf][tid]["strand"]) for exon in
                    #      shelf_stacks[shelf][tid]["features"]["exon"]],

                    exons[exon_set].append((tid, shelf))
                tids = []
                logger.debug("%d intron chains for pos %s",
                             len(exons), "{}:{}-{}".format(chrom, key[0], key[1]))
                for tid_list in exons.values():
                    if len(tid_list) > 1 and keep_redundant is False:
                        logger.debug("The following transcripts are redundant: %s",
                                     ",".join([_[0] for _ in tid_list]))
                        to_keep = random.choice(tid_list)
                        logger.debug("Keeping only %s out of the list",
                                     to_keep)
                        tids.append(to_keep)
                    else:
                        tids.extend(tid_list)

            # seq = chrom_seq[key[0]-1:key[1]]
            for tid in tids:
                # counter += 1
                yield [tid, chrom, key]
            # keys.extend([tid, seq] for tid in tids)

    # logger.info("Finished to sort %d transcripts, %d remain", len(exon_lines), counter)

    # return keys


def perform_check(keys, shelve_stacks, args, logger):

    """
    This is the most important method. After preparing the data structure,
    this function creates the real transcript instances and checks that
    they are correct when looking at the underlying genome sequence.
    This is also the point at which we start using multithreading, if
    so requested.
    :param keys: sorted list of [tid, sequence]
    :param shelve_stacks: dictionary containing the name and the handles of the shelf DBs
    :param args: the namespace
    :param logger: logger
    :return:
    """

    counter = 0

    # FASTA extraction *has* to be done at the main process level, it's too slow
    # to create an index in each process.

    if args.json_conf["prepare"]["single"] is True or args.json_conf["prepare"]["procs"] == 1:

        # Use functools to pre-configure the function
        # with all necessary arguments aside for the lines
        partial_checker = functools.partial(
            create_transcript,
            lenient=args.json_conf["prepare"]["lenient"],
            # strand_specific=args.json_conf["prepare"]["strand_specific"],
            canonical_splices=args.json_conf["prepare"]["canonical"],
            logger=logger)

        for tid, chrom, key in keys:
            tid, shelf_name = tid
            try:
                tobj = json.loads(next(shelve_stacks[shelf_name]["cursor"].execute(
                    "SELECT features FROM dump WHERE tid = ?", (tid,)))[0])
            except sqlite3.ProgrammingError as exc:
                raise sqlite3.ProgrammingError("{}. Tids: {}".format(exc, tid))

            transcript_object = partial_checker(
                tobj,
                str(args.json_conf["reference"]["genome"][chrom][key[0]-1:key[1]]),
                key[0], key[1],
                strand_specific=tobj["strand_specific"])
            if transcript_object is None:
                continue
            counter += 1
            if counter >= 10**4 and counter % (10**4) == 0:
                logger.info("Retrieved %d transcript positions", counter)
            elif counter >= 10**3 and counter % (10**3) == 0:
                logger.debug("Retrieved %d transcript positions", counter)
            print(transcript_object.format("gtf"),
                  file=args.json_conf["prepare"]["files"]["out"])
            print(transcript_object.fasta,
                  file=args.json_conf["prepare"]["files"]["out_fasta"])
    else:
        # pylint: disable=no-member

        submission_queue = multiprocessing.Queue(-1)

        working_processes = [CheckingProcess(
            submission_queue,
            args.logging_queue,
            args.json_conf["reference"]["genome"].filename,
            _ + 1,
            os.path.basename(args.json_conf["prepare"]["files"]["out_fasta"].name),
            os.path.basename(args.json_conf["prepare"]["files"]["out"].name),
            args.tempdir.name,
            lenient=args.json_conf["prepare"]["lenient"],
            # strand_specific=args.json_conf["prepare"]["strand_specific"],
            canonical_splices=args.json_conf["prepare"]["canonical"],
            log_level=args.level) for _ in range(args.json_conf["prepare"]["procs"])]

        [_.start() for _ in working_processes]

        for counter, keys in enumerate(keys):
            tid, chrom, (pos) = keys
            tid, shelf_name = tid
            tobj = json.loads(next(shelve_stacks[shelf_name]["cursor"].execute(
                "SELECT features FROM dump WHERE tid = ?", (tid,)))[0])
            submission_queue.put((tobj, pos[0], pos[1], counter + 1))

        submission_queue.put(tuple(["EXIT"]*4))

        [_.join() for _ in working_processes]

        partial_gtf = [os.path.join(args.tempdir.name,
                                    "{0}-{1}".format(
                                        os.path.basename(args.json_conf["prepare"]["files"]["out"].name),
                                        _ + 1)) for _ in range(args.json_conf["prepare"]["procs"])]
        merge_partial(partial_gtf, args.json_conf["prepare"]["files"]["out"])

        partial_fasta = [os.path.join(
            args.tempdir.name,
            "{0}-{1}".format(os.path.basename(args.json_conf["prepare"]["files"]["out_fasta"].name), _ + 1))
                         for _ in range(args.json_conf["prepare"]["procs"])]
        merge_partial(partial_fasta, args.json_conf["prepare"]["files"]["out_fasta"])

    args.json_conf["prepare"]["files"]["out_fasta"].close()
    args.json_conf["prepare"]["files"]["out"].close()

    logger.setLevel(logging.INFO)
    # logger.info("Finished to analyse %d transcripts (%d retained)",
    #             len(exon_lines), counter)
    logger.setLevel(args.level)
    return


def load_exon_lines(args, shelve_names, logger, min_length=0):

    """This function loads all exon lines from the GFF inputs into a
     defaultdict instance.
    :param args: the Namespace from the command line.
    :param shelve_names: list of names of the shelf DB files.
    :param logger: the logger instance.
    :type logger: logging.Logger
    :param min_length: minimal length of the transcript.
    If it is not met, the transcript will be discarded.
    :type min_length: int

    :return: exon_lines
    :rtype: collections.defaultdict[list]
    """

    threads = min([len(args.json_conf["prepare"]["files"]["gff"]),
                   args.json_conf["prepare"]["procs"]])
    strip_cds = args.json_conf["prepare"]["strip_cds"]

    if args.json_conf["prepare"]["single"] is True or threads == 1:

        logger.info("Starting to load lines from %d files (single-threaded)",
                    len(args.json_conf["prepare"]["files"]["gff"]))
        previous_file_ids = collections.defaultdict(set)
        for new_shelf, label, strand_specific, gff_name in zip(
                shelve_names,
                args.json_conf["prepare"]["files"]["labels"],
                args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
                args.json_conf["prepare"]["files"]["gff"]):
            logger.info("Starting with %s", gff_name)
            gff_handle = to_gff(gff_name)
            found_ids = set.union(set(), *previous_file_ids.values())
            if gff_handle.__annot_type__ == "gff3":
                new_ids = load_from_gff(new_shelf,
                                        gff_handle,
                                        label,
                                        found_ids,
                                        logger,
                                        min_length=min_length,
                                        strip_cds=strip_cds,
                                        strand_specific=strand_specific)
            else:
                new_ids = load_from_gtf(new_shelf,
                                        gff_handle,
                                        label,
                                        found_ids,
                                        logger,
                                        min_length=min_length,
                                        strip_cds=strip_cds,
                                        strand_specific=strand_specific)

            previous_file_ids[gff_handle.name] = new_ids
    else:
        logger.info("Starting to load lines from %d files (using %d processes)",
                    len(args.json_conf["prepare"]["files"]["gff"]), threads)
        submission_queue = multiprocessing.Queue(-1)

        working_processes = [AnnotationParser(
            submission_queue,
            args.logging_queue,
            _ + 1,
            log_level=args.level,
            min_length=min_length,
            strip_cds=strip_cds) for _ in range(threads)]

        [_.start() for _ in working_processes]
        for new_shelf, label, strand_specific, gff_name in zip(
                shelve_names,
                args.json_conf["prepare"]["files"]["labels"],
                args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
                args.json_conf["prepare"]["files"]["gff"]):

            submission_queue.put((label, gff_name, strand_specific, new_shelf))

        submission_queue.put(("EXIT", "EXIT", "EXIT", "EXIT"))

        [_.join() for _ in working_processes]

        tid_counter = Counter()
        for shelf in shelve_names:
            conn = sqlite3.connect(shelf)
            cursor = conn.cursor()
            tid_counter.update([_[0] for _ in cursor.execute("SELECT tid FROM dump")])
            if tid_counter.most_common()[0][1] > 1:
                if set(args.json_conf["prepare"]["files"]["labels"]) == {""}:
                    exception = exceptions.RedundantNames(
                        """Found redundant names during multiprocessed file analysis.
                        Please repeat using distinct labels for your input files. Aborting.""")
                else:
                    exception = exceptions.RedundantNames(
                        """Found redundant names during multiprocessed file analysis, even if
                        unique labels were provided. Please try to repeat with a different and
                        more unique set of labels. Aborting.""")
                logger.exception(exception)
                raise exception

        del working_processes
        gc.collect()

    logger.info("Finished loading lines from %d files",
                len(args.json_conf["prepare"]["files"]["gff"]))

    return


def prepare(args, logger):
    """Main script function.

    :param args: the ArgumentParser-derived namespace.
    :param logger: a logging instance
    :type logger: logging.Logger
    """

    if not isinstance(args.json_conf["reference"]["genome"], (io.TextIOWrapper, pyfaidx.Fasta)):
        if not (isinstance(args.json_conf["reference"]["genome"], str) and
                os.path.exists(args.json_conf["reference"]["genome"])):
            logger.critical("Invalid FASTA file: %s",
                            args.json_conf["reference"]["genome"])
            # sys.exit(1)
        else:
            pass
    else:
        args.json_conf["reference"]["genome"].close()
        if isinstance(args.json_conf["reference"]["genome"], pyfaidx.Fasta):
            args.json_conf["reference"]["genome"] = args.json_conf["reference"]["genome"].filename
        else:
            args.json_conf["reference"]["genome"] = args.json_conf["reference"]["genome"].name

    assert len(args.json_conf["prepare"]["files"]["gff"]) > 0
    assert len(args.json_conf["prepare"]["files"]["gff"]) == len(args.json_conf["prepare"]["files"]["labels"]), (
        args.json_conf["prepare"]["files"]["gff"],
        args.json_conf["prepare"]["files"]["labels"]
    )

    if args.json_conf["prepare"]["strand_specific"] is True:
        args.json_conf["prepare"]["files"]["strand_specific_assemblies"] = [True] * len(
            args.json_conf["prepare"]["files"]["gff"])
    else:
        args.json_conf["prepare"]["files"]["strand_specific_assemblies"] = [
            (member in args.json_conf["prepare"]["files"]["strand_specific_assemblies"])
            for member in args.json_conf["prepare"]["files"]["gff"]]

    shelve_names = [path_join(args.json_conf["prepare"]["files"]["output_dir"],
                              "mikado_shelf_{}.db".format(str(_).zfill(5))) for _ in
                    range(len(args.json_conf["prepare"]["files"]["gff"]))]

    logger.propagate = False
    if args.json_conf["prepare"]["single"] is False and args.json_conf["prepare"]["procs"] > 1:
        multiprocessing.set_start_method(args.json_conf["multiprocessing_method"],
                                         force=True)
        args.logging_queue = multiprocessing.Queue(-1)
        log_queue_handler = logging.handlers.QueueHandler(args.logging_queue)
        log_queue_handler.setLevel(logging.DEBUG)
        # logger.addHandler(log_queue_handler)
        args.tempdir = tempfile.TemporaryDirectory(dir=args.json_conf["prepare"]["files"]["output_dir"])
        args.listener = logging.handlers.QueueListener(args.logging_queue, logger)
        args.listener.propagate = False
        args.listener.start()

    args.json_conf["prepare"]["files"]["out_fasta"] = open(
        path_join(args.json_conf["prepare"]["files"]["output_dir"],
                  args.json_conf["prepare"]["files"]["out_fasta"]), 'w')
    args.json_conf["prepare"]["files"]["out"] = open(path_join(
        args.json_conf["prepare"]["files"]["output_dir"],
        args.json_conf["prepare"]["files"]["out"]), 'w')


    logger.info("Loading reference file")
    args.json_conf["reference"]["genome"] = pyfaidx.Fasta(args.json_conf["reference"]["genome"])

    logger.info("Finished loading genome file")
    logger.info("Started loading exon lines")

    shelf_stacks = dict()
    try:
        load_exon_lines(args,
                        shelve_names,
                        logger,
                        min_length=args.json_conf["prepare"]["minimum_length"])

        logger.info("Finished loading exon lines")

        # Prepare the sorted data structure
        sorter = functools.partial(
            store_transcripts,
            logger=logger,
            keep_redundant=args.json_conf["prepare"]["keep_redundant"]
            # min_length=args.json_conf["prepare"]["minimum_length"]
        )

        try:
            for shelf in shelve_names:
                conn = sqlite3.connect(shelf)
                shelf_stacks[shelf] = {"conn": conn, "cursor": conn.cursor()}
            # shelf_stacks = dict((_, shelve.open(_, flag="r")) for _ in shelve_names)
        except Exception as exc:
            raise TypeError((shelve_names, exc))
        perform_check(sorter(shelf_stacks), shelf_stacks, args, logger)
    except Exception as exc:
        logger.exception(exc)
        __cleanup(args, shelve_names)
        logger.error("Mikado has encountered an error, exiting")
        # sys.exit(1)

    if args.json_conf["prepare"]["single"] is False and args.json_conf["prepare"]["procs"] > 1:
        args.tempdir.cleanup()
        args.listener.enqueue_sentinel()

    logger.setLevel(logging.INFO)
    __cleanup(args, shelve_names)

    logger.info("Finished")
    logging.shutdown()
    # sys.exit(0)
    # for handler in logger.handlers:
    #     handler.close()
