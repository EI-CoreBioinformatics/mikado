import os
import tempfile
import gc
from .checking import create_transcript, CheckingProcess
from .annotation_parser import AnnotationParser, load_from_gtf, load_from_gff, load_from_bed12
from ..parsers import to_gff
import operator
import collections
from .. import exceptions
import logging.handlers
import numpy
import functools
import msgpack
import multiprocessing
import multiprocessing.connection
import multiprocessing.sharedctypes
import logging
from ..utilities import path_join, merge_partial
from collections import Counter
import sqlite3
import pysam
import numpy as np
try:
    import rapidjson as json
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

    for fname in shelves:
        [os.remove(fname + suff) for suff in ("", "-shm", "-wal", "-journal") if os.path.exists(fname + suff)]


def store_transcripts(shelf_stacks, logger, keep_redundant=False, seed=None):

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
        shelf_score = shelf_stacks[shelf_name]["score"]
        is_reference = shelf_stacks[shelf_name]["is_reference"]
        try:
            for values in shelf_stacks[shelf_name]["cursor"].execute("SELECT chrom, start, end, strand, tid FROM dump"):
                chrom, start, end, strand, tid = values
                if (start, end) not in transcripts[chrom]:
                    transcripts[chrom][(start, end)] = list()
                transcripts[chrom][(start, end)].append((tid, shelf_name, shelf_score, is_reference))
        except sqlite3.OperationalError:
            raise sqlite3.OperationalError("dump not found in {}".format(shelf_name))

    for chrom in sorted(transcripts.keys()):

        logger.debug("Starting with %s (%d positions)",
                     chrom,
                     len(transcripts[chrom]))

        for key in sorted(transcripts[chrom].keys(),
                          key=operator.itemgetter(0, 1)):
            tids = transcripts[chrom][key]
            if len(tids) > 1:
                exons = collections.defaultdict(list)
                di_features = dict()
                for tid, shelf, score, is_reference in tids:
                    strand, features = next(shelf_stacks[shelf]["cursor"].execute(
                        "select strand, features from dump where tid = ?", (tid,)))
                    features = json.loads(features)
                    try:
                        exon_set = tuple(sorted([(exon[0], exon[1], strand) for exon in
                                                features["features"]["exon"]],
                                                key=operator.itemgetter(0, 1)))
                    except (TypeError, IndexError, ValueError):
                        logger.error("Error in analysing %s. Skipping",
                                     tid)
                        continue
                    exons[exon_set].append((tid, shelf, score, is_reference))
                    di_features[tid] = features
                tids = []
                if len(exons) == 0:
                    continue
                logger.debug("%d exon chains for pos %s",
                             len(exons), "{}:{}-{}".format(chrom, key[0], key[1]))
                for tid_list in exons.values():
                    if len(tid_list) > 1 and keep_redundant is False:
                        cds = collections.defaultdict(list)
                        for tid, shelf, score, is_reference in tid_list:
                            cds_set = tuple(sorted([(exon[0], exon[1]) for exon in
                                                    di_features[tid]["features"].get("CDS", [])]))
                            cds[cds_set].append((tid, shelf, score, is_reference))
                        # Now checking the CDS
                        for cds_list in cds.values():
                            if len(cds_list) > 1:
                                logger.debug("The following transcripts are redundant: %s",
                                             ",".join([_[0] for _ in cds_list]))
                            # TODO: we have to select things so that we prioritise transcripts correctly
                            # First step: select things that are reference
                            to_keep = [_ for _ in cds_list if _[3] is True]
                            if len(to_keep) == 0:
                                max_score = max([_[2] for _ in cds_list])
                                cds_list = sorted([_ for _ in cds_list if _[2] == max_score])
                                np.random.seed(seed)
                                to_keep = [cds_list[numpy.random.choice(len(cds_list))]]

                            for tid in cds_list:
                                if tid not in to_keep:
                                    logger.info("Discarding %s as redundant", tid[0])
                                else:
                                    logger.info("Keeping %s amongst redundant transcripts", tid[0])
                            tids.extend(to_keep)
                    else:
                        tids.extend(tid_list)

            for tid in tids:
                # counter += 1
                tid = tid[:2]
                assert len(tid) == 2, tid
                yield [tid, chrom, key]


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

    if args.json_conf["prepare"]["single"] is True or args.json_conf["threads"] == 1:

        # Use functools to pre-configure the function
        # with all necessary arguments aside for the lines
        partial_checker = functools.partial(
            create_transcript,
            canonical_splices=args.json_conf["prepare"]["canonical"],
            logger=logger,
            force_keep_cds=not args.json_conf["prepare"]["strip_cds"])

        for tid, chrom, key in keys:
            tid, shelf_name = tid
            try:
                tobj = json.loads(next(shelve_stacks[shelf_name]["cursor"].execute(
                    "SELECT features FROM dump WHERE tid = ?", (tid,)))[0])
            except sqlite3.ProgrammingError as exc:
                raise sqlite3.ProgrammingError("{}. Tids: {}".format(exc, tid))

            transcript_object = partial_checker(
                tobj,
                str(args.json_conf["reference"]["genome"].fetch(chrom, key[0]-1, key[1])),
                key[0], key[1],
                lenient=args.json_conf["prepare"]["lenient"],
                is_reference=tobj["is_reference"],
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

        # submission_queue = multiprocessing.JoinableQueue(-1)

        batches = list(enumerate(keys, 1))
        np.random.shuffle(batches)
        kwargs = {
            "fasta_out": os.path.basename(args.json_conf["prepare"]["files"]["out_fasta"].name),
            "gtf_out": os.path.basename(args.json_conf["prepare"]["files"]["out"].name),
            "tmpdir": args.tempdir.name,
            "seed": args.json_conf["seed"],
            "lenient": args.json_conf["prepare"]["lenient"],
            "canonical_splices": args.json_conf["prepare"]["canonical"],
            "force_keep_cds": not args.json_conf["prepare"]["strip_cds"],
            "log_level": args.level
        }

        working_processes = []
        for idx, batch in enumerate(np.array_split(batches, args.json_conf["threads"]), 1):
            batch_file = tempfile.NamedTemporaryFile(delete=False, mode="wb")
            msgpack.dump(batch.tolist(), batch_file)
            batch_file.flush()
            batch_file.close()

            proc = CheckingProcess(
                batch_file.name,
                args.logging_queue,
                args.json_conf["reference"]["genome"].filename,
                idx,
                shelve_stacks.keys(),
                **kwargs)
            proc.start()
            working_processes.append(proc)

        [_.join() for _ in working_processes]

        partial_gtf = [os.path.join(args.tempdir.name,
                                    "{0}-{1}".format(
                                        os.path.basename(args.json_conf["prepare"]["files"]["out"].name),
                                        _ + 1)) for _ in range(args.json_conf["threads"])]
        merge_partial(partial_gtf, args.json_conf["prepare"]["files"]["out"])

        partial_fasta = [os.path.join(
            args.tempdir.name,
            "{0}-{1}".format(os.path.basename(args.json_conf["prepare"]["files"]["out_fasta"].name), _ + 1))
                         for _ in range(args.json_conf["threads"])]
        merge_partial(partial_fasta, args.json_conf["prepare"]["files"]["out_fasta"])

    args.json_conf["prepare"]["files"]["out_fasta"].close()
    args.json_conf["prepare"]["files"]["out"].close()

    logger.setLevel(logging.INFO)
    # logger.info("Finished to analyse %d transcripts (%d retained)",
    #             len(exon_lines), counter)
    logger.setLevel(args.level)
    return


def _load_exon_lines_single_thread(args, shelve_names, logger, min_length, strip_cds, max_intron):

    logger.info("Starting to load lines from %d files (single-threaded)",
                len(args.json_conf["prepare"]["files"]["gff"]))
    previous_file_ids = collections.defaultdict(set)
    to_do = list(zip(
            shelve_names,
            args.json_conf["prepare"]["files"]["labels"],
            args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
            args.json_conf["prepare"]["files"]["reference"],
            args.json_conf["prepare"]["files"]["gff"]))
    logger.info("To do: %d combinations", len(to_do))

    for new_shelf, label, strand_specific, is_reference, gff_name in to_do:
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
                                    max_intron=max_intron,
                                    strip_cds=strip_cds and not is_reference,
                                    is_reference=is_reference,
                                    strand_specific=strand_specific or is_reference)
        elif gff_handle.__annot_type__ == "bed12":
            new_ids = load_from_bed12(new_shelf,
                                      gff_handle,
                                      label,
                                      found_ids,
                                      logger,
                                      min_length=min_length,
                                      max_intron=max_intron,
                                      strip_cds=strip_cds and not is_reference,
                                      is_reference=is_reference,
                                      strand_specific=strand_specific or is_reference)
        else:
            new_ids = load_from_gtf(new_shelf,
                                    gff_handle,
                                    label,
                                    found_ids,
                                    logger,
                                    min_length=min_length,
                                    max_intron=max_intron,
                                    strip_cds=strip_cds and not is_reference,
                                    is_reference=is_reference,
                                    strand_specific=strand_specific or is_reference)

        previous_file_ids[gff_handle.name] = new_ids
    return


def _load_exon_lines_multi(args, shelve_names, logger, min_length, strip_cds, threads, max_intron=3*10**5):
    logger.info("Starting to load lines from %d files (using %d processes)",
                len(args.json_conf["prepare"]["files"]["gff"]), threads)
    submission_queue = multiprocessing.JoinableQueue(-1)

    working_processes = []
    # working_processes = [ for _ in range(threads)]

    for num in range(threads):
        proc = AnnotationParser(submission_queue,
                                args.logging_queue,
                                num + 1,
                                log_level=args.level,
                                min_length=min_length,
                                max_intron=max_intron,
                                strip_cds=strip_cds,
                                seed=args.json_conf["seed"])
        proc.start()
        working_processes.append(proc)

    # [_.start() for _ in working_processes]
    for new_shelf, label, strand_specific, is_reference, gff_name in zip(
            shelve_names,
            args.json_conf["prepare"]["files"]["labels"],
            args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
            args.json_conf["prepare"]["files"]["reference"],
            args.json_conf["prepare"]["files"]["gff"]):
        submission_queue.put((label, gff_name, strand_specific, is_reference, new_shelf))

    submission_queue.put(("EXIT", "EXIT", "EXIT", "EXIT", "EXIT"))

    [_.join() for _ in working_processes]

    tid_counter = Counter()
    for shelf in shelve_names:
        conn = sqlite3.connect("file:{}?mode=ro".format(shelf),
                               uri=True,  # Necessary to use the Read-only mode from file string
                               isolation_level="DEFERRED",
                               timeout=60,
                               check_same_thread=False  # Necessary for SQLite3 to function in multiprocessing
                               )
        cursor = conn.cursor()
        tid_counter.update([_[0] for _ in cursor.execute("SELECT tid FROM dump")])
        if tid_counter.most_common()[0][1] > 1:
            if set(args.json_conf["prepare"]["files"]["labels"]) == {""}:
                exception = exceptions.RedundantNames(
                    """Found redundant names during multiprocessed file analysis.\
Please repeat using distinct labels for your input files. Aborting. Redundant names:\n\
{}""".format("\n".join(tid_counter.most_common())))
            else:
                exception = exceptions.RedundantNames(
                    """Found redundant names during multiprocessed file analysis, even if \
unique labels were provided. Please try to repeat with a different and more unique set of labels. Aborting. Redundant names:\n\
{}""".format("\n".join([_[0] for _ in tid_counter.most_common() if _[1] > 1])))
            logger.exception(exception)
            raise exception

    del working_processes
    gc.collect()


def load_exon_lines(args, shelve_names, logger, min_length=0, max_intron=3*10**5):

    """This function loads all exon lines from the GFF inputs into a
     defaultdict instance.
    :param args: the Namespace from the command line.
    :param shelve_names: list of names of the shelf DB files.
    :param logger: the logger instance.
    :type logger: logging.Logger
    :param min_length: minimal length of the transcript. If it is not met, the transcript will be discarded.
    :param max_intron: maximum length for an intron. If it is not met, the transcript will be discarded.
    :type min_length: int
f
    :return: exon_lines
    :rtype: collections.defaultdict[list]
    """

    threads = min([len(args.json_conf["prepare"]["files"]["gff"]),
                   args.json_conf["threads"]])
    strip_cds = args.json_conf["prepare"]["strip_cds"]

    if args.json_conf["prepare"]["single"] is True or threads == 1:
        _load_exon_lines_single_thread(args, shelve_names, logger, min_length, strip_cds, max_intron)
    else:
        _load_exon_lines_multi(args, shelve_names, logger, min_length, strip_cds, threads, max_intron)

    logger.info("Finished loading lines from %d files",
                len(args.json_conf["prepare"]["files"]["gff"]))

    return


def prepare(args, logger):
    """Main script function.

    :param args: the ArgumentParser-derived namespace.
    :param logger: a logging instance
    :type logger: logging.Logger
    """

    if hasattr(args.json_conf["reference"]["genome"], "close"):
        args.json_conf["reference"]["genome"].close()
        if hasattr(args.json_conf["reference"]["genome"], "filename"):
            args.json_conf["reference"]["genome"] = getattr(args.json_conf["reference"]["genome"], "filename")
        elif hasattr(args.json_conf["reference"]["genome"], "name"):
            args.json_conf["reference"]["genome"] = getattr(args.json_conf["reference"]["genome"], "name")
        else:
            logger.critical("Invalid FASTA file: %s",
                            args.json_conf["reference"]["genome"])
            raise AttributeError
    elif not isinstance(args.json_conf["reference"]["genome"], (str, bytes)):
        logger.critical("Invalid FASTA file: %s",
                        args.json_conf["reference"]["genome"])
        raise AttributeError

    if not os.path.exists(args.json_conf["reference"]["genome"]):
        logger.critical("Invalid FASTA file: %s",
                        args.json_conf["reference"]["genome"])
        raise AttributeError

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

    args.json_conf["prepare"]["files"]["reference"] = [
        (member in args.json_conf["prepare"]["files"]["reference"] or
         label in args.json_conf["prepare"]["files"]["reference"])
        for member, label in zip(args.json_conf["prepare"]["files"]["gff"],
                                 args.json_conf["prepare"]["files"]["labels"])
    ]

    shelve_names = [path_join(args.json_conf["prepare"]["files"]["output_dir"],
                              "mikado_shelf_{}.db".format(str(_).zfill(5))) for _ in
                    range(len(args.json_conf["prepare"]["files"]["gff"]))]

    logger.propagate = False
    if args.json_conf["prepare"]["single"] is False and args.json_conf["threads"] > 1:
        multiprocessing.set_start_method(args.json_conf["multiprocessing_method"],
                                         force=True)
        args.logging_queue = multiprocessing.JoinableQueue(-1)
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

    logger.info("Output dir: %s. Output GTF: %s. Output Fasta: %s",
                args.json_conf["prepare"]["files"]["output_dir"],
                args.json_conf["prepare"]["files"]["out"].name,
                args.json_conf["prepare"]["files"]["out_fasta"].name)
    logger.info("Loading reference file")
    args.json_conf["reference"]["genome"] = pysam.FastaFile(args.json_conf["reference"]["genome"])
    logger.info("Finished loading genome file")
    logger.info("Started loading exon lines")
    shelf_stacks = dict()
    try:
        load_exon_lines(args,
                        shelve_names,
                        logger,
                        min_length=args.json_conf["prepare"]["minimum_cdna_length"],
                        max_intron=args.json_conf["prepare"]["max_intron_length"],)

        logger.info("Finished loading exon lines")

        # Prepare the sorted data structure
        sorter = functools.partial(
            store_transcripts,
            logger=logger,
            seed=args.json_conf["seed"],
            keep_redundant=args.json_conf["prepare"]["keep_redundant"]
        )

        shelve_source_scores = []
        for label in args.json_conf["prepare"]["files"]["labels"]:
            shelve_source_scores.append(
                args.json_conf["prepare"]["files"]["source_score"].get(label, 0)
            )

        try:
            for shelf, score, is_reference in zip(shelve_names, shelve_source_scores,
                                    args.json_conf["prepare"]["files"]["reference"]):
                assert isinstance(is_reference, bool)
                conn_string = "file:{shelf}?mode=ro&immutable=1".format(shelf=shelf)
                conn = sqlite3.connect(conn_string, uri=True,
                                       isolation_level="DEFERRED",
                                       check_same_thread=False)
                shelf_stacks[shelf] = {"conn": conn, "cursor": conn.cursor(), "score": score,
                                       "conn_string": conn_string,
                                       "is_reference": is_reference}
            # shelf_stacks = dict((_, shelve.open(_, flag="r")) for _ in shelve_names)
        except Exception as exc:
            raise TypeError((shelve_names, exc))
        perform_check(sorter(shelf_stacks), shelf_stacks, args, logger)
    except Exception as exc:
        logger.exception(exc)
        __cleanup(args, shelve_names)
        logger.error("Mikado has encountered an error, exiting")
        # sys.exit(1)

    if args.json_conf["prepare"]["single"] is False and args.json_conf["threads"] > 1:
        args.tempdir.cleanup()
        args.listener.enqueue_sentinel()

    logger.setLevel(logging.INFO)
    __cleanup(args, shelve_names)

    logger.addHandler(logging.StreamHandler())
    logger.info("""Mikado prepare has finished correctly. The output %s FASTA file can now be used for BLASTX \
and/or ORF calling before the next step in the pipeline, `mikado serialise`.""",
                args.json_conf["prepare"]["files"]["out_fasta"])
    logging.shutdown()
