import os
import tempfile
import gc
from .checking import create_transcript, CheckingProcess
from .annotation_parser import AnnotationParser, load_from_gtf, load_from_gff, load_from_bed12
from ..utilities import Interval, IntervalTree
from ..parsers import to_gff
import operator
import collections
from .. import exceptions
import logging.handlers
import functools
import msgpack
import multiprocessing
import multiprocessing.connection
import multiprocessing.sharedctypes
from collections import defaultdict
import logging
from ..utilities import path_join, merge_partial, overlap
from collections import Counter
import sqlite3
import pysam
import numpy as np
import rapidjson as json


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

    return


def _retrieve_data(shelf_name, shelve_stacks, tid, chrom, key, score, logger,
                   merged_transcripts, chains, monoexonic_tree):
    shelf = shelve_stacks[shelf_name]
    strand, dumped = next(shelf["cursor"].execute("select strand, features from dump where tid = ?", (tid,)))
    dumped = json.loads(dumped)
    try:
        features = dumped["features"]
        exon_set = tuple(sorted([(exon[0], exon[1], strand) for exon in features["exon"]],
                                key=operator.itemgetter(0, 1)))
        if len(exon_set) > 1:
            introns = tuple([(_[0] + 1, _[1] - 1) for _ in zip([_[1] for _ in exon_set][:-1],
                                                         [_[0] for _ in exon_set][1:])])
        else:
            introns = tuple([tuple([exon_set[0][0], exon_set[0][1]])])
        cds_set = tuple(sorted([(exon[0], exon[1]) for exon in features.get("CDS", [])],
                               key=operator.itemgetter(0, 1)))
        monoexonic = not (len(exon_set) > 1)
        data = dict()
        data["introns"], data["strand"], data["score"] = introns, strand, score
        data["monoexonic"] = monoexonic
        data["is_reference"], data["keep_redundant"] = dumped["is_reference"], dumped["keep_redundant"]
        data["start"], data["end"], data["cds_set"] = key[0], key[1], cds_set
        data["key"] = (tuple([tid, shelf_name]), chrom, (data["start"], data["end"]))
        if data["monoexonic"] is True:
            # Additional check at the end because the intervaltree class does not support item removal yet.
            caught = dict((i.value, merged_transcripts[i.value])
                          for i in monoexonic_tree.find(data["start"], data["end"])
                          if i.value in merged_transcripts)
        else:
            caught = dict((i, merged_transcripts[i]) for i in chains.get(data["introns"], []))
        return data, caught
    except (TypeError, IndexError, ValueError, KeyError) as exc:
        logger.error("Error in analysing %s. Skipping. Error: %s", tid, exc)
        return None, None


def _select_transcript(is_reference, score, other, start, end):
    # First off: check they are *actually* redundant
    ovl = overlap((start, end), (other["start"], other["end"]), positive=True)

    if is_reference is True and other["is_reference"] is True:
        return True, True
    length, olength = end - start, other["end"] - other["start"]
    if ovl < length and ovl < olength:
        # Insufficient overlap: the models are staggered!. Continue
        return True, True
    to_keep, other_to_keep = False, False
    # (26612190, 26614205), (26612165, 26614294)
    if other_to_keep is False and (is_reference is False and other["is_reference"] is True) or (score < other["score"]):
        other_to_keep = True
    elif to_keep is False and (is_reference is True and other["is_reference"] is False) or (score > other["score"]):
        to_keep = True
    elif ovl == length == olength:
        to_keep, other_to_keep = np.random.permutation([True, False])
    elif ovl == length:
        other_to_keep = True
    elif ovl == olength:
        to_keep = True
    if other["keep_redundant"] is True:
        other_to_keep = True  # Always keep transcripts marked like this
    if not any([other_to_keep, to_keep]):
        raise ValueError((
            (is_reference, other["is_reference"]),
            (score, other["score"]),
            ((start, end), (other["start"], other["end"])))
        )
    return to_keep, other_to_keep


def _check_correspondence(data: dict, other: dict):
    if data["monoexonic"] is True and other["monoexonic"] is True:
        # raise ValueError((other["introns"][0], key))
        check = ((other["strand"] == data["strand"]) and
                 (other["cds_set"] == data["cds_set"]) and
                 (overlap((other["start"], other["end"]), (data["start"], data["end"])) > 0))
    elif data["monoexonic"] is False and other["monoexonic"] is False:
        check = ((other["strand"] == data["strand"]) and
                 (other["cds_set"] == data["cds_set"]))
    else:
        check = False
    return check


def _analyse_chrom(chrom: str, keys: dict, shelve_stacks: dict, logger, keep_redundant=True):

    merged_transcripts, chains, monoexonic = dict(), defaultdict(set), IntervalTree()
    current = None
    for key in sorted(keys.keys(),
                      key=operator.itemgetter(0, 1)):
        tids = keys[key]
        if keep_redundant is True:
            for tid in tids:
                yield tid[:2], chrom, key
            continue
        if current is not None and overlap(current, key) < 0:
            if not merged_transcripts:
                logger.debug("No transcript kept for %s:%s-%s", chrom, key[0], key[1])
            for tid in merged_transcripts:
                logger.debug("Keeping %s as not redundant", tid)
                yield merged_transcripts[tid]["key"]
            current = key
            merged_transcripts, chains, monoexonic = dict(), defaultdict(set), IntervalTree()
        elif current is not None:
            current = tuple([min(key[0], current[0]), max(key[1], current[1])])
        else:
            current = key

        for tid, shelf_name, score, is_reference, _ in tids:
            to_keep, others_to_remove = True, set()
            data, caught = _retrieve_data(shelf_name, shelve_stacks, tid, chrom, key, score, logger,
                                          merged_transcripts, chains, monoexonic)
            if data is None:
                continue
            to_keep = True
            logger.debug("Checking %s (introns: %s; monoexonic: %s)", tid, data["introns"], data["monoexonic"])
            for otid, other in caught.items():
                check = _check_correspondence(data, other)
                if check:
                    # Redundancy!
                    to_keep, other_to_keep = _select_transcript(
                        is_reference, score, other, start=key[0], end=key[1])
                    if to_keep is True and other_to_keep is True:
                        logger.debug("Keeping both %s and %s.", tid, otid)
                    elif to_keep is True and other_to_keep is False:
                        others_to_remove.add(otid)
                    elif to_keep is False and other_to_keep is True:
                        to_keep = False
                        logger.info("Excluding %s as redundant with %s", tid, otid)
                        break
            if others_to_remove:
                logger.info("Excluding %s as redundant with %s", ",".join(others_to_remove), tid)
                [merged_transcripts.__delitem__(otid) for otid in others_to_remove]
                if data["monoexonic"] is False:
                    [chains[data["introns"]].remove(otid) for otid in others_to_remove]
            if to_keep is True or to_keep is None:
                merged_transcripts[tid] = data
                if data["monoexonic"] is False:
                    chains[data["introns"]].add(tid)
                else:
                    monoexonic.add_interval(Interval(data["start"], data["end"], value=tid))
                logger.debug("Keeping %s in the dataset", tid)

    if not merged_transcripts and current is not None:
        logger.debug("No transcript kept for %s:%s-%s", chrom, current[0], current[1])
    for tid in merged_transcripts:
        logger.debug("Keeping %s as not redundant", tid)
        yield merged_transcripts[tid]["key"]


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
    # Read from the temporary databases the indices (chrom, start, end, strand and transcript ID).
    # Then store
    for shelf_name in shelf_stacks:
        shelf_score = shelf_stacks[shelf_name]["score"]
        is_reference = shelf_stacks[shelf_name]["is_reference"]
        redundants_to_keep = shelf_stacks[shelf_name]["keep_redundant"]
        # redundants_to_keep = shelf_stacks[shelf_name]["keep_redundant"]
        try:
            for values in shelf_stacks[shelf_name]["cursor"].execute("SELECT chrom, start, end, strand, tid FROM dump"):
                chrom, start, end, strand, tid = values
                if (start, end) not in transcripts[chrom]:
                    transcripts[chrom][(start, end)] = list()
                transcripts[chrom][(start, end)].append((tid, shelf_name, shelf_score, is_reference,
                                                         redundants_to_keep))
        except sqlite3.OperationalError as exc:
            raise sqlite3.OperationalError("dump not found in {}; excecption: {}".format(shelf_name, exc))

    np.random.seed(seed)
    for chrom in sorted(transcripts.keys()):
        logger.debug("Starting with %s (%d positions)",
                     chrom,
                     len(transcripts[chrom]))
        yield from _analyse_chrom(chrom, transcripts[chrom], shelve_stacks=shelf_stacks,
                                  logger=logger, keep_redundant=keep_redundant)


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

            if chrom not in args.json_conf["reference"]["genome"].references:
                raise KeyError("Invalid chromosome name! {}, {}, {}, {}".format(tid, shelf_name, chrom, key))

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
    if args.json_conf["prepare"]["files"]["keep_redundant"] == []:
        args.json_conf["prepare"]["files"]["keep_redundant"] = [False] * len(args.json_conf["prepare"]["files"]["gff"])
    if not len(args.json_conf["prepare"]["files"]["keep_redundant"]) == len(args.json_conf["prepare"]["files"]["gff"]):
        raise AssertionError

    to_do = list(zip(
            shelve_names,
            args.json_conf["prepare"]["files"]["labels"],
            args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
            args.json_conf["prepare"]["files"]["reference"],
            args.json_conf["prepare"]["files"]["keep_redundant"],
            args.json_conf["prepare"]["files"]["gff"]
    ))
    logger.info("To do: %d combinations", len(to_do))

    for new_shelf, label, strand_specific, is_reference, keep_redundant, gff_name in to_do:
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
                                    keep_redundant=keep_redundant,
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
                                      keep_redundant=keep_redundant,
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
                                    keep_redundant=keep_redundant,
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

    if args.json_conf["prepare"]["files"]["keep_redundant"] == []:
        args.json_conf["prepare"]["files"]["keep_redundant"] = [False] * len(args.json_conf["prepare"]["files"]["gff"])
    for new_shelf, label, strand_specific, is_reference, keep_redundant, gff_name in zip(
            shelve_names,
            args.json_conf["prepare"]["files"]["labels"],
            args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
            args.json_conf["prepare"]["files"]["reference"],
            args.json_conf["prepare"]["files"]["keep_redundant"],
            args.json_conf["prepare"]["files"]["gff"]):
        submission_queue.put((label, gff_name, strand_specific, is_reference, keep_redundant, new_shelf))

    submission_queue.put(("EXIT", "EXIT", "EXIT", "EXIT", "EXIT", "EXIT"))

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
    errored = False
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
            for shelf, score, is_reference, keep_redundant in zip(
                    shelve_names, shelve_source_scores,
                    args.json_conf["prepare"]["files"]["reference"],
                    args.json_conf["prepare"]["files"]["keep_redundant"]):
                assert isinstance(is_reference, bool)
                conn_string = "file:{shelf}?mode=ro&immutable=1".format(shelf=shelf)
                conn = sqlite3.connect(conn_string, uri=True,
                                       isolation_level="DEFERRED",
                                       check_same_thread=False)
                shelf_stacks[shelf] = {"conn": conn, "cursor": conn.cursor(), "score": score,
                                       "is_reference": is_reference, "keep_redundant": keep_redundant,
                                       "conn_string": conn_string}
            # shelf_stacks = dict((_, shelve.open(_, flag="r")) for _ in shelve_names)
        except Exception as exc:
            raise TypeError((shelve_names, exc))
        perform_check(sorter(shelf_stacks), shelf_stacks, args, logger)
    except Exception as exc:
        logger.exception(exc)
        __cleanup(args, shelve_names)
        errored = True
        logger.error("Mikado has encountered an error, exiting")
        # sys.exit(1)

    if args.json_conf["prepare"]["single"] is False and args.json_conf["threads"] > 1:
        args.tempdir.cleanup()
        args.listener.enqueue_sentinel()

    logger.setLevel(logging.INFO)
    __cleanup(args, shelve_names)

    logger.addHandler(logging.StreamHandler())
    if errored is False:
        logger.info("""Mikado prepare has finished correctly. The output %s FASTA file can now be used for BLASTX \
    and/or ORF calling before the next step in the pipeline, `mikado serialise`.""",
                    args.json_conf["prepare"]["files"]["out_fasta"])
    else:
        logger.error("Mikado prepared has encountered a fatal error. Please check the logs and, if there is a bug,"\
                     "report it to https://github.com/EI-CoreBioinformatics/mikado/issues")
    logging.shutdown()
