import os
import tempfile
import gc
from .checking import create_transcript, CheckingProcess
from .annotation_parser import AnnotationParser, loaders, row_struct
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
import sqlite3
import pysam
import numpy as np
import random
import pandas as pd
import zlib


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


def _retrieve_data(shelf, shelf_name, tid, chrom, key, strand, score, write_start, write_length, logger,
                   merged_transcripts, chains):
    shelf.seek(write_start)
    dumped = shelf.read(write_length)
    dumped = msgpack.loads(zlib.decompress(dumped), raw=False)
    assert isinstance(dumped, dict), dumped
    try:
        features = dumped["features"]
        exon_set = tuple(sorted([(exon[0], exon[1], strand) for exon in features["exon"]],
                                key=operator.itemgetter(0, 1)))
        if len(exon_set) > 1:
            introns = tuple([(_[0] + 1, _[1] - 1) for _ in zip([_[1] for _ in exon_set][:-1],
                                                         [_[0] for _ in exon_set][1:])])
        else:
            introns = None
        cds_set = tuple(sorted([(exon[0], exon[1]) for exon in features.get("CDS", [])],
                               key=operator.itemgetter(0, 1)))
        data = dict()
        data["introns"], data["strand"], data["score"] = introns, strand, score
        data["monoexonic"] = (len(exon_set) == 1)
        data["is_reference"], data["exclude_redundant"] = dumped["is_reference"], dumped["exclude_redundant"]
        data["start"], data["end"], data["cds_set"] = key[0], key[1], cds_set
        data["key"] = (tuple([tid, shelf_name, write_start, write_length]),
                       chrom, (int(data["start"]), int(data["end"])))
        # Sorted by default
        try:
            caught = [(i.value, merged_transcripts[i.value])
                      for i in chains[data["introns"]].find(data["start"], data["end"], contained_check=True,
                                                            num_intervals=10**5)
                      if i.value in merged_transcripts]
        except AttributeError as exc:
            raise AttributeError("{}\n{}\n{}".format(exc, type(chains), type(chains[data["introns"]])))

        return data, caught
    except (TypeError, IndexError, ValueError, KeyError) as exc:
        logger.error("Error in analysing %s. Skipping. Error: %s", tid, exc)
        return None, None


def _select_transcript(is_reference, exclude_redundant, score, other, start, end):
    # First off: check they are *actually* redundant
    ovl = overlap((start, end), (other["start"], other["end"]), positive=True)

    if is_reference is True and other["is_reference"] is True:
        return True, True
    length, olength = end - start, other["end"] - other["start"]
    if ovl < min(length, olength):
        # Insufficient overlap: the models are staggered!. Continue
        return True, True

    to_keep, other_to_keep = False, False
    # (26612190, 26614205), (26612165, 26614294)
    # Special case: the two transcripts are identical.
    if ovl == length == olength:
        if (is_reference is False and other["is_reference"] is True) or (score < other["score"]):
            other_to_keep = True
        elif (is_reference is True and other["is_reference"] is False) or (score > other["score"]):
            to_keep = True
        else:
            # to_keep, other_to_keep = np.random.permutation([True, False])
            to_keep = random.choice([True, False])
            other_to_keep = not to_keep
    else:
        if (is_reference is False and other["is_reference"] is True) or (score < other["score"]):
            other_to_keep = True
        elif (is_reference is True and other["is_reference"] is False) or (score > other["score"]):
            to_keep = True
        elif ovl == length:
            other_to_keep = True
        elif ovl == olength:
            to_keep = True
        if other["exclude_redundant"] is False:
            other_to_keep = True  # Always keep transcripts marked like this
        if exclude_redundant is False:
            to_keep = True  # Always keep transcripts marked like this
    if not any([other_to_keep, to_keep]):
        raise ValueError((
            (is_reference, other["is_reference"]),
            (score, other["score"]),
            ((start, end), (other["start"], other["end"])))
        )
    return to_keep, other_to_keep


def _check_correspondence(data: dict, other: dict):
    """Verify that two transcripts might be redundant. They must have the same CDS for that to happen."""
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


def _analyse_chrom(chrom: str, keys: pd.DataFrame, shelves, logger):

    merged_transcripts, chains = dict(), defaultdict(IntervalTree)
    current = None
    keys = keys.sort_values(["start", "end"])

    for (start, end), tids in keys.groupby(["start", "end"]):
        if current is not None and overlap(current, (start, end)) < 0:
            if not merged_transcripts:
                logger.debug("No transcript kept for %s:%s-%s", chrom, start, end)
            for tid in merged_transcripts:
                logger.debug("Keeping %s as not redundant", tid)
                yield merged_transcripts[tid]["key"]
            current = (start, end)
            merged_transcripts, chains = dict(), defaultdict(IntervalTree)
        elif current is not None:
            current = tuple([min(start, current[0]), max(end, current[1])])
        else:
            current = (start, end)

        assert (tids.columns == ["start", "end", "strand", "tid", "write_start", "write_length", "shelf"] + [
            "score", "is_reference", "exclude_redundant"]).all(), tids.columns
        for row in tids.values:
            strand, tid, write_start, write_length, shelf_name, score, is_reference = row[2:-1]
            write_start, write_length = int(write_start), int(write_length)
            is_reference = bool(is_reference)
            shelf = shelves[shelf_name]
            to_keep, others_to_remove = True, set()
            data, caught = _retrieve_data(shelf, shelf_name, tid, chrom, (start, end), strand, score,
                                          write_start, write_length, logger,
                                          merged_transcripts, chains)
            if data is None:
                continue
            to_keep = True
            logger.debug("Checking %s (introns: %s; monoexonic: %s)", tid, data["introns"], data["monoexonic"])
            for otid, other in caught:
                check = _check_correspondence(data, other)
                if check:
                    # Redundancy!
                    to_keep, other_to_keep = _select_transcript(
                        is_reference, data["exclude_redundant"],
                        score, other, start=start, end=end)
                    logger.debug("Checking %s vs %s; initial check: %s; KR: %s; OKR: %s",
                                 tid, otid, check, data["exclude_redundant"], other["exclude_redundant"])
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
            if to_keep is True or to_keep is None:
                merged_transcripts[tid] = data
                chains[data["introns"]].add_interval(Interval(data["start"], data["end"], value=tid))
                logger.debug("Keeping %s in the dataset", tid)

    if not merged_transcripts and current is not None:
        logger.debug("No transcript kept for %s:%s-%s", chrom, current[0], current[1])
    for tid in merged_transcripts:
        logger.debug("Keeping %s as not redundant", tid)
        yield merged_transcripts[tid]["key"]


def perform_check(keys, shelve_names, args, logger):

    """
    This is the most important method. After preparing the data structure,
    this function creates the real transcript instances and checks that
    they are correct when looking at the underlying genome sequence.
    This is also the point at which we start using multithreading, if
    so requested.
    :param keys: sorted list of [tid, sequence]
    :param shelve_names: list of the temporary files.
    :param args: the namespace
    :param logger: logger
    :return:
    """

    counter = 0

    # FASTA extraction *has* to be done at the main process level, it's too slow
    # to create an index in each process.

    if args.json_conf["prepare"]["single"] is True or args.json_conf["threads"] == 1:

        shelve_stacks = dict((shelf, open(shelf, "rb")) for shelf in shelve_names)
        # Use functools to pre-configure the function
        # with all necessary arguments aside for the lines
        partial_checker = functools.partial(
            create_transcript,
            canonical_splices=args.json_conf["prepare"]["canonical"],
            logger=logger,
            force_keep_cds=not args.json_conf["prepare"]["strip_cds"])

        for tid, chrom, key in keys:
            tid, shelf_name, write_start, write_length = tid
            try:
                shelf = shelve_stacks[shelf_name]
                shelf.seek(write_start)
                tobj = msgpack.loads(zlib.decompress((shelf.read(write_length))), raw=False)
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
        # np.random.shuffle(batches)
        random.shuffle(batches)
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
        for idx, batch in enumerate(np.array_split(np.array(batches,
                                                            dtype=object), args.json_conf["threads"]), 1):
            batch_file = tempfile.NamedTemporaryFile(delete=False, mode="wb")
            msgpack.dump(batch.tolist(), batch_file)
            batch_file.flush()
            batch_file.close()

            proc = CheckingProcess(
                batch_file.name,
                args.logging_queue,
                args.json_conf["reference"]["genome"].filename,
                idx,
                shelve_names,
                **kwargs)
            try:
                proc.start()
            except TypeError as exc:
                logger.critical("Failed arguments: %s", (batch_file.name,
                args.logging_queue,
                args.json_conf["reference"]["genome"].filename,
                idx,
                shelve_names))
                logger.critical("Failed kwargs: %s", kwargs)
                logger.critical(exc)
                raise
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


row_columns = ["chrom", "start", "end", "strand", "tid", "write_start", "write_length", "shelf"]


def _load_exon_lines_single_thread(args, shelve_names, logger, min_length, strip_cds, max_intron):

    logger.info("Starting to load lines from %d files (single-threaded)",
                len(args.json_conf["prepare"]["files"]["gff"]))
    previous_file_ids = collections.defaultdict(set)
    if args.json_conf["prepare"]["files"]["exclude_redundant"] == []:
        args.json_conf["prepare"]["files"]["exclude_redundant"] = [False] * len(args.json_conf["prepare"]["files"]["gff"])
    if not len(args.json_conf["prepare"]["files"]["exclude_redundant"]) == len(args.json_conf["prepare"]["files"]["gff"]):
        raise AssertionError

    to_do = list(zip(
            shelve_names,
            args.json_conf["prepare"]["files"]["labels"],
            args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
            args.json_conf["prepare"]["files"]["reference"],
            args.json_conf["prepare"]["files"]["exclude_redundant"],
            args.json_conf["prepare"]["files"]["strip_cds"],
            args.json_conf["prepare"]["files"]["gff"],
    ))
    if len(to_do) == 0 and len(shelve_names) > 0:
        raise OSError(
            (
                shelve_names,
                args.json_conf["prepare"]["files"]["labels"],
                args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
                args.json_conf["prepare"]["files"]["reference"],
                args.json_conf["prepare"]["files"]["exclude_redundant"],
                args.json_conf["prepare"]["files"]["strip_cds"],
                args.json_conf["prepare"]["files"]["gff"]
            )
        )

    logger.debug("To do: %d combinations", len(to_do))

    rows = pd.DataFrame([], columns=row_columns)

    for new_shelf, label, strand_specific, is_reference, exclude_redundant, file_strip_cds, gff_name in to_do:
        if file_strip_cds is True:
            file_strip_cds = True
        else:
            file_strip_cds = strip_cds

        logger.info("Starting with %s", gff_name)
        gff_handle = to_gff(gff_name)
        found_ids = set.union(set(), *previous_file_ids.values())
        loader = loaders.get(gff_handle.__annot_type__, None)
        if loader is None:
            raise ValueError("Invalid file type: {}".format(gff_handle.name))
        new_ids, new_rows = loader(new_shelf, gff_handle, label, found_ids, logger,
                                   min_length=min_length, max_intron=max_intron,
                                   strip_cds=file_strip_cds and not is_reference,
                                   exclude_redundant=exclude_redundant, is_reference=is_reference,
                                   strand_specific=strand_specific or is_reference)
        previous_file_ids[gff_handle.name] = new_ids
        new_rows = pd.DataFrame(new_rows, columns=row_columns[:-1])
        new_rows[row_columns[-1]] = new_shelf
        rows = pd.concat([rows, new_rows])

    for column in ["chrom", "strand", "tid"]:
        rows[column] = rows[column].str.decode("utf-8")

    return rows


def _load_exon_lines_multi(args, shelve_names, logger, min_length, strip_cds, threads, max_intron=3*10**5):
    logger.info("Starting to load lines from %d files (using %d processes)",
                len(args.json_conf["prepare"]["files"]["gff"]), threads)
    manager = multiprocessing.Manager()
    submission_queue = manager.JoinableQueue(-1)
    return_queue = manager.JoinableQueue(-1)
    working_processes = []
    # working_processes = [ for _ in range(threads)]

    for num in range(threads):
        proc = AnnotationParser(submission_queue,
                                return_queue,
                                args.logging_queue,
                                num + 1,
                                log_level=args.level,
                                min_length=min_length,
                                max_intron=max_intron,
                                strip_cds=strip_cds,
                                seed=args.json_conf["seed"])
        proc.start()
        working_processes.append(proc)

    shelve_df = []
    for shelf_index, (new_shelf, label, strand_specific, is_reference,
                      exclude_redundant, file_strip_cds, gff_name) in enumerate(zip(
            shelve_names,
            args.json_conf["prepare"]["files"]["labels"],
            args.json_conf["prepare"]["files"]["strand_specific_assemblies"],
            args.json_conf["prepare"]["files"]["reference"],
            args.json_conf["prepare"]["files"]["exclude_redundant"],
            args.json_conf["prepare"]["files"]["strip_cds"],
            args.json_conf["prepare"]["files"]["gff"])):
        submission_queue.put((label, gff_name, strand_specific, is_reference, exclude_redundant, file_strip_cds,
                              new_shelf, shelf_index))
        shelve_df.append((shelf_index, new_shelf))

    shelve_df = pd.DataFrame(shelve_df, columns=["shelf_index", "shelf"])
    submission_queue.put(tuple(["EXIT"] * 8))
    
    rows = []

    retrieved = 0
    import sys
    while retrieved < len(working_processes):
        if return_queue.empty():
            continue
        row = return_queue.get(block=False)
        if row == "FINISHED":
            retrieved += 1
        else:
            rows.append(row)
        continue

    [_.join() for _ in working_processes]
    
    rows = pd.DataFrame(rows, columns=row_columns[:-1] + ["shelf_index"])
    rows = rows.merge(shelve_df, on="shelf_index", how="left").drop(["shelf_index"], axis=1)
    for key in ["chrom", "tid", "strand"]:
         rows[key] = rows[key].str.decode("utf-8")

    del working_processes
    gc.collect()
    logger.info("Finished parsing all input files")
    manager.shutdown()
    return rows


def load_exon_lines(args, shelve_names, logger, min_length=0, max_intron=3*10**5) -> pd.DataFrame:

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
        rows = _load_exon_lines_single_thread(args, shelve_names, logger, min_length, strip_cds, max_intron)
    else:
        rows = _load_exon_lines_multi(args, shelve_names, logger, min_length, strip_cds, threads, max_intron)

    logger.info("Finished loading lines from %d files",
                len(args.json_conf["prepare"]["files"]["gff"]))

    if rows["tid"].duplicated().any():
        if set(args.json_conf["prepare"]["files"]["labels"]) == {""}:
            exception = exceptions.RedundantNames(
                """Found redundant names during multiprocessed file analysis.\
Please repeat using distinct labels for your input files. Aborting.""")
        else:
            with open("rows.tsv", "wt") as out:
                rows.to_csv(out, sep="\t", header=True, index=False)

            exception = exceptions.RedundantNames(
                """Found redundant names during multiprocessed file analysis, even if \
unique labels were provided. Please try to repeat with a different and more unique set of labels. Aborting.""")
        logger.exception(exception)
        raise exception

    rows["strand"] = rows["strand"].mask(rows["strand"] == ".", None)
    return rows


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
        error = "Invalid FASTA file: {}".format(args.json_conf["reference"]["genome"])
        logger.critical(error)
        raise AttributeError(error)

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

    ref_len = len(args.json_conf["prepare"]["files"]["reference"])
    file_len = len(args.json_conf["prepare"]["files"]["gff"])
    if ref_len == 0:
        args.json_conf["prepare"]["files"]["reference"] = ([False] * file_len)
    elif (ref_len != file_len) or (args.json_conf["prepare"]["files"]["reference"][0] not in (True, False)):
        ref_set = set(args.json_conf["prepare"]["files"]["reference"])
        args.json_conf["prepare"]["files"]["reference"] = [
            (_ in ref_set) for _ in args.json_conf["prepare"]["files"]["gff"]
        ]

    if not args.json_conf["prepare"]["files"]["exclude_redundant"]:
        args.json_conf["prepare"]["files"]["exclude_redundant"] = (
            [getattr(args, "exclude_redundant", False)] * len(args.json_conf["prepare"]["files"]["gff"]))

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
    errored = False
    try:
        # chrom, start, end, strand, tid, write_start, write_length, shelf
        rows = load_exon_lines(
            args, shelve_names, logger,
            min_length=args.json_conf["prepare"]["minimum_cdna_length"],
            max_intron=args.json_conf["prepare"]["max_intron_length"],)

        logger.info("Finished loading exon lines")

        shelve_source_scores = []
        for label in args.json_conf["prepare"]["files"]["labels"]:
            shelve_source_scores.append(
                args.json_conf["prepare"]["files"]["source_score"].get(label, 0)
            )

        shelve_table = []

        for shelf, score, is_reference, exclude_redundant in zip(
                shelve_names, shelve_source_scores,
                args.json_conf["prepare"]["files"]["reference"],
                args.json_conf["prepare"]["files"]["exclude_redundant"]):
            assert isinstance(is_reference, bool), \
                (is_reference, args.json_conf["prepare"]["files"]["reference"])
            shelve_table.append((shelf, score, is_reference, exclude_redundant))

        shelve_table = pd.DataFrame(shelve_table, columns=["shelf", "score", "is_reference", "exclude_redundant"])

        rows = rows.merge(shelve_table, on="shelf", how="left")
        random.seed(args.json_conf["seed"])

        shelves = dict((shelf_name, open(shelf_name, "rb")) for shelf_name in shelve_table["shelf"].unique())

        def divide_by_chrom():
            # chrom, start, end, strand, tid, write_start, write_length, shelf
            transcripts = rows.groupby(["chrom"])
            columns = rows.columns[1:]
            for chrom in sorted(transcripts.groups.keys()):
                logger.debug("Starting with %s (%d positions)",
                             chrom, transcripts.size()[chrom])

                yield from _analyse_chrom(chrom, rows.loc[transcripts.groups[chrom], columns],
                                          shelves, logger=logger)

        perform_check(divide_by_chrom(), shelve_names, args, logger)
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
        logging.shutdown()
    else:
        logger.error("Mikado prepare has encountered a fatal error. Please check the logs and, if there is a bug,"\
                     "report it to https://github.com/EI-CoreBioinformatics/mikado/issues")
        logging.shutdown()
        exit(1)
