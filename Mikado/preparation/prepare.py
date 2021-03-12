import os
import tempfile
import gc
from .checking import create_transcript, CheckingProcess
from .annotation_parser import AnnotationParser, loaders, load_into_storage
from ..configuration import MikadoConfiguration
from ..exceptions import InvalidConfiguration, InvalidParsingFormat
from ..utilities import Interval, IntervalTree
from ..parsers import parser_factory
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
import sys
import pysam
import numpy as np
import random
import pandas as pd
import zlib


__author__ = 'Luca Venturini'


def __cleanup(mikado_config, shelves):
    """Private function to close the opened handles."""

    if hasattr(mikado_config.reference.genome, "close"):
        mikado_config.reference.genome.close()
        mikado_config.reference.genome = mikado_config.reference.genome.filename

    if hasattr(mikado_config.prepare.files.out, "close"):
        mikado_config.prepare.files.out.close()
        mikado_config.prepare.files.out = mikado_config.prepare.files.out.name

    if hasattr(mikado_config.prepare.files.out_fasta, "close"):
        mikado_config.prepare.files.out_fasta.close()
        mikado_config.prepare.files.out_fasta = mikado_config.prepare.files.out_fasta.name

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
                                                            n=10**5)
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
                chains[data["introns"]].add(Interval(data["start"], data["end"], value=tid))
                logger.debug("Keeping %s in the dataset", tid)

    if not merged_transcripts and current is not None:
        logger.debug("No transcript kept for %s:%s-%s", chrom, current[0], current[1])
    for tid in merged_transcripts:
        logger.debug("Keeping %s as not redundant", tid)
        yield merged_transcripts[tid]["key"]


def perform_check(keys, shelve_names, mikado_config: MikadoConfiguration, logger):

    """
    This is the most important method. After preparing the data structure,
    this function creates the real transcript instances and checks that
    they are correct when looking at the underlying genome sequence.
    This is also the point at which we start using multithreading, if
    so requested.
    :param keys: sorted list of [tid, sequence]
    :param shelve_names: list of the temporary files.
    :param mikado_config: MikadoConfiguration
    :param logger: logger
    :return:
    """

    counter = 0

    # FASTA extraction *has* to be done at the main process level, it's too slow
    # to create an index in each process.

    if mikado_config.prepare.single is True or mikado_config.threads == 1:

        shelve_stacks = dict((shelf, open(shelf, "rb")) for shelf in shelve_names)
        # Use functools to pre-configure the function
        # with all necessary arguments aside for the lines
        partial_checker = functools.partial(
            create_transcript,
            canonical_splices=mikado_config.prepare.canonical,
            logger=logger,
            strip_faulty_cds=mikado_config.prepare.strip_faulty_cds)

        for tid, chrom, key in keys:
            tid, shelf_name, write_start, write_length = tid
            try:
                shelf = shelve_stacks[shelf_name]
                shelf.seek(write_start)
                tobj = msgpack.loads(zlib.decompress((shelf.read(write_length))), raw=False)
            except sqlite3.ProgrammingError as exc:
                raise sqlite3.ProgrammingError("{}. Tids: {}".format(exc, tid))

            if chrom not in mikado_config.reference.genome.references:
                raise KeyError("Invalid chromosome name! {}, {}, {}, {}".format(tid, shelf_name, chrom, key))

            try:
                seq = str(mikado_config.reference.genome.fetch(chrom, key[0] - 1, key[1]))
            except ValueError:
                raise ValueError(tobj)

            transcript_object = partial_checker(
                tobj,
                seq,
                key[0], key[1],
                lenient=mikado_config.prepare.lenient,
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
                  file=mikado_config.prepare.files.out)
            print(transcript_object.fasta,
                  file=mikado_config.prepare.files.out_fasta)
    else:
        # pylint: disable=no-member

        # submission_queue = multiprocessing.JoinableQueue(-1)

        batches = list(enumerate(keys, 1))
        # np.random.shuffle(batches)
        random.shuffle(batches)
        kwargs = {
            "fasta_out": os.path.basename(mikado_config.prepare.files.out_fasta.name),
            "gtf_out": os.path.basename(mikado_config.prepare.files.out.name),
            "tmpdir": mikado_config.tempdir.name,
            "seed": mikado_config.seed,
            "lenient": mikado_config.prepare.lenient,
            "canonical_splices": mikado_config.prepare.canonical,
            "strip_faulty_cds": mikado_config.prepare.strip_faulty_cds,
            "log_level": mikado_config.log_settings.log_level
        }

        working_processes = []
        batch_files = []
        for idx, batch in enumerate(np.array_split(np.array(batches,
                                                            dtype=object), mikado_config.threads), 1):
            batch_file = tempfile.NamedTemporaryFile(delete=True, mode="wb")
            msgpack.dump(batch.tolist(), batch_file)
            batch_file.flush()
            batch_files.append(batch_file)

            proc = CheckingProcess(
                batch_file.name,
                mikado_config.logging_queue,
                mikado_config.reference.genome.filename,
                idx,
                shelve_names,
                **kwargs)
            try:
                proc.start()
            except TypeError as exc:
                logger.critical("Failed arguments: %s", (batch_file.name,
                                                         mikado_config.logging_queue,
                                                         mikado_config.reference.genome.filename,
                                                         idx,
                                                         shelve_names))
                logger.critical("Failed kwargs: %s", kwargs)
                logger.critical(exc)
                raise
            working_processes.append(proc)

        [_.join() for _ in working_processes]

        partial_gtf = [os.path.join(mikado_config.tempdir.name,
                                    "{0}-{1}".format(
                                        os.path.basename(mikado_config.prepare.files.out.name),
                                        _ + 1)) for _ in range(mikado_config.threads)]
        merge_partial(partial_gtf, mikado_config.prepare.files.out)

        partial_fasta = [os.path.join(
            mikado_config.tempdir.name,
            "{0}-{1}".format(os.path.basename(mikado_config.prepare.files.out_fasta.name), _ + 1))
                         for _ in range(mikado_config.threads)]
        merge_partial(partial_fasta, mikado_config.prepare.files.out_fasta)
        [batch_file.close() for batch_file in batch_files]

    mikado_config.prepare.files.out_fasta.close()
    mikado_config.prepare.files.out.close()

    logger.setLevel(logging.INFO)
    # logger.info("Finished to analyse %d transcripts (%d retained)",
    #             len(exon_lines), counter)
    logger.setLevel(mikado_config.log_settings.log_level)
    return


row_columns = ["chrom", "start", "end", "strand", "tid", "write_start", "write_length", "shelf"]


def _get_strand_specific_assemblies_boolean_vector(mikado_config):
    return [(member in mikado_config.prepare.files.strand_specific_assemblies)
            for member in mikado_config.prepare.files.gff]


def _load_exon_lines_single_thread(mikado_config, shelve_names, logger, min_length, strip_cds, max_intron):

    logger.info("Starting to load lines from %d files (single-threaded)",
                len(mikado_config.prepare.files.gff))
    previous_file_ids = collections.defaultdict(set)
    if mikado_config.prepare.files.exclude_redundant == []:
        mikado_config.prepare.files.exclude_redundant = [False] * len(mikado_config.prepare.files.gff)
    if not len(mikado_config.prepare.files.exclude_redundant) == len(mikado_config.prepare.files.gff):
        raise AssertionError

    to_do = list(zip(
            shelve_names,
            mikado_config.prepare.files.labels,
            _get_strand_specific_assemblies_boolean_vector(mikado_config),
            mikado_config.prepare.files.reference,
            mikado_config.prepare.files.exclude_redundant,
            mikado_config.prepare.files.strip_cds,
            mikado_config.prepare.files.gff,
    ))
    if len(to_do) == 0 and len(shelve_names) > 0:
        raise OSError(
            (
                shelve_names,
                mikado_config.prepare.files.labels,
                _get_strand_specific_assemblies_boolean_vector(mikado_config),
                mikado_config.prepare.files.reference,
                mikado_config.prepare.files.exclude_redundant,
                mikado_config.prepare.files.strip_cds,
                mikado_config.prepare.files.gff
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
        try:
            gff_handle = parser_factory(gff_name)
        except InvalidParsingFormat as exc:
            logger.exception("Invalid file: %s. Skipping it", gff_name)
            logger.exception(exc)
            previous_file_ids[gff_name] = set()
            load_into_storage(new_shelf, [], min_length, logger, strip_cds=True,
                              max_intron=3 * 10 ** 5)
            continue
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


def _load_exon_lines_multi(mikado_config, shelve_names, logger, min_length, strip_cds, threads, max_intron=3 * 10 ** 5):
    logger.info("Starting to load lines from %d files (using %d processes)",
                len(mikado_config.prepare.files.gff), threads)
    manager = multiprocessing.Manager()
    submission_queue = manager.JoinableQueue(-1)
    return_queue = manager.JoinableQueue(-1)
    working_processes = []
    # working_processes = [ for _ in range(threads)]

    for num in range(threads):
        proc = AnnotationParser(submission_queue,
                                return_queue,
                                mikado_config.logging_queue,
                                num + 1,
                                log_level=mikado_config.log_settings.log_level,
                                min_length=min_length,
                                max_intron=max_intron,
                                strip_cds=strip_cds,
                                seed=mikado_config.seed)
        proc.start()
        working_processes.append(proc)

    shelve_df = []
    for shelf_index, (new_shelf, label, strand_specific, is_reference,
                      exclude_redundant, file_strip_cds, gff_name) in enumerate(zip(
            shelve_names,
            mikado_config.prepare.files.labels,
            _get_strand_specific_assemblies_boolean_vector(mikado_config),
            mikado_config.prepare.files.reference,
            mikado_config.prepare.files.exclude_redundant,
            mikado_config.prepare.files.strip_cds,
            mikado_config.prepare.files.gff)):
        submission_queue.put((label, gff_name, strand_specific, is_reference, exclude_redundant, file_strip_cds,
                              new_shelf, shelf_index))
        shelve_df.append((shelf_index, new_shelf))

    shelve_df = pd.DataFrame(shelve_df, columns=["shelf_index", "shelf"])
    submission_queue.put(tuple(["EXIT"] * 8))
    
    rows = []

    retrieved = 0
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


def load_exon_lines(mikado_config, shelve_names, logger, min_length=0, max_intron=3 * 10 ** 5) -> pd.DataFrame:

    """This function loads all exon lines from the GFF inputs into a
     defaultdict instance.
    :param mikado_config: the Namespace from the command line.
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

    threads = min([len(mikado_config.prepare.files.gff),
                   mikado_config.threads])
    strip_cds = mikado_config.prepare.strip_cds

    if mikado_config.prepare.single is True or threads == 1:
        rows = _load_exon_lines_single_thread(mikado_config, shelve_names, logger, min_length, strip_cds, max_intron)
    else:
        rows = _load_exon_lines_multi(mikado_config, shelve_names, logger, min_length, strip_cds, threads, max_intron)

    logger.info("Finished loading lines from %d files",
                len(mikado_config.prepare.files.gff))

    if rows["tid"].duplicated().any():
        if set(mikado_config.prepare.files.labels) == {""}:
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


def prepare(mikado_config: MikadoConfiguration, logger):
    """Main script function.

    :param mikado_config: the ArgumentParser-derived namespace.
    :param logger: a logging instance
    :type logger: logging.Logger
    """

    if not hasattr(mikado_config.reference, "genome"):
        raise InvalidConfiguration("Invalid configuration; reference: {}".format(mikado_config))

    if hasattr(mikado_config.reference.genome, "close"):
        mikado_config.reference.genome.close()
        if hasattr(mikado_config.reference.genome, "filename"):
            mikado_config.reference.genome = getattr(mikado_config.reference.genome, "filename")
        elif hasattr(mikado_config.reference.genome, "name"):
            mikado_config.reference.genome = getattr(mikado_config.reference.genome, "name")
        else:
            logger.critical("Invalid FASTA file: %s", mikado_config.reference.genome)
            raise AttributeError
    elif not isinstance(mikado_config.reference.genome, (str, bytes)):
        logger.critical("Invalid FASTA file: %s", mikado_config.reference.genome)
        raise AttributeError

    if not os.path.exists(mikado_config.reference.genome):
        error = "Invalid FASTA file: {}".format(mikado_config.reference.genome)
        logger.critical(error)
        raise AttributeError(error)

    assert len(mikado_config.prepare.files.gff) > 0
    assert len(mikado_config.prepare.files.gff) == len(mikado_config.prepare.files.labels), (
        mikado_config.prepare.files.gff,
        mikado_config.prepare.files.labels
    )

    if mikado_config.prepare.strand_specific is True:
        mikado_config.prepare.files.strand_specific_assemblies = mikado_config.prepare.files.gff[:]

    ref_len = len(mikado_config.prepare.files.reference)
    file_len = len(mikado_config.prepare.files.gff)
    if ref_len == 0:
        mikado_config.prepare.files.reference = ([False] * file_len)
    elif (ref_len != file_len) or (mikado_config.prepare.files.reference[0] not in (True, False)):
        ref_set = set(mikado_config.prepare.files.reference)
        mikado_config.prepare.files.reference = [
            (_ in ref_set) for _ in mikado_config.prepare.files.gff
        ]

    if not mikado_config.prepare.files.exclude_redundant:
        mikado_config.prepare.files.exclude_redundant = (
                [getattr(mikado_config, "exclude_redundant", False)] * len(mikado_config.prepare.files.gff))

    shelve_names = [path_join(mikado_config.prepare.files.output_dir,
                              "mikado_shelf_{}.db".format(str(_).zfill(5))) for _ in
                    range(len(mikado_config.prepare.files.gff))]

    logger.propagate = False
    if mikado_config.prepare.single is False and mikado_config.threads > 1:
        multiprocessing.set_start_method(mikado_config.multiprocessing_method,
                                         force=True)

        mikado_config.logging_queue = multiprocessing.JoinableQueue(-1)
        log_queue_handler = logging.handlers.QueueHandler(mikado_config.logging_queue)
        log_queue_handler.setLevel(logging.DEBUG)
        # logger.addHandler(log_queue_handler)
        mikado_config.tempdir = tempfile.TemporaryDirectory(dir=mikado_config.prepare.files.output_dir)
        mikado_config.listener = logging.handlers.QueueListener(mikado_config.logging_queue, logger)
        mikado_config.listener.propagate = False
        mikado_config.listener.start()

    mikado_config.prepare.files.out_fasta = open(
        path_join(mikado_config.prepare.files.output_dir,
                  mikado_config.prepare.files.out_fasta), 'w')
    mikado_config.prepare.files.out = open(path_join(
        mikado_config.prepare.files.output_dir,
        mikado_config.prepare.files.out), 'w')

    logger.info("Output dir: %s. Output GTF: %s. Output Fasta: %s",
                mikado_config.prepare.files.output_dir,
                mikado_config.prepare.files.out.name,
                mikado_config.prepare.files.out_fasta.name)
    logger.info("Loading reference file")
    mikado_config.reference.genome = pysam.FastaFile(mikado_config.reference.genome)
    logger.info("Finished loading genome file")
    logger.info("Started loading exon lines")
    errored = False
    try:
        # chrom, start, end, strand, tid, write_start, write_length, shelf
        rows = load_exon_lines(
            mikado_config, shelve_names, logger,
            min_length=mikado_config.prepare.minimum_cdna_length,
            max_intron=mikado_config.prepare.max_intron_length,)

        logger.info("Finished loading exon lines")

        shelve_source_scores = []
        for label in mikado_config.prepare.files.labels:
            shelve_source_scores.append(
                mikado_config.prepare.files.source_score.get(label, 0)
            )

        shelve_table = []

        for shelf, score, is_reference, exclude_redundant in zip(
                shelve_names, shelve_source_scores,
                mikado_config.prepare.files.reference,
                mikado_config.prepare.files.exclude_redundant):
            assert isinstance(is_reference, bool), \
                (is_reference, mikado_config.prepare.files.reference)
            shelve_table.append((shelf, score, is_reference, exclude_redundant))

        shelve_table = pd.DataFrame(shelve_table, columns=["shelf", "score", "is_reference", "exclude_redundant"])

        rows = rows.merge(shelve_table, on="shelf", how="left")
        random.seed(mikado_config.seed)

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

        perform_check(divide_by_chrom(), shelve_names, mikado_config, logger)
    except Exception as exc:
        # TODO: Consider using stderr to signal errors here too?
        logger.exception(exc)
        __cleanup(mikado_config, shelve_names)
        errored = True
        logger.error("Mikado has encountered an error, exiting")
        # sys.exit(1)

    if mikado_config.prepare.single is False and mikado_config.threads > 1:
        mikado_config.tempdir.cleanup()
        mikado_config.listener.enqueue_sentinel()

    logger.setLevel(logging.INFO)
    __cleanup(mikado_config, shelve_names)

    logger.addHandler(logging.StreamHandler())
    if errored is False:
        logger.info("Mikado prepare has finished correctly with seed %s. The output %s FASTA file can now be "
                    "used for BLASTX and/or ORF calling before the next step in the pipeline, `mikado serialise`.",
                    mikado_config.seed, mikado_config.prepare.files.out_fasta)
        logging.shutdown()
    else:
        logger.error("Mikado prepare has encountered a fatal error. Please check the logs and, if there is a bug,"\
                     "report it to https://github.com/EI-CoreBioinformatics/mikado/issues")
        logging.shutdown()
        exit(1)
