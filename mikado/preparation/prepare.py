import os
import sys
from ..loci_objects.transcriptchecker import TranscriptChecker
from ..loci_objects import Transcript
import operator
import collections
import io
from .. import exceptions
import logging
import logging.handlers
import random
import functools
import multiprocessing
import multiprocessing.connection
import multiprocessing.sharedctypes
import pyfaidx
from ..utilities import path_join, to_gff

__author__ = 'Luca Venturini'


def grouper(iterable, num=100):

    """
    Simple function to prechunk an iterable into groups of fixed size.
    :param iterable: the iterable to chunk.
    :param num: Maximum number of elements per chunk. Default: 100.
    :return:
    """

    cache = []

    for el in iterable:
        cache.append(el)
        if len(cache) >= num:
            yield cache
            cache = []
    yield cache


def create_transcript(lines,
                      fasta_seq,
                      start,
                      end,
                      lenient=False,
                      strand_specific=False,
                      canonical_splices=(("GT", "AG"),
                                         ("GC", "AG"),
                                         ("AT", "AC")),
                      log_queue=None):
    """Function to create the checker.

    :param lines: all the exon lines for an object
    :type lines: dict

    :param fasta_seq: genomic sequence of the transcript

    :param start: start position for the transcript
    :type start: int
    :param end: end position for the transcript
    :type end: int

    :type lenient: bool
    :type strand_specific: bool

    :param canonical_splices: the splices considered as canonical for the species.
    :type canonical_splices: list[tuple]

    :param log_queue: the optional Queue to use for logging (during multiprocessing)

    :rtype: (None|TranscriptChecker)
    """
    logger = logging.getLogger("main")
    if log_queue is not None:
        handler = logging.handlers.QueueHandler(log_queue)
        logger.addHandler(handler)

    logger.debug("Starting with %s", lines["tid"])

    try:
        transcript_line = Transcript()
        transcript_line.chrom = lines["chrom"]
        transcript_line.strand = lines["strand"]
        transcript_line.attributes.update(lines["attributes"])
        transcript_line.feature = "transcript"
        transcript_line.start, transcript_line.end = sorted([start, end])
        transcript_line.logger = logger
        transcript_line.id = lines["tid"]
        transcript_line.parent = lines["parent"]

        for feature in lines["features"]:
            transcript_line.add_exons(lines["features"][feature],
                                      features=feature)
        transcript_object = TranscriptChecker(transcript_line,
                                              fasta_seq,
                                              lenient=lenient,
                                              strand_specific=strand_specific,
                                              canonical_splices=canonical_splices,
                                              logger=logger)
        logger.debug("Finished adding exon lines to %s", lines["tid"])
        transcript_object.finalize()
        transcript_object.check_strand()
    except exceptions.IncorrectStrandError:
        logger.info("Discarded %s because of incorrect fusions of splice junctions",
                    lines["tid"])
        # logger.exception(exc)
        transcript_object = None
    except exceptions.InvalidTranscript:
        logger.info("Discarded generically invalid transcript %s",
                    lines["tid"])
        transcript_object = None
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except Exception as exc:
        logger.exception(exc)
        transcript_object = None

    logger.debug("Finished with %s", lines["tid"])

    return transcript_object


def store_transcripts(exon_lines, logger, min_length=0):

    """
    Function that analyses the exon lines from the original file
    and organises the data into a proper dictionary.
    :param exon_lines: dictionary of exon lines, ordered by TID
    :type exon_lines: dict

    :param logger: logger instance.
    :type logger: logging.Logger

    :param min_length: minimal length of the transcript.
    If it is not met, the transcript will be discarded.
    :type min_length: int

    :return: transcripts: dictionary which will be the final output
    :rtype: transcripts
    """

    logger.info("Starting to organise %d transcripts", len(exon_lines))
    transcripts = collections.defaultdict(dict)
    for tid in exon_lines:

        if ("exon" not in exon_lines[tid]["features"] or
                len(exon_lines[tid]["features"]["exon"]) == 0):
            logger.warning("No valid exon feature for %s, continuing", tid)

        chrom = exon_lines[tid]["chrom"]
        start = min((_[0] for _ in exon_lines[tid]["features"]["exon"]))
        end = max((_[1] for _ in exon_lines[tid]["features"]["exon"]))
        tlength = sum(exon[1] + 1 - exon[0] for exon in exon_lines[tid]["features"]["exon"])
        # Discard transcript under a certain size
        if tlength < min_length:
            logger.debug("Discarding %s because its size (%d) is under the minimum of %d",
                         tid, tlength, min_length)
            continue

        if (start, end) not in transcripts[chrom]:
            transcripts[chrom][(start, end)] = []
        transcripts[chrom][(start, end)].append(tid)

    logger.info("Starting to sort %d transcripts", len(exon_lines))
    # keys = []
    # counter = 0
    for chrom in sorted(transcripts.keys()):

        logger.debug("Starting with %s (%d positions)", chrom, len(transcripts[chrom]))

        for key in sorted(transcripts[chrom].keys(),
                          key=operator.itemgetter(0, 1)):
            tids = transcripts[chrom][key]
            if len(tids) > 1:
                exons = collections.defaultdict(list)
                for tid in tids:
                    exon_set = tuple(sorted(
                        [(exon[0], exon[1], exon_lines[tid]["strand"]) for exon in
                         exon_lines[tid]["features"]["exon"]],
                        key=operator.itemgetter(0, 1)))
                    exons[exon_set].append(tid)
                tids = []
                logger.debug("%d intron chains for pos %s",
                             len(exons), "{}:{}-{}".format(chrom, key[0], key[1]))
                for tid_list in exons.values():
                    if len(tid_list) > 1:
                        logger.debug("The following transcripts are redundant: %s",
                                     ",".join(tid_list))
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

    # FASTA extraction *has* to be done at the main process level, it's too slow
    # to create an index in each process.

    if args.json_conf["prepare"]["single"] is True:

        # Use functools to pre-configure the function
        # with all necessary arguments aside for the lines
        partial_checker = functools.partial(
            create_transcript,
            lenient=args.json_conf["prepare"]["lenient"],
            strand_specific=args.json_conf["prepare"]["strand_specific"],
            log_queue=None)

        for tid, chrom, key in keys:
            transcript_object = partial_checker(
                exon_lines[tid],
                key[0], key[1],
                str(args.json_conf["prepare"]["fasta"][chrom][key[0]-1:key[1]]))
            if transcript_object is None:
                continue
            counter += 1
            if counter >= 10**4 and counter % (10**4) == 0:
                logger.info("Retrieved %d transcript positions", counter)
            elif counter >= 10**3 and counter % (10**3) == 0:
                logger.debug("Retrieved %d transcript positions", counter)
            print(transcript_object.format("gtf"),
                  file=args.json_conf["prepare"]["out"])
            print(transcript_object.fasta,
                  file=args.json_conf["prepare"]["out_fasta"])
    else:
        # pylint: disable=no-member
        context = multiprocessing.get_context()
        manager = context.Manager()
        logging_queue = manager.Queue(-1)
        queue_handler = logging.handlers.QueueHandler(logging_queue)
        logger.addHandler(queue_handler)

        partial_checker = functools.partial(
            create_transcript,
            lenient=args.json_conf["prepare"]["lenient"],
            strand_specific=args.json_conf["prepare"]["strand_specific"],
            log_queue=logging_queue)

        pool = context.Pool(args.threads)

        # pylint: enable=no-member

        for group in grouper(keys, args.threads):

            results = pool.starmap_async(
                partial_checker,
                [(exon_lines[tid],
                  str(args.json_conf["prepare"]["fasta"][chrom][key[0]-1:key[1]]),
                  key[0], key[1])
                 for (tid, chrom, key) in group]
            )

            for transcript_object in results.get():

                if transcript_object is None:
                    continue
                counter += 1
                if counter >= 10**4 and counter % (10**4) == 0:
                    logger.info("Retrieved %d transcript positions", counter)
                elif counter >= 10**3 and counter % (10**3) == 0:
                    logger.debug("Retrieved %d transcript positions", counter)
                print(transcript_object.format("gtf"),
                      file=args.json_conf["prepare"]["out"])
                print(transcript_object.fasta,
                      file=args.json_conf["prepare"]["out_fasta"])

        pool.close()
        pool.join()

    logger.info("Finished to analyse %d transcripts (%d retained)",
                len(exon_lines), counter)
    return


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
            if row.id in found_ids:
                assert label == ''
                raise ValueError("""{0} has already been found in another file,
                this will cause unsolvable collisions. Please rerun preparation using
                labels to tag each file.""")
            assert row.id not in exon_lines
            exon_lines[row.id]["attributes"] = row.attributes.copy()
            exon_lines[row.id]["chrom"] = row.chrom
            exon_lines[row.id]["strand"] = row.strand
            exon_lines[row.id]["tid"] = row.transcript
            exon_lines[row.id]["parent"] = row.parent
            exon_lines[row.id]["features"] = dict()
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
                    if tid not in exon_lines:
                        exon_lines[tid]["attributes"] = row.attributes.copy()
                        exon_lines[tid]["chrom"] = row.chrom
                        exon_lines[tid]["strand"] = row.strand
                        exon_lines[tid]["features"] = dict()
                        exon_lines[row.id]["tid"] = tid
                        exon_lines[row.id]["parent"] = transcript2genes[tid]
                    else:
                        if "exon_number" in row.attributes:
                            del row.attributes["exon_number"]
                        assert exon_lines[tid]["chrom"] == row.chrom
                        assert exon_lines[tid]["strand"] == row.strand
                        exon_lines[tid]["attributes"].update(row.attributes)

                    if row.feature not in exon_lines[tid]["features"]:
                        exon_lines[tid]["features"][row.feature] = []
                    exon_lines[tid]["features"][row.feature].append((row.start, row.end))
                    new_ids.add(tid)
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
        if row.is_transcript is True:
            if label != '':
                row.transcript = "{0}_{1}".format(label, row.transcript)
            if row.transcript in found_ids:
                assert label != ''
                raise ValueError("""{0} has already been found in another file,
                    this will cause unsolvable collisions. Please rerun preparation using
                    labels to tag each file.""")
            assert row.id not in exon_lines
            exon_lines[row.id]["features"] = dict()
            exon_lines[row.id]["chrom"] = row.chrom
            exon_lines[row.id]["strand"] = row.strand
            exon_lines[row.id]["attributes"] = row.attributes.copy()
            exon_lines[row.id]["tid"] = row.id
            exon_lines[row.id]["parent"] = row.gene
            if "exon_number" in exon_lines[row.id]["attributes"]:
                del exon_lines[row.id]["attributes"]["exon_number"]
            continue
        if row.is_exon is False or (row.is_cds is True and strip_cds is True):
            continue
        if label != '':
            row.transcript = "{0}_{1}".format(label, row.transcript)
        if row.transcript in found_ids:
            assert label != ''
            raise ValueError("""{0} has already been found in another file,
                this will cause unsolvable collisions. Please rerun preparation using
                labels to tag each file.""")
        assert row.transcript is not None
        if row.transcript not in exon_lines:
            exon_lines[row.transcript]["chrom"] = row.chrom
            exon_lines[row.transcript]["strand"] = row.strand
            exon_lines[row.transcript]["exon"] = []
            exon_lines[row.transcript]["attributes"] = row.attributes.copy()
            exon_lines[row.id]["tid"] = row.transcript
            exon_lines[row.id]["parent"] = row.gene
        else:
            if "exon_number" in row.attributes:
                del row.attributes["exon_number"]
            assert "chrom" in exon_lines[row.transcript], (row.transcript, exon_lines)
            assert exon_lines[row.transcript]["chrom"] == row.chrom, exon_lines[row.transcript]
            assert exon_lines[row.transcript]["strand"] == row.strand
            exon_lines[row.transcript]["attributes"].update(row.attributes)
        if row.feature not in exon_lines[row.transcript]["features"]:
            exon_lines[row.transcript]["features"][row.feature] = []
        exon_lines[row.transcript]["features"][row.feature].append((row.start, row.end))
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

    exon_lines = collections.defaultdict(dict)
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


def prepare(args, logger):
    """Main script function.

    :param args: the ArgumentParser-derived namespace.
    :param logger: a logging instance
    """

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

    args.json_conf["prepare"]["out_fasta"] = open(
        path_join(args.json_conf["prepare"]["output_dir"],
                  args.json_conf["prepare"]["out_fasta"]), 'w')
    args.json_conf["prepare"]["out"] = open(path_join(
        args.json_conf["prepare"]["output_dir"],
        args.json_conf["prepare"]["out"]), 'w')

    logger.info("Loading reference file")
    args.json_conf["prepare"]["fasta"] = pyfaidx.Fasta(args.json_conf["prepare"]["fasta"])

    logger.info("Finished loading reference file")
    logger.info("Started loading exon lines")

    exon_lines = load_exon_lines(args, logger)

    logger.info("Finished loading exon lines")

    # Prepare the sorted data structure
    sorter = functools.partial(
        store_transcripts,
        logger=logger,
        min_length=args.json_conf["prepare"]["minimum_length"]
    )

    perform_check(sorter(exon_lines), exon_lines, args, logger)

    args.json_conf["prepare"]["out"].close()
    args.json_conf["prepare"]["out_fasta"].close()

    logger.info("Finished")
