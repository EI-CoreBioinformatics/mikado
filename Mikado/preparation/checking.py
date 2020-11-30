import functools
import multiprocessing
import multiprocessing.queues
import os
import zlib
import pysam
import msgpack
from ..transcripts.transcriptchecker import TranscriptChecker
from .. import exceptions
from ..loci import Transcript
from ..utilities.log_utils import create_null_logger, create_queue_logger
import logging
import queue
import time
import sys
import rapidjson as json
import operator
import random
import zlib


__author__ = 'Luca Venturini'


def create_transcript(lines,
                      fasta_seq,
                      start,
                      end,
                      lenient=False,
                      is_reference=False,
                      strand_specific=False,
                      canonical_splices=(("GT", "AG"),
                                         ("GC", "AG"),
                                         ("AT", "AC")),
                      force_keep_cds=False,
                      logger=None):
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

    :param force_keep_cds: boolean. If set to true, coding transcripts that would be flipped are instead excluded.
                           The intention is that this flag will mirror strip_cds.
    :type force_keep_cds: bool

    :param logger: optional logger to use during processing.

    :param is_reference: boolean. If set, the transcript's strand will not be checked.


    :rtype: (None|TranscriptChecker)
    """

    if logger is None:
        logger = create_null_logger()

    if "tid" not in lines:
        logger.error("Lines datastore lacks the transcript ID. Exiting.")
        return None

    try:
        logger.debug("Starting with %s", lines["tid"])
        transcript_line = Transcript()
        transcript_line.chrom = lines["chrom"]
        if "source" in lines:
            transcript_line.source = lines["source"]
        transcript_line.strand = lines["strand"]
        transcript_line.attributes.update(lines["attributes"])
        transcript_line.feature = "transcript"
        transcript_line.start, transcript_line.end = sorted([start, end])
        transcript_line.logger = logger
        assert lines["tid"] is not None, lines
        transcript_line.id = lines["tid"]
        transcript_line.parent = lines["parent"]

        for feature in lines["features"]:
            coords, phases = [], []
            for feat in lines["features"][feature]:
                try:
                    assert isinstance(feat, (list, tuple)) and 2 <= len(feat) <= 3, feat
                except AssertionError:
                    raise exceptions.InvalidTranscript("Invalid feature")
                coords.append((feat[0], feat[1]))
                if len(feat) == 3 and feat[2] in (0, 1, 2, None):
                    phases.append(feat[2])
                else:
                    phases.append(None)
            try:
                assert len(phases) == len(coords)
            except AssertionError:
                raise exceptions.InvalidTranscript("Invalid phases/coords")
            transcript_line.add_exons(coords, features=feature, phases=phases)

        transcript_object = TranscriptChecker(transcript_line,
                                              fasta_seq,
                                              lenient=lenient,
                                              strand_specific=strand_specific,
                                              canonical_splices=canonical_splices,
                                              force_keep_cds=force_keep_cds,
                                              is_reference=is_reference,
                                              logger=logger)
        logger.debug("Finished adding exon lines to %s", lines["tid"])
        transcript_object.finalize()
        transcript_object.check_strand()
        transcript_object.check_orf()
    except exceptions.IncorrectStrandError:
        logger.info("Discarded %s because of incorrect fusions of splice junctions",
                    lines["tid"])
        # logger.exception(exc)
        transcript_object = None
    except exceptions.InvalidTranscript as exc:
        logger.info("Discarded generically invalid transcript %s, exception: %s",
                    lines["tid"], exc)
        transcript_object = None
    except AssertionError as exc:
        logger.info("Validation failed on %s, assertion failure: %s",
                    lines["tid"], exc)
        transcript_object = None
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except Exception as exc:
        logger.exception(exc)
        transcript_object = None

    return transcript_object


class CheckingProcess(multiprocessing.Process):

    def __init__(self,
                 batch_file,
                 logging_queue,
                 fasta,
                 identifier,
                 shelve_stacks,
                 fasta_out,
                 gtf_out,
                 tmpdir,
                 lenient=False,
                 seed=None,
                 force_keep_cds=False,
                 canonical_splices=(("GT", "AG"),
                                    ("GC", "AG"),
                                    ("AT", "AC")),
                 log_level="WARNING"
                 ):

        super().__init__()
        self.__identifier = ""
        self.__set_identifier(identifier)
        if seed is not None:
            # numpy.random.seed(seed % (2 ** 32 - 1))
            random.seed(seed % (2 ** 32 - 1))
        else:
            # numpy.random.seed(None)
            random.seed(None)
        # self.strand_specific = strand_specific
        self.__canonical = []
        self.__set_canonical(canonical_splices)
        self.__log_level = "DEBUG"
        self.log_level = log_level
        self.logger = None  # This gets populated by the create_queue_logger function below
        self.__logging_queue = None
        self.__set_logging_queue(logging_queue)
        self.name = "Checker-{0}".format(self.identifier)
        try:
            create_queue_logger(self)
        except AttributeError as exc:
            raise AttributeError(exc)
        if batch_file is None:
            raise ValueError("Invalid Batch file!")
        self.__lenient = False
        self.lenient = lenient
        self.__fasta = fasta
        self.__submission_queue = None
        self.fasta = pysam.FastaFile(self.__fasta)
        self.fasta_out = os.path.join(tmpdir, "{0}-{1}".format(
            fasta_out, self.identifier
        ))
        self.gtf_out = os.path.join(tmpdir, "{0}-{1}".format(
            gtf_out, self.identifier
        ))
        self.force_keep_cds = force_keep_cds
        self.shelve_stacks = shelve_stacks
        self.batch_file = batch_file

    def _get_stacks(self):
        shelve_stacks = dict()
        for shelf in self.shelve_stacks:
            shelve_stacks[shelf] = open(shelf, "rb")

        return shelve_stacks

    def _get_keys(self):
        keys = msgpack.loads(open(self.batch_file, "rb").read(), raw=False, strict_map_key=False)
        keys = sorted(keys, key=operator.itemgetter(0))
        return keys

    def run(self):
        checker = functools.partial(create_transcript,
                                    # lenient=self.lenient,
                                    force_keep_cds=self.force_keep_cds,
                                    canonical_splices=self.canonical,
                                    logger=self.logger)

        fasta_out = open(self.fasta_out, "w")
        gtf_out = open(self.gtf_out, "w")
        self.logger.debug("Starting %s", self.name)
        self.logger.debug("Created output FASTA {self.fasta_out} and GTF {self.gtf_out}".format(**locals()))
        time.sleep(0.1)
        self.logger.debug(self.canonical)

        __printed = 0
        shelve_stacks = self._get_stacks()
        file_keys = self._get_keys()

        try:
            for key in file_keys:
                counter, keys = key
                # lines, start, end, counter = self.submission_queue.get()
                tid, chrom, (pos) = keys
                try:
                    tid, shelf_name, write_start, write_length = tid
                except ValueError as exc:
                    raise ValueError(f"{exc}\t{tid}")
                start, end = pos
                try:
                    shelf = shelve_stacks[shelf_name]
                except KeyError:
                    exception = f"{shelf_name} not found in shelves, available: {shelve_stacks.keys()}"
                    self.logger.error(exception)
                    raise KeyError(exception)

                shelf.seek(write_start)
                lines = msgpack.loads(zlib.decompress(shelf.read(write_length)))

                self.logger.debug("Checking %s", lines["tid"])
                if "is_reference" not in lines:
                    raise KeyError(lines)

                transcript = checker(lines,
                                     str(self.fasta.fetch(lines["chrom"], start-1, end)),
                                     start,
                                     end,
                                     lenient=self.lenient,
                                     is_reference=lines["is_reference"],
                                     strand_specific=lines["strand_specific"])

                if transcript is None:
                    self.logger.debug("%s failed the check", lines["tid"])
                    continue
                else:
                    self.logger.debug("Printing %s", lines["tid"])
                    __printed += 1
                    print("\n".join(["{0}/{1}".format(counter, line) for line in
                                     transcript.format("gtf").split("\n")]), file=gtf_out)
                    print("\n".join(["{0}/{1}".format(counter, line) for line in
                                     transcript.fasta.split("\n")]), file=fasta_out)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except Exception as exc:
            self.logger.error(exc)
            # self.logging_queue.close()
            raise

        time.sleep(0.1)
        fasta_out.flush()
        fasta_out.close()
        gtf_out.flush()
        gtf_out.close()
        if __printed > 0:
            self.logger.debug("Size of FASTA out and GTF out: %s, %s",
                              os.stat(fasta_out.name).st_size, os.stat(gtf_out.name).st_size)
            assert os.stat(gtf_out.name).st_size > 0
            assert os.stat(fasta_out.name).st_size > 0
        time.sleep(0.1)
        return

    def __getstate__(self):
        state = self.__dict__.copy()
        for key in ("fasta", "logger", "_log_handler"):
            if key in state:
                del state[key]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        create_queue_logger(self)
        self.fasta = pysam.FastaFile(self.__fasta)

    @property
    def identifier(self):
        return self.__identifier

    def __set_identifier(self, identifier):

        if identifier is None:
            raise ValueError("The identifier must be defined!")
        self.__identifier = str(identifier)

    @property
    def log_level(self):
        return self.__log_level

    @log_level.setter
    def log_level(self, log_level):
        _ = logging._checkLevel(log_level)
        self.__log_level = log_level

    @property
    def lenient(self):
        return self.__lenient

    @lenient.setter
    def lenient(self, lenient):
        if lenient not in (False, True):
            raise ValueError("Invalid lenient value: {}".format(lenient))
        self.__lenient = lenient

    # @property
    # def submission_queue(self):
    #     return self.__submission_queue
    #
    # def __set_submission_queue(self, submission):
    #     if isinstance(submission, multiprocessing.queues.SimpleQueue):
    #         if sys.version_info.minor < 6:
    #             raise TypeError("Invalid queue object for Python 3.5 and earlier!")
    #         submission.put_nowait = submission.put
    #     elif not isinstance(submission, (multiprocessing.queues.Queue,
    #                                    queue.Queue)):
    #         raise TypeError("Invalid queue object: {}".format(type(submission)))
    #     self.__submission_queue = submission
    #
    @property
    def logging_queue(self):
        return self.__logging_queue

    def __set_logging_queue(self, logging_queue):
        if isinstance(logging_queue, multiprocessing.queues.SimpleQueue):
            if sys.version_info.minor < 6:
                raise TypeError("Invalid queue object for Python 3.5 and earlier!")
            logging_queue.put_nowait = logging_queue.put
        elif not isinstance(logging_queue, (multiprocessing.queues.Queue, queue.Queue)):
            raise TypeError("Invalid queue object: {}".format(type(logging_queue)))
        self.__logging_queue = logging_queue

    @property
    def canonical(self):
        return self.__canonical

    def __set_canonical(self, canonical):
        if not isinstance(canonical, (tuple, list)):
            raise TypeError("Canonical splices should be lists or tuples")

        if len(canonical) == 0:
            raise ValueError("The list of canonical splices cannot be empty!")

        for el in canonical:
            if (len(el) != 2 or (not (isinstance(el[0], str) and len(el[0]) == 2) or
                                not (isinstance(el[1], str) and len(el[1]) == 2  ))):
                raise ValueError("Invalid splicing pattern!")

        self.__canonical = canonical
