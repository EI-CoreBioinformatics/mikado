import functools
import multiprocessing
import os

import pyfaidx

from Mikado.transcripts.transcriptchecker import TranscriptChecker
from .. import exceptions
from ..loci import Transcript
from ..utilities.log_utils import create_null_logger, create_queue_logger

__author__ = 'Luca Venturini'


def create_transcript(lines,
                      fasta_seq,
                      start,
                      end,
                      lenient=False,
                      strand_specific=False,
                      canonical_splices=(("GT", "AG"),
                                         ("GC", "AG"),
                                         ("AT", "AC")),
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

    :param logger: optional logger to use during processing.

    :rtype: (None|TranscriptChecker)
    """

    if logger is None:
        logger = create_null_logger("checker")

    logger.debug("Starting with %s", lines["tid"])

    try:
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


class CheckingProcess(multiprocessing.Process):

    def __init__(self,
                 submission_queue,
                 logging_queue,
                 fasta,
                 identifier,
                 fasta_out,
                 gtf_out,
                 tmpdir,
                 lenient=False,
                 # strand_specific=False,
                 canonical_splices=(("GT", "AG"),
                                    ("GC", "AG"),
                                    ("AT", "AC")),
                 log_level="WARNING"
                 ):

        super().__init__()
        self.__identifier = identifier
        # self.strand_specific = strand_specific
        self.canonical = canonical_splices
        self.log_level = log_level
        self.logger = None
        self.logging_queue = logging_queue
        self.name = "Checker-{0}".format(self.identifier)
        create_queue_logger(self)
        self.lenient = lenient
        self.__fasta = fasta
        self.submission_queue = submission_queue
        self.fasta = pyfaidx.Fasta(self.__fasta)
        self.fasta_out = os.path.join(tmpdir, "{0}-{1}".format(
            fasta_out, self.identifier
        ))
        self.gtf_out = os.path.join(tmpdir, "{0}-{1}".format(
            gtf_out, self.identifier
        ))
        self.logger.debug(self.canonical)

    def run(self):

        checker = functools.partial(create_transcript,
                                    lenient=self.lenient,
                                    # strand_specific=self.strand_specific,
                                    canonical_splices=self.canonical,
                                    logger=self.logger)

        fasta_out = open(self.fasta_out, "w")
        gtf_out = open(self.gtf_out, "w")

        while True:
            lines, start, end, counter = self.submission_queue.get()
            if lines == "EXIT":
                self.logger.debug("Finished for %s", self.name)
                self.submission_queue.put((lines,
                                           start,
                                           end,
                                           counter))
                break
            self.logger.debug("Checking %s", lines["tid"])
            transcript = checker(lines,
                                 str(self.fasta[lines["chrom"]][start-1:end]),
                                 start,
                                 end,
                                 strand_specific=lines["strand_specific"])

            if transcript is None:
                self.logger.debug("%s failed the check", lines["tid"])
                continue
            else:
                self.logger.debug("Printing %s", lines["tid"])
                print("\n".join(["{0}/{1}".format(counter, line) for line in
                                 transcript.format("gtf").split("\n")]), file=gtf_out)
                print("\n".join(["{0}/{1}".format(counter, line) for line in
                                 transcript.fasta.split("\n")]), file=fasta_out)

        fasta_out.close()
        gtf_out.close()
        return

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["fasta"]
        # del state["handler"]
        del state["logger"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        create_queue_logger(self)
        self.fasta = pyfaidx.Fasta(self.__fasta)

    @property
    def identifier(self):
        return self.__identifier
