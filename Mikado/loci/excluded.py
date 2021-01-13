# coding: utf-8

"""
This module defines a containers that hold transcripts excluded from further consideration.
It is invoked when all transcripts in a locus have a score of 0 and the "purge"
option has been enabled.
"""

from .abstractlocus import Abstractlocus
from ..transcripts import Transcript


class Excluded(Abstractlocus):
    """This is a container of discarded transcripts. It is used only for completeness purposes -
    i.e. printing out the discarded transcripts to a separate file.
    """

    __name__ = "excluded_transcripts"

    def __init__(self, monosublocus_instance=None, json_conf=None, logger=None):
        """
        Constructor method

        :param monosublocus_instance:
        :type monosublocus_instance: Mikado.loci_objects.monosublocus.Monosublocus

        :param json_conf: configuration file
        :type json_conf: dict

        :param logger: logger instance
        :type logger: logging.Logger | None
        """

        Abstractlocus.__init__(self)
        self.splitted = False
        self.metrics_calculated = False
        self.json_conf = json_conf
        self.logger = logger
        if isinstance(monosublocus_instance, Transcript):
            Abstractlocus.__init__(self, transcript_instance=monosublocus_instance)
        elif isinstance(monosublocus_instance, Abstractlocus):
            # Add the transcript to the Locus
            self.add_monosublocus(monosublocus_instance)

    def add_transcript_to_locus(self, transcript, **kwargs):
        """Override of the sublocus method, and reversal to the original
        method in the Abstractlocus class.
        :param transcript: a transcript to add
        :type transcript: Mikado.loci_objects.transcript.Transcript

        :param kwargs: optional arguments are completely ignored by this method.
        """

        # Notice that check_in_locus is always set to False.
        _ = kwargs

        Abstractlocus.add_transcript_to_locus(self, transcript, check_in_locus=False)

    def add_monosublocus(self, monosublocus_instance):
        """Wrapper to extract the transcript from the monosubloci and pass it
        to the constructor.
        :param monosublocus_instance
        :type monosublocus_instance: Mikado.loci_objects.monosublocus.Monosublocus
        """
        assert len(monosublocus_instance.transcripts) == 1
        for tid in monosublocus_instance.transcripts:
            self.add_transcript_to_locus(monosublocus_instance.transcripts[tid])

    def __str__(self):
        """This special method is explicitly *not* implemented;
        this Locus object is not meant for printing, only for computation!"""
        message = """This is a container used for computational purposes only,
        it should not be printed out directly!"""
        raise NotImplementedError(message)

    def filter_and_calculate_scores(self):
        """
        Suppress the method from the base class
        """

        raise NotImplementedError("Scores are not calculated by this class!")

    def define_monosubloci(self):
        """
        Suppress the method from the base class
        """

        raise NotImplementedError("Monosubloci are not calculated by this class!!")

    @classmethod
    def is_intersecting(cls):
        """Present to fulfill the contract with Abstractlocus, but
        it only raises a NotImplementedError"""
        raise NotImplementedError()
