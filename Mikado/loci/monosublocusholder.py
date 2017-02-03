# coding: utf-8

"""
This module defines a holder of Monosublocus instances. It is the last step
before the definition of real loci.
"""

import itertools
import logging
from sys import version_info

from ..transcripts.transcript import Transcript
from .abstractlocus import Abstractlocus
from .locus import Locus
from .monosublocus import Monosublocus
from .sublocus import Sublocus
from ..parsers.GFF import GffLine
from ..scales.contrast import compare as c_compare
from ..utilities import overlap
from ..utilities.log_utils import create_null_logger

if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict

# Resolution order is important here!
# pylint: disable=too-many-instance-attributes
class MonosublocusHolder(Sublocus, Abstractlocus):
    """This is a container that groups together the transcripts
    surviving the selection for the Monosublocus.
    The class inherits from both sublocus and Abstractlocus
    (the main abstract class) in order to be able to reuse
    some of the code present in the former.
    Internally, the most important method is define_loci -
    which will select the best transcript(s) and remove all the overlapping ones.
    The intersection function for this object is quite laxer than in previous stages,
    and so are the requirements for the inclusion.
    """

    __name__ = "monosubloci_holder"

    # pylint: disable=super-init-not-called
    def __init__(self, monosublocus_instance: Monosublocus, json_conf=None, logger=None, verified_introns=None):

        # I know what I am doing by NOT calling the Sublocus super but rather
        # Abstractlocus
        Abstractlocus.__init__(self, verified_introns=verified_introns)
        self.logger = logger
        self.splitted = False
        self.metrics_calculated = False
        self.json_conf = json_conf
        self.excluded = None
        self.purge = self.json_conf["pick"]["run_options"]["purge"]
        self.feature = "MonosublocusHolder"
        self.score = monosublocus_instance.score
        self.scores_calculated = False
        # Add the transcript to the Locus
        self.locus_verified_introns = set()
        self.add_monosublocus(monosublocus_instance)
        self.loci = SortedDict()
        self.attributes = dict()

    # Overriding is correct here
    # pylint: disable=arguments-differ
    def add_transcript_to_locus(self, transcript, check_in_locus=True):
        """Override of the sublocus method, and reversal to the original
        method in the Abstractlocus class.
        The check_in_locus boolean flag is used to decide
        whether to check if the transcript is in the Locus or not.
        This should be set to False for the first transcript, and True afterwards.

        :param transcript: a Transcript instance
        :type transcript: Transcript

        :param check_in_locus: optional flag to pass to the basic method.
        If set to False it disables checks on whether the transcript
        is really in the locus
        :type check_in_locus: bool
        """

        Abstractlocus.add_transcript_to_locus(self, transcript,
                                              check_in_locus=check_in_locus)
        self.locus_verified_introns = set.union(self.locus_verified_introns,
                                                transcript.verified_introns)
    # pylint: enable=arguments-differ

    def add_monosublocus(self, monosublocus_instance: Monosublocus):
        """Wrapper to extract the transcript from the monosubloci and pass it to the constructor.

        :param monosublocus_instance
        :type monosublocus_instance: Monosublocus
        """
        assert len(monosublocus_instance.transcripts) == 1
        if len(self.transcripts) == 0:
            check_in_locus = False
        else:
            check_in_locus = True
        for tid in monosublocus_instance.transcripts:
            self.add_transcript_to_locus(monosublocus_instance.transcripts[tid],
                                         check_in_locus=check_in_locus)

    def __str__(self, print_cds=False, source_in_name=True):
        """This special method is explicitly *not* implemented;
        this Locus object is not meant for printing, only for computation!

        :param print_cds: flag. Ignored.
        """

        lines = []

        self_line = GffLine('')
        for attr in ["chrom", 'feature', 'source', 'start', 'end', 'strand']:
            setattr(self_line, attr, getattr(self, attr))
        self.calculate_scores()
        self.score = max([_.score for _ in self.transcripts.values()])

        self_line.phase, self_line.score = None, self.score
        if source_in_name is True:
            self_line.id = "{0}_{1}".format(self.source, self.id)
        else:
            self_line.id = self.id
        self_line.name = self.name
        self_line.parent = self.parent
        self_line.attributes.update(self.attributes)
        self_line.attributes["multiexonic"] = (not self.monoexonic)
        lines.append(str(self_line))

        for tid in self.transcripts:
            transcript_instance = self.transcripts[tid]
            transcript_instance.source = self.source
            transcript_instance.parent = self_line.id
            self.logger.debug(self.attributes)
            for attribute in self.attributes:
                if attribute not in transcript_instance.attributes:
                    if attribute == "is_fragment" and self.attributes[attribute] is False:
                        continue
                    transcript_instance.attributes[attribute] = self.attributes[attribute]

            lines.append(transcript_instance.format(
                "gff",
                all_orfs=self.json_conf["pick"]["output_format"]["report_all_orfs"],
                with_cds=print_cds).rstrip())

        return "\n".join(lines)

    def define_monosubloci(self, **kwargs):
        """Overriden and set to NotImplemented to avoid cross-calling it when inappropriate."""
        raise NotImplementedError("Monosubloci are the input of this object, not the output.")

    def define_loci(self, purge=False, excluded=None):
        """This is the main function of the class. It is analogous
        to the define_subloci class defined for sublocus objects,
        but it returns "Locus" objects (not "Monosublocus").

        Optional parameters:

        :param purge: flag. If set to True, all loci whose transcripts
        have scores of 0 will be thrown out
        into an Excluded holder.
        :type purge: bool

        :param excluded
        :type excluded: Excluded

        """
        if self.splitted is True:
            return

        self.excluded = excluded

        self.calculate_scores()

        graph = self.define_graph(
            self.transcripts,
            inters=self.is_intersecting,
            logger=self.logger,
            cds_only=self.json_conf["pick"]["clustering"]["cds_only"],
            min_cdna_overlap=self.json_conf["pick"]["clustering"]["min_cdna_overlap"],
            min_cds_overlap=self.json_conf["pick"]["clustering"]["min_cds_overlap"])

        loci = []
        while len(graph) > 0:
            cliques = self.find_cliques(graph)
            communities = self.find_communities(graph)
            to_remove = set()
            for locus_comm in communities:
                locus_comm = dict((x, self.transcripts[x]) for x in locus_comm)
                selected_tid = self.choose_best(locus_comm)
                selected_transcript = self.transcripts[selected_tid]
                to_remove.add(selected_tid)
                for clique in cliques:
                    if selected_tid in clique:
                        to_remove.update(clique)

                if purge is False or selected_transcript.score > 0:
                    new_locus = Locus(selected_transcript, logger=self.logger)
                    loci.append(new_locus)
            self.logger.debug("Removing {0} transcripts from {1}".format(len(to_remove), self.id))
            graph.remove_nodes_from(to_remove)  # Remove nodes from graph, iterate

        for locus in sorted(loci):
            self.loci[locus.id] = locus
            if self.regressor is not None:
                self.loci[locus.id].regressor = self.regressor
        self.splitted = True
        return

    @classmethod
    def is_intersecting(cls,
                        transcript,
                        other,
                        cds_only=False,
                        logger=None,
                        min_cdna_overlap=0.2,
                        min_cds_overlap=0.2) -> bool:
        """
        Implementation of the is_intersecting method. Now that we are comparing transcripts that
        by definition span multiple subloci, we have to be less strict in our definition of what
        counts as an intersection.
        Criteria:
        - the cDNA and CDS overlap is over a user-specified threshold
        OR
        - there is some intronic overlap
        OR
        - one intron of either transcript is completely contained within an exon of the other.

         :param transcript
         :type transcript; Transcript

         :param other:
         :type other: Transcript

         :param cds_only: boolean flag. If set to True, only
         the CDS component of the transcripts will be considered to determine
         whether they are intersecting or not.
         :type cds_only: bool

        :param min_cdna_overlap: float. This is the minimum cDNA overlap for two transcripts to be considered as intersecting,
         even when all other conditions fail.
        :type min_cdna_overlap: float

        :param min_cds_overlap: float. This is the minimum CDS overlap for two transcripts to be considered as intersecting,
         even when all other conditions fail.
        :type min_cds_overlap: float

        :param logger: either None or a logger instance. If None, a null logger will be created.

         :rtype : bool
        """

        if logger is None or not isinstance(logger, logging.Logger):
            logger = create_null_logger("MSH")

        if transcript == other or transcript.id == other.id or transcript.strand != other.strand:
            logger.debug("Cannot intersect with itself (%s vs %s) or a transcript on the other strand (%s and %s)",
                         transcript.id, other.id, transcript.strand, other.strand)
            return False

        logger.debug("Consider only the CDS: %s", cds_only)
        if cds_only is True:
            logger.debug("%s %s, %s %s, so %s",
                         transcript.id, transcript.is_coding,
                         other.id, other.is_coding,
                         "removing the UTRs" if other.is_coding and transcript.is_coding else "not removing the UTRs"
                         )
            if transcript.is_coding and other.is_coding:
                transcript, other = cls.__strip_utr(transcript, other, logger)

        # Calculate the relationship between the transcripts
        comparison, _ = c_compare(other, transcript)

        logger.debug("Starting to check %s and %s" ,transcript.id, other.id)
        if comparison.n_f1[0] == 0:
            # No overlap. Return False
            logger.debug("No genomic overlap between %s and %s. Comparison: %s",
                         transcript.id, other.id, comparison)
            return False  # We do not want intersection with oneself
        elif comparison.j_f1[0] > 0 or comparison.ccode[0] == "h":
            # Simple case: they do intersect!
            logger.debug("%s and %s intersect; class code: %s", transcript.id, other.id, comparison.ccode)
            return True
        else:
            # Is at least one intron completely contained?
            if cls._intron_contained_in_exon(transcript, other) or cls._intron_contained_in_exon(other, transcript):
                logger.debug("Intronic containment within an exon for the comparison %s and %s; intersecting",
                             transcript.id, other.id)
                return True

        cdna_overlap, cds_overlap = cls.__calculate_overlap(transcript, other, comparison, cds_only=cds_only)
        if cdna_overlap >= min_cdna_overlap and cds_overlap >= min_cds_overlap:
            logger.debug("%s and %s have enough of CDS (%s) and cDNA (%s) overlap, intersecting.",
                         transcript.id, other.id, cds_overlap, cdna_overlap)
            return True
        else:
            logger.debug("%s and %s do not have enough of CDS (%s) and cDNA (%s) overlap, not intersecting.",
                         transcript.id, other.id, cds_overlap, cdna_overlap)
            return False

    @staticmethod
    def __strip_utr(transcript, other, logger):
        """Private method to remove the UTRs from both transcripts. Creates deep copies of the original objects,
        to avoid bugs down the line."""
        logger.debug("Considering only the CDS of %s and %s, as they are both coding; stripping the UTR",
                     transcript.id, other.id)
        transcript = transcript.deepcopy()
        transcript.remove_utrs()
        other = other.deepcopy()
        other.remove_utrs()
        logger.debug("New coordinates: %s (%d-%d), %s (%d-%d)",
                     transcript.id, transcript.start, transcript.end,
                     other.id, other.start, other.end)
        return transcript, other

    @staticmethod
    def __calculate_overlap(transcript, other, comparison, cds_only=False) -> (float, float):
        """Private method to return the cDNA overlap and the CDS overlap of two transcripts."""

        cdna_overlap = max(comparison.n_recall[0], comparison.n_prec[0]) / 100
        if cds_only is True and (transcript.is_coding and other.is_coding):
            return cdna_overlap, cdna_overlap
        elif cds_only is False and (transcript.is_coding and other.is_coding):
            cds_transcript = transcript.deepcopy()
            cds_transcript.remove_utrs()
            cds_other = other.deepcopy()
            cds_other.remove_utrs()
            cds_comparison, _ = c_compare(cds_other, cds_transcript)
            cds_overlap = max(cds_comparison.n_recall[0], cds_comparison.n_prec[0]) / 100
            return cdna_overlap, cds_overlap
        elif not (transcript.is_coding and other.is_coding):
            return cdna_overlap, cdna_overlap
        else:
            raise SyntaxError("Unhandled behaviour!")

    @staticmethod
    def _intron_contained_in_exon(transcript: Transcript, other: Transcript) -> bool:

        """Mini-method to assess whether at least one intron of "transcript" is **completely** contained
        within an exon of "other"."""

        return any((overlap(*_) == (_[0][1] - _[0][0])) for _ in itertools.product(transcript.introns, other.exons))

    @classmethod
    def in_locus(cls, monosublocus: Abstractlocus,
                 transcript: Transcript,
                 flank=0,
                 logger=None,
                 cds_only=False,
                 min_cdna_overlap=0.2,
                 min_cds_overlap=0.2) -> bool:

        """This method checks whether a transcript / monosbulocus
        falls inside the Locus coordinates.

        :rtype: bool

        :param monosublocus: Monosublocus instance
        :type monosublocus: Monosublocus

        :param transcript: the transcript to be compared
        :type transcript: Transcript

        :param flank: optional flank argument
        :type flank: int
        """

        if hasattr(transcript, "transcripts"):
            assert len(transcript.transcripts) == 1
            transcript = transcript.transcripts[list(transcript.transcripts.keys())[0]]
            assert hasattr(transcript, "finalize")
        is_in_locus = Abstractlocus.in_locus(monosublocus, transcript, flank=flank)
        if is_in_locus is True:
            is_in_locus = False
            for tran in monosublocus.transcripts:
                tran = monosublocus.transcripts[tran]
                is_in_locus = cls.is_intersecting(tran,
                                                  transcript,
                                                  logger=logger,
                                                  cds_only=cds_only,
                                                  min_cds_overlap=min_cds_overlap,
                                                  min_cdna_overlap=min_cdna_overlap)
                if is_in_locus is True:
                    break
        return is_in_locus

    @property
    def id(self):
        """
        Wrapper for the id method of abstractlocus. Necessary to redefine the name.
        """

        return Abstractlocus.id.fget(self)  # @UndefinedVariable
