# coding: utf-8

"""
The Sublocus class is the first to be invoked during the Mikado pick analysis.
Each of these containers holds transcripts which either are monoexonic and overlapping,
or multiexonic and with at least one intron in common.
"""

import itertools
import logging
from ..transcripts.transcript import Transcript
from .abstractlocus import Abstractlocus
from .excluded import Excluded
from .monosublocus import Monosublocus
from ..parsers.GFF import GffLine
from ..scales.contrast import compare as c_compare
from ..utilities.log_utils import create_null_logger
from ..utilities import overlap


# pylint: disable=too-many-instance-attributes
class Sublocus(Abstractlocus):
    """
    The sublocus class is created either by the superlocus class during
    the subloci definition, or directly using a G(T|F)line-like object.
    It is used to define the monosubloci.
    """

    __name__ = "sublocus"
    available_metrics = Transcript.get_available_metrics()

    # ############### Class special methods ##############

    def __init__(self,
                 transcript_instance=None,
                 json_conf=None,
                 logger=None,
                 verified_introns=None,
                 **kwargs):

        """
        :param transcript_instance: an instance which describes a genomic interval
        :type transcript_instance: Transcript | Abstractlocus

        :param json_conf: a configuration dictionary
        :type json_conf: dict

        :param logger: a logger instance from the logging library
        :type logger: logging.Logger | None

        This class takes as input a "span" feature - e.g. a GffLine or a transcript_instance.
        The span instance should therefore have such attributes as
        chrom, strand, start, end, attributes.
        """

        self.counter = 0  # simple tag for avoiding collisions in the GFF output
        self.__splitted = False
        Abstractlocus.__init__(self, verified_introns=verified_introns,
                               json_conf=json_conf,
                               logger=logger,
                               **kwargs)
        self.feature = self.__name__
        self.logger.debug("Verified introns for %s: %s", self.id, verified_introns)
        self.fixed_size = False
        self.source = self.json_conf["pick"]["output_format"]["source"]

        self.excluded = None
        self._not_passing = set()
        self.splitted = False
        # Flag to indicate that we have not calculated the metrics for the transcripts
        # Flag to indicate that we have not calculated the scores for the transcripts
        setattr(self, "monoexonic", getattr(transcript_instance, "monoexonic", None))
        if transcript_instance is not None:
            if transcript_instance.__name__ == "transcript":
                transcript_instance.finalize()
            if json_conf is None and transcript_instance.json_conf is not None:
                json_conf = transcript_instance.json_conf
            if transcript_instance and transcript_instance.feature == "sublocus":
                self.fixed_size = True

        if json_conf is None or not isinstance(json_conf, dict):
            raise ValueError("I am missing the configuration for prioritizing transcripts!")

        # This part is necessary to import modules
        if "modules" in self.json_conf:
            import importlib
            for mod in self.json_conf["modules"]:
                globals()[mod] = importlib.import_module(mod)

        if isinstance(transcript_instance, Transcript):
            self.add_transcript_to_locus(transcript_instance)
            self.attributes = transcript_instance.attributes.copy()
            self.parent = transcript_instance.parent
        else:
            self.parent = getattr(transcript_instance, "parent", None)
            self.chrom = getattr(transcript_instance, "chrom", None)
            self.start = getattr(transcript_instance, "start", None)
            self.end = getattr(transcript_instance, "end", None)
            self.strand = getattr(transcript_instance, "strand", None)
            self.attributes = getattr(transcript_instance, "attributes", dict())

        self.monosubloci = []
        self.logger.debug("Initialized {0}".format(self.id))
        self.metric_lines_store = []  # This list will contain the lines to be printed in the metrics file
        self.excluded = Excluded(json_conf=json_conf)
        self.scores = dict()

    # pylint: disable=arguments-differ
    def __str__(self, print_cds=True):

        lines = []

        self_line = GffLine('')
        for attr in ["chrom", 'feature', 'source', 'start', 'end', 'strand']:
            setattr(self_line, attr, getattr(self, attr))
        self_line.phase, self_line.score = None, None
        self_line.id = "{0}_{1}".format(self.source, self.id)
        self_line.name = self.name
        self_line.parent = self.parent
        self_line.attributes["multiexonic"] = (not self.monoexonic)
        lines.append(str(self_line))

        for tid in sorted(self.transcripts, key=lambda ttid: self.transcripts[ttid]):
            self.transcripts[tid].source = self.source
            self.transcripts[tid].parent = self_line.id
            lines.append(self.transcripts[tid].format(
                "gff3",
                all_orfs=self.json_conf["pick"]["output_format"]["report_all_orfs"],
                with_cds=print_cds).rstrip())

        return "\n".join(lines)
    # pylint: enable=arguments-differ

    # ########## Class instance methods #####################

    def as_dict(self):
        state = super().as_dict()
        state["monosubloci"] = [_.as_dict() for _ in self.monosubloci]
        state["excluded"] = self.excluded.as_dict()
        return state

    def load_dict(self, state):
        super().load_dict(state)
        self.monosubloci = []
        for stat in state["monosubloci"]:
            s = Monosublocus()
            s.load_dict(stat)
            self.monosubloci.append(s)
        self.excluded = Excluded()
        self.excluded.load_dict(state["excluded"])

    def add_transcript_to_locus(self, transcript: Transcript, **kwargs):

        """
        :param transcript: the transcript which might be putatively added to the Locus.
        :type transcript: Transcript

        :param kwargs: eventual keyword arguments are ignored.

        This is an override of the original method, as at the sublocus stage
        we need to accomplish a couple of things more:
        - check that transcripts added to the sublocus are either all monoexonic or all multiexonic
        - change the id of the transcripts to the new ID
        """

        _ = kwargs  # Ignore any keyword arguments

        if len(self.transcripts) > 0:
            self.logger.debug("Adding %s to %s",
                              transcript.id, self.id)
        else:
            self.logger.debug("Initializing sublocus with %s at %s%s:%d-%d",
                              transcript.id,
                              transcript.chrom,
                              transcript.strand,
                              transcript.start,
                              transcript.end)

        if transcript is None:
            return

        cds_only = self.json_conf["pick"]["clustering"]["cds_only"]

        monoexonic = self._is_transcript_monoexonic(transcript)

        if self.initialized is False:
            self.monoexonic = monoexonic
            self.logger.debug("Locus %s is %s because of %s that is monoexonic: %s",
                              self.id, self.monoexonic, transcript.id, monoexonic)

        elif self.monoexonic != monoexonic:
            raise ValueError("""Sublocus and transcript are not compatible!
                            {0}\t{1}\t{2}\t{3}\t
                            Locus monoexonic: {4}; CDS only: {5}
                            Wrong transcript:
                            {6}
                            Others:
                            {7}
                            """.format(self.chrom, self.start, self.end, self.strand, self.monoexonic, cds_only,
                                       transcript.format("bed12") + "\t{}".format(
                                           self._is_transcript_monoexonic(transcript)),
                                       "\n".join([_.format("bed12") + "\t{}".format(
                                           self._is_transcript_monoexonic(_)) for _ in self.transcripts.values()])))
        super().add_transcript_to_locus(transcript)
        # add the verified introns from the outside
        self.logger.debug("Added %s to %s", transcript.id, self.id)
        # Update the id

    def define_monosubloci(self, purge=False):
        """
        :param purge: a flag which indicates whether loci whose
        best transcript has a score of 0 should be excluded (True) or retained (False)
        :type purge: bool

        :param excluded: the excluded Locus to which transcripts from purged loci will be added to
        :type excluded: None
        :type excluded: Mikado.loci_objects.excluded.Excluded

        This function retrieves the best non-overlapping transcripts inside
        the sublocus, according to the score calculated by
        calculate_scores (explicitly called inside the method).
        The "excluded" keyword must contain either None or
        a MonosublocusHolder object. It is used to contain
        transcripts that must be excluded from the Locus due to unmet requirements.
        """

        self.monosubloci = []
        # self.excluded = excluded
        self.logger.debug("Launching calculate scores for {0}".format(self.id))
        self.filter_and_calculate_scores()

        if len(self._excluded_transcripts) > 0 and self.purge:
            excluded_tids = list(self._excluded_transcripts.keys())
            for excluded_tid in excluded_tids:
                self.excluded.add_transcript_to_locus(self._excluded_transcripts[excluded_tid],
                                                      check_in_locus=False)
                self.remove_transcript_from_locus(excluded_tid)
                del self._excluded_transcripts[excluded_tid]

        self.logger.debug("Defining monosubloci for {0}".format(self.id))

        transcript_graph = self.define_graph(self.transcripts,
                                             inters=self.is_intersecting,
                                             logger=self.logger)

        while len(transcript_graph) > 0:
            # cliques = self.find_cliques(transcript_graph)
            communities = self.find_communities(transcript_graph)
            # self.logger.debug("Cliques: {0}".format(cliques))
            self.logger.debug("Communities: {0}".format(communities))
            to_remove = set()
            for msbl in communities:
                msbl = dict((x, self.transcripts[x]) for x in msbl)
                selected_tid = self.choose_best(msbl)
                selected_transcript = self.transcripts[selected_tid]
                to_remove.add(selected_tid)
                self.logger.debug("Selected: %s (score: %f)",
                                  selected_tid, selected_transcript.score)
                self.logger.debug("Removing as intersecting {0}: {1}".format(
                            selected_tid,
                            ",".join(set(transcript_graph.neighbors(selected_tid)))
                        ))
                to_remove.update(set(transcript_graph.neighbors(selected_tid)))
                if purge is False or selected_transcript.score > 0:
                    new_locus = Monosublocus(selected_transcript,
                                             logger=self.logger,
                                             json_conf=self.json_conf,
                                             use_transcript_scores=self._use_transcript_scores)
                    new_locus.json_conf = self.json_conf
                    self.monosubloci.append(new_locus)
            if len(to_remove) < 1:
                message = "No transcripts to remove from the pool for {0}\n".format(self.id)
                message += "Transcripts remaining: {0}".format(communities)
                exc = ValueError(message)
                self.logger.exception(exc)
                raise exc
            self.logger.debug("Removing %d from transcripts for %s",
                              len(to_remove), self.id)
            transcript_graph.remove_nodes_from(to_remove)
        self.logger.debug("Defined monosubloci for %s", self.id)
        self.splitted = True
        self.logger.debug("Defined monosubloci for %s", self.id)
        return self.excluded

    def load_scores(self, scores):
        """
        :param scores: an external dictionary with scores
        :type scores: dict

        Simple mock function to load scores for the transcripts from an external dictionary.
        """

        for tid in self.transcripts:
            if tid in scores:
                self.transcripts[tid].score = scores[tid]
            else:
                self.transcripts[tid].score = 0

    def prepare_metrics(self):

        """This method prepares the dictionary "rows"
        that will be given to a csv.DictWriter class."""

        for row in Abstractlocus.print_metrics(self):
            yield row
        if self.excluded is not None:
            for row in self.excluded.print_metrics():
                yield row
        return

    def print_metrics(self):

        """
        THis method yields the dictionary "rows" that will be given to a csv.DictWriter class.
        :return:
        """

        self.filter_and_calculate_scores()
        self.metric_lines_store = [_ for _ in self.prepare_metrics()]
        for row in self.metric_lines_store:
            yield row

    def print_scores(self):
        """This method yields dictionary rows that are given to a csv.DictWriter class."""
        self.filter_and_calculate_scores()
        if self.regressor is None:
            score_keys = sorted(list(self.json_conf["scoring"].keys()) + ["source_score"])
        else:
            score_keys = sorted(self.regressor.metrics + ["source_score"])
        keys = ["tid", "alias", "parent", "score"] + sorted(score_keys)

        for tid in self.scores:
            row = dict().fromkeys(keys)
            row["tid"] = tid
            row["alias"] = self.transcripts[tid].alias
            row["parent"] = self.id
            row["score"] = round(self.scores[tid]["score"], 2)
            calculate_total = (self.regressor is None)
            for key in score_keys:
                if calculate_total:
                    assert key in self.scores[tid] and self.scores[tid][key] != "NA" and self.scores[tid][key] is not None, (key, self.scores[tid].keys())
                    row[key] = round(self.scores[tid][key], 2)

            if calculate_total is True and tid not in self._not_passing:
                score_sum = sum(row[key] for key in score_keys)

                if round(score_sum, 2) != round(self.scores[tid]["score"], 2):
                    try:
                        scores = dict(_ for _ in self.scores[tid].items())
                    except KeyError as exc:
                        raise KeyError((exc, row.keys()))
                    recalc = dict(_ for _ in row.items() if _[0] not in ["tid", "parent", "score"])
                    error = AssertionError("Tid: {}; Sum: {}; Calculated: {}\nScores: {}\nRecalculated: {}".format(
                        tid,
                        score_sum,
                        self.transcripts[tid].score,
                        scores,
                        recalc
                    ))
                    self.logger.exception(error)
                    raise error

            yield row

    # ############## Class methods ################

    # Class specific implementation differs from abstract blueprint
    # pylint: disable=arguments-differ
    @classmethod
    def is_intersecting(cls,
                        transcript,
                        other,
                        cds_only=False,
                        logger=None,
                        min_cdna_overlap=0.2,
                        min_cds_overlap=0.2,
                        simple_overlap_for_monoexonic=True) -> bool:
        """
        Implementation of the is_intersecting method. Now that we are comparing transcripts that
        by definition span multiple subloci, we have to be less strict in our definition of what
        counts as an intersection.
        Criteria:
        - the cDNA and CDS overlap is over a user-specified threshold
        OR
        - either transcript is monoexonic and simple_overlap_for_monoexonic is True, if there is exonic overlap
        OR
        - one intron of either transcript is completely contained within an exon of the other.

        The user can specify whether she prefers to consider the whole transcript (default) or whether to consider
        instead the **selected ORF** of the transcripts for the comparison. Please note that intersection in secondary
        ORFs will not be valid under this scenario.

         :param transcript
         :type transcript; Transcript

         :param other:
         :type other: Transcript

         :param cds_only: boolean flag. If set to True, only the CDS component of the transcripts will be
         considered to determine whether they are intersecting or not.
         :type cds_only: bool

        :param min_cdna_overlap: float. This is the minimum cDNA overlap for two transcripts to be considered as intersecting,
         even when all other conditions fail.
        :type min_cdna_overlap: float

        :param min_cds_overlap: float. This is the minimum CDS overlap for two transcripts to be considered as intersecting,
         even when all other conditions fail.
        :type min_cds_overlap: float

        :param simple_overlap_for_monoexonic: boolean flag. If set to true, any overlap for monoexonic transcripts
        will be enough to trigger incorporation in the locus.
        :type simple_overlap_for_monoexonic: bool

        :param logger: either None or a logger instance. If None, a null logger will be created.

         :rtype : bool
        """

        if logger is None or not isinstance(logger, logging.Logger):
            logger = create_null_logger("MSH")
        transcript.finalize()
        other.finalize()

        logger.debug("Comparing %s vs. %s", transcript.id, other.id)

        if transcript.id == other.id or transcript.strand != other.strand:
            logger.debug("Cannot intersect with itself (%s vs %s) or a transcript on the other strand (%s and %s)",
                         transcript.id, other.id, transcript.strand, other.strand)
            return False

        if cds_only is True and transcript.is_coding and other.is_coding:
            logger.debug("Consider only the CDS: %s", cds_only)
            if overlap((transcript._selected_orf_transcript.start, transcript._selected_orf_transcript.end),
                       (other._selected_orf_transcript.start, other._selected_orf_transcript.end)) <= 0:
                intersecting, reason = False, "No genomic overlap between {} and {}".format(
                    transcript.id, other.id
                )
            else:
                intersecting, reason = cls._transcripts_are_intersecting(
                    transcript._selected_orf_transcript,
                    other._selected_orf_transcript,
                    min_cdna_overlap=min_cdna_overlap,
                    min_cds_overlap=min_cds_overlap,
                    simple_overlap_for_monoexonic=simple_overlap_for_monoexonic)
        else:
            if overlap((transcript.start, transcript.end),
                       (other.start, other.end)) <= 0:
                intersecting, reason = False, "No genomic overlap between {} and {}".format(
                    transcript.id, other.id
                )
            else:
                intersecting, reason = cls._transcripts_are_intersecting(
                    transcript,
                    other,
                    min_cdna_overlap=min_cdna_overlap,
                    min_cds_overlap=min_cds_overlap,
                    simple_overlap_for_monoexonic=simple_overlap_for_monoexonic)

        logger.debug(reason)
        return intersecting

    @classmethod
    def _transcripts_are_intersecting(cls,
                                      transcript: Transcript,
                                      other: Transcript,
                                      min_cdna_overlap=0.2,
                                      min_cds_overlap=0.2,
                                      simple_overlap_for_monoexonic=True):
        """Private method which is called by is_intersecting. It decouples the determination of whether two transcripts
        intersect from the public interface of the method.
        :param transcript
         :type transcript; Transcript

         :param other:
         :type other: Transcript

        :param min_cdna_overlap: float. This is the minimum cDNA overlap for
        two transcripts to be considered as intersecting, even when all other conditions fail.
        :type min_cdna_overlap: float

        :param min_cds_overlap: float. This is the minimum CDS overlap for two transcripts to be considered as intersecting,
         even when all other conditions fail.
        :type min_cds_overlap: float

        :param simple_overlap_for_monoexonic: boolean flag. If set to true, any overlap for monoexonic transcripts
        will be enough to trigger incorporation in the locus.
        :type simple_overlap_for_monoexonic: bool

        :param is_internal_orf: boolean. Set to True if we are considering only the CDS for this run.
        :type is_internal_orf: bool
        """

        comparison, _ = c_compare(other, transcript)
        if comparison.n_f1[0] == 0:
            reason = "No genomic overlap between {} and {}".format(transcript.id, other.id)
            intersecting = False
            return intersecting, reason

        if comparison.j_f1[0] > 0 or comparison.ccode[0] == "h":
            reason = "{} and {} intersect; class code: {}".format(transcript.id, other.id, comparison.ccode[0])
            intersecting = True
        elif simple_overlap_for_monoexonic is True and any(_.monoexonic is True for _ in (transcript, other)):
            reason = "Simple overlap for monoexonic transcripts, for {} and {}".format(transcript.id, other.id)
            intersecting = True
        elif cls._intron_contained_in_exon(transcript, other) or cls._intron_contained_in_exon(other, transcript):
            reason = "Intronic containment within an exon for the comparison {} and {}; intersecting".format(
                transcript.id, other.id)
            intersecting = True
        else:
            intersecting, reason = cls._evaluate_transcript_overlap(
                transcript, other,
                min_cdna_overlap=min_cdna_overlap,
                min_cds_overlap=min_cds_overlap,
                comparison=comparison,
                fixed_perspective=False)

        return intersecting, reason

    @staticmethod
    def _intron_contained_in_exon(transcript: Transcript, other: Transcript) -> bool:

        """Mini-method to assess whether at least one intron of "transcript" is **completely** contained
        within an exon of "other"."""

        return any((overlap(*_) == (_[0][1] - _[0][0])) for _ in itertools.product(transcript.introns, other.exons))

    @property
    def splitted(self):
        """The splitted flag indicates whether a sublocus
        has already been processed to produce the necessary monosubloci.
        It must be set as a boolean flag (hence why it is coded as a property)

        :rtype bool
        """
        return self.__splitted

    @splitted.setter
    def splitted(self, verified):
        """
        :param verified: boolean flag to set.
        :type verified: bool
        """

        if not isinstance(verified, bool):
            raise TypeError()
        self.__splitted = verified

    @property
    def id(self):
        """
        :return: The name of the Locus.
        :rtype str
        """
        if self.monoexonic is True:
            addendum = "mono"
        else:
            addendum = "multi"
        if self.counter > 0:
            addendum = "{0}.{1}".format(addendum, self.counter)

        return "{0}.{1}".format(super().id, addendum)
