# coding: utf-8

"""
The Sublocus class is the first to be invoked during the Mikado pick analysis.
Each of these containers holds transcripts which either are monoexonic and overlapping,
or multiexonic and with at least one intron in common.
"""

import itertools
from .abstractlocus import Abstractlocus
from .excluded import Excluded
from .monosublocus import Monosublocus
from .transcript import Transcript
from ..parsers.GFF import GffLine
from sys import version_info
if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict
import operator


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

    def __init__(self, span, json_conf=None, logger=None):

        """
        :param span: an instance which describes a genomic interval
        :type span: Transcript | Abstractlocus

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
        Abstractlocus.__init__(self)
        self.feature = self.__name__
        self.logger = logger
        self.json_conf = json_conf
        self.fixed_size = True if span.feature == "sublocus" else False
        if span.__name__ == "transcript":
            span.finalize()
        self.purge = self.json_conf["pick"]["run_options"]["purge"]
        self.source = self.json_conf["pick"]["output_format"]["source"]

        self.excluded = None
        self.splitted = False
        # Flag to indicate that we have not calculated the metrics for the transcripts
        self.metrics_calculated = False
        # Flag to indicate that we have not calculated the scores for the transcripts
        self.scores_calculated = False
        setattr(self, "monoexonic", getattr(span, "monoexonic", None))
        if json_conf is None or not isinstance(json_conf, dict):
            raise ValueError("I am missing the configuration for prioritizing transcripts!")
        self.locus_verified_introns = set()

        # This part is necessary to import modules
        if "modules" in self.json_conf:
            import importlib
            for mod in self.json_conf["modules"]:
                globals()[mod] = importlib.import_module(mod)

        if isinstance(span, Transcript):
            self.add_transcript_to_locus(span)
        else:
            self.parent = getattr(span, "parent")
            self.chrom = getattr(span, "chrom")
            self.start = getattr(span, "start")
            self.end = getattr(span, "end")
            self.strand = getattr(span, "strand")
            self.attributes = getattr(span, "attributes")

        self.monosubloci = []
        self.logger.debug("Initialized {0}".format(self.id))
        self.metric_lines_store = []  # This list will contain the lines to be printed in the metrics file

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
        if self.initialized is False:
            self.monoexonic = transcript.monoexonic
        elif self.monoexonic != transcript.monoexonic:
            raise ValueError(
                """Sublocus and transcript are not compatible!
                {0}\t{1}\t{2}\t{3}\t{4}
                {5}""".format(self.chrom,
                              self.start,
                              self.end,
                              self.strand,
                              self.monoexonic,
                              Transcript))

        super().add_transcript_to_locus(transcript)
        # add the verified introns from the outside
        self.locus_verified_introns = set.union(self.locus_verified_introns,
                                                transcript.verified_introns)
        self.logger.debug("Added %s to %s", transcript.id, self.id)
        # Update the id

    def define_monosubloci(self, purge=False, excluded=None):
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
        self.excluded = excluded
        self.logger.debug("Launching calculate scores for {0}".format(self.id))
        self.calculate_scores()

        self.logger.debug("Defining monosubloci for {0}".format(self.id))

        transcript_graph = self.define_graph(self.transcripts,
                                             inters=self.is_intersecting,
                                             logger=self.logger)

        while len(transcript_graph) > 0:
            cliques = self.find_cliques(transcript_graph)
            communities = self.find_communities(transcript_graph)
            self.logger.debug("Cliques: {0}".format(cliques))
            self.logger.debug("Communities: {0}".format(communities))
            to_remove = set()
            for msbl in communities:
                msbl = dict((x, self.transcripts[x]) for x in msbl)
                selected_tid = self.choose_best(msbl)
                selected_transcript = self.transcripts[selected_tid]
                to_remove.add(selected_tid)
                self.logger.debug("Selected: %s (score: %f)",
                                  selected_tid, selected_transcript.score)
                for clique in cliques:
                    if selected_tid in clique:
                        self.logger.debug("Removing as intersecting {0}: {1}".format(
                            selected_tid,
                            ",".join(list(clique))
                        ))
                        to_remove.update(clique)
                if purge is False or selected_transcript.score > 0:
                    new_locus = Monosublocus(selected_transcript, logger=self.logger)
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
        return

    def calculate_metrics(self, tid: str):
        """
        :param tid: the name of the transcript to be analysed
        :type tid: str

        This function will calculate the metrics for a transcript which are relative in nature
        i.e. that depend on the other transcripts in the sublocus. Examples include the fraction
        of introns or exons in the sublocus, or the number/fraction of retained introns.
        """

        self.logger.debug("Calculating metrics for %s", tid)
        # The transcript must be finalized before we can calculate the score.
        self.transcripts[tid].finalize()

        _ = len(set.intersection(self.exons, self.transcripts[tid].exons))
        fraction = _ / len(self.exons)

        self.transcripts[tid].exon_fraction = fraction

        if len(self.introns) > 0:
            _ = len(set.intersection(self.transcripts[tid].introns, self.introns))
            fraction = _ / len(self.introns)
            self.transcripts[tid].intron_fraction = fraction
        else:
            self.transcripts[tid].intron_fraction = 0
        if len(self.selected_cds_introns) > 0:
            intersecting_introns = len(set.intersection(
                self.transcripts[tid].selected_cds_introns,
                set(self.selected_cds_introns)))
            fraction = intersecting_introns / len(self.selected_cds_introns)
            self.transcripts[tid].selected_cds_intron_fraction = fraction
        else:
            self.transcripts[tid].selected_cds_intron_fraction = 0

        if len(self.combined_cds_introns) > 0:
            intersecting_introns = len(
                set.intersection(
                    set(self.transcripts[tid].combined_cds_introns),
                    set(self.combined_cds_introns)))
            fraction = intersecting_introns / len(self.combined_cds_introns)
            self.transcripts[tid].combined_cds_intron_fraction = fraction
        else:
            self.transcripts[tid].combined_cds_intron_fraction = 0

        self.find_retained_introns(self.transcripts[tid])
        assert isinstance(self.transcripts[tid], Transcript)
        retained_bases = sum(e[1] - e[0] + 1
                             for e in self.transcripts[tid].retained_introns)
        fraction = retained_bases / self.transcripts[tid].cdna_length
        self.transcripts[tid].retained_fraction = fraction
        if len(self.locus_verified_introns) > 0:
            verified = len(
                set.intersection(self.transcripts[tid].verified_introns,
                                 self.locus_verified_introns))
            fraction = verified / len(self.locus_verified_introns)
            self.transcripts[tid].proportion_verified_introns_inlocus = fraction
        else:
            self.transcripts[tid].proportion_verified_introns_inlocus = 0
        self.logger.debug("Calculated metrics for {0}".format(tid))

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

    def __check_requirements(self):
        """
        This private method will identify and delete all transcripts which do not pass
        the minimum muster specified in the configuration.
        :return:
        """

        self.json_conf["requirements"]["compiled"] = compile(
            self.json_conf["requirements"]["expression"], "<json>",
            "eval")
        previous_not_passing = set()
        while True:
            not_passing = set()
            for tid in iter(tid for tid in self.transcripts if
                            tid not in previous_not_passing):
                evaluated = dict()
                for key in self.json_conf["requirements"]["parameters"]:
                    value = getattr(self.transcripts[tid],
                                    self.json_conf["requirements"]["parameters"][key]["name"])
                    evaluated[key] = self.evaluate(
                        value,
                        self.json_conf["requirements"]["parameters"][key])
                # pylint: disable=eval-used
                if eval(self.json_conf["requirements"]["compiled"]) is False:
                    not_passing.add(tid)
                # pylint: enable=eval-used
            if len(not_passing) == 0:
                return
            for tid in not_passing:
                self.transcripts[tid].score = 0
                if self.purge is True:
                    self.metrics_calculated = False
                    if self.excluded is None:
                        excluded = Monosublocus(self.transcripts[tid], logger=self.logger)
                        excluded.json_conf = self.json_conf
                        self.excluded = Excluded(excluded)
                    else:
                        self.excluded.add_transcript_to_locus(self.transcripts[tid])
                    self.remove_transcript_from_locus(tid)
            if len(self.transcripts) == 0:
                return
            else:
                # Recalculate the metrics
                self.get_metrics()

    def calculate_scores(self):
        """
        Function to calculate a score for each transcript, given the metrics derived
        with the calculate_metrics method and the scoring scheme provided in the JSON configuration.
        If any requirements have been specified, all transcripts which do not pass them
        will be assigned a score of 0 and subsequently ignored.
        Scores are rounded to the nearest integer.
        """

        if self.scores_calculated is True:
            return

        self.get_metrics()
        not_passing = set()
        if not hasattr(self, "logger"):
            self.logger = None
            self.logger.setLevel("DEBUG")
        self.logger.debug("Calculating scores for {0}".format(self.id))
        if "requirements" in self.json_conf:
            self.__check_requirements()

        if len(self.transcripts) == 0:
            self.logger.warning("No transcripts pass the muster for {0}".format(self.id))
            self.scores_calculated = True
            return
        self.scores = dict()

        for tid in self.transcripts:
            self.scores[tid] = dict()
            # Add the score for the transcript source
            self.scores[tid]["source_score"] = self.transcripts[tid].source_score

        if self.regressor is None:
            for param in self.json_conf["scoring"]:
                self._calculate_score(param)

            for tid in self.scores:
                self.transcripts[tid].scores = self.scores[tid].copy()

            for tid in self.transcripts:
                if tid in not_passing:
                    self.logger.debug("Excluding %s as it does not pass minimum requirements",
                                      tid)
                    self.transcripts[tid].score = 0
                else:
                    self.transcripts[tid].score = sum(self.scores[tid].values())
                    if self.transcripts[tid].score == 0:
                        self.logger.debug("Excluding %s as it has a score of 0", tid)

                if tid not in not_passing:
                    assert self.transcripts[tid].score == sum(self.scores[tid].values()), (
                        tid, self.transcripts[tid].score, sum(self.scores[tid].values())
                    )
                self.scores[tid]["score"] = self.transcripts[tid].score

        else:
            valid_metrics = self.regressor.metrics
            metric_rows = SortedDict()
            for tid, transcript in sorted(self.transcripts.items(), key=operator.itemgetter(0)):
                for param in valid_metrics:
                    self.scores[tid][param] = "NA"
                row = []
                for attr in valid_metrics:
                    val = getattr(transcript, attr)
                    if isinstance(val, bool):
                        if val:
                            val = 1
                        else:
                            val = 0
                    row.append(val)
                metric_rows[tid] = row
            # scores = SortedDict.fromkeys(metric_rows.keys())
            for pos, score in enumerate(self.regressor.predict(list(metric_rows.values()))):
                self.scores[list(metric_rows.keys())[pos]]["score"] = score
                self.transcripts[list(metric_rows.keys())[pos]].score = score

        self.metric_lines_store = [_ for _ in self.prepare_metrics()]
        self.scores_calculated = True

    def _calculate_score(self, param):
        """
        Private method that calculates a score for each transcript,
        given a target parameter.
        :param param:
        :return:
        """

        rescaling = self.json_conf["scoring"][param]["rescaling"]
        metrics = [getattr(self.transcripts[tid], param) for tid in self.transcripts]
        if rescaling == "target":
            target = self.json_conf["scoring"][param]["value"]
            denominator = max(abs(x - target) for x in metrics)
        else:
            target = None
            denominator = (max(metrics) - min(metrics))
        if denominator == 0:
            denominator = 1

        scores = []
        for tid in self.transcripts:
            tid_metric = getattr(self.transcripts[tid], param)
            score = 0
            check = True
            if ("filter" in self.json_conf["scoring"][param] and
                    self.json_conf["scoring"][param]["filter"] != {}):
                check = self.evaluate(tid_metric, self.json_conf["scoring"][param]["filter"])

            if check is True:
                if rescaling == "target":
                    score = 1 - abs(tid_metric - target) / denominator
                else:
                    if min(metrics) == max(metrics):
                        score = 1
                    elif rescaling == "max":
                        score = abs((tid_metric - min(metrics)) / denominator)
                    elif rescaling == "min":
                        score = abs(1 - (tid_metric - min(metrics)) / denominator)

            score *= self.json_conf["scoring"][param]["multiplier"]
            self.scores[tid][param] = round(score, 2)
            scores.append(score)

        # This MUST be true
        if "filter" not in self.json_conf["scoring"][param] and max(scores) == 0:
            self.logger.warning("All transcripts have a score of 0 for %s in %s",
                                param, self.id)

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

        self.calculate_scores()

        for row in self.metric_lines_store:
            yield row

    def print_scores(self):
        """This method yields dictionary rows that are given to a csv.DictWriter class."""
        self.calculate_scores()
        if self.regressor is None:
            score_keys = sorted(list(self.json_conf["scoring"].keys()) + ["source_score"])
        else:
            score_keys = sorted(self.regressor.metrics + ["source_score"])
        keys = ["tid", "parent", "score"] + sorted(score_keys)

        for tid in self.scores:
            row = dict().fromkeys(keys)
            row["tid"] = tid
            row["parent"] = self.id
            row["score"] = round(self.scores[tid]["score"], 2)
            calculate_total = (self.regressor is None)
            for key in score_keys:
                if calculate_total:
                    assert self.scores[tid][key] != "NA" and self.scores[tid][key] is not None
                    row[key] = round(self.scores[tid][key], 2)
            if calculate_total is True:
                score_sum = sum(row[key] for key in score_keys)
                #
                assert round(score_sum, 2) == round(self.scores[tid]["score"], 2), (
                    "Tid: {}; Sum: {}; Calculated: {}\nScores: {}\nRecalculated: {}".format(
                        tid,
                        score_sum,
                        self.transcripts[tid].score,
                        dict(_ for _ in self.scores[tid].items() if
                             self.scores[tid][_[0]] != row[_[0]]),
                        dict(_ for _ in row.items() if
                             self.scores[tid][_[0]] != row[_[0]])
                    ))
            yield row

    def get_metrics(self):

        """Quick wrapper to calculate the metrics for all the transcripts."""

        # TODO: Find an intelligent way ot restoring this check
        # I disabled it because otherwise the values for surviving transcripts would be wrong
        # But this effectively leads to a doubling of run time. A possibility would be to cache the results.
        if self.metrics_calculated is True:
            return

        # self.logger.info("Calculating the intron tree for %s", self.id)
        assert len(self._cds_introntree) == len(self.combined_cds_introns)

        # self.logger.info("Calculated the intron tree for %s, length %d",
        #                  self.id, len(self._cds_introntree))

        for tid in sorted(self.transcripts):
            self.calculate_metrics(tid)

        self.logger.debug("Finished to calculate the metrics for %s", self.id)

        self.metrics_calculated = True
        return

    # ############## Class methods ################

    # Class specific implementation differs from abstract blueprint
    # pylint: disable=arguments-differ
    @classmethod
    def is_intersecting(cls, transcript, other, logger=None):
        """
        :param transcript: the first transcript to check
        :type transcript: Transcript

        :param other: the second transcript to check
        :type other: Transcript

        :param logger: the logger to be used.
        :type logger: (None | logging.Logger)

        Implementation of the is_intersecting method. Here at the level of the sublocus,
        the intersection is seen as overlap between exons.
        """

        if transcript.id == other.id:
            # We do not want intersection with oneself
            if logger is not None:
                logger.debug("Self-comparison for {0}".format(transcript.id))
            return False
        if logger is not None:
            logger.debug("Comparing {0} and {1}".format(transcript.id, other.id))
        if any(True for comb in itertools.product(transcript.exons, other.exons) if
               cls.overlap(*comb) >= 0):
            if logger is not None:
                logger.debug("{0} and {1} are intersecting".format(transcript.id, other.id))

            return True
        if logger is not None:
            logger.debug("{0} and {1} are not intersecting".format(transcript.id, other.id))
        return False
    # pylint: enable=arguments-differ

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
