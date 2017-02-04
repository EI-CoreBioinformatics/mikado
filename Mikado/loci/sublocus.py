# coding: utf-8

"""
The Sublocus class is the first to be invoked during the Mikado pick analysis.
Each of these containers holds transcripts which either are monoexonic and overlapping,
or multiexonic and with at least one intron in common.
"""

import itertools
from sys import version_info
import numpy
from sklearn.ensemble import RandomForestClassifier
from ..transcripts.transcript import Transcript
from .abstractlocus import Abstractlocus
from .excluded import Excluded
from .monosublocus import Monosublocus
from ..parsers.GFF import GffLine
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

    def __init__(self, span, json_conf=None, logger=None, verified_introns=None):

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
        Abstractlocus.__init__(self, verified_introns=verified_introns)
        self.feature = self.__name__
        self.logger = logger
        self.logger.debug("Verified introns for %s: %s", self.id, verified_introns)
        self.json_conf = json_conf
        self.fixed_size = True if span.feature == "sublocus" else False
        if span.__name__ == "transcript":
            span.finalize()
        self.purge = self.json_conf["pick"]["clustering"]["purge"]
        self.source = self.json_conf["pick"]["output_format"]["source"]

        self.excluded = None
        self.splitted = False
        # Flag to indicate that we have not calculated the metrics for the transcripts
        # Flag to indicate that we have not calculated the scores for the transcripts
        self.scores_calculated = False
        setattr(self, "monoexonic", getattr(span, "monoexonic", None))
        if json_conf is None or not isinstance(json_conf, dict):
            raise ValueError("I am missing the configuration for prioritizing transcripts!")

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

        self.get_metrics()
        if self.purge is False:
            self.logger.debug("No purging for %s, returning", self.id)
            return

        previous_not_passing = set()
        while True:
            not_passing = self._check_not_passing(
                previous_not_passing=previous_not_passing)
            if len(not_passing) == 0:
                return
            for tid in not_passing:
                self.transcripts[tid].score = 0

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
            self.logger.debug("Scores calculation already effectuated for %s",
                              self.id)
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
                # if self.json_conf["pick"]["external_scores"]:
                #     assert any("external" in _ for _ in self.scores[tid].keys()), self.scores[tid].keys()

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
                # Necessary for sklearn ..
                row = numpy.array(row)
                # row = row.reshape(1, -1)
                metric_rows[tid] = row
            # scores = SortedDict.fromkeys(metric_rows.keys())
            if isinstance(self.regressor, RandomForestClassifier):
                # We have to pick the second probability (correct)
                for tid in metric_rows:
                    score = self.regressor.predict_proba(metric_rows[tid])[0][1]
                    self.scores[tid]["score"] = score
                    self.transcripts[tid].score = score
            else:
                pred_scores = self.regressor.predict(list(metric_rows.values()))
                for pos, score in enumerate(pred_scores):
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
        use_raw = self.json_conf["scoring"][param]["use_raw"]

        metrics = dict((tid, getattr(self.transcripts[tid], param)) for tid in self.transcripts)

        if use_raw is True and not param.startswith("external") and getattr(Transcript, param).usable_raw is False:
            self.logger.warning("The \"%s\" metric cannot be used as a raw score for %s, switching to False",
                                param, self.id)
            use_raw = False
        if use_raw is True and rescaling == "target":
            self.logger.warning("I cannot use a raw score for %s in %s when looking for a target. Switching to False",
                                param, self.id)
            use_raw = False

        if rescaling == "target":
            target = self.json_conf["scoring"][param]["value"]
            denominator = max(abs(x - target) for x in metrics.values())
        else:
            target = None
            if use_raw is True and rescaling == "max":
                denominator = 1
            elif use_raw is True and rescaling == "min":
                denominator = -1
            else:
                denominator = (max(metrics.values()) - min(metrics.values()))
        if denominator == 0:
            denominator = 1

        scores = []
        for tid in metrics:
            tid_metric = metrics[tid]
            score = 0
            check = True
            if ("filter" in self.json_conf["scoring"][param] and
                    self.json_conf["scoring"][param]["filter"] != {}):
                check = self.evaluate(tid_metric, self.json_conf["scoring"][param]["filter"])

            if check is True:
                if use_raw is True:
                    if not isinstance(tid_metric, (float, int)) and 0 <= tid_metric <= 1:
                        error = ValueError(
                            "Only scores with values between 0 and 1 can be used raw. Please recheck your values.")
                        self.logger.exception(error)
                        raise error
                    score = tid_metric / denominator
                elif rescaling == "target":
                    score = 1 - abs(tid_metric - target) / denominator
                else:
                    if min(metrics.values()) == max(metrics.values()):
                        score = 1
                    elif rescaling == "max":
                        score = abs((tid_metric - min(metrics.values())) / denominator)
                    elif rescaling == "min":
                        score = abs(1 - (tid_metric - min(metrics.values())) / denominator)

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
                    assert key in self.scores[tid] and self.scores[tid][key] != "NA" and self.scores[tid][key] is not None, (key, self.scores[tid].keys())
                    row[key] = round(self.scores[tid][key], 2)

            if calculate_total is True:
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
