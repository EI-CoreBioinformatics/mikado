# coding: utf-8

"""
The Sublocus class is the first to be invoked during the Mikado pick analysis.
Each of these containers holds transcripts which either are monoexonic and overlapping,
or multiexonic and with at least one intron in common.
"""

import itertools


from ..transcripts.transcript import Transcript
from .abstractlocus import Abstractlocus
from .excluded import Excluded
from .monosublocus import Monosublocus
from ..parsers.GFF import GffLine


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
                 transcript_instance,
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
        self.fixed_size = True if transcript_instance.feature == "sublocus" else False
        if transcript_instance.__name__ == "transcript":
            transcript_instance.finalize()
        self.source = self.json_conf["pick"]["output_format"]["source"]

        self.excluded = None
        self._not_passing = set()
        self.splitted = False
        # Flag to indicate that we have not calculated the metrics for the transcripts
        # Flag to indicate that we have not calculated the scores for the transcripts
        self.scores_calculated = False
        setattr(self, "monoexonic", getattr(transcript_instance, "monoexonic", None))
        if json_conf is None or not isinstance(json_conf, dict):
            raise ValueError("I am missing the configuration for prioritizing transcripts!")

        # This part is necessary to import modules
        if "modules" in self.json_conf:
            import importlib
            for mod in self.json_conf["modules"]:
                globals()[mod] = importlib.import_module(mod)

        if isinstance(transcript_instance, Transcript):
            self.add_transcript_to_locus(transcript_instance)
        else:
            self.parent = getattr(transcript_instance, "parent")
            self.chrom = getattr(transcript_instance, "chrom")
            self.start = getattr(transcript_instance, "start")
            self.end = getattr(transcript_instance, "end")
            self.strand = getattr(transcript_instance, "strand")
            self.attributes = getattr(transcript_instance, "attributes")

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

        if self._excluded_transcripts and self.purge:
            self.excluded = Excluded(self._excluded_transcripts.pop())
            while self._excluded_transcripts:
                self.excluded.add_transcript_to_locus(self._excluded_transcripts.pop(),
                                                      check_in_locus=False)

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
                            ",".join(transcript_graph.neighbors(selected_tid))
                        ))
                to_remove.update(set(transcript_graph.neighbors(selected_tid)))
                # for tid in transcript_graph.neighbors(selected_tid)
                # for clique in cliques:
                #     if selected_tid in clique:
                #         self.logger.debug("Removing as intersecting {0}: {1}".format(
                #             selected_tid,
                #             ",".join(list(clique))
                #         ))
                #         to_remove.update(clique)
                if purge is False or selected_transcript.score > 0:
                    new_locus = Monosublocus(selected_transcript,
                                             logger=self.logger,
                                             json_conf=self.json_conf)
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

    def calculate_scores(self):

        self.metric_lines_store = []
        super().calculate_scores()
        self.metric_lines_store = [_ for _ in self.prepare_metrics()]

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
            return False
        if any(True for comb in itertools.product(transcript.exons, other.exons) if
               cls.overlap(*comb) >= 0):

            return True
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
