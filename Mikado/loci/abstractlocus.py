# coding: utf-8

"""
Module that defines the blueprint for all loci classes.
"""

import abc
import itertools
import logging
import random
from sys import maxsize
import networkx
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
import numpy
from ..transcripts.clique_methods import find_cliques, find_communities, define_graph
from ..transcripts.transcript import Transcript
from ..configuration.configurator import to_json, check_json
from ..exceptions import NotInLocusError
from ..utilities import overlap, merge_ranges
import operator
from ..utilities.intervaltree import Interval, IntervalTree
from ..utilities.log_utils import create_null_logger
from sys import version_info
from ..scales.contrast import compare as c_compare
if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict


# I do not care that there are too many attributes: this IS a massive class!
# pylint: disable=too-many-instance-attributes,too-many-public-methods
json_conf = to_json(None)


class Abstractlocus(metaclass=abc.ABCMeta):
    """This abstract class defines the basic features of any Locus-like object.
    It also defines methods/properties that are needed throughout the program,
    e.g. the Bron-Kerbosch algorithm for defining cliques, or the find_retained_introns method.
    """

    __name__ = "Abstractlocus"
    available_metrics = Transcript.get_available_metrics()

    # ##### Special methods #########

    __json_conf = json_conf.copy()

    # This dictionary contains the correspondence between transcript attributes
    # and the relevant containers in Locus classes.
    __locus_to_transcript_attrs = {"splices": "splices",
                                   "introns": "introns",
                                   "combined_cds_introns": "combined_cds_introns",
                                   "combined_cds_exons": "combined_cds",
                                   "selected_cds_introns": "selected_cds_introns",
                                   "selected_cds_exons": "selected_cds",
                                   "exons": "exons",
                                   "locus_verified_introns": "verified_introns"}

    @abc.abstractmethod
    def __init__(self,
                 transcript_instance=None,
                 logger=None,
                 source="",
                 verified_introns=None,
                 json_conf=None):

        # Mock values
        self.__source = source

        self.__logger = None
        self.logger = logger
        self.__stranded = False
        self._not_passing = set()
        self._excluded_transcripts = set()
        self.transcripts = dict()
        self.start, self.end, self.strand = maxsize, -maxsize, None
        # self.introns, self.exons, self.splices = set(), set(), set()
        # Consider only the CDS part
        # self.combined_cds_introns, self.selected_cds_introns = set(), set()
        # self.combined_cds_exons, self.selected_cds_exons = set(), set()
        self.stranded = True
        self.initialized = False
        self.monoexonic = True
        self.chrom = None
        self.__locus_verified_introns = set()
        for locattr in self.__locus_to_transcript_attrs:
            setattr(self, locattr, set())
        if verified_introns is not None:
            self.locus_verified_introns = verified_introns

        self.scores_calculated = False
        self.scores = dict()
        self.__cds_introntree = IntervalTree()
        self.__segmenttree = IntervalTree()
        self.__cds_segmenttree = IntervalTree()
        self.__regressor = None
        self.session = None
        self.metrics_calculated = False
        self.__internal_graph = networkx.DiGraph()
        self.json_conf = json_conf
        if transcript_instance is not None and isinstance(transcript_instance, Transcript):
            self.add_transcript_to_locus(transcript_instance)

    @abc.abstractmethod
    def __str__(self, *args, **kwargs):
        raise NotImplementedError("This is an abstract class and it cannot be printed directly!")

    def __repr__(self):

        if len(self.transcripts) > 0:
            transcript_list = ",".join(list(self.transcripts.keys()))
        else:
            transcript_list = "NA"

        return "\t".join([self.__name__,
                          self.chrom,
                          str(self.start),
                          str(self.end),
                          self.strand,
                          transcript_list])

    def __eq__(self, other):
        if not isinstance(self, type(other)):
            return False
        for feature in ["chrom", "strand", "start",
                        "end", "exons", "introns",
                        "splices", "stranded"]:
            if getattr(self, feature) != getattr(other, feature):
                return False
        return True

    def __hash__(self):
        """This has to be defined, otherwise abstractloci objects won't be hashable
        (and therefore operations like adding to sets will be forbidden)"""
        return super().__hash__()

    def __len__(self):
        return self.end - self.start + 1

    def __lt__(self, other):
        if self.strand != other.strand or self.chrom != other.chrom:
            return False
        if self == other:
            return False
        if self.start < other.start:
            return True
        elif self.start == other.start and self.end < other.end:
            return True
        return False

    def __gt__(self, other):
        return not self < other

    def __le__(self, other):
        return (self == other) or (self < other)

    def __ge__(self, other):
        return (self == other) or (self > other)

    def __getstate__(self):
        """Method to allow serialisation - we remove the byte-compiled eval expression."""

        logger = self.logger
        del self.logger
        state = self.__dict__.copy()
        self.logger = logger

        if hasattr(self, "json_conf"):
            # This removes unpicklable compiled attributes, eg in "requirements" or "as_requirements"
            for key in self.json_conf:
                if (isinstance(self.json_conf[key], dict) and
                        self.json_conf[key].get("compiled", None) is not None):
                    del state["json_conf"][key]["compiled"]

        if hasattr(self, "session"):
            if self.session is not None:
                self.session.expunge_all()
                state["session"].expunge_all()
            state["sessionmaker"] = None
            state["session"] = None

        if hasattr(self, "engine"):
            del state["engine"]

        return state

    def __setstate__(self, state):
        """Method to recreate the object after serialisation."""
        self.__dict__.update(state)

        if hasattr(self, "json_conf"):
            if "requirements" in self.json_conf and "expression" in self.json_conf["requirements"]:
                self.json_conf["requirements"]["compiled"] = compile(
                    self.json_conf["requirements"]["expression"],
                    "<json>", "eval")
        # Set the logger to NullHandler
        self.logger = None

    def __iter__(self):
        return iter(self.transcripts.keys())

    def __getitem__(self, item):

        return self.transcripts[item]

    # #### Static methods #######

    @staticmethod
    def overlap(first_interval: (int, int),
                second_interval: (int, int),
                flank=0,
                positive=False) -> int:

        """:param first_interval: a tuple of integers
        :type first_interval: [int,int]

        :param second_interval: a tuple of integers
        :type second_interval: [int,int | intervaltree.Interval]

        :param flank: an optional extending parameter to check for neighbours
        :type flank: int

        :param positive: if True, negative overlaps will return 0. Otherwise, the negative overlap is returned.
        :type positive: bool

        This static method returns the overlap between two intervals.

        Values<=0 indicate no overlap.

        The optional "flank" argument (default 0) allows to expand a locus
        upstream and downstream.
        As a static method, it can be used also outside of any instance -
        "abstractlocus.overlap()" will function.
        Input: two 2-tuples of integers.
        """

        return overlap(first_interval, second_interval, flank, positive=positive)

    @staticmethod
    def evaluate(param: str, conf: dict) -> bool:

        """
        :param param: string to be checked according to the expression in the configuration
        :type param: str

        :param conf: a dictionary containing the expressions to evaluate
        :type conf: dict

        This static method evaluates a single parameter using the requested
        operation from the JSON dict file.
        """

        if conf["operator"] == "eq":
            comparison = (float(param) == float(conf["value"]))
        elif conf["operator"] == "ne":
            comparison = (float(param) != float(conf["value"]))
        elif conf["operator"] == "gt":
            comparison = (float(param) > float(conf["value"]))
        elif conf["operator"] == "lt":
            comparison = (float(param) < float(conf["value"]))
        elif conf["operator"] == "ge":
            comparison = (float(param) >= float(conf["value"]))
        elif conf["operator"] == "le":
            comparison = (float(param) <= float(conf["value"]))
        elif conf["operator"] == "in":
            comparison = (param in conf["value"])
        elif conf["operator"] == "not in":
            comparison = (param not in conf["value"])
        elif conf["operator"] == "within":
            comparison = (param in range(*sorted([conf["value"][0], conf["value"][1] + 1])))
        elif conf["operator"] == "not within":
            comparison = (param not in range(*sorted([conf["value"][0], conf["value"][1] + 1])))
        else:
            raise ValueError("Unknown operator: {0}".format(conf["operator"]))
        return comparison

    # #### Class methods ########

    @classmethod
    def in_locus(cls, locus_instance, transcript, flank=0) -> bool:
        """
        :param locus_instance: an inheritor of this class
        :param transcript: a transcript instance

        :param flank: an optional extending parameter to check for neighbours
        :type flank: int

        Function to determine whether a transcript should be added or not to the locus_instance.
        This is a class method, i.e. it can be used also unbound from any
        specific instance of the class.
        It will be possible therefore to use it to compare any locus_instance to any transcript.
        Arguments:
        - a "locus_instance" object
        - a "transcript" object (it must possess the "finalize" method)
        - flank - optional keyword"""

        transcript.finalize()
        # We want to check for the strand only if we are considering the strand
        if locus_instance is None:
            return False

        if locus_instance.chrom == transcript.chrom:
            if locus_instance.stranded is False or locus_instance.strand == transcript.strand:
                lbound = (locus_instance.start, locus_instance.end)
                tbound = (transcript.start, transcript.end)
                if cls.overlap(lbound, tbound, flank=flank) > 0:
                    return True
        return False

    def define_graph(self, objects: dict, inters=None, **kwargs) -> networkx.Graph:
        """
        :param objects: a dictionary of objects to be grouped into a graph
        :type objects: dict

        :param inters: the intersecting function to be used to define the graph
        :type inters: callable

        :param kwargs: optional arguments to be passed to the inters function
        :type kwargs: dict

        This function will compute the graph which will later be used by find_communities.
        The method takes as mandatory inputs the following:
            - "objects" a dictionary of objects that form the graph
            - "inters" a function/method that determines whether two objects are connected or not.

        It will then return a graph.
        The method accepts also kwargs that can be passed to the inters function.
        WARNING: the kwargs option is really stupid and does not check
        for correctness of the arguments!
        """

        if inters is None:
            inters = self.is_intersecting

        return define_graph(objects, inters, **kwargs)

    def find_communities(self, graph: networkx.Graph) -> list:
        """

        :param graph: a Graph instance from networkx
        :type graph: networkx.Graph

        This function is a wrapper around the networkX methods to find
        cliques and communities inside a graph.
        The method takes as input a precomputed graph and returns
        two lists:
            - cliques
            - communities
        """

        return find_communities(graph, self.logger)

    def find_cliques(self, graph: networkx.Graph) -> (networkx.Graph, list):
        """

        :param graph: graph to which it is necessary to call the cliques for.

        Wrapper for the BronKerbosch algorithm, which returns the maximal cliques in the graph.
        It is the new interface for the BronKerbosch function, which is not called directly
        from outside this class any longer.
        The "inters" keyword provides the function used to determine
        whether two vertices are connected or not in the graph.
        """

        return find_cliques(graph, self.logger)

    @classmethod
    def choose_best(cls, transcripts: dict) -> str:
        """
        :param transcripts: the dictionary of transcripts of the instance
        :type transcripts: dict

        Given a transcript dictionary, this function will choose the one with the highest score.
        If multiple transcripts have exactly the same score, one will be chosen randomly.

        """

        # Choose one transcript randomly between those that have the maximum score
        max_score = max(transcripts.values(),
                        key=operator.attrgetter("score")).score
        return random.choice(
            [transc for transc in transcripts if transcripts[transc].score == max_score])

    # ###### Class instance methods  #######

    def add_transcript_to_locus(self, transcript, check_in_locus=True, **kwargs):
        """
        :param transcript
        :type transcript: Mikado.loci_objects.transcript.Transcript

        :param check_in_locus: flag to indicate whether the function
        should check the transcript before adding it
        or instead whether to trust the assignment to be correct
        :type check_in_locus: bool

        This method checks that a transcript is contained within the superlocus
        (using the "in_superlocus" class method)
        and upon a successful check extends the superlocus with the new transcript.
        More precisely, it updates the boundaries (start and end) it adds the transcript
        to the internal "transcripts" store, and finally it extends
        the splices and introns with those found inside the transcript.
        """

        transcript.finalize()
        self.monoexonic = self.monoexonic and transcript.monoexonic

        if self.initialized is True:
            if check_in_locus is False:
                pass
            elif not self.in_locus(self, transcript, **kwargs):
                raise NotInLocusError("""Trying to merge a Locus with an incompatible transcript!
                Locus: {lchrom}:{lstart}-{lend} {lstrand} [{stids}]
                Transcript: {tchrom}:{tstart}-{tend} {tstrand} {tid}
                """.format(
                    lchrom=self.chrom, lstart=self.start, lend=self.end, lstrand=self.strand,
                    tchrom=transcript.chrom,
                    tstart=transcript.start,
                    tend=transcript.end,
                    tstrand=transcript.strand,
                    tid=transcript.id,
                    stids=", ".join(list(self.transcripts.keys())),

                ))
        else:
            self.strand = transcript.strand
            self.chrom = transcript.chrom

        self.start = min(self.start, transcript.start)
        self.end = max(self.end, transcript.end)

        self.transcripts[transcript.id] = transcript

        for locattr, tranattr in self.__locus_to_transcript_attrs.items():
            getattr(self, locattr).update(set(getattr(transcript, tranattr)))

        if transcript.monoexonic is False:
            assert len(self.introns) > 0

        self.add_path_to_graph(transcript, self._internal_graph)
        assert len(transcript.combined_cds) <= len(self.combined_cds_exons)
        assert len(self.locus_verified_introns) >= transcript.verified_introns_num

        self.initialized = True
        self.metrics_calculated = False
        self.scores_calculated = False

        return

    def remove_transcript_from_locus(self, tid: str):
        """
        :param tid: name of the transcript to remove
        :type tid: str

         This method will remove a transcript from an Abstractlocus-like
         instance and reset appropriately all derived attributes
         (e.g. introns, start, end, etc.).
        """

        if tid not in self.transcripts:
            self.logger.warning("Transcript %s is not present in the Locus. Ignoring it.", tid)
            return

        self.remove_path_from_graph(self.transcripts[tid], self._internal_graph)
        del self.transcripts[tid]
        for locattr, tranattr in self.__locus_to_transcript_attrs.items():
            setattr(self, locattr,
                    set.union(*[set()] + [  # The [set()] is necessary to prevent a crash from set.union for empty lists
                        set(getattr(self.transcripts[_], tranattr)) for _ in self.transcripts]
                              ))

        if self.transcripts:
            for tid in self.transcripts:
                self.transcripts[tid].parent = self.id
            self.end = max(self.transcripts[_].end for _ in self.transcripts)
            self.start = min(self.transcripts[_].start for _ in self.transcripts)
        else:
            self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
            self.stranded = False
            self.initialized = False

        self.logger.debug("Deleted %s from %s", tid, self.id)
        self.metrics_calculated = False
        self.scores_calculated = False

    @staticmethod
    def _exon_to_be_considered(exon,
                               transcript,
                               consider_truncated=False):
        """Private static method to evaluate whether an exon should be considered for being a retained intron.

        :param exon: the exon to be considered.
        :type exon: (tuple|Interval)

        :param transcript: the candidate transcript from which the exon comes from.
        :type transcript: Transcript

        :param consider_truncated: boolean flag. If set, also terminal exons can be considered for retained intron
        events.
        :type consider_truncated: bool

        :returns: boolean flag (True if it has to be considered, False otherwise), a list of sections of the exons
        which are non-coding, and a boolean flag indicating whether the exon is terminal or not
        """

        cds_segments = sorted(transcript.cds_tree.search(*exon))
        terminal = bool(set.intersection(
            set(exon),
            {transcript.start, transcript.end, transcript.combined_cds_end, transcript.combined_cds_start}))
        if cds_segments == [Interval(*exon)]:
            # It is completely coding
            if terminal is False or (not consider_truncated):
                to_consider = False
                frags = []
            else:
                to_consider = True
                frags = cds_segments
        else:
            frags = []
            to_consider = True
            if cds_segments:
                if cds_segments[0].start > exon[0]:
                    frags.append((exon[0], cds_segments[0].start - 1))
                for before, after in zip(cds_segments[:-1], cds_segments[1:]):
                    frags.append((before.end + 1, max(after.start - 1, before.end + 1)))
                if cds_segments[-1].end < exon[1]:
                    frags.append((cds_segments[-1].end + 1, exon[1]))
            else:
                frags = [Interval(*exon)]

        return to_consider, frags, terminal

    @staticmethod
    def _is_exon_retained(exon: tuple,
                          segmenttree: IntervalTree,
                          digraph: networkx.DiGraph,
                          frags: list,
                          consider_truncated=False,
                          terminal=False,
                          logger=create_null_logger()):

        """Private static method to verify whether a given exon is a retained intron in the current locus.
        The method is as follows:
        - Find all introns which are intersecting the candidate exon, within the exon/intron segment tree.
          Use the CDS graph for coding transcripts, the global graph otherwise.
        - For each candidate intron, find its parent exons that have a positive overlap with the candidate exon.
        - Cases:
           - If no overlapping parent exon is found, the exon is not retained.
           - If parent exons are found only on one side, the exon can be considered as "retained" only if
             if it is a *terminal* exons and we are considering truncations as indications of retained introns.
           - If overlapping exons are found on both sides, consider it as retained only if the non-coding part(s)
             of the query exon fall within the intron.

        :param exon: the exon to be considered.
        :type exon: (tuple|Interval)

        :param frags: a list of intervals that are non-coding within the exon.
        :type frags: list[(tuple|Interval)]

        :param consider_truncated: boolean flag. If set, also terminal exons can be considered for retained intron
        events.
        :type consider_truncated: bool

        :param terminal: whether the exon is at the 3' end.
        :type terminal: bool

        :rtype: bool
        """

        is_retained = False

        found_introns = set(
            [_._as_tuple() for _ in segmenttree.find(exon[0], exon[1], strict=not consider_truncated, value="intron")]
        )

        if not found_introns:
            return is_retained

        logger.debug("Found introns for %s: %s", exon, found_introns)

        for intron in found_introns:
            if is_retained:
                break
            # Only consider exons for which there is an overlap.
            before = {_ for _ in networkx.ancestors(digraph, intron) if overlap(_, exon) > 0}

            after = {_ for _ in networkx.descendants(digraph, intron) if overlap(_, exon) > 0}

            # Now we have to check whether the matched introns contain both coding and non-coding parts
            # Let us exclude any intron which is outside of the exonic span of interest.
            logger.debug("Exon: %s; Frags: %s; Intron: %s; Before: %s; After: %s", exon, frags, intron, before, after)
            if len(before) == 0 and len(after) == 0:
                # A retained intron must be overlapping some other exons!
                logger.debug("No intersecting exon before or after, False")
                is_retained = False
            elif len(before) == 0 or len(after) == 0:
                is_retained = (consider_truncated and terminal)
                logger.debug("Only one exon intersecting, %s", is_retained)
            else:
                # The *only* case where we do *not* consider it a retained intron is
                # if the CDS spans completely the intron.
                for frag, intron in itertools.product(frags, [intron]):
                    is_retained = (overlap(frag, intron, positive=True) > 0)
                    if is_retained:
                        logger.debug("Frag %s intersecting intron %s", frag, intron)
                        break

        return is_retained

    def find_retained_introns(self, transcript: Transcript):

        """This method checks the number of exons that are possibly retained
        introns for a given transcript.
        A retained intron is defined as an exon which:
        - spans completely an intron of another model *between coding exons*
        - is not completely coding itself
        - has at least part of its non-coding sections *within the intron*

        If the "pick/run_options/consider_truncated_for_retained" flag in the configuration is set to true,
         an exon will be considered as a retained intron event also if:
         - it is the last exon of the transcript
         - it ends *within* an intron of another model *between coding exons*
         - is not completely coding itself
         - if the model is coding, the exon has *part* of the non-coding section lying inside the intron
         (ie the non-coding section must not be starting in the exonic part).

        The results are stored inside the transcript instance, in the "retained_introns" tuple.

        If the transcript is non-coding, then all intron retentions are considered as valid: even those spanning introns
        across non-coding exons of a coding model.

        :param transcript: a Transcript instance
        :type transcript: Transcript
        :returns : None
        :rtype : None"""

        self.logger.debug("Starting to calculate retained introns for %s", transcript.id)
        if len(self.introns) == 0:
            transcript.retained_introns = tuple()
            self.logger.debug("No introns in the locus to check against. Exiting.")
            return
        transcript.logger = self.logger
        transcript.finalize()

        # A retained intron is defined as an exon which
        # - is not completely coding
        # - EITHER spans completely the intron of another transcript.
        # - OR is the last exon of the transcript and it ends within the intron of another transcript

        retained_introns = []
        consider_truncated = self.json_conf["pick"]["run_options"]["consider_truncated_for_retained"]

        for exon in transcript.exons:

            self.logger.debug("Checking exon %s of %s", exon, transcript.id)
            # is_retained = False
            to_consider, frags, terminal = self._exon_to_be_considered(
                exon, transcript, consider_truncated=consider_truncated)
            if not to_consider:
                self.logger.debug("Exon %s of %s is not to be considered", exon, transcript.id)
                continue

            self.logger.debug("Number of exons, introns, intervals in segmenttree: %d, %d, %d",
                              len(self.exons), len(self.introns), len(self.segmenttree))

            is_retained = self._is_exon_retained(
                exon,
                self.segmenttree,
                self._internal_graph,
                frags,
                consider_truncated=consider_truncated,
                terminal=terminal,
                logger=self.logger)

            if is_retained:
                self.logger.debug("Exon %s of %s is a retained intron",
                                  exon, transcript.id)
                retained_introns.append(exon)

        self.logger.debug("%s has %d retained introns%s",
                          transcript.id,
                          len(retained_introns),
                          " ({})".format(retained_introns) if retained_introns else "")

        transcript.retained_introns = tuple(sorted(retained_introns))
        return

    @staticmethod
    def _evaluate_transcript_overlap(
            transcript,
            other,
            min_cdna_overlap=0.2,
            min_cds_overlap=0.2,
            comparison=None,
            strict_cds_overlap=False,
            is_internal_orf=False):

        """This private static method evaluates whether the cDNA and CDS overlap of two transcripts
        is enough to consider them as intersecting.

         :param transcript
         :type transcript; Transcript

         :param other:
         :type other: Transcript

        :param min_cdna_overlap: float. This is the minimum cDNA overlap for two transcripts to be considered as intersecting,
         even when all other conditions fail.
        :type min_cdna_overlap: float

        :param min_cds_overlap: float. This is the minimum CDS overlap for two transcripts to be considered as intersecting,
         even when all other conditions fail.
        :type min_cds_overlap: float

        :param is_internal_orf: boolean. Set to True if we are considering only the CDS for this run.
        :type is_internal_orf: bool
        """

        if comparison is None:
            comparison, _ = c_compare(other, transcript)

        cdna_overlap = max(comparison.n_prec[0], comparison.n_recall[0]) / 100
        if strict_cds_overlap is False and (
                        is_internal_orf is True or not (transcript.is_coding and other.is_coding)):
            cds_overlap = cdna_overlap
        else:
            cds_overlap = 0
            for segment in transcript.selected_cds:
                for o_segment in other.selected_cds:
                    cds_overlap += Abstractlocus.overlap(segment, o_segment, positive=True, flank=0)
            cds_overlap /= min(transcript.selected_cds_length, other.selected_cds_length)
            assert cds_overlap <= 1
        intersecting = (cdna_overlap >= min_cdna_overlap and cds_overlap >= min_cds_overlap)
        reason = "{} and {} {}share enough cDNA ({}%, min. {}%) and CDS ({}%, min. {}%), {}intersecting".format(
            transcript.id, other.id,
            "do not " if not intersecting else "",
            cdna_overlap * 100, min_cdna_overlap * 100,
            cds_overlap * 100, min_cds_overlap * 100,
            "not " if not intersecting else "")
        return intersecting, reason

    def print_metrics(self):

        """This method yields dictionary "rows" that will be given to a csv.DictWriter class."""

        # Check that rower is an instance of the csv.DictWriter class

        # The rower is an instance of the DictWriter class from the standard CSV module

        for tid in sorted(self.transcripts.keys(), key=lambda ttid: self.transcripts[ttid]):
            row = {}
            assert self.available_metrics != []
            for key in self.available_metrics:
                if key.lower() in ("id", "tid"):
                    row[key] = tid
                elif key.lower() == "parent":
                    row[key] = self.id
                else:
                    row[key] = getattr(self.transcripts[tid], key, "NA")
                if isinstance(row[key], float):
                    row[key] = round(row[key], 2)
                elif row[key] is None or row[key] == "":
                    row[key] = "NA"
            for source in self.transcripts[tid].external_scores:
                # Each score from external files also contains a multiplier.
                row["external.{}".format(source)] = self.transcripts[tid].external_scores.get(source)

            assert row != {}
            yield row

        return

    def get_metrics(self):

        """Quick wrapper to calculate the metrics for all the transcripts."""

        # if self.metrics_calculated is True:
        #     return

        cds_bases = sum(_[1] - _[0] + 1 for _ in merge_ranges(
            itertools.chain(*[
                self.transcripts[_].combined_cds for _ in self.transcripts
                if self.transcripts[_].combined_cds])))

        selected_bases = sum(_[1] - _[0] + 1 for _ in merge_ranges(
            itertools.chain(*[
                self.transcripts[_].selected_cds for _ in self.transcripts
                if self.transcripts[_].selected_cds])))

        for tid in self.transcripts:
            if cds_bases == 0:
                self.transcripts[tid].combined_cds_locus_fraction = 0
                self.transcripts[tid].selected_cds_locus_fraction = 0
            else:
                selected_length = self.transcripts[tid].selected_cds_length
                combined_length = self.transcripts[tid].combined_cds_length

                self.transcripts[tid].combined_cds_locus_fraction = combined_length / cds_bases
                self.transcripts[tid].selected_cds_locus_fraction = selected_length / selected_bases

        for tid in sorted(self.transcripts):
            self.calculate_metrics(tid)

        self.logger.debug("Finished to calculate the metrics for %s", self.id)

        self.metrics_calculated = True
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
        # self.transcripts[tid].finalize()

        if len(self.locus_verified_introns) == 0 and self.transcripts[tid].verified_introns_num > 0:
            raise ValueError("Locus {} has 0 verified introns, but its transcript {} has {}!".format(
                self.id, tid, len(self.transcripts[tid].verified_introns)))

        if len(self.locus_verified_introns) > 0:
            verified = len(
                set.intersection(self.transcripts[tid].verified_introns,
                                 self.locus_verified_introns))
            fraction = verified / len(self.locus_verified_introns)

            self.transcripts[tid].proportion_verified_introns_inlocus = fraction
        else:
            self.transcripts[tid].proportion_verified_introns_inlocus = 1

        assert (not (len(self.locus_verified_introns) and self.transcripts[tid].verified_introns_num > 0) or
                (self.transcripts[tid].proportion_verified_introns_inlocus > 0))

        # if len(self.locus_verified_introns) and self.transcripts[tid].verified_introns_num > 0:
        #     assert self.transcripts[tid].proportion_verified_introns_inlocus > 0

        _ = len(set.intersection(self.exons, self.transcripts[tid].exons))
        fraction = _ / len(self.exons)

        self.transcripts[tid].exon_fraction = fraction
        self.logger.debug("Calculated exon fraction for %s", tid)

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

        self.logger.debug("Calculating CDS intron fractions for %s", tid)

        if len(self.combined_cds_introns) > 0:
            intersecting_introns = len(
                set.intersection(
                    set(self.transcripts[tid].combined_cds_introns),
                    set(self.combined_cds_introns)))
            fraction = intersecting_introns / len(self.combined_cds_introns)
            self.transcripts[tid].combined_cds_intron_fraction = fraction
        else:
            self.transcripts[tid].combined_cds_intron_fraction = 0

        self.logger.debug("Starting to calculate retained introns for %s", tid)
        self.find_retained_introns(self.transcripts[tid])
        assert isinstance(self.transcripts[tid], Transcript)
        retained_bases = sum(e[1] - e[0] + 1
                             for e in self.transcripts[tid].retained_introns)
        fraction = retained_bases / self.transcripts[tid].cdna_length
        self.transcripts[tid].retained_fraction = fraction

        self.logger.debug("Calculated metrics for {0}".format(tid))

    def _check_not_passing(self, previous_not_passing=set()):
        """
        This private method will identify all transcripts which do not pass
        the minimum muster specified in the configuration. It will *not* delete them;
        how to deal with them is left to the specifics of the subclass.
        :return:
        """

        self.get_metrics()
        # self.logger.debug("Expression: %s", self.json_conf["requirements"]["expression"])

        if ("compiled" not in self.json_conf["requirements"] or
                self.json_conf["requirements"]["compiled"] is None):
            self.json_conf["requirements"]["compiled"] = compile(
                self.json_conf["requirements"]["expression"], "<json>",
                "eval")

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
        self.logger.debug("The following transcripts in %s did not pass the minimum check for requirements: %s",
                          self.id, ", ".join(list(not_passing)))

        return not_passing

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
        # not_passing = set()
        if not hasattr(self, "logger"):
            self.logger = None
            self.logger.setLevel("DEBUG")
        self.logger.debug("Calculating scores for {0}".format(self.id))
        if "requirements" in self.json_conf:
            self._check_requirements()

        if len(self.transcripts) == 0:
            self.logger.warning("No transcripts pass the muster for %s (requirements:\n%s)",
                                self.id,
                                self.json_conf["requirements"])
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
                if tid in self._not_passing:
                    self.logger.debug("Excluding %s as it does not pass minimum requirements",
                                      tid)
                    self.transcripts[tid].score = 0
                else:
                    self.transcripts[tid].score = sum(self.scores[tid].values())
                    if self.transcripts[tid].score <= 0:
                        self.logger.debug("Excluding %s as it has a score <= 0", tid)
                        self.transcripts[tid].score = 0
                        self._not_passing.add(tid)

                if tid in self._not_passing:
                    pass
                else:
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

        self.scores_calculated = True

    def _check_requirements(self):
        """
        This private method will identify and delete all transcripts which do not pass
        the minimum muster specified in the configuration.
        :return:
        """

        self.get_metrics()

        previous_not_passing = set()
        beginning = len(self.transcripts)
        while True:
            not_passing = self._check_not_passing(
                previous_not_passing=previous_not_passing)

            if len(not_passing) == 0:
                self.metrics_calculated = True
                return
            self.metrics_calculated = not ((len(not_passing) > 0) and self.purge)
            self._not_passing.update(not_passing)
            for tid in not_passing:
                if self.purge is False:
                    self.logger.debug("%s has been assigned a score of 0 because it fails basic requirements",
                                      self.id)
                    self.transcripts[tid].score = 0
                else:
                    self.logger.debug("Excluding %s from %s because of failed requirements",
                                      tid, self.id)
                    self.remove_transcript_from_locus(tid)

            if not self.purge:
                assert len(self.transcripts) == beginning

            if len(self.transcripts) == 0 or self.metrics_calculated is True:
                return
            elif self.purge and len(not_passing) > 0:
                assert self._not_passing
            else:
                # Recalculate the metrics
                self.get_metrics()

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
        if "filter" not in self.json_conf["scoring"][param] and max(scores) <= 0:
            self.logger.warning("All transcripts have a score of 0 for %s in %s",
                                param, self.id)

    @classmethod
    def _calculate_graph(cls, transcripts):

        """Private method to calculate the internal graph of exons/introns in a locus."""

        graph = networkx.DiGraph()

        [cls.add_path_to_graph(transcript, graph) for transcript in transcripts]

        return graph

    @staticmethod
    def add_path_to_graph(transcript, graph):
        weights = networkx.get_node_attributes(graph, "weight")

        segments = sorted(list(transcript.exons) + list(transcript.introns), reverse=(transcript.strand == "-"))
        # Add path FAILS if the transcript is monoexonic!
        graph.add_nodes_from(segments)
        graph.add_path(segments)

        # assert len(graph.nodes()) >= len(segments), (len(graph.nodes()), len(graph.edges()), len(segments))

        for segment in segments:
            weights[segment] = weights.get(segment, 0) + 1

        networkx.set_node_attributes(graph, name="weight", values=weights)
        return

    @staticmethod
    def remove_path_from_graph(transcript: Transcript, graph: networkx.DiGraph):

        weights = networkx.get_node_attributes(graph, "weight")
        segments = sorted(list(transcript.exons) + list(transcript.introns), reverse=(transcript.strand == "-"))
        for segment in segments:
            weights[segment] -= 1

        nodes_to_remove = [interval for interval, weight in weights if weight == 0]
        graph.remove_nodes_from(nodes_to_remove)
        for node in nodes_to_remove:
            assert node not in graph.nodes()

        networkx.set_node_attributes(graph,
                                     name="weight",
                                     values=dict((k, v) for k, v in weights.items() if v > 0))
        return

    @classmethod
    @abc.abstractmethod
    def is_intersecting(cls, *args, **kwargs):
        """

        :param args: positional arguments
        :param kwargs: keyword arguments

        This class method defines how two transcript objects will be considered as overlapping.
        It is used by the BronKerbosch method, and must be implemented
        at the class level for each child object.
        """
        raise NotImplementedError("The is_intersecting method should be defined for each child!")

    # ##### Properties #######

    @property
    def json_conf(self):
        return self.__json_conf

    @json_conf.setter
    def json_conf(self, conf):
        if conf is None:
            conf = json_conf.copy()
        elif isinstance(conf, str):
            conf = to_json(conf)
        elif not isinstance(conf, dict):
            raise TypeError("Invalid configuration!")
        self.__json_conf = conf

    def check_configuration(self):
        """Method to be invoked to verify that the configuration is correct.
        Quite expensive to run, especially if done multiple times."""

        self.json_conf = check_json(self.json_conf)

    @property
    def stranded(self):
        """This property determines whether a Monosublocus will consider
        the strand for e.g. the in_locus method.
        By default, the parameter is set to True (i.e. the loci are strand-specific).
        At the moment, the only class which modifies the parameter is the superlocus class."""
        return self.__stranded

    @stranded.setter
    def stranded(self, flag):
        """
        :param flag: boolean value
        :type flag: bool
        """

        if not isinstance(flag, bool):
            raise ValueError("The stranded attribute must be boolean!")
        self.__stranded = flag

    # pylint: disable=invalid-name
    @property
    def id(self) -> str:
        """
        This is a generic string generator for all inherited children.
        :rtype : str
        """
        return "{0}:{1}{2}:{3}-{4}".format(
            self.__name__,
            self.chrom,
            self.strand,
            self.start,
            self.end)
    # pylint: enable=invalid-name

    @property
    def name(self) -> str:
        """
        Alias for id.
        :rtype : str
        """
        return self.id

    @property
    def logger(self):
        """
        Logger instance for the class.
        :rtype : logging.Logger
        """
        return self.__logger

    @logger.setter
    def logger(self, logger):
        """Set a logger for the instance.
        :param logger
        :type logger: logging.Logger | Nonell
        """

        # if self.__logger is not None:
        #     self.__logger.disabled = True

        if logger is None:
            self.__logger = create_null_logger(self)
            self.__logger.propagate = False
        elif not isinstance(logger, logging.Logger):
            raise TypeError("Invalid logger: {0}".format(type(logger)))
        else:
            # while len(logger.handlers) > 1:
            #     logger.handlers.pop()
            # logger.propagate = False
            self.__logger = logger

        # self.__logger.setLevel("DEBUG")

    @logger.deleter
    def logger(self):
        """
        Deleter method. It sets the logger to None. Used specifically for pickling.
        """
        self.__logger = None

    @property
    def source(self):
        """
        Property. Returns the source field.
        :rtype : str
        """
        return self.__source

    @property
    def _internal_graph(self):
        return self.__internal_graph

    @source.setter
    def source(self, value):
        """
        Setter for source. It accepts only strings.
        :param value:
        :type value: str

        """
        if not value:
            value = "Mikado"
        assert isinstance(value, str)
        self.__source = value

    @property
    def score(self):

        if len(self.transcripts):
            return max(_.score for _ in self.transcripts.values())
        else:
            return None

    @property
    def segmenttree(self):
        if len(self.__segmenttree) != len(self.exons) + len(self.introns):
            self.__segmenttree = self._calculate_segment_tree(self.exons, self.introns)

        return self.__segmenttree

    @property
    def cds_segmenttree(self):

        if len(self.__cds_segmenttree) != len(self.combined_cds_exons) + len(self.combined_cds_introns):
            self.__cds_segmenttree = self._calculate_segment_tree(self.combined_cds_exons, self.combined_cds_introns)

        return self.__cds_segmenttree

    @staticmethod
    def _calculate_segment_tree(exons, introns):

        return IntervalTree.from_intervals(
                [Interval(*_, value="exon") for _ in exons] + [Interval(*_, value="intron") for _ in introns]
            )

    @property
    def _cds_introntree(self):

        """
        :rtype: IntervalTree
        """

        if len(self.__cds_introntree) != len(self.combined_cds_introns):
            self.__cds_introntree = IntervalTree.from_tuples(
                [(_[0], _[1] + 1) for _ in self.combined_cds_introns])
        return self.__cds_introntree

    @property
    def longest_transcript(self):
        return max([len(_) for _ in self.transcripts.values()])

    @property
    def regressor(self):
        return self.__regressor

    @regressor.setter
    def regressor(self, regr):

        if isinstance(regr, dict) and isinstance(regr["scoring"], (RandomForestRegressor, RandomForestClassifier)):
            self.__regressor = regr["scoring"]
        elif regr is None or isinstance(regr, (RandomForestRegressor, RandomForestClassifier)):
            self.__regressor = regr
        else:
            raise TypeError("Invalid regressor provided, type: %s", type(regr))

        self.logger.debug("Set regressor")

    @property
    def locus_verified_introns(self):
        return self.__locus_verified_introns

    @locus_verified_introns.setter
    def locus_verified_introns(self, *args):

        if not isinstance(args[0], set):
            raise ValueError("Invalid value for verified introns: %s",
                             type(args[0]))

        self.__locus_verified_introns = args[0]

    @property
    def purge(self):
        """This property relates to pick/clustering/purge."""

        return self.json_conf.get("pick", dict()).get("clustering", {}).get("purge", True)
