# coding: utf-8

"""
Module that defines the blueprint for all loci classes.
"""

import abc
import itertools
import numpy as np
import logging
from sys import maxsize
import networkx
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
import numpy
from ..transcripts.clique_methods import find_communities, define_graph
from ..transcripts.transcript import Transcript
from ..configuration.configurator import to_json, check_json
from ..exceptions import NotInLocusError
from ..utilities import overlap, merge_ranges, rhasattr, rgetattr
import operator
from ..utilities.intervaltree import Interval, IntervalTree
from ..utilities.log_utils import create_null_logger
from sys import version_info
from ..scales.contrast import compare as c_compare
if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict
import random
import rapidjson as json
from typing import Union

# I do not care that there are too many attributes: this IS a massive class!
# pylint: disable=too-many-instance-attributes,too-many-public-methods
json_conf = to_json(None)


def to_bool(param: Union[str,bool,int,float]):
    if isinstance(param, bool):
        return param
    elif isinstance(param, (int, float)):
        if param == 1:
            return True
        elif param == 0:
            return False
    else:
        lparam = param.lower()
        if lparam == 'true':
            return True
        elif lparam == 'false':
            return False

    raise ValueError


class Abstractlocus(metaclass=abc.ABCMeta):
    """This abstract class defines the basic features of any Locus-like object.
    It also defines methods/properties that are needed throughout the program,
    e.g. the Bron-Kerbosch algorithm for defining cliques, or the find_retained_introns method.
    """

    __name__ = "Abstractlocus"
    cast_to = {'int': int,
               'float': float,
               'bool': to_bool}
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
                 json_conf=None,
                 use_transcript_scores=False,
                 flank=None):

        # Mock values
        self.__source = source
        self._attribute_metrics = dict()
        self.__logger = None
        self.logger = logger
        self.__stranded = False
        self._not_passing = set()
        self._excluded_transcripts = dict()
        self.transcripts = dict()
        self.start, self.end, self.strand = maxsize, -maxsize, None
        # Consider only the CDS part
        self.stranded = True
        self.initialized = False
        self.monoexonic = True
        self.chrom = None
        self.__locus_verified_introns = set()
        # This will create the .splices, .introns, etc. stores for the locus
        for locattr in self.__locus_to_transcript_attrs:
            setattr(self, locattr, set())
        if verified_introns is not None:
            self.locus_verified_introns = verified_introns

        self.scores_calculated = False
        self.scores = dict()
        self.__segmenttree = IntervalTree()
        self.__regressor = None
        self.session = None
        self.metrics_calculated = False
        self._metrics = dict()
        self.__scores = dict()
        self.__internal_graph = networkx.DiGraph()
        self.json_conf = json_conf
        if transcript_instance is not None and isinstance(transcript_instance, Transcript):
            self.add_transcript_to_locus(transcript_instance)
        self.__use_transcript_scores = use_transcript_scores
        self.__flank = 0
        if flank is not None:
            self.flank = flank
        else:
            self.flank = self.json_conf["pick"]["clustering"]["flank"]

    @abc.abstractmethod
    def __str__(self, *args, **kwargs):
        """Printing method for the locus class. Each child class must have its own defined method."""

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

        return (
            self.start, self.end, self.chrom, self.strand, self.stranded,
            self.introns, self.splices, self.exons
        ) == (other.start, other.end, other.chrom, other.strand, other.stranded,
              other.introns, other.splices, other.exons)

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
            if "json_conf" not in state:
                state["json_conf"] = self.json_conf.copy()
            for key in self.json_conf:
                if (isinstance(self.json_conf[key], dict) and
                        self.json_conf[key].get("compiled", None) is not None):
                    assert key in state["json_conf"]
                    del state["json_conf"][key]["compiled"]

        if hasattr(self, "session"):
            if self.session is not None:
                self.session.expunge_all()
                state["session"].expunge_all()
            state["sessionmaker"] = None
            state["session"] = None

        if self.__internal_graph.nodes():
            try:
                nodes = json.dumps(list(self.__internal_graph.nodes())[0], number_mode=json.NM_NATIVE)
            except ValueError:
                nodes = json.dumps(list(self.__internal_graph.nodes())[0])
        else:
            nodes = "[]"
        state["_Abstractlocus__internal_nodes"] = nodes

        edges = [edge for edge in self.__internal_graph.edges()]
        # Remember that the graph is in form [((start, end), (start, end)), etc.]
        # So that each edge is composed by a couple of tuples.
        try:
            state["_Abstractlocus__internal_edges"] = json.dumps(edges, number_mode=json.NM_NATIVE)
        except ValueError:
            state["_Abstractlocus__internal_edges"] = json.dumps(edges)
        if hasattr(self, "engine"):
            del state["engine"]

        del state["_Abstractlocus__segmenttree"]
        return state

    def __setstate__(self, state):
        """Method to recreate the object after serialisation."""
        self.__dict__.update(state)
        self.__segmenttree = IntervalTree()
        self.__internal_graph = networkx.DiGraph()
        try:
            nodes = json.loads(state["_Abstractlocus__internal_nodes"])
        except json.decoder.JSONDecodeError:
            raise json.decoder.JSONDecodeError(state["_Abstractlocus__internal_nodes"])
        edges = []
        for edge in json.loads(state["_Abstractlocus__internal_edges"]):
            edges.append((tuple(edge[0]), tuple(edge[1])))
        try:
            self.__internal_graph.add_nodes_from(nodes)
        except TypeError:
            raise TypeError(nodes)
        try:
            self.__internal_graph.add_edges_from(edges)
        except TypeError:
            raise TypeError(edges)

        if hasattr(self, "json_conf"):
            if "requirements" in self.json_conf and "expression" in self.json_conf["requirements"]:
                self.json_conf["requirements"]["compiled"] = compile(
                    self.json_conf["requirements"]["expression"],
                    "<json>", "eval")
        # Set the logger to NullHandler
        _ = self.__segmenttree
        self.logger = None

    def as_dict(self):
        self.get_metrics()
        state = self.__getstate__()
        state["transcripts"] = dict((tid, state["transcripts"][tid].as_dict()) for tid in state["transcripts"])
        assert "metrics_calculated" in state
        return state

    def load_dict(self, state, load_transcripts=True):
        assert isinstance(state, dict)
        self.__setstate__(state)

        assert self.metrics_calculated is True
        if load_transcripts is True:
            self.transcripts = dict((tid, Transcript()) for tid in state["transcripts"])
            [self[tid].load_dict(state["transcripts"][tid], trust_orf=True) for tid in state["transcripts"]]

        for attr in ["locus_verified_introns", "introns", "exons",
                     "selected_cds_introns", "combined_cds_introns"]:
            setattr(self, attr, set([tuple(_) for _ in getattr(self, attr)]))

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

        if not isinstance(transcript, Transcript):
            raise TypeError("I can only perform this operation on transcript classes, not {}".format(
                type(transcript)))

        transcript.finalize()
        # We want to check for the strand only if we are considering the strand
        if not isinstance(locus_instance, cls):
            raise TypeError("I cannot perform this operation on non-locus classes, this is a {}".format(
                type(locus_instance)))

        if not hasattr(locus_instance, "chrom"):
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

    def choose_best(self, transcripts: dict) -> str:
        """
        :param transcripts: the dictionary of transcripts of the instance
        :type transcripts: dict

        Given a transcript dictionary, this function will choose the one with the highest score.
        If multiple transcripts have exactly the same score, one will be chosen randomly.

        """

        # Choose one transcript randomly between those that have the maximum score
        if len(transcripts) == 1:
            return list(transcripts.keys())[0]
        # np.random.seed(self.json_conf["seed"])
        random.seed(self.json_conf["seed"])
        if self.reference_update is True:
            # We need to select only amongst the reference transcripts
            reference_sources = {source for source, is_reference in
                                 zip(self.json_conf["prepare"]["files"]["labels"],
                                     self.json_conf["prepare"]["files"]["reference"]) if is_reference is True}
            reference_tids = set()
            for tid, transcript in transcripts.items():
                is_reference = transcript.original_source in reference_sources or transcript.is_reference is True
                if is_reference:
                    reference_tids.add(tid)
            if len(reference_tids) > 0:
                transcripts = dict((tid, transcripts[tid]) for tid in reference_tids)

        max_score = max(transcripts.values(),
                        key=operator.attrgetter("score")).score
        valid = sorted([transc for transc in transcripts if transcripts[transc].score == max_score])
        # chosen = valid[numpy.random.choice(len(valid))]
        chosen = valid[random.choice(range(len(valid)))]
        self.logger.debug("Chosen {chosen} out of {}".format(", ".join(valid), chosen=chosen))
        return chosen

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
        self.monoexonic = self.monoexonic and self._is_transcript_monoexonic(transcript)

        if "flank" in kwargs:
            pass
        else:
            kwargs["flank"] = self.flank

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
        if self._use_transcript_scores is False:
            self.scores_calculated = False
        else:
            self.scores[transcript.id] = transcript.score
            self.scores_calculated = True

        return

    def _swap_transcript(self,
                         original_transcript: Transcript,
                         transcript: Transcript):

        """This method is needed to exchange transcripts that might have been modified by padding."""
        if original_transcript.tid != transcript.tid:
            raise KeyError("I cannot hot swap two transcripts with two different IDs!")
        # if hash(original_transcript) == hash(transcript):  # Expensive operation, let us try to avoid it if possible
        #     return
        if original_transcript == transcript:
            return

        self.logger.debug("Swapping %s with a new transcript", original_transcript.id)
        self.transcripts[original_transcript.id] = original_transcript
        Abstractlocus.remove_transcript_from_locus(self, original_transcript.id)
        Abstractlocus.add_transcript_to_locus(self, transcript, check_in_locus=False)
        return

    def _is_transcript_monoexonic(self, transcript: Transcript):
        if self._cds_only and transcript.is_coding is True:
            return len(transcript.selected_cds) == 1
        else:
            return transcript.monoexonic

    def remove_transcript_from_locus(self, tid: str):
        """
        :param tid: name of the transcript to remove
        :type tid: str

         This method will remove a transcript from an Abstractlocus-like
         instance and reset appropriately all derived attributes
         (e.g. introns, start, end, etc.).
        """

        if tid not in self.transcripts:
            self.logger.debug("Transcript %s is not present in the Locus. Ignoring it.", tid)
            return

        self.logger.debug("Deleting %s from %s", tid, self.id)
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
            # self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
            import sys
            self.start, self.end, self.strand = sys.maxsize, -sys.maxsize, None
            self.stranded = False
            self.initialized = False

        self.logger.debug("Deleted %s from %s", tid, self.id)
        if tid in self._metrics:
            del self._metrics[tid]
        if tid in self.scores:
            del self.scores[tid]

        self.metrics_calculated = False
        self.scores_calculated = False

    def _remove_all(self):
        """This method will remove all transcripts from the locus."""
        self.logger.warning("Removing all transcripts from %s", self.id)
        self.__internal_graph = networkx.DiGraph()
        self.transcripts = dict()
        self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
        self.stranded = False
        self.initialized = False
        self.metrics_calculated = False
        self.scores_calculated = False

    @staticmethod
    def _exon_to_be_considered(exon,
                               transcript,
                               logger=create_null_logger()):
        """Private static method to evaluate whether an exon should be considered for being a retained intron.

        :param exon: the exon to be considered.
        :type exon: (tuple|Interval)

        :param transcript: the candidate transcript from which the exon comes from.
        :type transcript: Transcript

        :returns: boolean flag (True if it has to be considered, False otherwise), a list of sections of the exons
        which are non-coding, and a set of the boundaries that are splicing sites (0, 1 or 2).
        """

        cds_segments = sorted(transcript.cds_tree.search(*exon))

        internal_splices = set.difference(set(exon), {transcript.start, transcript.end})
        before_met = False
        if transcript.is_coding:
            # Avoid considering exons that are full 3' UTR exons in a coding transcript.
            if transcript.strand == "-":
                if max(transcript.combined_cds_start, transcript.combined_cds_end) < exon[1]:
                    before_met = True
            elif transcript.strand != "-":
                if min(transcript.combined_cds_start, transcript.combined_cds_end) > exon[0]:
                    before_met = True

        if cds_segments == [Interval(*exon)]:
            # It is completely coding
            if len(internal_splices) == 2:
                logger.debug("%s is internal and completely coding. Not considering it.", exon)
                to_consider = False
                frags = []
            else:
                logger.debug("%s is terminal and completely coding. Considering it in its entirety.", exon)
                to_consider = True
                frags = cds_segments
        else:
            frags = []
            to_consider = True
            if cds_segments:
                logger.debug("Calculating fragments for %s (CDS segments: %s)", exon, cds_segments)
                if cds_segments[0].start > exon[0]:
                    if before_met and transcript.strand == "+":
                        frags.append((exon[0], cds_segments[0].start - 1))
                        logger.debug("New fragment for %s: %s", exon, frags[-1])
                    elif transcript.strand == "-":  # Negative UTR
                        frags.append((exon[0], cds_segments[0].start - 1))
                        logger.debug("New fragment for %s: %s", exon, frags[-1])
                    else:
                        frags.append((cds_segments[0].start - 1, cds_segments[0].start))
                        logger.debug("New fragment for %s: %s", exon, frags[-1])
                for before, after in zip(cds_segments[:-1], cds_segments[1:]):
                    frags.append((before.end, before.end + 1))
                    frags.append((after.start - 1, after.start))
                    logger.debug("New fragments for %s: %s, %s", exon, frags[-2], frags[-1])
                    # frags.append((before.end + 1, max(after.start - 1, before.end + 1)))
                if cds_segments[-1].end < exon[1]:
                    if before_met and transcript.strand == "-":
                        frags.append((cds_segments[-1].end + 1, exon[1]))
                        logger.debug("New fragment for %s: %s", exon, frags[-1])
                    elif transcript.strand == "+":
                        frags.append((cds_segments[-1].end, exon[1]))
                        logger.debug("New fragment for %s: %s", exon, frags[-1])
                    else:
                        frags.append((cds_segments[-1].end, cds_segments[-1].end + 1))
                        logger.debug("New fragment for %s: %s", exon, frags[-1])
                logger.debug("%s is partially non-coding. Considering it, with frags: %s.", exon, frags)
            elif before_met is True:
                logger.debug("%s is completely non-coding in the 5'UTR. Considering it as it might disrupt the CDS.",
                             exon)
                frags.append((exon[0], exon[1]))
            else:
                logger.debug("%s is completely non-coding. Considering it, but with no frags.", exon)
                frags = []  # [Interval(*exon)]

        return to_consider, frags, internal_splices

    @staticmethod
    def _is_exon_retained(exon: tuple,
                          strand: str,
                          segmenttree: IntervalTree,
                          digraph: networkx.DiGraph,
                          frags: list,
                          introns: set,
                          internal_splices: set,
                          cds_introns: set,
                          coding=True,
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

        :param internal_splices: which of the two boundaries of the exon are actually internal (if any).
        :type internal_splices: set

        :param terminal: whether the exon is at the 3' end.
        :type terminal: bool

        :rtype: bool
        """

        is_retained = False
        cds_broken = False

        found_introns = set(
            [_._as_tuple() for _ in segmenttree.find(exon[0], exon[1], strict=False, value="intron")]
        )

        logger.debug("Analysing exon %s with frags %s", exon, frags)

        if not found_introns:
            logger.debug("No intron found for %s, returning False.", exon)
            return is_retained, cds_broken

        # logger.debug("Found introns for %s: %s", exon, found_introns)

        for intron in found_introns:
            # Only consider exons for which there is an overlap.
            if is_retained:
                if intron in cds_introns and cds_broken is True:
                    continue
                elif intron not in cds_introns:
                    continue
                elif coding is False:
                    break

            before = {_ for _ in networkx.ancestors(digraph, intron) if
                      _ not in introns and overlap(_, exon) > 0}

            after = {_ for _ in networkx.descendants(digraph, intron) if
                     _ not in introns and overlap(_, exon) > 0}

            # Now we have to check whether the matched introns contain both coding and non-coding parts
            # Let us exclude any intron which is outside of the exonic span of interest.
            # logger.debug("Exon: %s; Frags: %s; Intron: %s; Before: %s; After: %s", exon, frags, intron, before, after)
            if len(before) == 0 and len(after) == 0:
                # A retained intron must be overlapping some other exons!
                logger.debug("No before/after exonic overlap found for exon %s vs intron %s. Skipping", exon, intron)
                continue
            else:
                if len(internal_splices) == 0:
                    start_found = True
                    end_found = True
                elif strand == "-":
                    # Negative strand
                    end_found = (exon[0] in set([e[0] for e in after - introns]))
                    start_found = (exon[1] in set([e[1] for e in before - introns]))
                    logger.debug("Exon %s vs intron %s: strand %s, start found %s, end found %s (I.S. %s)",
                                 exon, intron, strand, start_found, end_found, internal_splices)
                    if len(internal_splices) != 2:
                        if exon[1] in internal_splices:  # This means that the end is dangling
                            end_found = True
                        elif exon[0] in internal_splices:  # This means that the start is dangling
                            start_found = True
                else:
                    start_found = (exon[0] in set([e[0] for e in before - introns]))
                    end_found = (exon[1] in set([e[1] for e in after - introns]))
                    logger.debug("Exon %s vs intron %s: strand %s, start found %s, end found %s (I.S. %s)",
                                 exon, intron, strand, start_found, end_found, internal_splices)
                    if len(internal_splices) == 1:
                        if exon[0] in internal_splices:  # This means that the end is dangling
                            end_found = True
                        elif exon[1] in internal_splices:  # This means that the start is dangling
                            start_found = True

                logger.debug(
                    "Exon %s vs intron %s: strand %s, retained %s, start found %s, end found %s (I.S. %s); frags: %s",
                    exon, intron, strand, start_found and end_found, start_found, end_found, internal_splices, frags)
                if start_found and end_found:
                    logger.debug("Exon %s is retained (for intron %s)", exon, intron)
                    # Now we have to check whether the CDS breaks within the intron
                    if intron in cds_introns:
                        for frag, intron in itertools.product(frags, [intron]):
                            cds_broken = cds_broken or (overlap(frag, intron, positive=True) > 0)
                            if cds_broken is True:
                                logger.debug("Frag %s intersecting intron %s: CDS interrupted", frag, intron)
                                break
                            else:
                                logger.debug("Frag %s of exon %s does not intersect intron %s.",
                                             frag, exon, intron)

                is_retained = is_retained or (start_found and end_found)

        return is_retained, cds_broken

    def _load_scores(self, scores: dict):
        """This private method is present *strictly for testing purposes only*.
        Its aim is to load some pre-calculated scores for the transcripts in the locus,
        *completely bypassing the normal method of calculating scores*."""

        if not isinstance(scores, dict):
            raise ValueError("This private method takes strictly a dictionary as input")

        if set.difference(set(self.transcripts.keys()), set(scores.keys())):
            raise KeyError("I am missing transcripts from the scores dictionary. Aborting")

        for tid in self.transcripts:
            self.scores[tid] = scores[tid]
            self[tid].score = scores[tid]

        self.scores_calculated = True
        self.metrics_calculated = True

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

        # self.logger.debug("Starting to calculate retained introns for %s", transcript.id)
        if self.stranded is False:
            self.logger.error("Trying to find retained introns in a non-stranded locus (%s) is invalid. Aborting.",
                              self.id)
            return

        transcript.logger = self.logger
        transcript.finalize()
        # Reset this flag
        transcript.cds_disrupted_by_ri = False
        if len(self.introns) == 0:
            transcript.retained_introns = tuple()
            self.logger.debug("No introns in the locus to check against. Exiting.")
            return

        # A retained intron is defined as an exon which
        # - is not completely coding
        # - EITHER spans completely the intron of another transcript.
        # - OR is the last exon of the transcript and it ends within the intron of another transcript

        retained_introns = []
        # consider_truncated = self.json_conf["pick"]["run_options"]["consider_truncated_for_retained"]

        cds_broken = False

        for exon in transcript.exons:
            # self.logger.debug("Checking exon %s of %s", exon, transcript.id)
            # is_retained = False
            to_consider, frags, internal_splices = self._exon_to_be_considered(exon, transcript, logger=self.logger)
            if not to_consider:
                # self.logger.debug("Exon %s of %s is not to be considered", exon, transcript.id)
                continue

            result = self._is_exon_retained(
                exon,
                self.strand,
                self.segmenttree,
                self._internal_graph,
                frags,
                internal_splices=internal_splices,
                introns=self.introns,
                cds_introns=self.combined_cds_introns,
                coding=transcript.is_coding,
                logger=self.logger)

            is_retained = result[0]
            if is_retained:
                self.logger.debug("Exon %s of %s is a retained intron", exon, transcript.id)
                retained_introns.append(exon)

            cds_broken = transcript.is_coding and (cds_broken or result[1])

            self.logger.debug("After exon %s, CDS broken: %s", exon, cds_broken)

        self.logger.debug("%s has %s retained intron%s%s",
                          transcript.id,
                          len(retained_introns) if retained_introns else "no",
                          "s" if len(retained_introns) != 1 else "",
                          " ({})".format(retained_introns) if retained_introns else "")
        self.logger.debug("%s has its CDS interrupted by a retained intron: %s",
                          transcript.id, cds_broken)

        transcript.cds_disrupted_by_ri = cds_broken
        transcript.retained_introns = tuple(sorted(retained_introns))
        return

    @staticmethod
    def _evaluate_transcript_overlap(
            other: Transcript,
            transcript: Transcript,
            min_cdna_overlap=0.2,
            min_cds_overlap=0.2,
            comparison=None,
            check_references=True,
            fixed_perspective=True):

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

        :param comparison: default None. If one is provided, it should be a pre-calculated output of Mikado compare.

        :param fixed_perspective: boolean. If True, the method will consider the second transcript as the "reference" and
        therefore the cDNA and CDS overlap percentages will be calculated using its ratio. Otherwise, the method will
        consider the shortest cDNA and shortest CDS for calculating the ratio. Mikado uses the **former** method when
        evaluating potential alternative splicing events (as we are interested in the relationship with the primary)
        and the **latter** during the monosublocus stage.
        :type fixed_perspective: bool
        """

        if comparison is None:
            comparison, _ = c_compare(other, transcript)

        cds_overlap = 0
        if fixed_perspective is False:
            cdna_overlap = max(comparison.n_prec[0], comparison.n_recall[0]) / 100
        else:
            cdna_overlap = comparison.n_recall[0] / 100

        if not (transcript.is_coding and other.is_coding):
            cds_overlap = cdna_overlap
        else:
            # Only consider the selected internal ORF
            t_orfs = [(_[1][0], _[1][1], _[2]) for _ in transcript.selected_internal_orf if _[0] == "CDS"]
            o_orfs = [(_[1][0], _[1][1], _[2]) for _ in other.selected_internal_orf if _[0] == "CDS"]

            for start, end, phase in t_orfs:
                for ostart, oend, ophase in o_orfs:
                    seg_overlap = overlap((start, end), (ostart, oend))
                    if seg_overlap <= 0:
                        continue
                    # For this, we do have to check whether the start + phase of the upstream transcript
                    # (considering the strand) is in frame with the downstream start + phase
                    if transcript.strand == "-":
                        first, last = sorted([(start, end, phase), (ostart, oend, ophase)],
                                             key=operator.itemgetter(1), reverse=True)
                        in_frame = (0 == ((last[1] - last[2]) - (first[1] - first[2])) % 3)
                    else:
                        first, last = sorted([(start, end, phase), (ostart, oend, ophase)],
                                             key=operator.itemgetter(0), reverse=False)
                        in_frame = (0 == ((last[0] + last[2]) - (first[0] + first[2])) % 3)
                    if in_frame:
                        cds_overlap += seg_overlap

            if fixed_perspective:
                cds_overlap /= transcript.combined_cds_length
            else:
                cds_overlap /= min(transcript.selected_cds_length, other.selected_cds_length)
            assert cds_overlap <= 1

        if other.is_reference is True and check_references is False and transcript.is_reference is True:
            intersecting = True
            reason = "{} is a reference transcript being added to a reference locus. Keeping it.".format(other.id)
        elif transcript.is_coding and other.is_coding:
            intersecting = (cdna_overlap >= min_cdna_overlap and cds_overlap >= min_cds_overlap)

            reason = "{} and {} {}share enough cDNA ({}%, min. {}%) and CDS ({}%, min. {}%), {}intersecting".format(
                transcript.id, other.id,
                "do not " if not intersecting else "",
                cdna_overlap * 100, min_cdna_overlap * 100,
                cds_overlap * 100, min_cds_overlap * 100,
                "not " if not intersecting else "")
        else:
            intersecting = (cdna_overlap >= min_cdna_overlap)
            reason = "{} and {} {}share enough cDNA ({}%, min. {}%), {}intersecting".format(
                transcript.id, other.id,
                "do not " if not intersecting else "",
                cdna_overlap * 100, min_cdna_overlap * 100,
                "not " if not intersecting else "")

        return intersecting, reason

    def print_metrics(self):

        """This method yields dictionary "rows" that will be given to a csv.DictWriter class."""

        # Check that rower is an instance of the csv.DictWriter class

        # The rower is an instance of the DictWriter class from the standard CSV module

        self.get_metrics()

        for tid, transcript in sorted(self.transcripts.items(), key=operator.itemgetter(1)):
            row = {}
            if tid not in self._metrics and transcript.alias in self._metrics:
                metrics = self._metrics[transcript.alias]
            else:
                metrics = self._metrics[tid]

            for num, key in enumerate(self.available_metrics):

                if num == 0:  # transcript id
                    value = tid
                elif num == 2:  # Parent
                    value = self.id
                else:
                    value = metrics.get(key, "NA")
                    # value = getattr(transcript, key, "NA")
                if isinstance(value, float):
                    value = round(value, 2)
                elif value is None or value == "":
                    if key == "score":
                        value = self.scores.get(tid, dict()).get("score", None)
                        self.transcripts[tid].score = value
                        if isinstance(value, float):
                            value = round(value, 2)
                        elif value is None:
                            value = "NA"
                    else:
                        value = "NA"
                row[key] = value

            for source in transcript.external_scores:
                # Each score from external files also contains a multiplier.
                key = "external.{}".format(source)
                value = transcript.external_scores.get(source)[0]
                if isinstance(value, float):
                    value = round(value, 2)
                elif value is None or value == "":
                    value = "NA"
                row[key] = value

            # assert row != {}
            yield row

        return

    def get_metrics(self):

        """Quick wrapper to calculate the metrics for all the transcripts."""

        # if self.metrics_calculated is True:
        #     return

        if self.metrics_calculated is True:
            return
        self._metrics = dict()
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

        if self.metrics_calculated is True:
            return

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

        try:
            self.transcripts[tid].exon_fraction = fraction
        except ValueError:
            raise ValueError("Invalid fraction. Exons:\n{}\n{}".format(
                self.transcripts[tid].exons,
                self.exons
            ))

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
        if len(self.transcripts) > 1 and any(len(self[_].introns) > 0 for _ in self if _ != tid):
            self.find_retained_introns(self.transcripts[tid])
        assert isinstance(self.transcripts[tid], Transcript)

        retained_bases = sum(e[1] - e[0] + 1
                             for e in self.transcripts[tid].retained_introns)
        fraction = retained_bases / self.transcripts[tid].cdna_length
        self.transcripts[tid].retained_fraction = fraction

        self._metrics[tid] = dict((metric, rgetattr(self.transcripts[tid], metric))
                                   for metric in self.available_metrics)

        for metric, values in self._attribute_metrics.items():
            # 11 == len('attributes.') removes 'attributes.' to keep the metric name same as in the file attributes
            rtype = self.cast_to[values['rtype']]
            attribute_metric_value = self.transcripts[tid].attributes.get(metric[11:], values['default'])
            try:
                attribute_metric_value = rtype(attribute_metric_value)
                if values['percentage']:
                    attribute_metric_value = attribute_metric_value / 100.
            except ValueError:
                message = f"Error encountered when processing Transcript {tid}. " \
                          f"The 'attributes' based metric {metric}, with value {attribute_metric_value} " \
                          f"could not be created as type {values['rtype']}"
                self.logger.error(message)
                raise ValueError(message)
            if values['use_raw'] is True and not (0 <= attribute_metric_value <= 1):
                message = f"Error encountered when processing Transcript {tid}. " \
                          f"The 'attributes' based metric {metric}, with value {attribute_metric_value} " \
                          f"defined as 'use_raw' is not between 0 and 1 inclusive."
                self.logger.error(message)
                raise ValueError(message)

            self._metrics[tid].update([(metric, attribute_metric_value)])

        self.logger.debug("Calculated metrics for {0}".format(tid))

    def _check_not_passing(self, previous_not_passing=set(), section="requirements"):
        """
        This private method will identify all transcripts which do not pass
        the minimum muster specified in the configuration. It will *not* delete them;
        how to deal with them is left to the specifics of the subclass.
        :return:
        """

        self.get_metrics()

        if ("compiled" not in self.json_conf[section] or
                self.json_conf[section]["compiled"] is None):
            if "expression" not in self.json_conf[section]:
                raise KeyError(self.json_conf[section])

            self.json_conf[section]["compiled"] = compile(
                self.json_conf[section]["expression"], "<json>",
                "eval")

        not_passing = set()
        reference_sources = {source for source, is_reference in
                             zip(self.json_conf["prepare"]["files"]["labels"],
                                 self.json_conf["prepare"]["files"]["reference"]) if is_reference}

        for tid in iter(tid for tid in self.transcripts if
                        tid not in previous_not_passing):
            if self.transcripts[tid].json_conf is None:
                self.transcripts[tid].json_conf = self.json_conf
            elif self.transcripts[tid].json_conf["prepare"]["files"]["reference"] != self.json_conf["prepare"]["files"]["reference"]:
                self.transcripts[tid].json_conf = self.json_conf

            is_reference = ((self.transcripts[tid].is_reference is True) or
                             self.transcripts[tid].original_source in reference_sources)

            if is_reference is False:
                self.logger.debug("Transcript %s (source %s) is not a reference transcript (references: %s; in it: %s)",
                                  tid, self.transcripts[tid].original_source,
                                  self.json_conf["prepare"]["files"]["reference"],
                                  self.transcripts[tid].original_source in self.json_conf["prepare"]["files"][
                                      "reference"])
            elif is_reference is True and self.json_conf["pick"]["run_options"]["check_references"] is False:
                self.logger.debug("Skipping %s from the requirement check as it is a reference transcript", tid)
                continue
            elif is_reference is True and self.json_conf["pick"]["run_options"]["check_references"] is True:
                self.logger.debug("Performing the requirement check for %s even if it is a reference transcript", tid)

            evaluated = dict()
            for key in self.json_conf[section]["parameters"]:
                value = rgetattr(self.transcripts[tid],
                                 self.json_conf[section]["parameters"][key]["name"])
                if "external" in key:
                    value = value[0]

                evaluated[key] = self.evaluate(
                    value,
                    self.json_conf[section]["parameters"][key])
            # pylint: disable=eval-used
            if eval(self.json_conf[section]["compiled"]) is False:

                not_passing.add(tid)
        self.logger.debug("The following transcripts in %s did not pass the minimum check for requirements: %s",
                          self.id, ", ".join(list(not_passing)))

        return not_passing

    def filter_and_calculate_scores(self, check_requirements=True):
        """
        Function to calculate a score for each transcript, given the metrics derived
        with the calculate_metrics method and the scoring scheme provided in the JSON configuration.
        If any requirements have been specified, all transcripts which do not pass them
        will be either purged according to the configuration file ('pick.clustering.purge')
        or assigned a score of 0 and subsequently ignored.
        Scores are rounded to the nearest integer.
        """

        if self.scores_calculated is True:
            # assert self.transcripts[list(self.transcript.keys())[0]].score is not None
            test = set.difference(set(self.transcripts.keys()),
                                  set(self._excluded_transcripts.keys()))
            if test and self.transcripts[test.pop()].score is None:
                for tid in self.transcripts:
                    if tid in self._excluded_transcripts:
                        self.transcripts[tid] = 0
                    else:
                        self.transcripts[tid] = self.scores[tid]["score"]
                return
            else:
                self.logger.debug("Scores calculation already effectuated for %s", self.id)
                return

        self.get_metrics()
        self.logger.debug("Calculating scores for {0}".format(self.id))
        if "requirements" in self.json_conf and check_requirements:
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
                assert self.transcripts[tid].score is not None

                if tid in self._not_passing:
                    pass
                else:
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
                    self.logger.debug("Excluding %s from %s because of failed requirements", tid, self.id)
                    self._excluded_transcripts[tid] = self.transcripts[tid]
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
        multiplier = self.json_conf["scoring"][param]["multiplier"]

        metrics = dict()
        for tid, transcript in self.transcripts.items():
            try:
                # metric = rgetattr(self.transcripts[tid], param)
                if tid not in self._metrics and transcript.alias in self._metrics:
                    if param in self._metrics[transcript.alias]:
                        metric = self._metrics[transcript.alias][param]
                    else:
                        metric = rgetattr(self.transcripts[tid], param)
                        self._metrics[transcript.alias][param] = metric
                else:
                    if tid not in self._metrics:
                        self._metrics[tid] = dict()
                    if param in self._metrics[tid]:
                        metric = self._metrics[tid][param]
                    else:
                        metric = rgetattr(self.transcripts[tid], param)
                        self._metrics[tid][param] = metric
                if isinstance(metric, (tuple, list)):
                    metric = metric[0]
                metrics[tid] = metric
            except TypeError:
                raise TypeError(param)
            except KeyError:
                metric = rgetattr(self.transcripts[tid], param)
                raise KeyError((tid, param, metric))
            except AttributeError:
                raise AttributeError(param)

        for tid in self.transcripts.keys():
            tid_metric = metrics[tid]

            if ("filter" in self.json_conf["scoring"][param] and
                    self.json_conf["scoring"][param]["filter"] != {}):
                if "metric" not in self.json_conf["scoring"][param]["filter"]:
                    metric_to_evaluate = tid_metric
                else:
                    metric_key = self.json_conf["scoring"][param]["filter"]["metric"]
                    if not rhasattr(self.transcripts[tid], metric_key):
                        raise KeyError("Asked for an invalid metric in filter: {}".format(metric_key))
                    if tid not in self._metrics and self.transcripts[tid].alias in self._metrics:
                        metric_to_evaluate = self._metrics[self.transcripts[tid].alias][metric_key]
                    else:
                        metric_to_evaluate = self._metrics[tid][metric_key]
                    # metric_to_evaluate = rgetattr(self.transcripts[tid], metric_key)
                    if "external" in metric_key:
                        metric_to_evaluate = metric_to_evaluate[0]

                check = self.evaluate(metric_to_evaluate, self.json_conf["scoring"][param]["filter"])
                if not check:
                    del metrics[tid]
            else:
                continue

        if len(metrics) == 0:
            for tid in self.transcripts:
                self.scores[tid][param] = 0
        else:
            if param.startswith("external"):
                # Take any transcript and verify
                try:
                    transcript = self.transcripts[list(self.transcripts.keys())[0]]
                except (IndexError, TypeError, KeyError):
                    raise TypeError("No transcripts left!")
                try:
                    metric = rgetattr(transcript, param)
                except (IndexError, TypeError, KeyError):
                    raise TypeError("{param} not found in transcripts of {self.id}".format(**locals()))
                try:
                    usable_raw = metric[1]
                    if usable_raw not in (True, False):
                        raise TypeError
                except (IndexError, TypeError, KeyError):
                    raise TypeError(
                        "Value of {param} is {metric}. It should be a tuple with a boolean second element".format(
                            **locals()))
            elif param.startswith('attributes'):
                usable_raw = use_raw

            else:
                usable_raw = getattr(Transcript, param).usable_raw

            assert usable_raw in (False, True)
            if use_raw is True and usable_raw is False:
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
                    try:
                        denominator = (max(metrics.values()) - min(metrics.values()))
                    except TypeError:
                        raise TypeError([param, metrics])
            if denominator == 0:
                denominator = 1

            for tid in self.transcripts.keys():
                score = 0
                if tid in metrics:
                    tid_metric = metrics[tid]
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

                score *= multiplier
                self.scores[tid][param] = round(score, 2)

        # This MUST be true
        if "filter" not in self.json_conf["scoring"][param] and max(
                [self.scores[tid][param] for tid in self.transcripts.keys()]) == 0:
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
        networkx.add_path(graph, segments)

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
            if segment in weights:
                weights[segment] -= 1
            else:
                weights[segment] = 0

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
        # Get the value for each attribute defined metric
        self._attribute_metrics = dict((param,
                                  {
                                      'default': self.json_conf["scoring"][param]["default"],
                                      'rtype': self.json_conf['scoring'][param]['rtype'],
                                      'use_raw': self.json_conf['scoring'][param]['use_raw'],
                                      'percentage': self.json_conf['scoring'][param]['percentage']
                                  }
                                  )
                                 for param in self.__json_conf["scoring"] if param.startswith("attributes."))

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

    @property
    def flank(self):
        return self.__flank

    @flank.setter
    def flank(self, flank):
        if isinstance(flank, int) and flank >= 0:
            self.__flank = flank
        else:
            raise TypeError("Flank must be either null or an integer greater than 0")

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
            self.__logger = create_null_logger()
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
    def _cds_only(self):
        return self.json_conf["pick"]["clustering"]["cds_only"]

    @_cds_only.setter
    def _cds_only(self, value):
        if value not in (True, False):
            raise ValueError(value)
        self.json_conf["pick"]["clustering"]["cds_only"] = value

    @property
    def segmenttree(self):
        if len(self.__segmenttree) != len(self.exons) + len(self.introns):
            self.__segmenttree = self._calculate_segment_tree(self.exons, self.introns)

        return self.__segmenttree

    @staticmethod
    def _calculate_segment_tree(exons, introns):

        return IntervalTree.from_intervals(
                [Interval(*_, value="exon") for _ in exons] + [Interval(*_, value="intron") for _ in introns]
            )

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

    @property
    def _use_transcript_scores(self):
        return self.__use_transcript_scores

    @_use_transcript_scores.setter
    def _use_transcript_scores(self, val):
        if not isinstance(val, bool):
            raise ValueError("This method takes only boolean values")
        self.__use_transcript_scores = val
        self.scores_calculated = val

    @property
    def perform_padding(self):
        return self.json_conf["pick"]["alternative_splicing"]["pad"]

    def _set_padding(self, value):
        if value not in (False, True):
            raise ValueError
        self.json_conf["pick"]["alternative_splicing"]["pad"] = value

    @property
    def only_reference_update(self):
        return self.json_conf.get("pick", dict()).get("run_options", dict()).get("only_reference_update", False)

    @property
    def reference_update(self):
        ref_update = self.json_conf.get("pick", dict()).get("run_options", dict()).get("reference_update", False)
        return ref_update or self.only_reference_update
