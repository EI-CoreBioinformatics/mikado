# coding: utf-8

"""
Module that defines the blueprint for all loci classes.
"""

import abc
import itertools
import logging
import operator
import random
from sys import maxsize

import networkx
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

from ..transcripts.clique_methods import find_cliques, find_communities, define_graph
from ..transcripts.transcript import Transcript
from ..configuration.configurator import to_json, check_json
from ..exceptions import NotInLocusError
from ..utilities import overlap, merge_ranges
from ..utilities.intervaltree import Interval, IntervalTree
from ..utilities.log_utils import create_null_logger


# I do not care that there are too many attributes: this IS a massive class!
# pylint: disable=too-many-instance-attributes,too-many-public-methods
class Abstractlocus(metaclass=abc.ABCMeta):
    """This abstract class defines the basic features of any Locus-like object.
    It also defines methods/properties that are needed throughout the program,
    e.g. the Bron-Kerbosch algorithm for defining cliques, or the find_retained_introns method.
    """

    __name__ = "Abstractlocus"
    available_metrics = Transcript.get_available_metrics()

    # ##### Special methods #########

    __json_conf = to_json(None)

    @abc.abstractmethod
    def __init__(self, source="", verified_introns=None):

        # Mock values
        self.__source = source

        self.__logger = None
        self.__stranded = False

        self.transcripts = dict()
        self.introns, self.exons, self.splices = set(), set(), set()
        # Consider only the CDS part
        self.combined_cds_introns, self.selected_cds_introns = set(), set()
        self.start, self.end, self.strand = maxsize, -maxsize, None
        self.stranded = True
        self.initialized = False
        self.monoexonic = True
        self.chrom = None
        self.cds_introns = set()
        self.__locus_verified_introns = set()
        if verified_introns is not None:
            self.locus_verified_introns = verified_introns

        self.__cds_introntree = IntervalTree()
        self.__regressor = None
        self.session = None
        self.metrics_calculated = False

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
                if "compiled" in self.json_conf[key]:
                    del state["json_conf"][key]["compiled"]

        if hasattr(self, "session"):
            if self.session is not None:
                self.session.expunge_all()
                state["session"].expunge_all()
            state["sessionmaker"] = None
            state["session"] = None

        if hasattr(self, "engine"):
            del state["engine"]

        # import pickle
        # try:
        #     _ = pickle.dumps(state)
        # except pickle.PicklingError:
        #     raise pickle.PicklingError("Failed to serialise {}".format(self.id))
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
            comparison = (param in range(conf["value"][0], conf["value"][1]))
        elif conf["operator"] == "not within":
            comparison = (param not in range(conf["value"][0], conf["value"][1]))
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

    def add_transcript_to_locus(self, transcript, check_in_locus=True, logger=None, **kwargs):
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
        self.splices.update(transcript.splices)
        self.introns.update(transcript.introns)
        if transcript.monoexonic is False:
            assert len(self.introns) > 0

        self.combined_cds_introns = set.union(
            self.combined_cds_introns, transcript.combined_cds_introns)
        assert len(transcript.combined_cds_introns) <= len(self.combined_cds_introns)

        self.selected_cds_introns.update(transcript.selected_cds_introns)

        self.exons.update(set(transcript.exons))
        assert isinstance(self.locus_verified_introns, set)
        assert isinstance(transcript.verified_introns, set)
        self.locus_verified_introns = set.union(self.locus_verified_introns,
                                                transcript.verified_introns)
        if transcript.verified_introns_num > 0:
            assert len(self.locus_verified_introns) > 0

        if self.initialized is False:
            self.initialized = True
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
            raise KeyError("Transcript {0} is not present in the Locus.".format(tid))

        if len(self.transcripts) == 1:
            self.transcripts = dict()
            self.introns, self.exons, self.splices = set(), set(), set()
            self.cds_introns, self.selected_cds_introns = set(), set()
            self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
            self.stranded = True
            self.initialized = False

        else:
            keys = list(key for key in self.transcripts if key != tid)
            self.end = max(self.transcripts[t].end for t in self.transcripts if t != tid)
            self.start = min(self.transcripts[t].start for t in self.transcripts if t != tid)

            # Remove excess exons
            other_exons = [
                set(self.transcripts[otid].exons if
                    otid in self.transcripts else []) for otid in keys]
            other_exons = set.union(*other_exons)
            exons_to_remove = set.difference(set(self.transcripts[tid].exons), other_exons)
            self.exons.difference_update(exons_to_remove)

            # Remove excess introns
            other_introns = set.union(
                *[set(self.transcripts[otid].introns if otid in self.transcripts else [])
                  for otid in keys])
            introns_to_remove = set.difference(set(self.transcripts[tid].introns), other_introns)
            self.introns.difference_update(introns_to_remove)

            # Remove excess cds introns
            other_cds_introns = set.union(
                *[set(
                    self.transcripts[otid].combined_cds_introns
                    if otid in self.transcripts else [])
                  for otid in keys])
            for otid in keys:
                if otid in self.transcripts:
                    other_cds_introns.update(set(self.transcripts[otid].combined_cds_introns))

            cds_introns_to_remove = set.difference(
                set(self.transcripts[tid].combined_cds_introns),
                other_cds_introns)
            self.combined_cds_introns.difference_update(cds_introns_to_remove)

            # Remove excess selected_cds_introns
            other_selected_cds_introns = set.union(
                *[set(
                    self.transcripts[otid].selected_cds_introns if otid in self.transcripts else []
                ) for otid in keys])
            selected_cds_introns_to_remove = set.difference(
                set(self.transcripts[tid].selected_cds_introns),
                other_selected_cds_introns)
            self.selected_cds_introns.difference_update(selected_cds_introns_to_remove)

            # Remove excess splices
            other_splices = set.union(
                *[self.transcripts[otid].splices
                  if otid in self.transcripts else set() for otid in keys])
            splices_to_remove = set.difference(self.transcripts[tid].splices, other_splices)
            self.splices.difference_update(splices_to_remove)

            del self.transcripts[tid]
            self.logger.debug("Deleted %s from %s", tid, self.id)
            for tid in self.transcripts:
                self.transcripts[tid].parent = self.id

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

        :returns: boolean flag (True if it has to be considered, False otherwise) and a list of sections of the exons
        which are non-coding.
        """

        cds_segments = sorted(transcript.cds_tree.search(*exon))
        terminal = bool(set.intersection(
            set(exon),
            {transcript.start, transcript.end, transcript.combined_cds_end, transcript.combined_cds_start}))
        if cds_segments == [Interval(*exon)]:
            # It is completely coding
            if terminal is False:
                return False, []
            elif not consider_truncated:
                return False, []
            else:
                return True, cds_segments
        else:
            frags = []
            if cds_segments:
                if cds_segments[0].start > exon[0]:
                    frags.append((exon[0], cds_segments[0].start - 1))
                for before, after in zip(cds_segments[:-1], cds_segments[1:]):
                    frags.append((before.end + 1, max(after.start - 1, before.end + 1)))
                if cds_segments[-1].end < exon[1]:
                    frags.append((cds_segments[-1].end + 1, exon[1]))
            else:
                frags = [Interval(*exon)]
            return True, frags

    @staticmethod
    def _is_exon_retained_in_transcript(exon: tuple,
                                        frags: list,
                                        candidate: Transcript,
                                        consider_truncated=False,
                                        terminal=False):

        """Private static method to verify whether a given exon is a retained intron of the candidate Transcript.
        :param exon: the exon to be considered.
        :type exon: (tuple|Interval)

        :param frags: a list of intervals that are non-coding within the exon.
        :type frags: list[(tuple|Interval)]

        :param candidate: a transcript to be evaluated to verify whether the exon is a retained intron event.
        :type candidate: Transcript

        :param consider_truncated: boolean flag. If set, also terminal exons can be considered for retained intron
        events.
        :type consider_truncated: bool

        :param terminal: whether the exon is at the 3' end.
        :type terminal: bool

        :rtype: bool
        """

        found_exons = sorted(
            candidate.segmenttree.find(exon[0], exon[1], strict=False, value="exon"),
            reverse=(candidate.strand == "-"))
        found_introns = sorted(
            candidate.segmenttree.find(exon[0], exon[1], strict=not consider_truncated, value="intron"),
            reverse=(candidate.strand == "-"))

        if len(found_exons) == 0 or len(found_introns) == 0:
            return False
        elif len(found_exons) == 1 and len(found_introns) == 1:
            found_exons = found_exons.pop()
            found_introns = found_introns.pop()
            if candidate.strand != "-" and found_exons[1] + 1 == found_introns[0]:
                return consider_truncated and terminal
            elif candidate.strand == "-" and found_exons[0] - 1 == found_introns[1]:
                return consider_truncated and terminal
            else:
                return False
        else:
            # Now we have to check whether the matched introns contain both coding and non-coding parts
            assert len(found_exons) >= 2, (found_exons, found_introns)
            for index, exon in enumerate(found_exons[:-1]):
                intron = found_introns[index]
                if candidate.strand == "-":
                    assert intron[1] == exon[0] - 1
                else:
                    assert exon[1] == intron[0] - 1
                for frag in frags:
                    # The fragment is just a sub-section of the exon
                    if (overlap(frag, exon) < exon[1] - exon[0] and
                            overlap(frag, exon, positive=True) == 0 and
                            overlap(frag, intron, positive=True)):
                        return True
                    elif overlap(frag, exon) == exon[1] - exon[0]:
                        return True

        return False

    def find_retained_introns(self, transcript: Transcript):

        """This method checks the number of exons that are possibly retained
        introns for a given transcript.
        A retained intron is defined as an exon which:
        - spans completely an intron of another model *between coding exons*
        - is not completely coding itself
        - if the model is coding, the exon has *part* of the non-coding section lying inside the intron (ie the non-coding section must not be starting in the exonic part).

        If the "pick/run_options/consider_truncated_for_retained" flag in the configuration is set to true,
         an exon will be considered as a retained intron event also if:
         - it is the last exon of the transcript
         - it ends *within* an intron of another model *between coding exons*
         - is not completely coding itself
         - if the model is coding, the exon has *part* of the non-coding section lying inside the intron (ie the non-coding section must not be starting in the exonic part).

        The results are stored inside the transcript instance, in the "retained_introns" tuple.
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
            is_retained = False
            to_consider, frags = self._exon_to_be_considered(
                exon, transcript, consider_truncated=consider_truncated)
            if not to_consider:
                continue

            if exon[0] == transcript.start and transcript.strand == "-":
                terminal = True
            elif exon[1] == transcript.end:
                terminal = True
            else:
                terminal = False

            for candidate in (_ for _ in self.transcripts.values()):
                if candidate == transcript:
                    continue
                elif candidate.strand != transcript.strand and None not in (transcript.strand, candidate.strand):
                    continue

                is_retained = self._is_exon_retained_in_transcript(exon,
                                                                   frags,
                                                                   # transcript,
                                                                   candidate,
                                                                   terminal=terminal,
                                                                   consider_truncated=consider_truncated)
                if is_retained:
                    self.logger.debug("Exon %s of %s is a retained intron of %s",
                                      exon, transcript.id, candidate.id)
                    retained_introns.append(exon)
                    break

        transcript.retained_introns = tuple(sorted(retained_introns))
        return

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

        # TODO: Find an intelligent way ot restoring this check
        # I disabled it because otherwise the values for surviving transcripts would be wrong
        # But this effectively leads to a doubling of run time. A possibility would be to cache the results.
        if self.metrics_calculated is True:
            return

        assert len(self._cds_introntree) == len(self.combined_cds_introns)

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

        if len(self.locus_verified_introns) and self.transcripts[tid].verified_introns_num > 0:
            assert self.transcripts[tid].proportion_verified_introns_inlocus > 0

        _ = len(set.intersection(self.exons, self.transcripts[tid].exons))
        fraction = _ / len(self.exons)

        self.transcripts[tid].exon_fraction = fraction

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
                self.transcripts[tid].combined_cds_locus_fraction = self.transcripts[tid].combined_cds_length / cds_bases
                self.transcripts[tid].selected_cds_locus_fraction = self.transcripts[tid].selected_cds_length / selected_bases

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
        self.__json_conf = conf

    def _check_json(self):
        """Private method to be invoked to verify that the configuration is correct.
        Quite expensive to run, especially if done multiple times."""

        conf = self.__json_conf
        if conf is None or isinstance(conf, (str, bytes)):
            conf = to_json(conf)
        elif isinstance(conf, dict):
            conf = check_json(conf)
        else:
            raise TypeError("Unrecognized type for configuration: {}".format(type(conf)))
        self.__json_conf = conf

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
