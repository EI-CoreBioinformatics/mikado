# coding: utf-8

"""
Module that defines the blueprint for all loci classes.
"""

import operator
import abc
import random
import logging
from sys import maxsize
from .clique_methods import find_cliques, find_communities, define_graph
import intervaltree
import networkx
from ..exceptions import NotInLocusError
from ..utilities import overlap
from ..utilities.log_utils import create_null_logger
from .transcript import Transcript
from sklearn.ensemble import RandomForestRegressor


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

    @abc.abstractmethod
    def __init__(self, source=""):

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
        self.json_conf = dict()
        self.__cds_introntree = intervaltree.IntervalTree()
        self.__regressor = None
        self.session = None

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
            if "requirements" in self.json_conf and "compiled" in self.json_conf["requirements"]:
                del state["json_conf"]["requirements"]["compiled"]

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
        self.logger = create_null_logger(self)

    # #### Static methods #######
    @staticmethod
    def overlap(first_interval: tuple([int, int]),
                second_interval: tuple([int, int]), flank=0) -> int:
        """

        :param first_interval: a tuple of integers
        :type first_interval: (int,int)

        :param second_interval: a tuple of integers
        :type second_interval: (int,int | intervaltree.Interval)

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

        return overlap(first_interval, second_interval, flank)

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

        return define_graph(objects, inters)

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

    def add_transcript_to_locus(self, transcript, check_in_locus=True):
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
            elif not self.in_locus(self, transcript):
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
            for tid in self.transcripts:
                self.transcripts[tid].parent = self.id

    def find_retained_introns(self, transcript):

        """This method checks the number of exons that are possibly retained
        introns for a given transcript.
        A retained intron is defined as an exon which:

         - spans completely an intron of another model *between coding exons*
         - is not completely coding itself
         - has *part* of the non-coding section lying inside the intron

        The results are stored inside the transcript instance,
        in the "retained_introns" tuple.

        :param transcript: a Transcript instance
        :type transcript: Transcript

        :returns : None
        :rtype : None
        """

        # introns = intervaltree.IntervalTree([
        #     intervaltree.Interval(*intron) for intron in self.combined_cds_introns
        # ])

        self.logger.debug("Starting to calculate retained introns for %s", transcript.id)
        if len(self._cds_introntree) == 0:
            transcript.retained_introns = tuple()
            self.logger.debug("No CDS intron found in %s, exiting for %s",
                              self.id,
                              transcript.id)
            return

        # self.logger.debug("Introns: %d (%d orig, %d (%d, %d CDS segs) transcript)",
        #                   len(self._cds_introntree), len(self.combined_cds_introns),
        #                   len(transcript.combined_cds_introns),
        #                   len(transcript.selected_cds_introns),
        #                   len(transcript.selected_cds))

        retained_introns = []
        if len(self.introns) == 0:
            self.logger.debug("No introns in the locus to check against. Exiting.")

        if transcript.cds_tree is None:
            # Enlarge the CDS segments if they are of length 1
            transcript.cds_tree = intervaltree.IntervalTree(
                [intervaltree.Interval(cds[0], max(cds[1], cds[0]+1))
                 for cds in transcript.combined_cds])

        for exon in iter(_ for _ in transcript.exons if _ not in transcript.combined_cds):
            # Monobase exons are a problem
            if exon[0] == exon[1]:
                self.logger.debug("Monobase exon found in %s: %s:%d-%d",
                                  self.id, self.chrom, exon[0], exon[1])
                continue

            exon_interval = intervaltree.IntervalTree([
                intervaltree.Interval(exon[0], max(exon[1], exon[0] + 1))])

            for cds_segment in transcript.cds_tree.search(*exon):
                exon_interval.chop(cds_segment[0], cds_segment[1])

            # Exclude from consideration any exon which is fully coding
            for frag in exon_interval:
                self.logger.debug("Checking %s from exon %s for retained introns for %s",
                                  frag, exon, transcript.id)
                if self._cds_introntree.overlaps_range(frag[0], frag[1]):
                    self.logger.debug("Exon %s of %s is a retained intron",
                                      exon, transcript.id)
                    retained_introns.append(exon)
                    break

        # Sort the exons marked as retained introns
        # self.logger.info("Finished calculating retained introns for %s", transcript.id)
        transcript.retained_introns = tuple(sorted(retained_introns))
        # self.logger.info("Returning retained introns for %s", transcript.id)
        # return transcript

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
            assert row != {}
            yield row

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
        if logger is None:
            logger = create_null_logger(self)
        elif not isinstance(logger, logging.Logger):
            raise TypeError("Invalid logger: {0}".format(type(logger)))
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
        :rtype: intervaltree.IntervalTree
        """

        if len(self.__cds_introntree) != len(self.combined_cds_introns):
            self.__cds_introntree = intervaltree.IntervalTree(
                [intervaltree.Interval(_[0], _[1] + 1) for _ in self.combined_cds_introns])
        return self.__cds_introntree

    @property
    def longest_transcript(self):
        return max([len(_) for _ in self.transcripts.values()])

    @property
    def regressor(self):
        return self.__regressor

    @regressor.setter
    def regressor(self, regr):

        if not (isinstance(regr, RandomForestRegressor) or regr is None):
            raise TypeError("Invalid regressor provided, type: %s", type(regr))

        self.logger.debug("Set regressor")
        self.__regressor = regr
