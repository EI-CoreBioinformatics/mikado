import operator
import abc
import random
import logging
import networkx
from mikado_lib.exceptions import NotInLocusError


class Abstractlocus(metaclass=abc.ABCMeta):
    """This abstract class defines the basic features of any locus-like object.
    It also defines methods/properties that are needed throughout the program,
    e.g. the Bron-Kerbosch algorithm for defining cliques, or the find_retained_introns method."""

    __name__ = "Abstractlocus"

    ###### Special methods #########

    @abc.abstractmethod
    def __init__(self):
        self.transcripts = dict()
        self.introns, self.exons, self.splices = set(), set(), set()
        # Consider only the CDS part
        self.combined_cds_introns, self.selected_cds_introns = set(), set()
        self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
        self.stranded = True
        self.initialized = False
        self.monoexonic = False

    @abc.abstractmethod
    def __str__(self):
        raise NotImplementedError("This is an abstract class and it cannot be printed directly!")

    def __repr__(self):

        return "\t".join([self.__name__,
                          self.chrom,
                          str(self.start),
                          str(self.end),
                          self.strand,
                          ",".join(list(self.transcripts.keys())) if len(self.transcripts) > 0 else "NA"])

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        for feature in ["chrom", "strand", "start", "end", "exons", "introns", "splices", "stranded"]:
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

        state = self.__dict__.copy()
        if hasattr(self, "json_dict"):
            if "requirements" in self.json_dict and "compiled" in self.json_dict["requirements"]:
                del state["json_dict"]["requirements"]["compiled"]

        if hasattr(self, "session"):
            self.session.expunge_all()
            state["session"].expunge_all()
            del state["sessionmaker"]
            del state["session"]

        if hasattr(self, "engine"):
            del state["engine"]

        if hasattr(self, "logger"):
            del state["logger"]

        return state

    def __setstate__(self, state):
        """Method to recreate the object after serialisation."""
        self.__dict__.update(state)
        if hasattr(self, "json_dict"):
            if "requirements" in self.json_dict and "expression" in self.json_dict["requirements"]:
                self.json_dict["requirements"]["compiled"] = compile(self.json_dict["requirements"]["expression"],
                                                                     "<json>", "eval")
        self.set_logger(None)

    def set_logger(self, logger):
        """Set a logger for the instance."""
        if logger is None:
            logger = self.create_default_logger()
        self.logger = logger

    # #### Static methods #######
    @staticmethod
    def overlap(a: tuple([int, int]), b: tuple([int, int]), flank=0) -> int:
        """

        :param a: a tuple of integers
        :type a: (int,int)

        :param b: a tuple of integers
        :type b: (int,int)

        :param flank: an optional extending parameter to check for neighbours
        :type flank: int

        This static method returns the overlap between two intervals. Values<=0 indicate no overlap.
        The optional "flank" argument (default 0) allows to expand a superlocus upstream and downstream.
        As a static method, it can be used also outside of any instance - "superlocus.overlap()" will function.
        Input: two 2-tuples of integers.
        """
        a = sorted(a)
        b = sorted(b)

        left_boundary = max(a[0] - flank, b[0] - flank)
        right_boundary = min(a[1] + flank, b[1] + flank)

        return right_boundary - left_boundary

    @staticmethod
    def evaluate(param: str, conf: dict) -> bool:

        """
        :param param: string to be checked according to the expression in the configuration
        :type param: str

        :param conf: a dictionary containing the expressions to evaluate
        :type conf: dict

        This static method evaluates a single parameter using the requested operation from the JSON dict file.
        """

        if conf["operator"] == "eq":
            return float(param) == float(conf["value"])
        elif conf["operator"] == "ne":
            return float(param) != float(conf["value"])
        elif conf["operator"] == "gt":
            return float(param) > float(conf["value"])
        elif conf["operator"] == "lt":
            return float(param) < float(conf["value"])
        elif conf["operator"] == "ge":
            return float(param) >= float(conf["value"])
        elif conf["operator"] == "le":
            return float(param) <= float(conf["value"])
        elif conf["operator"] == "in":
            return param in conf["value"]
        elif conf["operator"] == "not in":
            return param not in conf["value"]

    # #### Class methods ########


    @classmethod
    def create_default_logger(cls):
        """Static method to create a default logging instance for the loci.
        The default is a null handler (no log)
        """

        formatter = logging.Formatter("{asctime} - {levelname} - {lineno} - {funcName} - {processName} - {message}",
                                      style="{"
                                      )

        logger = logging.getLogger("{0}_logger".format(cls.__name__))
        handler = logging.NullHandler()
        handler.setFormatter(formatter)
        logger.setLevel(logging.WARN)
        logger.addHandler(handler)
        return logger

    @classmethod
    def in_locus(cls, locus_instance, transcript, flank=0) -> bool:
        """
        :param locus_instance: an inheritor of this class
        :type locus_instance: Abstractlocus

        :param transcript: a transcript instance
        :type transcript: mikado_lib.loci_objects.transcript.Transcript

        :param flank: an optional extending parameter to check for neighbours
        :type flank: int

        Function to determine whether a transcript should be added or not to the locus_instance.
        This is a class method, i.e. it can be used also unbound from any specific instance of the class.
        It will be possible therefore to use it to compare any locus_instance to any transcript.
        Arguments: 
        - a "locus_instance" object
        - a "transcript" object (it must possess the "finalize" method)
        - flank - optional keyword"""

        transcript.finalize()
        # We want to check for the strand only if we are considering the strand
        if locus_instance is None:
            return False
        if locus_instance.chrom == transcript.chrom and \
                (locus_instance.stranded is False or locus_instance.strand == transcript.strand) and \
                        cls.overlap((locus_instance.start, locus_instance.end), (transcript.start, transcript.end),
                                    flank=flank) > 0:
            return True
        return False

    @classmethod
    def define_graph(cls, objects: dict, inters=None, **kwargs) -> networkx.Graph:
        """
        :param objects: a dictionary of objects to be grouped into a graph
        :type objects: dict[str]=object

        :param inters: the intersecting function to be used to define the graph
        :type inters: function

        :param kwargs: optional arguments to be passed to the inters function
        :type kwargs: dict

        This function will compute the graph which will later be used by find_communities.
        The method takes as mandatory inputs the following:
            - "objects" a dictionary of objects that form the graph
            - "inters" a function/method that determines whether two objects are connected or not.

        It will then return a graph.
        The method accepts also kwargs that can be passed to the inters function.
        WARNING: the kwargs option is really stupid and does not check for correctness of the arguments!        
        """

        if inters is None:
            inters = cls.is_intersecting

        graph = networkx.Graph()

        # As we are using intern for transcripts, this should prevent memory usage to increase too much
        graph.add_nodes_from( objects.keys() )

        for obj in objects:
            for other_obj in filter(lambda x: x != obj, objects):
                if inters(objects[obj], objects[other_obj], **kwargs):
                    graph.add_edge(*tuple(sorted([obj, other_obj])))  # Connections are not directional

        return graph

    @classmethod
    def find_communities(cls, graph: networkx.Graph) -> (list, list):
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
        cliques = [frozenset(x) for x in networkx.find_cliques(graph)]
        communities = list(networkx.k_clique_communities(graph, 2, cliques)) + [frozenset(x) for x in cliques if
                                                                                len(x) == 1]
        return cliques, communities

    @classmethod
    def find_cliques(cls, objects: list, inters=None) -> (networkx.Graph, list):
        """
        :param objects: list of objects to find the cliques of.
        :type objects: list

        :param inters: the intersecting function to be used
        :type inters: function

        Wrapper for the BronKerbosch algorithm, which returns the maximal cliques in the graph.
        It is the new interface for the BronKerbosch function, which is not called directly from outside this class any longer.
        The "inters" keyword provides the function used to determine whether two vertices are connected or not in the graph.
        """
        if inters is None:
            inters = cls.is_intersecting
        assert hasattr(inters, "__call__")

        graph = dict()
        for obj in objects:
            graph[obj] = [other_obj for other_obj in objects if (obj != other_obj) and inters(obj, other_obj) is True]

        ngraph = networkx.Graph()
        ngraph.add_nodes_from(list(graph.keys()))
        for node in graph:
            for other_node in graph[node]:
                ngraph.add_edge(node, other_node)
        graph = ngraph
        del ngraph

        final_cliques = list(networkx.find_cliques(graph))
        final_cliques = [set(x) for x in final_cliques]

        return graph, final_cliques

    @classmethod
    def choose_best(cls, transcripts: dict) -> str:
        """
        :param transcripts: the dictionary of transcripts of the instance
        :type transcripts: dict

        Given a transcript dictionary, this function will choose the one with the highest score.
        If multiple transcripts have exactly the same score, one will be chosen randomly.

        """

        # Choose one transcript randomly between those that have the maximum score
        return random.choice(
            list(filter(lambda t: t.score == max(transcripts.values(), key=operator.attrgetter("score")).score,
                        transcripts.values()))).id

    ####### Class instance methods  #######

    def add_transcript_to_locus(self, transcript, check_in_locus=True):
        """
        :param transcript
        :type transcript: mikado_lib.loci_objects.transcript.Transcript

        :param check_in_locus: flag to indicate whether the function should check the transcript before adding it
        or instead whether to trust the assignment to be correct
        :type check_in_locus: bool

        This method checks that a transcript is contained within the superlocus (using the "in_superlocus" class method) and
        upon a successful check extends the superlocus with the new transcript.
        More precisely, it updates the boundaries (start and end), adds the transcript to the internal "transcripts" store,
        and extends the splices and introns with those found inside the transcript.
        """

        transcript.finalize()
        self.monoexonic = self.monoexonic and transcript.monoexonic

        if self.initialized is True:
            if check_in_locus is False:
                pass
            elif not self.in_locus(self, transcript):
                raise NotInLocusError("""Trying to merge a locus with an incompatible transcript!
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

        self.combined_cds_introns = set.union(self.combined_cds_introns, transcript.combined_cds_introns)
        self.selected_cds_introns.update(transcript.selected_cds_introns)

        self.exons.update(set(transcript.exons))

        if self.initialized is False:
            self.initialized = True
        self.source = transcript.source
        #         self.source = "Mikado"
        return

    def remove_transcript_from_locus(self, tid: str):
        """
        :param tid: name of the transcript to remove
        :type tid: str

         This method will remove a transcript from an Abstractlocus-like instance and reset all derived
         attributes (e.g. introns, start, end, etc.) appropriately.
        """

        if tid not in self.transcripts:
            raise KeyError("Transcript {0} is not present in the locus.".format(tid))

        if len(self.transcripts) == 1:
            self.transcripts = dict()
            self.introns, self.exons, self.splices = set(), set(), set()
            self.cds_introns, self.selected_cds_introns = set(), set()
            self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
            self.stranded = True
            self.initialized = False

        else:
            keys = list(filter(lambda x: x != tid, self.transcripts.keys()))
            self.end = max(self.transcripts[t].end for t in self.transcripts if t != tid)
            self.start = min(self.transcripts[t].start for t in self.transcripts if t != tid)

            # Remove excess exons
            other_exons = [set(self.transcripts[otid].exons if otid in self.transcripts else []) for otid in keys]
            other_exons = set.union(*other_exons)
            exons_to_remove = set.difference(set(self.transcripts[tid].exons), other_exons)
            self.exons.difference_update(exons_to_remove)

            # Remove excess introns
            other_introns = set.union(
                *[set(self.transcripts[otid].introns if otid in self.transcripts else []) for otid in keys])
            introns_to_remove = set.difference(set(self.transcripts[tid].introns), other_introns)
            self.introns.difference_update(introns_to_remove)

            # Remove excess cds introns
            other_cds_introns = set.union(
                *[set(self.transcripts[otid].combined_cds_introns if otid in self.transcripts else []) for otid in
                  keys])
            for otid in keys:
                if otid in self.transcripts:
                    other_cds_introns.update(set(self.transcripts[otid].combined_cds_introns))

            cds_introns_to_remove = set.difference(set(self.transcripts[tid].combined_cds_introns), other_cds_introns)
            self.combined_cds_introns.difference_update(cds_introns_to_remove)

            # Remove excess selected_cds_introns
            other_selected_cds_introns = set.union(
                *[set(self.transcripts[otid].selected_cds_introns if otid in self.transcripts else []) for otid in
                  keys])
            selected_cds_introns_to_remove = set.difference(set(self.transcripts[tid].selected_cds_introns),
                                                            other_selected_cds_introns)
            self.selected_cds_introns.difference_update(selected_cds_introns_to_remove)

            # Remove excess splices
            other_splices = set.union(
                *[self.transcripts[otid].splices if otid in self.transcripts else set() for otid in keys])
            splices_to_remove = set.difference(self.transcripts[tid].splices, other_splices)
            self.splices.difference_update(splices_to_remove)

            del self.transcripts[tid]
            for tid in self.transcripts:
                self.transcripts[tid].parent = self.id

    def find_retained_introns(self, transcript_instance):

        """
        :param transcript_instance: the transcript to be searched for retained introns.
        :type transcript_instance: mikado_lib.loci_objects.transcript.Transcript

        This method checks the number of exons that are possibly retained introns for a given transcript.
        To perform this operation, it checks for each non-CDS exon whether it exists a sublocus intron that
        is *completely* contained within a transcript exon.
        CDS exons are ignored because their retention might be perfectly valid.
        The results are stored inside the transcript instance, in the "retained_introns" tuple.
        """

        transcript_instance.retained_introns = []
        for exon in filter(lambda e: e not in transcript_instance.combined_cds, transcript_instance.exons):
            # Check that the overlap is at least as long as the minimum between the exon and the intron.
            if any(filter(
                    lambda junction: self.overlap(exon, junction) >= junction[1] - junction[0],
                    self.introns
            )) is True:
                transcript_instance.retained_introns.append(exon)
        transcript_instance.retained_introns = tuple(transcript_instance.retained_introns)

    @classmethod
    @abc.abstractmethod
    def is_intersecting(self):
        """
        This class method defines how two transcript objects will be considered as overlapping.
        It is used by the BronKerbosch method, and must be implemented at the class level for each child object.
        """
        raise NotImplementedError("The is_intersecting method should be defined for each child!")

    ###### Properties #######

    @property
    def stranded(self):
        """This property determines whether a Monosublocus will consider the strand for e.g. the in_locus method.
        By default, the parameter is set to True (i.e. the loci are strand-specific).
        At the moment, the only class which modifies the parameter is the superlocus class."""
        return self.__stranded

    @stranded.setter
    def stranded(self, *args):
        """
        :param args: boolean value
        :type args: list(bool)
        :type args: bool
        """
        if len(args) == 0:
            args = [True]
        stranded = args[0]
        if type(stranded) != bool:
            raise ValueError("The stranded attribute must be boolean!")
        self.__stranded = stranded

    @property
    def id(self) -> str:
        return "{0}:{1}{2}:{3}-{4}".format(
            self.__name__,
            self.chrom,
            self.strand,
            self.start,
            self.end)

    @property
    def name(self) -> str:
        return self.id
