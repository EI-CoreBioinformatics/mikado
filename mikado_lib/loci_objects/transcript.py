# coding: utf-8

"""
This module defines the RNA objects. It also defines Metric, a property alias.
"""

# pylint: disable=bad-builtin, too-many-lines

import operator
import logging
import sys
import re
import functools
from collections import OrderedDict, Counter
import inspect
# import asyncio
from mikado_lib.log_utils import create_null_logger
from mikado_lib.exceptions import InvalidTranscript
# SQLAlchemy imports
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.sql.expression import desc, asc
from sqlalchemy import and_
from sqlalchemy.ext import baked
from sqlalchemy import bindparam
# mikado imports
from mikado_lib.serializers.junction import Junction
import mikado_lib.serializers.orf
from mikado_lib.serializers.blast_serializer import Query, Hit
from mikado_lib.serializers.orf import Orf
from mikado_lib.parsers import bed12
from mikado_lib.loci_objects.abstractlocus import Abstractlocus
from mikado_lib.parsers.GTF import GtfLine
from mikado_lib.parsers.GFF import GffLine
import mikado_lib.exceptions


# from memory_profiler import profile
# import logging
# if "line_profiler" not in dir():
#     def profile(function):
#         """
#         Mock wrapper to imitate the profile decorator
#         :param function: the function to be wrapped
#         :return:
#         """
#         def inner(*args, **kwargs):
#             """
#             Returns the wrapped function
#             :param args: arguments to be passed
#             :param kwargs: keyword arguments to be passed
#             :return:
#             """
#             return function(*args, **kwargs)
#         return inner


class Metric(property):
    """Simple aliasing of property. All transcript metrics
    should use this alias, not "property", as a decorator."""
    pass


# noinspection PyPropertyAccess
# I do not care that there are too many attributes: this IS a massive class!
# pylint: disable=too-many-instance-attributes,too-many-public-methods
class Transcript:
    """
    This class defines a transcript, down to its exon/CDS/UTR components.
    It is instantiated by a transcript GTF/GFF3 line.
    Key attributes:

    :param chrom: The chromosome of the transcript
    :type chrom: str
    :type source: str
    :param feature: mRNA if at least one CDS is defined, else use the one
    derived from input; default is "transcript"
    :type feature: str
    :param start: Start of the transcript. Checked against the exons.
    :type start: int
    :param end: End of the transcript. Checked against the exons.
    :type end: int
    :param score: The score assigned to the transcript. Modified inside Mikado.
    :type score: float
    :param strand: one of +,-,None
    :type strand: str
    :param id            the ID of the transcripts (or tid)
    :type id: str
    :param parent: The parent leaves of the transcript
    :type parent: list
    :param attributes: a dictionary with additional informations from the GFFline
    :type attributes: dict

    After all exons have been loaded into the instance (see "addExon"),
    the class must be finalized with the appropriate method.
    CDS locations can be uploaded from the external, using a dictionary of indexed BED12 entries.
    The database queries are baked at the *class* level in order to minimize overhead.
    """

    __name__ = "transcript"
    __logger = create_null_logger(__name__)

    # Query baking to minimize overhead
    bakery = baked.bakery()
    query_baked = bakery(lambda session: session.query(Query))
    query_baked += lambda q: q.filter(Query.query_name == bindparam("query_name"))

    blast_baked = bakery(lambda session: session.query(Hit))
    blast_baked += lambda q: q.filter(and_(Hit.query_id == bindparam("query_id"),
                                           Hit.evalue <= bindparam("evalue")),)

    blast_baked += lambda q: q.order_by(asc(Hit.evalue))
    blast_baked += lambda q: q.limit(bindparam("max_target_seqs"))

    orf_baked = bakery(lambda session: session.query(Orf))
    orf_baked += lambda q: q.filter(
        mikado_lib.serializers.orf.Orf.query_id == bindparam("query_id"))
    orf_baked += lambda q: q.filter(
        mikado_lib.serializers.orf.Orf.cds_len >= bindparam("cds_len"))
    orf_baked += lambda q: q.order_by(desc(Orf.cds_len))

    # ######## Class special methods ####################

    def __init__(self, *args,
                 source=None,
                 logger=None,
                 intron_range=(0, sys.maxsize)):

        """Initialise the transcript object, using a mRNA/transcript line.
        Note: I am assuming that the input line is an object from my own "GFF" class.
        The transcript instance must be initialised by a "(m|r|lnc|whatever)RNA" or
        "transcript" GffLine.

        :param intron_range: range of valid intron size. Any intron shorter
        or longer than this will be flagged.
        :type intron_range: list(int,int)

        """

        # Mock setting of base hidden variables
        self.__id = ""
        self.__strand = self.__score = None
        self.__has_start_codon, self.__has_stop_codon = False, False
        self.__max_internal_orf_index = None
        self.__max_internal_orf_length = self.__intron_fraction = self.__exon_fraction = 0
        # Metrics might have queer names
        # pylint: disable=invalid-name
        self.__proportion_verified_introns_inlocus = 0
        self.__retained_fraction = 0
        self.__combined_cds_intron_fraction = self.__selected_cds_intron_fraction = 0
        self.__non_overlapping_cds = set()
        self.__exons = set()
        self.__parent = []
        self.__combined_cds = []
        self.__selected_cds = []
        self.__combined_utr = []
        # pylint: enable=invalid-name
        self._selected_internal_orf_cds = []
        # This is used to set the phase if the CDS is loaded from the GFF
        self._first_phase = 0
        self.__phases = []  # will contain (start, phase) for each CDS exon

        # Starting settings for everything else
        self.chrom = None
        self.source = source
        self.feature = "transcript"
        self.start, self.end = None, None
        self.attributes = dict()
        self.exons, self.combined_cds, self.combined_utr = [], [], []
        self.logger = logger
        self.introns = []
        self.splices = []
        self.finalized = False  # Flag. We do not want to repeat the finalising more than once.
        self.selected_internal_orf_index = None
        self.non_overlapping_cds = None
        self.verified_introns = set()
        self.segments = []
        self.intron_range = intron_range
        self.internal_orfs = []
        self.blast_hits = []

        # Relative properties
        self.retained_introns = []
        self.retained_fraction = 0
        self.exon_fraction = self.intron_fraction = 1
        self.cds_intron_fraction = self.selected_cds_intron_fraction = 1

        # Json configuration
        self.json_conf = None

        # Things that will be populated by querying the database
        self.loaded_bed12 = []
        self.engine, self.session, self.sessionmaker = None, None, None
        self.query_id = None

        if len(args) == 0:
            return
        else:
            self.__initialize_with_line(args[0])

    def __initialize_with_line(self, transcript_row):
        """
        Private method to copy the necessary attributes from
        an external GTF/GFF3 row.
        :param transcript_row:
        :return:
        """

        if not isinstance(transcript_row, (GffLine, GtfLine)):
            raise TypeError("Invalid data type: {0}".format(type(transcript_row)))

        self.chrom = transcript_row.chrom
        assert transcript_row.is_transcript is True
        self.feature = transcript_row.feature
        # pylint: disable=invalid-name
        self.id = transcript_row.id
        # pylint: enable=invalid-name
        self.name = transcript_row.name
        if self.source is None:
            self.source = transcript_row.source
        self.start = transcript_row.start
        self.strand = transcript_row.strand
        self.end = transcript_row.end
        self.score = transcript_row.score
        self.parent = transcript_row.parent
        self.attributes = transcript_row.attributes
        self.blast_hits = []
        self.json_conf = None

    def __str__(self, to_gtf=False, print_cds=True):
        """
        :type to_gtf: bool
        :type print_cds: bool

        Each transcript will be printed out in the GFF style.
        This is pretty rudimentary, as the class does not hold
        any information on the original source,
        feature, score, etc.
        """

        self.finalize()  # Necessary to sort the exons
        if print_cds is True:
            lines = self.create_lines_cds(to_gtf=to_gtf)
        else:
            lines = self.create_lines_no_cds(to_gtf=to_gtf)

        return "\n".join(lines)

    def __eq__(self, other) -> bool:
        """
        :param other: another transcript instance to compare to
        :type other: mikado_lib.loci_objects.transcript.Transcript

        Two transcripts are considered identical if they have the same
        start, end, chromosome, strand and internal exons.
        IDs are not important for this comparison; two transcripts coming from different
        methods and having different IDs can still be identical."""

        if not isinstance(self, type(other)):
            return False
        self.finalize()
        other.finalize()

        if self.strand == other.strand and self.chrom == other.chrom:
            if other.start == self.start:
                if self.end == other.end:
                    if self.exons == other.exons:
                        return True

        return False

    def __hash__(self):
        """Returns the hash of the object (call to super().__hash__()).
        Necessary to be able to add these objects to hashes like sets.
        """

        return super().__hash__()

    def __len__(self) -> int:
        """Returns the length occupied by the unspliced transcript on the genome."""
        return self.end - self.start + 1

    def __lt__(self, other) -> bool:
        """A transcript is lesser than another if it is on a lexicographic inferior chromosome,
        or if it begins before the other, or (in the case where they begin at the same location)
        it ends earlier than the other.
        """
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        if self == other:
            return False
        if self.start < other.start:
            return True
        elif self.start == other.start and self.end < other.end:
            return True
        return False

    def __gt__(self, other) -> bool:
        return not self < other

    def __le__(self, other) -> bool:
        return (self == other) or (self < other)

    def __ge__(self, other) -> bool:
        return (self == other) or (self > other)

    def __getstate__(self):

        logger = self.logger
        del self.logger
        state = self.__dict__.copy()
        self.logger = logger

        if hasattr(self, "json_conf") and self.json_conf is not None:
            if "requirements" in self.json_conf and "compiled" in self.json_conf["requirements"]:
                del state["json_conf"]["requirements"]["compiled"]

        if hasattr(self, "session"):
            if state["session"] is not None:
                state["session"].expunge_all()
                state["session"].close()

            del state["session"]
        if hasattr(self, "sessionmaker"):
            del state["sessionmaker"]
            del state["engine"]

        if "blast_baked" in state:
            del state["blast_baked"]
            del state["query_baked"]

        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # Set the logger to NullHandler
        self.logger = None

    # ######## Class instance methods ####################

    def add_exon(self, gffline):
        """This function will append an exon/CDS feature to the object.
        :param gffline: an annotation line
        :type gffline: mikado_lib.parsers.GFF.GffLine, mikado_lib.parsers.GTF.GtfLine
        """

        if self.finalized is True:
            raise mikado_lib.exceptions.ModificationError(
                "You cannot add exons to a finalized transcript!")

        if self.id not in gffline.parent:
            raise mikado_lib.exceptions.InvalidTranscript(
                """Mismatch between transcript and exon:
                {0}
                {1}
                {2}""".format(self.id, gffline.parent, gffline))
        assert gffline.is_exon is True, str(gffline)

        if gffline.feature.upper().endswith("CDS"):
            store = self.combined_cds
            self.__phases.append((gffline.start, gffline.phase))
        elif "combined_utr" in gffline.feature or "UTR" in gffline.feature.upper():
            store = self.combined_utr
        elif gffline.feature.endswith("exon"):
            store = self.exons
        elif gffline.feature == "start_codon":
            self.has_start_codon = True
            return
        elif gffline.feature == "stop_codon":
            self.has_stop_codon = True
            return
        elif gffline.feature == "intron":
            store = self.introns
        else:
            raise mikado_lib.exceptions.InvalidTranscript(
                "Unknown feature: {0}".format(gffline.feature))

        start, end = sorted([gffline.start, gffline.end])
        store.append((start, end))

    def __create_cds_lines(self, cds_run, tid, to_gtf=False):

        """
        Private method to create the exon/UTR/CDS lines for printing
        out in GTF/GFF format.
        :param cds_run: the internal orf run we are preparing
        :param tid: name of the transcript
        :param to_gtf: boolean, indicates whether
        we want GTF or GFF output
        :return:
        """

        exon_lines = []
        cds_begin = False
        counter = Counter()

        line_creator = functools.partial(self.__create_exon_line,
                                         **{"to_gtf": to_gtf,
                                            "tid": tid})
        for segment in cds_run:
            exon_line, counter, cds_begin = line_creator(segment,
                                                         counter,
                                                         cds_begin)
            exon_lines.append(exon_line)
        exon_lines = [str(exon_line) for exon_line in
                      self.__add_phase(exon_lines)]

        return exon_lines

    # pylint: disable=too-many-arguments
    def __create_exon_line(self, segment, counter, cds_begin,
                           tid="", to_gtf=False):
        """
        Private method that creates an exon line for printing.
        :param segment: a segment of the form (feature, start, end)
        :type segment: list(str, int, int)
        :param counter: a Counter object that keeps track of how many exons,
        CDS, UTR segments we have already seen
        :type counter: Counter
        :param cds_begin: boolean flag that indicates whether the CDS has already begun
        :type cds_begin: bool
        :param tid: name of the transcript
        :param to_gtf: boolean flag
        :return: exon_line, counter, cds_begin
        :rtype: ((GtfLine | GffLine), Counter, bool)
        """

        if to_gtf is False:
            constructor = GffLine
            utr3_feature = "three_prime_UTR"
            utr5_feature = "five_prime_UTR"
        else:
            constructor = GtfLine
            utr3_feature = "3UTR"
            utr5_feature = "5UTR"

        if cds_begin is False and segment[0] == "CDS":
            cds_begin = True
        if segment[0] == "UTR":
            if (cds_begin is True and self.strand == "-") or \
                    (self.strand == "+" and cds_begin is False):
                feature = utr5_feature
                counter.update(["five"])
                index = counter["five"]
            else:
                feature = utr3_feature
                counter.update(["three"])
                index = counter["three"]
        elif segment[0] == "CDS":
            counter.update(["CDS"])
            index = counter["CDS"]
            feature = "CDS"
        else:
            counter.update(["exon"])
            index = counter["exon"]
            feature = segment[0]
        exon_line = constructor(None)

        for attr in ["chrom", "source", "strand"]:
            setattr(exon_line, attr, getattr(self, attr))

        exon_line.feature = feature
        exon_line.start, exon_line.end = segment[1], segment[2]
        exon_line.phase = None
        exon_line.score = None
        if to_gtf is True:
            # noinspection PyPropertyAccess
            exon_line.gene = self.parent
            exon_line.transcript = tid
        else:
            exon_line.id = "{0}.{1}{2}".format(tid, feature, index)
            exon_line.parent = tid
        return exon_line, counter, cds_begin
    # pylint: enable=too-many-arguments

    def create_lines_cds(self, to_gtf=False):

        """
        Method to create the GTF/GFF lines for printing in the presence of CDS information.
        WARNING: at the moment, the phase support is disabled.

        :param to_gtf:
        :return:
        """

        if to_gtf is False:
            constructor = GffLine
        else:
            constructor = GtfLine

        lines = []
        transcript_counter = 0

        parent_line = constructor(None)
        for index in range(len(self.internal_orfs)):
            if self.number_internal_orfs > 1:
                transcript_counter += 1
                tid = "{0}.orf{1}".format(self.id, transcript_counter)

                if index == self.selected_internal_orf_index:
                    self.attributes["maximal"] = True
                else:
                    self.attributes["maximal"] = False
            else:
                tid = self.id
            cds_run = self.internal_orfs[index]

            for attr in ["chrom", "source", "feature", "start", "end",
                         "score", "strand", "attributes", "parent"]:
                setattr(parent_line, attr, getattr(self, attr))

            parent_line.phase = '.'

            parent_line.id = tid
            parent_line.name = self.id

            exon_lines = self.__create_cds_lines(cds_run,
                                                 tid,
                                                 to_gtf=to_gtf)

            lines.append(str(parent_line))
            lines.extend(exon_lines)
        return lines

    def __add_phase(self, exon_lines):

        """
        Private method to add the phase to a transcript.
        :param exon_lines:
        :return:
        """

        # We start by 0 if no CDS loaded, else
        # we use the first phase

        previous = ((self._first_phase - 1) % 3 + 1) % 3

        new_lines = []
        for line in sorted(exon_lines, reverse=(self.strand == "-")):
            if line.feature == "CDS":
                line.phase = ((previous % 3 - 1) % 3 + 1) % 3
                previous += len(line)
            new_lines.append(line)
        return sorted(new_lines)

    def create_lines_no_cds(self, to_gtf=False):

        """
        Method to create the GTF/GFF lines for printing in the absence of CDS information.
        """

        if to_gtf is True:
            constructor = GtfLine
        else:
            constructor = GffLine

        parent_line = constructor(None)

        for attr in ["chrom", "source", "feature", "start", "end",
                     "score", "strand", "attributes", "parent"]:
            setattr(parent_line, attr, getattr(self, attr))

        parent_line.phase = '.'
        parent_line.id = self.id
        parent_line.name = self.id

        lines = [str(parent_line)]
        exon_lines = []

        exon_count = 0
        for exon in self.exons:
            exon_count += 1
            exon_line = constructor(None)
            for attr in ["chrom", "source", "strand", "attributes"]:
                setattr(exon_line, attr, getattr(self, attr))
            exon_line.feature = "exon"
            exon_line.start, exon_line.end = exon[0], exon[1]
            exon_line.score = None
            exon_line.phase = None

            exon_line.id = "{0}.{1}{2}".format(self.id, "exon", exon_count)
            exon_line.parent = self.id
            exon_lines.append(str(exon_line))

        lines.extend(exon_lines)
        return lines

    @staticmethod
    def orf_sorter(orf):
        """Sorting function for the ORFs."""
        return ((orf.has_start_codon and orf.has_stop_codon),
                (orf.has_start_codon or orf.has_stop_codon),
                orf.cds_len)

    def find_candidate_orfs(self, graph, orf_dictionary) -> list:

        """
        Function that returns the best non-overlapping ORFs
        :param graph: The NetworkX graph to be analysed
        :return:
        """

        candidate_orfs = []

        while len(graph) > 0:
            cliques, communities = Abstractlocus.find_communities(graph)
            clique_str = []
            for clique in cliques:
                clique_str.append(str([(orf_dictionary[x].thick_start,
                                        orf_dictionary[x].thick_end) for x in clique]))
            comm_str = []
            for comm in communities:
                comm_str.append(str([(orf_dictionary[x].thick_start,
                                      orf_dictionary[x].thick_end) for x in comm]))
            self.logger.debug("{0} communities for {1}:\n\t{2}".format(len(communities),
                                                                       self.id,
                                                                       "\n\t".join(
                                                                           comm_str
                                                                       )))
            self.logger.debug("{0} cliques for {1}:\n\t{2}".format(len(cliques),
                                                                   self.id,
                                                                   "\n\t".join(clique_str)))
            to_remove = set()
            for comm in communities:
                comm = [orf_dictionary[x] for x in comm]
                best_orf = sorted(comm, key=self.orf_sorter, reverse=True)[0]
                candidate_orfs.append(best_orf)
                for clique in iter(cl for cl in cliques if best_orf.name in cl):
                    to_remove.update(clique)
            graph.remove_nodes_from(to_remove)

        candidate_orfs = sorted(candidate_orfs,
                                key=self.orf_sorter, reverse=True)
        return candidate_orfs

    def find_overlapping_cds(self, candidates: list) -> list:
        """
        :param candidates: candidate ORFs to analyse
        :type candidates: list(mikado_lib.serializers.orf.Orf)

        Wrapper for the Abstractlocus method, used for finding overlapping ORFs.
        It will pass to the function the class's "is_overlapping_cds" method
        (which would be otherwise be inaccessible from the Abstractlocus class method).
        As we are interested only in the communities, not the cliques,
        this wrapper discards the cliques
        (first element of the Abstractlocus.find_communities results)
        """

        # If we are looking at a multiexonic transcript
        if not (self.monoexonic is True and self.strand is None):
            candidates = list(filter(lambda co: co.strand == "+", candidates))

        # Prepare the minimal secondary length parameter
        if self.json_conf is not None:
            minimal_secondary_orf_length = \
                self.json_conf["orf_loading"]["minimal_secondary_orf_length"]
        else:
            minimal_secondary_orf_length = 0

        self.logger.debug("{0} input ORFs for {1}".format(len(candidates), self.id))
        candidates = list(filter(lambda x: x.invalid is False, candidates))
        self.logger.debug("{0} filtered ORFs for {1}".format(len(candidates), self.id))
        if len(candidates) == 0:
            return []

        orf_dictionary = dict((x.name, x) for x in candidates)

        # First define the graph
        graph = Abstractlocus.define_graph(orf_dictionary, inters=self.is_overlapping_cds)
        candidate_orfs = self.find_candidate_orfs(graph, orf_dictionary)

        self.logger.debug("{0} candidate retained ORFs for {1}: {2}".format(
            len(candidate_orfs),
            self.id,
            [x.name for x in candidate_orfs]))
        final_orfs = [candidate_orfs[0]]
        if len(candidate_orfs) > 1:
            others = list(filter(lambda o: o.cds_len >= minimal_secondary_orf_length,
                                 candidate_orfs[1:]))
            self.logger.debug("Found {0} secondary ORFs for {1} of length >= {2}".format(
                len(others), self.id,
                minimal_secondary_orf_length
            ))
            final_orfs.extend(others)

        self.logger.debug("Retained %d ORFs for %s: %s",
                          len(final_orfs),
                          self.id,
                          [orf.name for orf in final_orfs])
        return final_orfs

    def check_split_by_blast(self, cds_boundaries):

        """
        This method verifies if a transcript with multiple ORFs has support by BLAST to
        NOT split it into its different components.
        :param cds_boundaries:
        :return:
        """

        # Establish the minimum overlap between an ORF and a BLAST hit to consider it
        # to establish belongingness

        minimal_overlap = self.json_conf["chimera_split"]["blast_params"]["minimal_hsp_overlap"]

        cds_hit_dict = OrderedDict().fromkeys(cds_boundaries.keys())
        for key in cds_hit_dict:
            cds_hit_dict[key] = set()

        # BUG, this is a hacky fix
        if not hasattr(self, "blast_hits"):
            self.logger.warning(
                "BLAST hits store lost for %s! Creating a mock one to avoid a crash",

                self.id)
            self.blast_hits = []

        # Determine for each CDS which are the hits available
        for hit in self.blast_hits:
            for hsp in filter(lambda lambda_hsp:
                              lambda_hsp["hsp_evalue"] <=
                              self.json_conf['chimera_split']['blast_params']['hsp_evalue'],
                              hit["hsps"]):

                for cds_run in cds_boundaries:
                    # If I have a valid hit b/w the CDS region and the hit,
                    # add the name to the set
                    overlap_threshold = minimal_overlap * (cds_run[1] + 1 - cds_run[0])
                    if Abstractlocus.overlap(cds_run, (
                            hsp['query_hsp_start'],
                            hsp['query_hsp_end'])) >= overlap_threshold:
                        cds_hit_dict[cds_run].add(hit["target"])

        final_boundaries = OrderedDict()
        for boundary in self.__get_boundaries_from_blast(cds_boundaries, cds_hit_dict):
            if len(boundary) == 1:
                assert len(boundary[0]) == 2
                boundary = boundary[0]
                final_boundaries[boundary] = cds_boundaries[boundary]
            else:
                nboun = (boundary[0][0], boundary[-1][1])
                final_boundaries[nboun] = []
                for boun in boundary:
                    final_boundaries[nboun].extend(cds_boundaries[boun])

        cds_boundaries = final_boundaries.copy()
        return cds_boundaries

    def __get_boundaries_from_blast(self, cds_boundaries, cds_hit_dict):

        """
        Private method that calculates the CDS boundaries to keep
        given the blast hits. Called by check_split_by_blast
        :param cds_boundaries:
        :return:
        """
        new_boundaries = []
        for cds_boundary in cds_boundaries:
            if not new_boundaries:
                new_boundaries.append([cds_boundary])
            else:
                old_boundary = new_boundaries[-1][-1]
                cds_hits = cds_hit_dict[cds_boundary]
                old_hits = cds_hit_dict[old_boundary]
                if cds_hits == set() and old_hits == set():  # No hit found for either CDS
                    # If we are stringent, we DO NOT split
                    if self.json_conf['chimera_split']['blast_params']['leniency'] == "STRINGENT":
                        new_boundaries[-1].append(cds_boundary)
                    else:  # Otherwise, we do split
                        new_boundaries.append([cds_boundary])
                elif cds_hits == set() or old_hits == set():  # We have hits for only one
                    # If we are permissive, we split
                    if self.json_conf["chimera_split"]["blast_params"]["leniency"] == "PERMISSIVE":
                        new_boundaries.append([cds_boundary])
                    else:
                        new_boundaries[-1].append(cds_boundary)
                else:
                    # We do not have any hit in common
                    if set.intersection(cds_hits, old_hits) == set():
                        new_boundaries.append([cds_boundary])
                    # We have hits in common
                    else:
                        new_boundaries[-1].append(cds_boundary)
            # } # Finish BLAST check
        return new_boundaries

    def __split_complex_exon(self, exon, texon, left, right, boundary):

        """
        Private method used to split an exon when it is only partially coding,
        :param exon: Exon to be analysed
        :param texon: Transcriptomic coordinates of the exon
        :param left: boolean flag, it indicates wheter there is another transcript
        on the left of the current one.
        :param right: boolean flag, it indicates wheter there is another transcript
        on the left of the current one.
        :param boundary: Transcriptomic coordinates of the ORF boundary.
        :return:
        """

        to_discard = None
        new_exon = list(exon)

        if texon[1] == boundary[0]:
            # In this case we have that the exon ends exactly at the end of the
            # UTR, so we have to keep a one-base exon
            if left is False:
                self.logger.debug("Appending mixed UTR/CDS 5' exon %s", exon)
            else:
                if self.strand == "+":
                    # Keep only the LAST base
                    to_discard = (exon[0], exon[1]-1)
                    new_exon = (exon[1]-1, exon[1])
                    texon = (texon[1]-1, texon[1])
                    self.logger.debug("Appending monobase CDS exon %s (Texon %s)",
                                      new_exon,
                                      texon)
                else:
                    # Keep only the FIRST base
                    to_discard = (exon[0]+1, exon[1])
                    new_exon = (exon[0], exon[0]+1)
                    texon = (texon[1]-1, texon[1])
                    self.logger.debug(
                        "Appending monobase CDS exon %s (Texon %s)",
                        new_exon,
                        texon)

        elif texon[0] == boundary[1]:
            # In this case we have that the exon ends exactly at the end of the
            # CDS, so we have to keep a one-base exon
            if right is False:
                self.logger.debug(
                    "Appending mixed UTR/CDS right exon %s",
                    exon)
            else:
                if self.strand == "+":
                    # In this case we have to keep only the FIRST base
                    to_discard = (exon[0]+1, exon[1])
                    new_exon = (exon[0], exon[0]+1)
                    texon = (texon[0], texon[0]+1)
                    self.logger.debug(
                        "Appending monobase CDS exon %s (Texon %s)",
                        new_exon,
                        texon)
                else:
                    # In this case we have to keep only the LAST base
                    to_discard = (exon[0], exon[1]-1)
                    new_exon = (exon[1]-1, exon[1])
                    texon = (texon[0], texon[0]+1)
                    self.logger.debug(
                        "Appending monobase CDS exon %s (Texon %s)",
                        new_exon,
                        texon)

        elif texon[0] <= boundary[0] <= boundary[1] <= texon[1]:
            # Monoexonic
            self.logger.debug("Exon %s, case 3.1", exon)
            if self.strand == "-":
                if left is True:
                    new_exon[1] = exon[0] + (texon[1] - boundary[0])
                if right is True:
                    new_exon[0] = exon[1] - (boundary[1] - texon[0])
            else:
                if left is True:
                    new_exon[0] = exon[1] - (texon[1] - boundary[0])
                if right is True:
                    new_exon[1] = exon[0] + (boundary[1] - texon[0])
            self.logger.debug(
                "[Monoexonic] Tstart shifted for %s, %d to %d",
                self.id, texon[0], boundary[0])
            self.logger.debug(
                "[Monoexonic] GStart shifted for %s, %d to %d",
                self.id, exon[0], new_exon[1])
            self.logger.debug(
                "[Monoexonic] Tend shifted for %s, %d to %d",
                self.id, texon[1], boundary[1])
            self.logger.debug(
                "[Monoexonic] Gend shifted for %s, %d to %d",
                self.id, exon[1], new_exon[1])

            if left is True:
                texon[0] = boundary[0]
            if right is True:
                texon[1] = boundary[1]

        elif texon[0] <= boundary[0] <= texon[1] <= boundary[1]:
            # In this case we have that exon is sitting halfway
            # i.e. there is a partial 5'UTR
            if left is True:
                if self.strand == "-":
                    new_exon[1] = exon[0] + (texon[1] - boundary[0])
                else:
                    new_exon[0] = exon[1] - (texon[1] - boundary[0])
                self.logger.debug(
                    "Tstart shifted for %s, %d to %d", self.id, texon[0], boundary[0])
                self.logger.debug(
                    "GStart shifted for %s, %d to %d", self.id, exon[0], new_exon[1])
                texon[0] = boundary[0]

        elif texon[1] >= boundary[1] >= texon[0] >= boundary[0]:
            # In this case we have that exon is sitting halfway
            # i.e. there is a partial 3'UTR
            if right is True:
                if self.strand == "-":
                    new_exon[0] = exon[1] - (boundary[1] - texon[0])
                else:
                    new_exon[1] = exon[0] + (boundary[1] - texon[0])
                self.logger.debug(
                    "Tend shifted for %s, %d to %d",
                    self.id, texon[1], boundary[1])
                self.logger.debug(
                    "Gend shifted for %s, %d to %d",
                    self.id, exon[1], new_exon[1])
                texon[1] = boundary[1]
            else:
                self.logger.debug("New exon: %s", new_exon)
                self.logger.debug("New texon: %s", texon)

        return new_exon, texon, to_discard

    def __create_splitted_exons(self, boundary, left, right):

        """
        Given a boundary in transcriptomic coordinates, this method will extract the
        exons retained in the splitted part of the model.

        :param boundary: the *transcriptomic* coordinates of start/end of the ORF(s)
        to be included in the new transcript
        :type boundary: (int,int)

        :param left: boolean flag indicating whether there is another sub-transcript
        to the left of the one we mean to create, irrespective of *genomic* strand
        :type left: bool

        :param left: boolean flag indicating whether there is another sub-transcript
        to the right of the one we mean to create, irrespective of *genomic* strand
        :type right: bool


        :return: my_exons (final exons), discarded_exons (eventual discarded exons),
        tstart (new transcript start), tend (new transcript end)
        :rtype: (list(int,int),list(int,int),int,int)
        """

        my_exons = []

        discarded_exons = []
        tlength = 0
        tstart = float("Inf")
        tend = float("-Inf")

        if self.strand == "-":
            reversal = True
        else:
            reversal = False

        for exon in sorted(self.exons, key=operator.itemgetter(0), reverse=reversal):
            # Translate into transcript coordinates
            elength = exon[1] - exon[0] + 1
            texon = [tlength + 1, tlength + elength]
            tlength += elength
            self.logger.debug("Analysing exon %s [%s] for %s",
                              exon, texon, self.id)

            # SIMPLE CASES
            # Exon completely contained in the ORF
            if boundary[0] <= texon[0] < texon[1] <= boundary[1]:
                self.logger.debug("Appending CDS exon %s", exon)
                my_exons.append(exon)
            # Exon on the left of the CDS
            elif texon[1] < boundary[0]:
                if left is False:
                    self.logger.debug("Appending 5'UTR exon %s", exon)
                    my_exons.append(exon)
                else:
                    self.logger.debug("Discarding 5'UTR exon %s", exon)
                    discarded_exons.append(exon)
                    continue
            elif texon[0] > boundary[1]:
                if right is False:
                    self.logger.debug("Appending 3'UTR exon %s", exon)
                    my_exons.append(exon)
                else:
                    self.logger.debug("Discarding 3'UTR exon %s", exon)
                    discarded_exons.append(exon)
                    continue
            else:
                # exon with partial UTR, go to the relative function
                # to handle these complex cases
                exon, texon, to_discard = self.__split_complex_exon(
                    exon, texon, left, right, boundary)
                my_exons.append(tuple(sorted(exon)))
                if to_discard is not None:
                    discarded_exons.append(to_discard)

            tstart = min(tstart, texon[0])
            tend = max(tend, texon[1])

        return my_exons, discarded_exons, tstart, tend

    @staticmethod
    def relocate_orfs(bed12_objects, tstart, tend):
        """
        Function to recalculate the coordinates of BED12 objects based on
        the new transcriptomic start/end
        :param bed12_objects: list of the BED12 ORFs to relocate
        :param tstart: New transcriptomic start
        :param tend: New transcriptomic end
        :return:
        """
        new_bed12s = []
        for obj in bed12_objects:
            assert isinstance(obj, bed12.BED12), (obj, bed12_objects)

            obj.start = 1
            obj.end = min(obj.end, tend) - tstart + 1
            obj.fasta_length = obj.end
            obj.thick_start = min(obj.thick_start, tend) - tstart + 1
            obj.thick_end = min(obj.thick_end, tend) - tstart + 1
            obj.blockSizes = [obj.end]
            assert obj.invalid is False, (len(obj), obj.cds_len, obj.fasta_length,
                                          obj.invalid_reason,
                                          str(obj))
            new_bed12s.append(obj)
        return new_bed12s

    def __check_collisions(self, nspan, spans):

        """
        This method checks whether a new transcript collides with a previously
        defined transcript.
        :param nspan:
        :param spans:
        :return:
        """

        if len(spans) == 0:
            return
        for span in spans:
            overl = Abstractlocus.overlap(span, nspan)

            self.logger.debug(
                "Comparing start-ends for split of %s. SpanA: %s SpanB: %s Overlap: %d",
                self.id, span,
                nspan, overl)

            if overl > 0:
                err_message = "Invalid overlap for {0}! T1: {1}. T2: {2}".format(
                    self.id, span, nspan)
                self.logger.error(err_message)
                raise InvalidTranscript(err_message)

    def __create_splitted_transcripts(self, cds_boundaries):

        """
        Private method called by split_by_cds to create the various (N>1) transcripts
        that are its output.
        :param cds_boundaries: a list of int tuples, containing the boundaries
         of the new transcripts.
        :return:
        """

        spans = []
        new_transcripts = []

        for counter, (boundary, bed12_objects) in enumerate(
                sorted(cds_boundaries.items(),
                       key=operator.itemgetter(0))):
            new_transcript = self.__class__()
            new_transcript.feature = "mRNA"
            for attribute in ["chrom", "source", "score", "strand", "attributes"]:
                setattr(new_transcript, attribute, getattr(self, attribute))
            # Determine which ORFs I have on my right and left
            new_transcript.parent = self.parent
            left = True
            right = True
            if counter == 0:  # leftmost
                left = False
            if 1 + counter == len(cds_boundaries):  # rightmost
                right = False
            counter += 1  # Otherwise they start from 0
            new_transcript.id = "{0}.split{1}".format(self.id, counter)
            new_transcript.logger = self.logger

            my_exons, discarded_exons, tstart, tend = self.__create_splitted_exons(
                boundary, left, right)

            self.logger.debug("""TID %s counter %d, boundary %s, left %s right %s""",
                              self.id,
                              counter,
                              boundary,
                              left,
                              right)

            if right is True:
                self.logger.debug("TID %s TEND %d Boun[1] %s",
                                  self.id, tend, boundary[1])
            if left is True:
                self.logger.debug("TID %s TSTART %d Boun[0] %s",
                                  self.id, tstart, boundary[0])

            assert len(my_exons) > 0, (discarded_exons, boundary)

            new_transcript.exons = my_exons

            new_transcript.start = min(exon[0] for exon in new_transcript.exons)
            new_transcript.end = max(exon[1] for exon in new_transcript.exons)
            new_transcript.json_conf = self.json_conf
            # Now we have to modify the BED12s to reflect
            # the fact that we are starting/ending earlier
            new_transcript.finalize()
            if new_transcript.monoexonic is True:
                new_transcript.strand = None

            new_bed12s = self.relocate_orfs(bed12_objects, tstart, tend)
            self.logger.debug("Loading %d ORFs into the new transcript",
                              len(new_bed12s))
            new_transcript.load_orfs(new_bed12s)

            if new_transcript.selected_cds_length <= 0:
                err_message = "No CDS information retained for {0} split {1}\n".format(
                    self.id, counter)
                err_message += "BED: {0}".format("\n\t".join([str(x) for x in new_bed12s]))
                raise InvalidTranscript(err_message)

            new_transcripts.append(new_transcript)
            nspan = (new_transcript.start, new_transcript.end)
            self.logger.debug(
                "Transcript {0} split {1}, discarded exons: {2}".format(
                    self.id, counter, discarded_exons))
            self.__check_collisions(nspan, spans)
            spans.append([new_transcript.start, new_transcript.end])

        return new_transcripts

    def split_by_cds(self):
        """This method is used for transcripts that have multiple ORFs.
        It will split them according to the CDS information into multiple transcripts.
        UTR information will be retained only if no ORF is down/upstream.
        The minimal overlap is defined inside the JSON at the key
            ["chimera_split"]["blast_params"]["minimal_hsp_overlap"]
        basically, we consider a HSP a hit only if the overlap is over a certain threshold
        and the HSP evalue under a certain threshold.

        The split by CDS can be executed in three different ways - CLEMENT, LENIENT, STRINGENT:

        - PERMISSIVE: split if two CDSs do not have hits in common,
        even when one or both do not have a hit at all.
        - STRINGENT: split only if two CDSs have hits and none
        of those is in common between them.
        - LENIENT: split if *both* lack hits, OR *both* have hits and none
        of those is in common.
        """
        self.finalize()

        # List of the transcript that will be retained

        if self.number_internal_orfs < 2:
            new_transcripts = [self]  # If we only have one ORF this is easy
        else:

            cds_boundaries = OrderedDict()
            for orf in sorted(self.loaded_bed12,
                              key=operator.attrgetter("thick_start", "thick_end")):
                cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

            # Check whether we have to split or not based on BLAST data
            if self.json_conf is not None:
                if self.json_conf["chimera_split"]["blast_check"] is True:
                    cds_boundaries = self.check_split_by_blast(cds_boundaries)

            if len(cds_boundaries) == 1:
                # Recheck how many boundaries we have - after the BLAST check
                # we might have determined that the transcript has not to be split
                new_transcripts = [self]
            else:
                new_transcripts = self.__create_splitted_transcripts(cds_boundaries)

        assert len(new_transcripts) > 0, str(self)
        for new_transc in new_transcripts:
            yield new_transc

        return

    def remove_utrs(self):
        """Method to strip a transcript from its UTRs.
        It will not execute anything if the transcript lacks a CDS or
        it has more than one ORF defined.
        """

        self.finalize()
        if self.selected_cds_length == 0:
            return
        elif self.three_utr_length + self.five_utr_length == 0:
            return  # No UTR to strip

        elif self.number_internal_orfs > 1:
            return
        elif re.search(r"\.orf[0-9]+$", self.id):
            return

        self.finalized = False
        exons = []
        cds_start, cds_end = self.combined_cds[0][0], self.combined_cds[-1][1]
        assert isinstance(cds_start, int)
        assert isinstance(cds_end, int)
        if len(self.selected_cds) == 1:
            self.exons = self.selected_cds
        else:
            for exon in self.exons:
                if exon in self.combined_utr:
                    continue
                elif exon in self.selected_cds:
                    exons.append(exon)
                elif exon[0] <= cds_start <= exon[1]:
                    exons.append((cds_start, exon[1]))
                elif exon[0] <= cds_end <= exon[1]:
                    exons.append((exon[0], cds_end))
                else:
                    raise InvalidTranscript(
                        "Exon: {0}; cds_start: {1}; cds_end: {2}; ID: {3}".format(
                            exon, self.selected_cds_start,
                            self.selected_cds_end, self.id))

            assert (len(exons) < len(self.exons) or
                    exons[0][0] > self.exons[0][0] or
                    exons[-1][1] < self.exons[-1][1]),\
                (exons, self.exons)
            self.exons = exons
        self.start = cds_start
        self.end = cds_end
        self.combined_utr = []
        self.finalize()

    def strip_cds(self):
        """Method to completely remove CDS information from a transcript.
        Necessary for those cases where the input is malformed."""

        self.logger.warning("Stripping CDS from {0}".format(self.id))
        self.finalized = False
        self.combined_cds = []
        self.combined_utr = []
        self.finalize()

    def __check_cdna_vs_utr(self):

        """
        Verify that cDNA + UTR in the transcript add up.
        :return:
        """

        if self.cdna_length > self.combined_utr_length + self.combined_cds_length:
            if self.combined_utr == [] and self.combined_cds != []:
                self.combined_cds = sorted(self.combined_cds,
                                           key=operator.itemgetter(0, 1))
                for exon in self.exons:
                    if exon in self.combined_cds:
                        continue
                    elif (exon[1] < self.combined_cds[0][0] or
                          exon[0] > self.combined_cds[-1][1]):
                        self.combined_utr.append(exon)
                    elif (exon[0] < self.combined_cds[0][0] and
                          exon[1] == self.combined_cds[0][1]):
                        self.combined_utr.append((exon[0], self.combined_cds[0][0] - 1))
                    elif (exon[1] > self.combined_cds[-1][1] and
                          exon[0] == self.combined_cds[-1][0]):
                        self.combined_utr.append((self.combined_cds[-1][1] + 1, exon[1]))
                    else:
                        if len(self.combined_cds) == 1:
                            self.combined_utr.append(
                                (exon[0], self.combined_cds[0][0] - 1))
                            self.combined_utr.append(
                                (self.combined_cds[-1][1] + 1, exon[1]))
                        else:
                            raise mikado_lib.exceptions.InvalidCDS(
                                "Error while inferring the UTR",
                                exon, self.id,
                                self.exons, self.combined_cds)

                equality_one = (self.combined_cds_length == self.combined_utr_length == 0)
                equality_two = (self.cdna_length ==
                                self.combined_utr_length + self.combined_cds_length)
                if not (equality_one or equality_two):
                    raise mikado_lib.exceptions.InvalidCDS(
                        "Failed to create the UTR",
                        self.id, self.exons,
                        self.combined_cds, self.combined_utr)
            else:
                pass

    def __basic_final_checks(self):

        """
        Function that verifies minimal criteria of a transcript before finalising.
        :return:
        """

        if len(self.exons) == 0:
            raise mikado_lib.exceptions.InvalidTranscript(
                "No exon defined for the transcript {0}. Aborting".format(self.id))

        if len(self.exons) > 1 and self.strand is None:
            raise mikado_lib.exceptions.InvalidTranscript(
                "Multiexonic transcripts must have a defined strand! Error for {0}".format(
                    self.id))

        if self.combined_utr != [] and self.combined_cds == []:
            raise mikado_lib.exceptions.InvalidTranscript(
                "Transcript {tid} has defined UTRs but no CDS feature!".format(
                    tid=self.id))

    def __verify_boundaries(self):

        """
        Method to verify that the start/end of the transcripts are exactly where they should.
        Called from finalise.
        :return:
        """

        try:
            if self.exons[0][0] != self.start or self.exons[-1][1] != self.end:
                if self.exons[0][0] > self.start and self.selected_cds[0][0] == self.start:
                    self.exons[0] = (self.start, self.exons[0][0])
                if self.exons[-1][1] < self.end and self.selected_cds[-1][1] == self.end:
                    self.exons[-1] = (self.exons[-1][0], self.end)

                if self.exons[0][0] != self.start or self.exons[-1][1] != self.end:
                    raise mikado_lib.exceptions.InvalidTranscript(
                        """The transcript {id} has coordinates {tstart}:{tend},
                    but its first and last exons define it up until {estart}:{eend}!
                    Exons: {exons}
                    """.format(id=self.id,
                               tstart=self.start,
                               tend=self.end,
                               estart=self.exons[0][0],
                               eend=self.exons[-1][1],
                               exons=self.exons))
        except IndexError as err:
            raise mikado_lib.exceptions.InvalidTranscript(
                err, self.id, str(self.exons))

    def finalize(self):
        """Function to calculate the internal introns from the exons.
        In the first step, it will sort the exons by their internal coordinates.
        """

        if self.finalized is True:
            return

        self.__basic_final_checks()
        # Sort the exons by start then stop
        self.exons = sorted(self.exons, key=operator.itemgetter(0, 1))

        self.__check_cdna_vs_utr()

        self.__calculate_introns()

        self.combined_cds = sorted(self.combined_cds,
                                   key=operator.itemgetter(0, 1))

        self.combined_utr = sorted(self.combined_utr,
                                   key=operator.itemgetter(0, 1))
        self.__check_completeness()

        # assert self.selected_internal_orf_index > -1
        self.segments = [("exon", e[0], e[1]) for e in self.exons] + \
                        [("CDS", c[0], c[1]) for c in self.combined_cds] + \
                        [("UTR", u[0], u[1]) for u in self.combined_utr]
        self.segments = sorted(self.segments, key=operator.itemgetter(1, 2, 0))

        self.internal_orfs = [self.segments]
        if self.combined_cds_length > 0:
            self.selected_internal_orf_index = 0
            if len(self.__phases) > 0:
                self._first_phase = sorted(self.__phases, key=operator.itemgetter(0),
                                           reverse=(self.strand == "-"))[0][1]
            else:
                self._first_phase = 0

        # Necessary to set it to the default value
        _ = self.selected_internal_orf

        if len(self.combined_cds) > 0:
            self.feature = "mRNA"

        self.__verify_boundaries()

        if len(self.combined_cds) == 0:
            self.selected_internal_orf_cds = tuple([])
        else:
            self.selected_internal_orf_cds = tuple(
                filter(lambda x: x[0] == "CDS",
                       self.internal_orfs[self.selected_internal_orf_index])
            )

        self.finalized = True
        return

    def __calculate_introns(self):

        """Private method to create the stores of intron
        and splice sites positions.
        """

        introns = []
        splices = []

        if len(self.exons) > 1:
            for index in range(len(self.exons) - 1):
                exona, exonb = self.exons[index:index + 2]
                if exona[1] >= exonb[0]:
                    raise mikado_lib.exceptions.InvalidTranscript(
                        "Overlapping exons found!\n{0} {1}/{2}\n{3}".format(
                            self.id, exona, exonb, self.exons))
                # Append the splice junction
                introns.append((exona[1] + 1, exonb[0] - 1))
                # Append the splice locations
                splices.extend([exona[1] + 1, exonb[0] - 1])
        self.introns = set(introns)
        self.splices = set(splices)

    def __check_completeness(self):

        """Private method that checks whether a transcript is complete
        or not based solely on the presence of CDS/UTR information."""

        if len(self.combined_utr) > 0:
            if self.combined_utr[0][0] < self.combined_cds[0][0]:
                if self.strand == "+":
                    self.has_start_codon = True
                elif self.strand == "-":
                    self.has_stop_codon = True
            if self.combined_utr[-1][1] > self.combined_cds[-1][1]:
                if self.strand == "+":
                    self.has_stop_codon = True
                elif self.strand == "-":
                    self.has_start_codon = True

    def reverse_strand(self):
        """Method to reverse the strand"""
        if self.strand == "+":
            self.strand = "-"
        elif self.strand == "-":
            self.strand = "+"
        elif self.strand is None:
            pass
        return

    def __connect_to_db(self):

        """This method will connect to the database using the information
        contained in the JSON configuration."""

        self.engine = mikado_lib.serializers.dbutils.connect(
            self.json_conf, self.logger)

        self.sessionmaker = sessionmaker()
        self.sessionmaker.configure(bind=self.engine)
        self.session = self.sessionmaker()

    def load_information_from_db(self, json_conf, introns=None, session=None,
                                 data_dict=None):
        """This method will invoke the check for:

        :param json_conf: Necessary configuration file
        :type json_conf: dict

        :param introns: the verified introns in the Locus
        :type introns: None,set

        :param session: an SQLAlchemy session
        :type session: sqlalchemy.orm.session

        :param data_dict: a dictionary containing the information directly
        :type data_dict: dict

        Verified introns can be provided from outside using the keyword.
        Otherwise, they will be extracted from the database directly.
        """

        self.logger.debug("Loading {0}".format(self.id))
        self.load_json(json_conf)

        if data_dict is not None:
            self.retrieve_from_dict(data_dict)
            # yield from self.retrieve_from_dict(data_dict)
        else:
            if session is None:
                self.__connect_to_db()
            else:
                self.session = session
            self.__load_verified_introns(introns)
            # yield from self.load_verified_introns(introns)
            self.query_id = self.query_baked(self.session).params(
                query_name=self.id).all()
            if len(self.query_id) == 0:
                self.logger.warning(
                    "Transcript not in database: %s", self.id)
            else:
                self.query_id = self.query_id[0].query_id
                self.load_orfs(list(self.retrieve_orfs()))
                self.__load_blast()
            self.logger.debug("Loaded %s", self.id)

    def retrieve_from_dict(self, data_dict):
        """
        Method to retrieve transcript data directly from a dictionary.
        :param data_dict: the dictionary with loaded data from DB
        """

        self.logger.debug(
            "Retrieving information from DB dictionary for %s",
            self.id)
        # Intron data
        for intron in self.introns:
            if (self.chrom, intron[0], intron[1], self.strand) in data_dict["junctions"]:
                self.verified_introns.add(intron)

        # ORF data
        trust_strand = self.json_conf["orf_loading"]["strand_specific"]
        min_cds_len = self.json_conf["orf_loading"]["minimal_orf_length"]

        self.logger.debug("Retrieving ORF information from DB dictionary for %s",
                          self.id)
        if self.id in data_dict["orfs"]:
            candidate_orfs = list(filter(
                lambda orf: orf.cds_len >= min_cds_len,
                data_dict["orfs"][self.id]))
        else:
            candidate_orfs = []

        # They must already be as ORFs
        if (self.monoexonic is False) or (self.monoexonic is True and trust_strand is True):
            # Remove negative strand ORFs for multiexonic transcripts,
            # or monoexonic strand-specific transcripts
            candidate_orfs = list(filter(lambda orf: orf.strand != "-",
                                         candidate_orfs))

        self.load_orfs(candidate_orfs)

        if self.json_conf["chimera_split"]["blast_check"] is True:
            self.logger.debug("Retrieving BLAST hits for %s",
                              self.id)
            maximum_evalue = self.json_conf["chimera_split"]["blast_params"]["evalue"]

            if self.id in data_dict["hits"]:
                # this is a dictionary full of lists of dictionary
                hits = data_dict["hits"][self.id]
            else:
                hits = list()

            self.logger.debug("Found %d potential BLAST hits for %s with evalue <= %f",
                              len(hits),
                              self.id,
                              maximum_evalue)

            self.blast_hits.extend(hits)
            self.logger.debug("Loaded %d BLAST hits for %s",
                              len(self.blast_hits), self.id)

        self.logger.debug("Retrieved information from DB dictionary for %s",
                          self.id)

    def load_json(self, json_conf):
        """
        Setter for the json configuration dictionary.
        :param json_conf: The configuration dictionary
        :type json_conf: dict
        """
        self.json_conf = json_conf

    def __load_verified_introns(self, introns=None):

        """This method will load verified junctions from the external
        (usually the superlocus class).

        :param introns: verified introns
        :type introns: set,None
        """

        if introns is None:
            for intron in self.introns:
                # Disable checks as the hybridproperties confuse
                # both pycharm and pylint
                # noinspection PyCallByClass,PyTypeChecker
                # pylint: disable=no-value-for-parameter
                if self.session.query(Junction).filter(
                        Junction.is_equal(self.chrom, intron[0],
                                          intron[1], self.strand)).count() == 1:
                    self.verified_introns.add(intron)

        else:
            for intron in introns:
                if intron in self.introns:
                    self.verified_introns.add(intron)
        return

    def retrieve_orfs(self):

        """This method will look up the ORFs loaded inside the database.
        During the selection, the function will also remove overlapping ORFs.
        """

        if self.query_id is None:
            return []

        trust_strand = self.json_conf["orf_loading"]["strand_specific"]
        min_cds_len = self.json_conf["orf_loading"]["minimal_orf_length"]

        orf_results = self.orf_baked(self.session).params(query_id=self.query_id,
                                                          cds_len=min_cds_len)

        if (self.monoexonic is False) or (self.monoexonic is True and trust_strand is True):
            # Remove negative strand ORFs for multiexonic transcripts,
            # or monoexonic strand-specific transcripts
            candidate_orfs = list(filter(lambda orf: orf.strand != "-", orf_results))
        else:
            candidate_orfs = orf_results.all()

        if len(candidate_orfs) == 0:
            return []
        else:
            return [orf.as_bed12() for orf in candidate_orfs]

    def __create_internal_orf(self, orf):

        """
        Private method that calculates the assignment of the exons given the
        coordinates of the transcriptomic ORF.
        """

        cds_exons = []
        current_start, current_end = 0, 0

        for exon in sorted(self.exons, key=operator.itemgetter(0, 1),
                           reverse=(self.strand == "-")):
            cds_exons.append(("exon", exon[0], exon[1]))
            current_start += 1
            current_end += exon[1] - exon[0] + 1
            # Whole UTR
            if current_end < orf.thick_start or current_start > orf.thick_end:
                cds_exons.append(("UTR", exon[0], exon[1]))
            else:
                if self.strand == "+":
                    c_start = exon[0] + max(0, orf.thick_start - current_start)
                    c_end = exon[1] - max(0, current_end - orf.thick_end)
                else:
                    c_start = exon[0] + max(0, current_end - orf.thick_end)
                    c_end = exon[1] - max(0, orf.thick_start - current_start)
                if c_start > exon[0]:
                    u_end = c_start - 1
                    cds_exons.append(("UTR", exon[0], u_end))
                if c_start <= c_end:
                    cds_exons.append(("CDS", c_start, c_end))
                if c_end < exon[1]:
                    cds_exons.append(("UTR", c_end + 1, exon[1]))
            current_start = current_end

        # if self.strand == "+":
        #     for exon in sorted(self.exons, key=operator.itemgetter(0, 1)):
        #         cds_exons.append(("exon", exon[0], exon[1]))
        #         current_start += 1
        #         current_end += exon[1] - exon[0] + 1
        #         # Whole UTR
        #         if current_end < orf.thick_start or current_start > orf.thick_end:
        #             cds_exons.append(("UTR", exon[0], exon[1]))
        #         else:
        #             c_start = exon[0] + max(0, orf.thick_start - current_start)
        #             c_end = exon[1] - max(0, current_end - orf.thick_end)
        #             if c_start > exon[0]:
        #                 u_end = c_start - 1
        #                 cds_exons.append(("UTR", exon[0], u_end))
        #             if c_start <= c_end:
        #                 cds_exons.append(("CDS", c_start, c_end))
        #             if c_end < exon[1]:
        #                 cds_exons.append(("UTR", c_end + 1, exon[1]))
        #         current_start = current_end
        #
        # elif self.strand == "-":
        #     for exon in sorted(self.exons, key=operator.itemgetter(0, 1), reverse=True):
        #         cds_exons.append(("exon", exon[0], exon[1]))
        #         current_start += 1
        #         current_end += exon[1] - exon[0] + 1
        #         if current_end < orf.thick_start or current_start > orf.thick_end:
        #             cds_exons.append(("UTR", exon[0], exon[1]))
        #         else:
        #             c_start = exon[0] + max(0, current_end - orf.thick_end)
        #             c_end = exon[1] - max(0, orf.thick_start - current_start)
        #             if c_start > exon[0]:
        #                 cds_exons.append(("UTR", exon[0], c_start - 1))
        #             if c_start <= c_end:
        #                 cds_exons.append(("CDS", c_start, c_end))
        #             if c_end < exon[1]:
        #                 cds_exons.append(("UTR", c_end + 1, exon[1]))
        #         current_start = current_end
        return cds_exons

    def load_orfs(self, candidate_orfs):

        """
        :param candidate_orfs: The ORFs to be inspected for loading.
        :type candidate_orfs: list[mikado_lib.parsers.serializers.orf.Orf

        This method replicates what is done internally by the
        "cdna_alignment_orf_to_genome_orf.pl"
        utility in the TransDecoder suite. It takes as argument "candidate_orfs"
        i.e. a list of BED12 serialised objects.
        The method expects as argument a dictionary containing BED entries,
        and indexed by the transcript name. The indexed name *must* equal the
        "id" property, otherwise the method returns immediately.
        If no entry is found for the transcript, the method exits immediately.
        Otherwise, any CDS information present in the original GFF/GTF
        file is completely overwritten.
        Briefly, it follows this logic:
        - Finalise the transcript
        - Retrieve from the dictionary (input) the CDS object
        - Sort in decreasing order the CDSs on the basis of:
            - Presence of start/stop codon
            - CDS length (useful for monoexonic transcripts where we might want to set the strand)
        - For each CDS:
            - If the ORF is on the + strand:
                - all good
            - If the ORF is on the - strand:
                - if the transcript is monoexonic: invert its genomic strand
                - if the transcript is multiexonic: skip
            - Start looking at the exons
        """

        # Prepare the transcript
        self.finalize()

        candidate_orfs = self.find_overlapping_cds(candidate_orfs)

        if candidate_orfs is None or len(candidate_orfs) == 0:
            self.logger.debug("No ORF for {0}".format(self.id))
            return

        self.combined_utr = []
        self.combined_cds = []
        self.internal_orfs = []
        self.finalized = False
        # Token to be set to False after the first CDS is exhausted
        primary_orf = True
        primary_strand = None
        # This will keep in memory the original BED12 objects
        self.loaded_bed12 = []

        for orf in candidate_orfs:
            # Minimal check
            if primary_orf is True:
                (self.has_start_codon, self.has_stop_codon) = (orf.has_start_codon,
                                                               orf.has_stop_codon)
                primary_orf = False
                primary_strand = orf.strand
            elif primary_orf is False and orf.strand != primary_strand:
                continue

            check_sanity = (orf.thick_start >= 1 and orf.thick_end <= self.cdna_length)
            if len(orf) != self.cdna_length or not check_sanity:
                message = "Wrong ORF for {0}: ".format(orf.id)
                message += "cDNA length: {0}; ".format(self.cdna_length)
                message += "orf length: {0}; ".format(len(orf))
                message += "CDS: {0}-{1}".format(orf.thick_start, orf.thick_end)
                self.logger.warning(message)
                continue

            if self.strand is None:
                self.strand = orf.strand

            self.loaded_bed12.append(orf)
            cds_exons = self.__create_internal_orf(orf)
            self.internal_orfs.append(sorted(
                cds_exons, key=operator.itemgetter(1, 2)))

        # Now verify the loaded content
        self.check_loaded_orfs()

    def check_loaded_orfs(self):

        """
        This function verifies the ORF status after
        loading from an external data structure/database.
        :return:
        """

        if len(self.internal_orf_lengths) == 0:
            self.logger.warning("No candidate ORF retained for %s",
                                self.id)

        if len(self.internal_orfs) == 1:
            self.logger.debug("Found 1 ORF for %s", self.id)
            self.combined_cds = sorted(
                [(a[1], a[2]) for a in filter(lambda x: x[0] == "CDS",
                                              self.internal_orfs[0])],
                key=operator.itemgetter(0, 1)

            )
            self.combined_utr = sorted(
                [(a[1], a[2]) for a in filter(lambda x: x[0] == "UTR",
                                              self.internal_orfs[0])],
                key=operator.itemgetter(0, 1)

            )

        elif len(self.internal_orfs) > 1:
            self.logger.debug("Found %d ORFs for %s",
                              len(self.internal_orfs),
                              self.id)
            cds_spans = []
            candidates = []
            for internal_cds in self.internal_orfs:
                candidates.extend(
                    [tuple([a[1], a[2]]) for a in filter(
                        lambda tup: tup[0] == "CDS", internal_cds)])

            for comm in self.find_communities(candidates):
                span = tuple([min(t[0] for t in comm), max(t[1] for t in comm)])
                cds_spans.append(span)

            self.combined_cds = sorted(cds_spans, key=operator.itemgetter(0, 1))

            # This method is probably OBSCENELY inefficient,
            # but I cannot think of a better one for the moment.
            curr_utr_segment = None

            utr_pos = set.difference(
                set.union(*[set(range(exon[0], exon[1] + 1)) for exon in self.exons]),
                set.union(*[set(range(cds[0], cds[1] + 1)) for cds in self.combined_cds])
            )
            for pos in sorted(list(utr_pos)):
                if curr_utr_segment is None:
                    curr_utr_segment = (pos, pos)
                else:
                    if pos == curr_utr_segment[1] + 1:
                        curr_utr_segment = (curr_utr_segment[0], pos)
                    else:
                        self.combined_utr.append(curr_utr_segment)
                        curr_utr_segment = (pos, pos)

            if curr_utr_segment is not None:
                self.combined_utr.append(curr_utr_segment)

            equality = (self.cdna_length ==
                        self.combined_cds_length + self.combined_utr_length)
            assert equality, (self.cdna_length, self.combined_cds, self.combined_utr)

        if not self.internal_orfs:
            self.finalize()
        else:
            self.feature = "mRNA"
            self.finalized = True
        return

    def __load_blast(self):

        """This method looks into the DB for hits corresponding to the desired requirements.
        Hits will be loaded into the "blast_hits" list;
        we will not store the SQLAlchemy query object,
        but rather its representation as a dictionary
        (using the Hit.as_dict() method).
        """

        if self.query_id is None:
            return

        if self.json_conf["chimera_split"]["blast_check"] is False:
            return

        max_target_seqs = self.json_conf["chimera_split"]["blast_params"]["max_target_seqs"]
        maximum_evalue = self.json_conf["chimera_split"]["blast_params"]["evalue"]

        blast_hits_query = self.blast_baked(self.session).params(
            query_id=self.query_id,
            evalue=maximum_evalue,
            max_target_seqs=max_target_seqs)
        counter = 0
        self.logger.debug("Starting to load BLAST data for %s",
                          self.id)
        for hit in blast_hits_query:
            counter += 1
            self.blast_hits.append(hit.as_dict())
        self.logger.debug("Loaded %d BLAST hits for %s",
                          counter, self.id)

    @property
    def logger(self):
        """
        Property. It returns the logger instance attached to the class.
        :rtype : logging.Logger | None
        """

        return self.__logger

    @logger.setter
    def logger(self, logger):
        """Set a logger for the instance.
        :param logger: a Logger instance
        :type logger: logging.Logger | None
        """
        if logger is None:
            if self.__logger is None:
                logger = create_null_logger(self)
                self.__logger = logger
            else:
                pass
        else:
            assert isinstance(logger, logging.Logger)
            self.__logger = logger

    @logger.deleter
    def logger(self):
        """
        Destroyer for the logger. It sets the internal __logger attribute to None.
        """
        self.__logger = None

    # ###################Class methods#####################################

    @classmethod
    def is_overlapping_cds(cls, first, second):
        """
        :param first: first ORF to check for overlap
        :param second: second ORF to check for overlap
        :rtype bool
        """
        if first == second or cls.overlap(
                (first.thick_start, first.thick_end),
                (second.thick_start, second.thick_end)) < 0:
            return False
        return True

    @classmethod
    def is_intersecting(cls, first, second):
        """
        :param first: first exon to check
        :type first: tuple([int, int])

        :param second: second exon to check
        :type second: tuple([int, int])

        :rtype bool

        Implementation of the is_intersecting method.
        It checks overlaps between exons.
        """

        if first == second or cls.overlap(first, second) < 0:
            return False
        return True

    @classmethod
    def overlap(cls, first, second):
        """
        :param first: first exon to check
        :type first: tuple([int, int])

        :param second: second exon to check
        :type second: tuple([int, int])
        :rtype: int

        This method checks the overlap between two int duplexes.
        """

        lend = max(first[0], second[0])
        rend = min(first[1], second[1])
        return rend - lend

    @classmethod
    def find_communities(cls, objects: list) -> list:
        """

        :param objects: a list of objects to analyse
        :type objects: list,set

        Wrapper for the Abstractlocus method.
        As we are interested only in the communities, not the cliques,
        this wrapper discards the cliques
        (first element of the Abstractlocus.find_communities results)
        """
        data = dict((obj, obj) for obj in objects)
        communities = Abstractlocus.find_communities(
            Abstractlocus.define_graph(data,
                                       inters=cls.is_intersecting))[1]

        return communities

    @classmethod
    def get_available_metrics(cls) -> list:
        """This function retrieves all metrics available for the class."""
        metrics = list(x[0] for x in filter(
            lambda y: "__" not in y[0] and isinstance(cls.__dict__[y[0]], Metric),
            inspect.getmembers(cls)))
        assert "tid" in metrics and "parent" in metrics and "score" in metrics
        final_metrics = ["tid", "parent", "score"] + sorted(
            list(filter(lambda x: x not in ["tid", "parent", "score"], metrics)))
        return final_metrics

    # ###################Class properties##################################

    # This will be id, no changes.
    # pylint: disable=invalid-name
    @property
    def id(self):
        """ID of the transcript - cannot be an undefined value."""
        return self.__id

    @id.setter
    def id(self, newid):
        """
        :param newid: a string which will become the ID of the instance.
        :type newid: str
        """

        if not isinstance(newid, str):
            raise ValueError("Invalid value for id: {0}, type {1}".format(
                newid, type(newid)))
        self.__id = sys.intern(newid)
    # pylint: enable=invalid-name

    @property
    def available_metrics(self) -> list:
        """Return the list of available metrics, using the "get_metrics" function."""
        return self.get_available_metrics()

    @property
    def strand(self):
        """
        Strand of the transcript. One of None, "-", "+"

        :rtype str | None
        """
        return self.__strand

    @strand.setter
    def strand(self, strand):
        """

        :param strand
        :type strand: None | str

        Setter for the strand of the transcript. It must be one of None, "-", "+"
        """
        if strand in ("+", "-"):
            self.__strand = strand
        elif strand in (None, ".", "?"):
            self.__strand = None
        else:
            raise ValueError("Invalid value for strand: {0}".format(strand))

    @property
    def selected_internal_orf(self):
        """This property will return the tuple of tuples of the ORF selected as "best".
        To avoid memory wasting, the tuple is accessed in real-time using
        a token (__max_internal_orf_index) which holds the position in the
        __internal_cds list of the longest CDS.
        """

        # Non-sense to calculate the maximum CDS for transcripts without it
        if len(self.combined_cds) == 0:
            self.__max_internal_orf_length = 0
            self.selected_internal_orf_index = 0
            return tuple([])
        else:
            return self.internal_orfs[self.selected_internal_orf_index]

    @property
    def selected_internal_orf_cds(self):
        """This property will return the tuple of tuples of the CDS segments of
        the selected ORF inside the transcript. To avoid memory wasting,
        the tuple is accessed in real-time using a token
        (__max_internal_orf_index) which holds the position
        in the __internal_cds list of the longest CDS.
        """

        # Non-sense to calculate the maximum CDS for transcripts without it
        return self._selected_internal_orf_cds

    @selected_internal_orf_cds.setter
    def selected_internal_orf_cds(self, internal_orf):
        """
        Setter for selected_internal_orf_cds
        :param internal_orf:
        :return:
        """

        if not isinstance(internal_orf, tuple):
            raise TypeError("Invalid internal ORF type ({0}): {1}".format(
                type(internal_orf),
                internal_orf
            ))

        self._selected_internal_orf_cds = internal_orf

    @property
    def five_utr(self):
        """Returns the exons in the 5' UTR of the selected ORF.
        If the start codon is absent, no UTR is given."""
        if len(self.combined_cds) == 0:
            return []
        if self.strand == "-":
            return list(
                filter(lambda exon: exon[0] == "UTR" and exon[1] > self.selected_cds_start,
                       self.selected_internal_orf))
        else:
            return list(
                filter(lambda exon: exon[0] == "UTR" and exon[2] < self.selected_cds_start,
                       self.selected_internal_orf))

    @property
    def three_utr(self):
        """Returns the exons in the 3' UTR of the selected ORF.
        If the end codon is absent, no UTR is given."""
        if len(self.combined_cds) == 0:
            return []
        if self.strand == "-":
            return list(
                filter(lambda exon: exon[0] == "UTR" and exon[2] < self.selected_cds_end,
                       self.selected_internal_orf))
        else:
            return list(
                filter(lambda exon: exon[0] == "UTR" and exon[1] > self.selected_cds_end,
                       self.selected_internal_orf))

    @property
    def selected_internal_orf_index(self):
        """Token which memorizes the position in the ORF list of the selected ORF.
        :rtype : None | int
        """
        return self.__max_internal_orf_index

    @selected_internal_orf_index.setter
    def selected_internal_orf_index(self, index):
        """Setter for selected_internal_orf_index.
        :param index:
        :type index: None,int
        """
        if index is None:
            self.__max_internal_orf_index = index
            return
        if not isinstance(index, int):
            raise TypeError()
        if index < 0 or index >= len(self.internal_orfs):
            raise IndexError("No ORF corresponding to this index: {0}".format(index))
        self.__max_internal_orf_index = index

    @property
    def internal_orf_lengths(self):
        """This property returns a list of the lengths of the internal ORFs.
        :rtype : list[int]
        """
        lengths = []
        for internal_cds in self.internal_orfs:
            lengths.append(sum(x[2] - x[1] + 1 for x in filter(
                lambda c: c[0] == "CDS", internal_cds)))
        lengths = sorted(lengths, reverse=True)
        return lengths

    @property
    def non_overlapping_cds(self):
        """This property returns a set containing the set union of all CDS segments
        inside the internal CDSs. In the case of a transcript with no CDS, this is empty.
        In the case where there is only one CDS, this returns the combined_cds holder.
        In the case instead where there are multiple CDSs, the property will calculate
        the set union of all CDS segments.
        """
        if self.__non_overlapping_cds is None:
            self.finalize()
            self.__non_overlapping_cds = set()
            for internal_cds in self.internal_orfs:
                segments = set([(x[1], x[2]) for x in filter(
                    lambda segment: segment[0] == "CDS", internal_cds)])
                self.__non_overlapping_cds.update(segments)
        return self.__non_overlapping_cds

    @non_overlapping_cds.setter
    def non_overlapping_cds(self, arg):
        """
        :param arg: the unioin of all non-overlapping CDS segments.
        :type arg: set
        Setter for the non_overlapping_cds property."""
        self.__non_overlapping_cds = arg

    @property
    def exons(self):
        """This property stores the exons of the transcript as (start,end) tuples.

        :rtype : list
        """
        return self.__exons

    @exons.setter
    def exons(self, *args):
        """
        :param args: a list/set of exons
        :type args: set | list

        """

        if not isinstance(args[0], (set, list)):
            raise TypeError(type(args[0]))
        self.__exons = list(args[0])

    @property
    def combined_cds_introns(self):
        """This property returns the introns which are located between CDS
        segments in the combined CDS."""
        if self.number_internal_orfs < 2:
            return self.selected_cds_introns
        if self.number_internal_orfs == 0 or len(self.combined_cds) < 2:
            return set()

        cintrons = []
        for position in range(len(self.combined_cds) - 1):
            former = self.combined_cds[position]
            latter = self.combined_cds[position + 1]
            junc = (former[1] + 1, latter[0] - 1)
            if junc in self.introns:
                cintrons.append(junc)
        cintrons = set(cintrons)
        return cintrons

    @property
    def selected_cds_introns(self):
        """This property returns the introns which are located between
        CDS segments in the selected ORF."""

        if len(self.selected_cds) < 2:
            return set()
        if self.number_internal_orfs == 0 or len(self.combined_cds) < 2:
            return set()

        cintrons = []
        for position in range(len(self.selected_internal_orf_cds) - 1):
            cintrons.append(
                (self.selected_internal_orf_cds[position][2] + 1,
                 self.selected_internal_orf_cds[position + 1][1] - 1)
            )
        cintrons = set(cintrons)
        return cintrons

    @property
    def combined_cds_start(self):
        """This property returns the location of the start of the combined
        CDS for the transcript. If no CDS is defined, it defaults
        to the transcript start."""

        if len(self.combined_cds) == 0:
            if self.strand == "+":
                return self.start
            else:
                return self.end
        if self.strand == "+":
            return self.combined_cds[0][0]
        else:
            return self.combined_cds[-1][1]

    @property
    def combined_cds(self):
        """This is a list which contains all the non-overlapping CDS
        segments inside the cDNA. The list comprises the segments
        as duples (start,end)."""
        return self.__combined_cds

    @combined_cds.setter
    def combined_cds(self, combined):
        """
        Setter for combined_cds. It performs some basic checks,
        e.g. that all the members of the list are integer duplexes.

        :param combined: list
        :type combined: list[(int,int)]
        """

        error = TypeError("Invalid value for combined CDS: {0}".format(combined))

        if not isinstance(combined, list):
            raise error
        elif any(self.__wrong_combined_entry(comb) for comb in combined):
            raise error

        self.__combined_cds = combined

    @staticmethod
    def __wrong_combined_entry(to_test):
        """
        Private method to test the correctness of entries for "combined"
        data
        :param to_test:
        :return:
        """
        if len(to_test) != 2:
            return True
        elif not isinstance(to_test[0], int):
            return True
        elif not isinstance(to_test[1], int):
            return True
        return False

    @property
    def combined_utr(self):
        """This is a list which contains all the non-overlapping UTR
        segments inside the cDNA.
        The list comprises the segments as duples (start,end)."""
        return self.__combined_utr

    @combined_utr.setter
    def combined_utr(self, combined):
        """Setter for combined UTR. It performs some basic checks,
        e.g. that all the members of the list
        are integer duplexes.

        :param combined: UTR list
        :type combined: list[(int,int)]

        """

        if not isinstance(combined, list):
            raise TypeError("Invalid value for combined UTR: {0}".format(combined))
        elif any(self.__wrong_combined_entry(comb) for comb in combined):
            raise TypeError("Invalid value for combined UTR: {0}".format(combined))

        self.__combined_utr = combined

    @property
    def combined_cds_end(self):
        """This property returns the location of the end of the combined CDS
        for the transcript. If no CDS is defined, it defaults
        to the transcript end."""
        if len(self.combined_cds) == 0:
            if self.strand == "+":
                return self.end
            else:
                return self.start
        if self.strand == "-":
            return self.combined_cds[0][0]
        else:
            return self.combined_cds[-1][1]

    @property
    def selected_cds(self):
        """This property return the CDS exons of the ORF selected as best
         inside the cDNA, in the form of duplices (start, end)"""
        if len(self.combined_cds) == 0:
            self.__selected_cds = []
        else:
            self.__selected_cds = list(
                (x[1], x[2]) for x in filter(lambda x: x[0] == "CDS",
                                             self.selected_internal_orf))
        return self.__selected_cds

    @property
    def selected_cds_start(self):
        """This property returns the location of the start
        of the best CDS for the transcript.
        If no CDS is defined, it defaults to the transcript start."""

        if len(self.combined_cds) == 0:
            return None

        if self.strand == "-":
            return self.selected_cds[-1][1]
        else:
            return self.selected_cds[0][0]

    @property
    def selected_cds_end(self):
        """This property returns the location of the end
        of the best CDS for the transcript.
        If no CDS is defined, it defaults to the transcript start."""

        if len(self.combined_cds) == 0:
            return None
        if self.strand == "-":
            return self.selected_cds[0][0]
        else:
            return self.selected_cds[-1][1]

    @property
    def monoexonic(self):
        """
        Property. True if the transcript has only one exon, False otherwise.
        :rtype bool
        """
        if len(self.exons) == 1:
            return True
        return False

    # ################### Class metrics ##################################

    # Disable normal checks on names and hidden methods, as
    # checkers get confused by the Metric method
    # pylint: disable=method-hidden,invalid-name
    @Metric
    def tid(self):
        """ID of the transcript - cannot be an undefined value. Alias of id.
        :rtype str
        """
        return self.id

    @tid.setter
    def tid(self, tid):
        """
        :param tid: ID of the transcript.
        :type tid: str
        """
        self.id = tid

    @Metric
    def parent(self):
        """Name of the parent feature of the transcript."""
        return self.__parent

    @parent.setter
    def parent(self, parent):
        """
        :param parent: the parent of the transcript.
        :type parent: list
        :type parent: str
        """
        if isinstance(parent, (list, type(None))):
            self.__parent = parent
        elif isinstance(parent, str):
            if "," in parent:
                self.__parent = parent.split(",")
            else:
                self.__parent = [parent]
        else:
            raise ValueError("Invalid value for parent: {0}, type {1}".format(
                parent, type(parent)))

    @Metric
    def score(self):
        """Numerical value which summarizes the reliability of the transcript."""
        return self.__score

    @score.setter
    def score(self, score):

        """Setter for the numerical value which summarizes the reliability
        of the transcript.
        :param score: the new score of the transcript
        :type score: None
        :type score: int
        :type score: float
        """

        if score is not None:
            if not isinstance(score, (float, int)):
                try:
                    score = float(score)
                except:
                    raise ValueError(
                        "Invalid value for score: {0}, type {1}".format(score, type(score)))
        self.__score = score

    @Metric
    def combined_cds_length(self):
        """This property return the length of the CDS part of the transcript."""
        return sum([c[1] - c[0] + 1 for c in self.combined_cds])

    @Metric
    def combined_cds_num(self):
        """This property returns the number of non-overlapping CDS segments
        in the transcript."""
        return len(self.combined_cds)

    @Metric
    def combined_cds_num_fraction(self):
        """This property returns the fraction of non-overlapping CDS segments
        in the transcript
        vs. the total number of exons"""
        return len(self.combined_cds) / len(self.exons)

    @Metric
    def combined_cds_fraction(self):
        """This property return the percentage of the CDS part of the transcript
        vs. the cDNA length"""
        return self.combined_cds_length / self.cdna_length

    @Metric
    def combined_utr_length(self):
        """This property return the length of the UTR part of the transcript."""
        return sum([e[1] - e[0] + 1 for e in self.combined_utr])

    @Metric
    def combined_utr_fraction(self):
        """This property returns the fraction of the cDNA which is not coding according
        to any ORF. Complement of combined_cds_fraction"""
        return 1 - self.combined_cds_fraction

    @Metric
    def cdna_length(self):
        """This property returns the length of the transcript."""
        return sum([e[1] - e[0] + 1 for e in self.exons])

    @Metric
    def number_internal_orfs(self):
        """This property returns the number of ORFs inside a transcript."""
        return len(self.internal_orfs)

    @Metric
    def selected_cds_length(self):
        """This property calculates the length of the CDS selected as best inside
        the cDNA."""
        if len(self.combined_cds) == 0:
            self.__max_internal_orf_length = 0
        else:
            self.__max_internal_orf_length = sum(
                x[2] - x[1] + 1 for x in filter(lambda x: x[0] == "CDS",
                                                self.selected_internal_orf))
        return self.__max_internal_orf_length

    @Metric
    def selected_cds_num(self):
        """This property calculates the number of CDS exons for the selected ORF"""
        return len(list(filter(lambda exon: exon[0] == "CDS",
                               self.selected_internal_orf)))

    @Metric
    def selected_cds_fraction(self):
        """This property calculates the fraction of the selected CDS vs. the cDNA length."""
        return self.selected_cds_length / self.cdna_length

    @Metric
    def highest_cds_exons_num(self):
        """Returns the number of CDS segments in the selected ORF
        (irrespective of the number of exons involved)"""
        return len(list(filter(lambda x: x[0] == "CDS", self.selected_internal_orf)))

    @Metric
    def selected_cds_exons_fraction(self):
        """Returns the fraction of CDS segments in the selected ORF
        (irrespective of the number of exons involved)"""
        return len(list(filter(lambda x: x[0] == "CDS",
                               self.selected_internal_orf))) / len(self.exons)

    @Metric
    def highest_cds_exon_number(self):
        """This property returns the maximum number of CDS segments
        among the ORFs; this number can refer to an ORF *DIFFERENT*
        from the maximal ORF."""
        cds_numbers = []
        for cds in self.internal_orfs:
            cds_numbers.append(len(list(filter(lambda x: x[0] == "CDS", cds))))
        return max(cds_numbers)

    @Metric
    def selected_cds_number_fraction(self):
        """This property returns the proportion of best possible CDS segments
        vs. the number of exons. See selected_cds_number."""
        return self.selected_cds_num / self.exon_num

    @Metric
    def cds_not_maximal(self):
        """This property returns the length of the CDS excluded from the selected ORF."""
        if len(self.internal_orfs) < 2:
            return 0
        return self.combined_cds_length - self.selected_cds_length

    @Metric
    def cds_not_maximal_fraction(self):
        """This property returns the fraction of bases not in the selected ORF compared to
        the total number of CDS bases in the cDNA."""
        if self.combined_cds_length == 0:
            return 0
        else:
            return self.cds_not_maximal / self.combined_cds_length

    @Metric
    def five_utr_length(self):
        """Returns the length of the 5' UTR of the selected ORF."""
        if len(self.combined_cds) == 0:
            return 0
        return sum(x[2] - x[1] + 1 for x in self.five_utr)

    @Metric
    def five_utr_num(self):
        """This property returns the number of 5' UTR segments for the selected ORF."""
        return len(self.five_utr)

    @Metric
    def five_utr_num_complete(self):
        """This property returns the number of 5' UTR segments for the selected ORF,
        considering only those which are complete exons."""
        return len(list(filter(lambda utr: (utr[1], utr[2]) in self.exons,
                               self.five_utr)))

    @Metric
    def three_utr_length(self):
        """Returns the length of the 5' UTR of the selected ORF."""
        if len(self.combined_cds) == 0:
            return 0
        return sum(x[2] - x[1] + 1 for x in self.three_utr)

    @Metric
    def three_utr_num(self):
        """This property returns the number of 3' UTR segments
        (referred to the selected ORF)."""
        return len(self.three_utr)

    @Metric
    def three_utr_num_complete(self):
        """This property returns the number of 3' UTR segments for the selected ORF,
        considering only those which are complete exons."""
        return len(list(filter(lambda utr: (utr[1], utr[2]) in self.exons,
                               self.three_utr)))

    @Metric
    def utr_num(self):
        """Returns the number of UTR segments (referred to the selected ORF)."""
        return len(self.three_utr + self.five_utr)

    @Metric
    def utr_num_complete(self):
        """Returns the number of UTR segments which are
        complete exons (referred to the selected ORF)."""
        return self.three_utr_num_complete + self.five_utr_num_complete

    @Metric
    def utr_fraction(self):
        """This property calculates the length of the UTR
        of the selected ORF vs. the cDNA length."""
        return 1 - self.selected_cds_fraction

    @Metric
    def utr_length(self):
        """Returns the sum of the 5'+3' UTR lengths"""
        return self.three_utr_length + self.five_utr_length

    @Metric
    def has_start_codon(self):
        """Boolean. True if the selected ORF has a start codon.
        :rtype: bool"""
        return self.__has_start_codon

    @has_start_codon.setter
    def has_start_codon(self, value):
        """Setter. Checks that the argument is boolean.
        :param value: boolean flag
        :type value: bool
        """

        if value not in (None, False, True):
            raise TypeError(
                "Invalid value for has_start_codon: {0}".format(type(value)))
        self.__has_start_codon = value

    @Metric
    def has_stop_codon(self):
        """Boolean. True if the selected ORF has a stop codon.
        :rtype bool
        """
        return self.__has_stop_codon

    @has_stop_codon.setter
    def has_stop_codon(self, value):
        """Setter. Checks that the argument is boolean.
        :param value: list
        :type value: bool
        """

        if value not in (None, False, True):
            raise TypeError(
                "Invalid value for has_stop_codon: {0}".format(type(value)))
        self.__has_stop_codon = value

    @Metric
    def is_complete(self):
        """Boolean. True if the selected ORF has both start and end."""
        return (self.__has_start_codon is True) and (self.__has_stop_codon is True)

    @Metric
    def exon_num(self):
        """This property returns the number of exons of the transcript."""
        return len(self.exons)

    @Metric
    def exon_fraction(self):
        """This property returns the fraction of exons of the transcript
        which are contained in the sublocus.
        If the transcript is by itself, it returns 1. Set from outside."""

        return self.__exon_fraction

    @exon_fraction.setter
    def exon_fraction(self, *args):
        """Setter for exon_fraction. Set from the Locus-type classes.
        :param args: list of values, only the first is retained
        :type args: list(float) | float
        """

        if not isinstance(args[0], (float, int)) or (args[0] <= 0 or args[0] > 1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__exon_fraction = args[0]

    @Metric
    def intron_fraction(self):
        """This property returns the fraction of introns of the transcript
        vs. the total number of introns in the Locus.
        If the transcript is by itself, it returns 1. Set from outside."""
        return self.__intron_fraction

    @intron_fraction.setter
    def intron_fraction(self, *args):
        """Setter for intron_fraction. Set from the Locus-type classes.
        :param args: list of values, only the first is retained
        :type args: list(float) | float
        """

        if not isinstance(args[0], (float, int)) or (args[0] < 0 or args[0] > 1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        if not self.monoexonic and args[0] == 0:
            raise ValueError(
                """It is impossible that the intron fraction is null
                when the transcript has at least one intron!""")
        self.__intron_fraction = args[0]

    @Metric
    def max_intron_length(self):
        """This property returns the greatest intron length for the transcript."""
        if len(self.introns) == 0:
            return 0
        return max(intron[1] + 1 - intron[0] for intron in self.introns)

    @Metric
    def start_distance_from_tss(self):
        """This property returns the distance of the start of the combined CDS
        from the transcript start site.
        If no CDS is defined, it defaults to 0."""
        if len(self.internal_orfs) < 2:
            return self.selected_start_distance_from_tss
        distance = 0
        if self.strand == "+" or self.strand is None:
            for exon in self.exons:
                distance += min(exon[1], self.combined_cds_start - 1) - exon[0] + 1
                if self.combined_cds_start <= exon[1]:
                    break
        elif self.strand == "-":
            exons = reversed(list(self.exons[:]))
            for exon in exons:
                distance += exon[1] + 1 - max(self.combined_cds_start + 1, exon[0])
                if self.combined_cds_start >= exon[0]:
                    break
        return distance

    # pylint: disable=invalid-name
    @Metric
    def selected_start_distance_from_tss(self):
        """This property returns the distance of the start of the best CDS
        from the transcript start site.
        If no CDS is defined, it defaults to 0."""
        if len(self.combined_cds) == 0:
            return 0
        distance = 0
        if self.strand == "+" or self.strand is None:
            for exon in self.exons:
                distance += min(exon[1], self.selected_cds_start - 1) - exon[0] + 1
                if self.selected_cds_start <= exon[1]:
                    break
        elif self.strand == "-":
            exons = reversed(list(self.exons[:]))
            for exon in exons:
                distance += exon[1] + 1 - max(self.selected_cds_start + 1, exon[0])
                if self.selected_cds_start >= exon[0]:
                    break
        return distance

    @Metric
    def selected_end_distance_from_tes(self):
        """This property returns the distance of the end of the best CDS
        from the transcript end site.
        If no CDS is defined, it defaults to 0."""
        if len(self.combined_cds) == 0:
            return 0
        distance = 0
        if self.strand == "-":
            for exon in self.exons:
                distance += min(exon[1], self.selected_cds_end - 1) - exon[0] + 1
                if self.selected_cds_end <= exon[1]:
                    break
        elif self.strand == "+" or self.strand is None:
            exons = reversed(list(self.exons[:]))
            for exon in exons:
                distance += exon[1] + 1 - max(self.selected_cds_end + 1, exon[0])
                if self.selected_cds_end >= exon[0]:
                    break
        return distance

    @Metric
    def selected_end_distance_from_junction(self):
        """This metric returns the distance between the stop codon and the
        nearest downstream junction. In many eukaryotes, this distance
        cannot exceed 50-55 bps, otherwise the transcript becomes a target of NMD.
        If the transcript is not coding or there is no junction downstream of
        the stop codon, the metric returns 0."""

        if len(self.combined_cds) == 0 or self.exon_num == 1:
            return 0
        if self.strand == "+":
            # Case 1: the stop is after the latest junction
            if self.selected_cds_end > max(self.splices):
                return 0
            else:
                return min(list(filter(lambda s: s > self.selected_cds_end,
                                       self.splices))) - self.selected_cds_end
        elif self.strand == "-":
            if self.selected_cds_end < min(self.splices):
                return 0
            else:
                return self.selected_cds_end - max(list(
                    filter(lambda s: s < self.selected_cds_end,
                           self.splices)))

    @Metric
    def end_distance_from_junction(self):
        """This metric returns the distance between the stop codon and
        the nearest downstream junction.
        In many eukaryotes, this distance cannot exceed 50-55 bps
        otherwise the transcript becomes a target of NMD.
        If the transcript is not coding or there is no junction downstream
        of the stop codon, the metric returns 0.
        This metric considers the combined CDS end."""

        if len(self.combined_cds) == 0 or self.exon_num == 1:
            return 0
        if self.strand == "+":
            # Case 1: the stop is after the latest junction
            if self.combined_cds_end > max(self.splices):
                return 0
            else:
                return min(list(filter(
                    lambda s: s > self.combined_cds_end, self.splices))) - self.combined_cds_end
        elif self.strand == "-":
            if self.combined_cds_end < min(self.splices):
                return 0
            else:
                return self.combined_cds_end - max(
                    list(filter(
                        lambda s: s < self.combined_cds_end, self.splices)))

    @Metric
    def end_distance_from_tes(self):
        """This property returns the distance of the end of the combined CDS
        from the transcript end site.
        If no CDS is defined, it defaults to 0."""
        if len(self.internal_orfs) < 2:
            return self.selected_end_distance_from_tes
        distance = 0
        if self.strand == "-":
            for exon in self.exons:
                distance += min(exon[1], self.combined_cds_end - 1) - exon[0] + 1
                if self.combined_cds_end <= exon[1]:
                    break
        elif self.strand == "+" or self.strand is None:
            exons = reversed(list(self.exons[:]))
            for exon in exons:
                distance += exon[1] + 1 - max(self.combined_cds_end + 1, exon[0])
                if self.combined_cds_end >= exon[0]:
                    break
        return distance

    @Metric
    def combined_cds_intron_fraction(self):
        """This property returns the fraction of CDS introns of the transcript
        vs. the total number of CDS introns in the Locus.
        If the transcript is by itself, it returns 1."""
        return self.__combined_cds_intron_fraction

    @combined_cds_intron_fraction.setter
    def combined_cds_intron_fraction(self, value):
        """
        This is the setter for combined_cds_intron_fraction. It checks that the value is
        a valid type, i.e. a float or integer between 0 and 1, before setting it.
        :param value
        :type value: int,float
        """

        if not isinstance(value, (float, int)) or (value < 0 or value > 1):
            raise TypeError(
                "Invalid value for the fraction: {0}".format(value))
        self.__combined_cds_intron_fraction = value

    @Metric
    def selected_cds_intron_fraction(self):
        """This property returns the fraction of CDS introns of
        the selected ORF of the transcript vs. the total number
        of CDS introns in the Locus (considering only the selected ORF).
        If the transcript is by itself, it should return 1.
        """
        return self.__selected_cds_intron_fraction

    @selected_cds_intron_fraction.setter
    def selected_cds_intron_fraction(self, *args):
        """Setter for selected_cds_intron_fraction.
        :param args: either a single float/int or a list (only the first value is retained)
        :type args: list(int) | list(float)
        """

        if not isinstance(args[0], (float, int)) or (args[0] < 0 or args[0] > 1):
            raise TypeError(
                "Invalid value for the fraction: {0}".format(args[0]))
        self.__selected_cds_intron_fraction = args[0]

    @Metric
    def retained_intron_num(self):
        """This property records the number of introns in the transcripts
        which are marked as being retained.
        See the corresponding method in the sublocus class."""
        return len(self.retained_introns)

    @Metric
    def retained_fraction(self):
        """This property returns the fraction of the cDNA which
        is contained in retained introns."""
        return self.__retained_fraction

    @retained_fraction.setter
    def retained_fraction(self, *args):
        """Setter for retained_intron_fraction.
        :param args: either a single float/int or a list (only the first value is retained)
        :type args: list(int) | list(float)
        """

        if not isinstance(args[0], (float, int)) or (args[0] < 0 or args[0] > 1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__retained_fraction = args[0]

    @Metric
    def proportion_verified_introns(self):
        """This metric returns, as a fraction, how many of the transcript introns
        are validated by external data."""
        if self.monoexonic is True:
            return 0
        else:
            return len(self.verified_introns) / len(self.introns)

    @Metric
    def non_verified_introns_num(self):
        """
        This metric returns the number of introns of the transcript which are not validated
        by external data.
        :rtype : int
        """
        return len(self.introns) - len(self.verified_introns)

    @Metric
    def verified_introns_num(self):
        """
        This metric returns the number of introns of the transcript which are validated
        by external data.
        :rtype : int
        """
        return len(self.verified_introns)

    @Metric
    def proportion_verified_introns_inlocus(self):
        """This metric returns, as a fraction, how many of the
        verified introns inside the Locus
        are contained inside the transcript."""
        return self.__proportion_verified_introns_inlocus

    @proportion_verified_introns_inlocus.setter
    def proportion_verified_introns_inlocus(self, *args):
        """Setter for retained_intron_fraction.
        :param args: either a single float/int or a list
        (only the first value is retained)
        :type args: list(int) | list(float)
        """

        if not isinstance(args[0], (float, int)) or (args[0] < 0 or args[0] > 1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))

        value = args[0]
        if value == 0:
            assert len(self.verified_introns) == 0
        assert 0 <= value <= 1
        self.__proportion_verified_introns_inlocus = value

    @Metric
    def num_introns_greater_than_max(self):
        """
        This metric returns the number of introns greater
        than the maximum acceptable intron size
        indicated in the constructor.
        :rtype : int
        """

        return len(list(filter(lambda x: x[1]-x[0]+1 > self.intron_range[1],
                               self.introns)))

    @Metric
    def num_introns_smaller_than_min(self):
        """
        This metric returns the number of introns smaller
        than the mininum acceptable intron size
        indicated in the constructor.
        :rtype : int
        """

        return len(list(filter(lambda x: x[1]-x[0]+1 < self.intron_range[0],
                               self.introns)))
