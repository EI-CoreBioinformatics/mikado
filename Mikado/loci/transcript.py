# coding: utf-8

"""
This module defines the RNA objects. It also defines Metric, a property alias.
"""

# pylint: disable=too-many-lines

import logging
import copy
from sys import intern, maxsize
import re
import inspect
import intervaltree
from ..utilities.log_utils import create_null_logger
from sqlalchemy.sql.expression import desc, asc  # SQLAlchemy imports
from sqlalchemy import and_
from sqlalchemy.ext import baked
from sqlalchemy import bindparam
from ..exceptions import ModificationError, InvalidTranscript, CorruptIndex
from ..serializers.blast_serializer import Query, Hit
from ..serializers.orf import Orf
from .clique_methods import find_communities, define_graph
from ..parsers.GTF import GtfLine
from ..parsers.GFF import GffLine
from .transcript_methods import splitting, retrieval
from .transcript_methods.printing import create_lines_cds
from .transcript_methods.printing import create_lines_no_cds, create_lines_bed, as_bed12
from .transcript_methods.finalizing import finalize
import functools


class Metric(property):
    """Simple aliasing of property. All transcript metrics
    should use this alias, not "property", as a decorator."""

    __category__ = None

    def category_getter(self):
        return self.__category__

    def category_setter(self, category):
        if category is not None and not isinstance(category, (str, bytes)):
            raise TypeError("Invalid category specified: {}".format(category))
        self.__category__ = category

    category = property(lambda self: self.category_getter(),
                        lambda self, v: self.category_setter(v),
                        lambda self: None)

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

    __name__ = intern("transcript")
    __logger = create_null_logger(__name__)

    # Query baking to minimize overhead
    bakery = baked.bakery()
    query_baked = bakery(lambda session: session.query(Query))
    query_baked += lambda q: q.filter(Query.query_name == bindparam("query_name"))

    blast_baked = bakery(lambda session: session.query(Hit))
    blast_baked += lambda q: q.filter(and_(Hit.query == bindparam("query"),
                                           Hit.evalue <= bindparam("evalue")),)

    blast_baked += lambda q: q.order_by(asc(Hit.evalue))
    # blast_baked += lambda q: q.limit(bindparam("max_target_seqs"))

    orf_baked = bakery(lambda session: session.query(Orf))
    orf_baked += lambda q: q.filter(
        Orf.query == bindparam("query"))
    orf_baked += lambda q: q.filter(
        Orf.cds_len >= bindparam("cds_len"))
    orf_baked += lambda q: q.order_by(desc(Orf.cds_len))

    # ######## Class special methods ####################

    def __init__(self, *args,
                 source=None,
                 logger=None,
                 intron_range=(0, maxsize),
                 trust_orf=False):

        """Initialise the transcript object, using a mRNA/transcript line.
        Note: I am assuming that the input line is an object from my own "GFF" class.
        The transcript instance must be initialised by a "(m|r|lnc|whatever)RNA" or
        "transcript" GffLine.

        :param intron_range: range of valid intron size. Any intron shorter
        or longer than this will be flagged.
        :type intron_range: list[int, int]

        :param source: optional source for the transcript
        :type logger: (str|None)
        :param logger: optional logger for the object. WARNING: creating a logger is
        an expensive operation, it is better to do it once per program!
        :type logger: logging.Logger

        :param trust_orf: optional boolean flag. If set, the ORF won't be checked for consistency.
        :type trust_orf: bool

        """

        # Mock setting of base hidden variables
        self.__id = ""
        self._first_phase = None
        self.__logger = None
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
        self.__cdna_length = None
        self._combined_cds_introns = set()
        self._selected_cds_introns = set()
        self._trust_orf = trust_orf  # This is used inside finalizing.__check_internal_orf
        self.__combined_utr = []
        # pylint: enable=invalid-name
        self._selected_internal_orf_cds = []
        # This is used to set the phase if the CDS is loaded from the GFF
        self.__phases = dict()  # will contain (start, phase) for each CDS exon
        self.__blast_score = 0  # Homology score
        self.__derived_children = set()

        # Starting settings for everything else
        self.chrom = None
        self.source = source
        self.feature = "transcript"
        self.start, self.end = None, None
        self.attributes = dict()
        self.exons, self.combined_cds, self.combined_utr = [], [], []
        self.logger = logger
        self.introns = set()
        self.splices = set()
        self.finalized = False  # Flag. We do not want to repeat the finalising more than once.
        self.selected_internal_orf_index = None
        self.non_overlapping_cds = None
        self.verified_introns = set()
        self.segments = []
        self.intron_range = intron_range
        self.internal_orfs = []
        self.blast_hits = []

        # Relative properties
        self.retained_introns = ()
        self.retained_fraction = 0
        self.exon_fraction = self.intron_fraction = 1
        self.cds_intron_fraction = self.selected_cds_intron_fraction = 1

        # Json configuration
        self.__json_conf = None

        # Things that will be populated by querying the database
        self.loaded_bed12 = []
        self.engine, self.session, self.sessionmaker = None, None, None
        # Initialisation of the CDS segments used for finding retained introns
        self.__cds_tree = None
        self.__expandable = False
        # self.query_id = None

        if len(args) == 0:
            return
        else:
            self.__initialize_with_line(args[0])

        self.feature = intern(self.feature)

    def __initialize_with_line(self, transcript_row):
        """
        Private method to copy the necessary attributes from
        an external GTF/GFF3 row.
        :param transcript_row:
        :return:
        """

        if not isinstance(transcript_row, (GffLine, GtfLine)):
            raise TypeError("Invalid data type: {0}".format(type(transcript_row)))
        self.chrom = intern(transcript_row.chrom)
        self.feature = intern(transcript_row.feature)
        # pylint: disable=invalid-name

        # pylint: enable=invalid-name
        self.name = transcript_row.name
        if self.source is None:
            self.source = intern(transcript_row.source)
        self.start = transcript_row.start
        self.strand = transcript_row.strand
        self.end = transcript_row.end
        self.score = transcript_row.score
        self.scores = dict()

        for key, val in transcript_row.attributes.items():
            self.attributes[intern(key)] = val
        self.blast_hits = []
        self.json_conf = None

        if transcript_row.is_transcript is False:
            if transcript_row.is_exon is False or isinstance(transcript_row, GffLine):
                raise TypeError(
                    "Invalid transcript line, the feature should be a transcript:\n{0}".format(
                        transcript_row))
            self.__expandable = True
            self.parent = transcript_row.gene
            self.id = transcript_row.transcript
            self.add_exon(transcript_row)
        else:
            self.parent = transcript_row.parent
            self.id = transcript_row.id

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
            lines = create_lines_cds(self, to_gtf=to_gtf)
        else:
            lines = create_lines_no_cds(self, to_gtf=to_gtf)

        return "\n".join(lines)

    def __eq__(self, other) -> bool:
        """
        :param other: another transcript instance to compare to
        :type other: Mikado.loci_objects.transcript.Transcript

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

        # if hasattr(self, "json_conf") and self.json_conf is not None:
        #     if "requirements" in self.json_conf and "compiled" in self.json_conf["requirements"]:
        #         del state["json_conf"]["requirements"]["compiled"]

        if hasattr(self, "session"):
            if state["session"] is not None:
                state["session"].expunge_all()
                state["session"].close()

            del state["session"]
        if hasattr(self, "sessionmaker"):
            del state["sessionmaker"]
            del state["engine"]

        # if "blast_baked" in state:
        #     del state["blast_baked"]
        #     del state["query_baked"]

        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # Set the logger to NullHandler
        self.logger = None

    # ######## Class instance methods ####################

    def add_exon(self, gffline, feature=None):
        """This function will append an exon/CDS feature to the object.
        :param gffline: an annotation line
        :type gffline: (Mikado.parsers.GFF.GffLine | Mikado.parsers.GTF.GtfLine | tuple | list)
        :type feature: flag to indicate what kind of feature we are adding
        """

        if isinstance(gffline, (tuple, list)):
            assert len(gffline) == 2
            start, end = sorted(gffline)
            phase = None
            if feature is None:
                feature = "exon"
        elif isinstance(gffline, intervaltree.Interval):
            start, end = gffline[0], gffline[1]
            phase = None
            if feature is None:
                feature = "exon"
        elif not isinstance(gffline, (GtfLine, GffLine)):
            raise InvalidTranscript("Unkwown feature type! %s",
                                    type(gffline))
        else:
            start, end = sorted([gffline.start, gffline.end])
            if feature is None:
                feature = gffline.feature
            if self.id not in gffline.parent:
                raise InvalidTranscript(
                    """Mismatch between transcript and exon:
                    {0}
                    {1}
                    {2}""".format(self.id, gffline.parent, gffline))
            assert gffline.is_exon is True, str(gffline)
            assert gffline.chrom == self.chrom, (gffline.chrom, self.chrom)
            phase = gffline.phase

        assert isinstance(start, int) and isinstance(end, int)
        if self.finalized is True:
            raise ModificationError("You cannot add exons to a finalized transcript!")

        if feature.upper().endswith("CDS"):
            store = self.combined_cds
            if phase is not None:
                self.phases[(start, end)] = phase

        elif "combined_utr" in feature or "UTR" in feature.upper():
            store = self.combined_utr
        elif feature.endswith("exon"):
            store = self.exons
        elif feature == "start_codon":
            self.has_start_codon = True
            return
        elif feature == "stop_codon":
            self.has_stop_codon = True
            return
        elif feature == "intron":
            store = self.introns
        else:
            raise InvalidTranscript("Unknown feature: {0}".format(gffline.feature))

        segment = tuple([start, end])
        if segment in store:
            # raise InvalidTranscript(
            #     "Attempt to add {} to {}, but it is already present!".format(
            #         segment, self.id))
            return
        # assert isinstance(segment[0], int) and isinstance(segment[1], int)
        if self.__expandable is True:
            self.start = min([self.start, start])
            self.end = max([self.end, end])
        store.append(segment)
        return

    def add_exons(self, exons, features=None):

        """
        Wrapper of the add_exon method for multiple lines.
        :param exons: An iterable of G(tf|ff)lines to iterate, or of tuples of values.
        :param features: Optional array with the feature types
        :return:
        """
        if features is None:
            features = [None] * len(exons)
        elif isinstance(features, str):
            features = [features] * len(exons)
        else:
            if len(features) != len(exons):
                raise InvalidTranscript("Mismatch between exons and features! %s,\t%s",
                                        exons,
                                        features)

        for exon, feature in zip(exons, features):
            self.add_exon(exon, feature)
        return

    def format(self, format_name,
               with_introns=False,
               with_cds=True,
               all_orfs=False):

        """
        Method to format the string representation of the object. Available formats:
        - GFF3 ("gff3", "gff")
        - GTF ("gtf")
        - BED12 ("bed", "bed12")
        :param format_name: the name of the format to use
        :param with_introns: if True, introns will be printed as well.
        :type with_introns: bool
        :param with_cds: if set to False, CDS lines will be omitted from the output
        :type with_cds: bool
        :param all_orfs: boolean flags that indicates whether all ORFs of a transcript
        should be printed, or only the primary one (default).
        :type: all_orfs: bool
        :return: the formatted string
        :rtype: str
        """

        if format_name not in ("gff", "gtf", "gff3", "bed", "bed12"):
            raise ValueError(
                "Invalid format: {0}. Accepted formats: gff/gff3 (equivalent), gtf".format(
                    format_name))

        self.finalize()  # Necessary to sort the exons
        if format_name in ("bed", "bed12"):
            lines = [create_lines_bed(self)]
        else:
            to_gtf = (format_name == "gtf")
            if with_cds is True:
                lines = create_lines_cds(self,
                                         to_gtf=to_gtf,
                                         with_introns=with_introns,
                                         all_orfs=all_orfs
                                         )
            else:
                lines = create_lines_no_cds(self, to_gtf=to_gtf, with_introns=with_introns)

        return "\n".join(lines)

    def split_by_cds(self):
        """This method is used for transcripts that have multiple ORFs.
        It will split them according to the CDS information into multiple transcripts.
        UTR information will be retained only if no ORF is down/upstream.
        """

        for new_transcript in splitting.split_by_cds(self):
            yield new_transcript

        return

    def as_bed12(self):

        """
        Method to return a BED12 representation of the transcript object.

        :return: bed12 representation of the instance
         :rtype: Mikado.parsers.bed12.BED12
        """

        return as_bed12(self)

    def remove_exon(self, exon):

        """
        Function to remove an exon properly from a Transcript instance.
        :param exon: remove an exon from a Transcript.
        :return:
        """

        if self.finalized is True:
            raise ValueError("Cannot remove a segment from a finalised transcript!")

        if hasattr(exon, "start"):
            exon = tuple([exon.start, exon.end])
        elif not isinstance(exon, tuple):
            exon = ([exon[0], exon[1]])
        else:
            pass

        if exon in self.exons:
            self.exons.remove(exon)
            self.segments.remove(("exon", exon))
            tr = set()
            for segment in self.segments:
                if self.overlap(segment[1], exon) >= 0:
                    tr.add(segment)
            for _ in tr:
                self.segments.remove(_)
            tr = set()
            for segment in self.combined_cds:
                if self.overlap(segment, exon) >= 0:
                    tr.add(segment)
            for _ in tr:
                self.combined_cds.remove(_)
            tr = set()
            for segment in self.combined_utr:
                if self.overlap(segment, exon) >= 0:
                    tr.add(segment)
            for _ in tr:
                self.combined_utr.remove(_)

    def remove_exons(self, exons):

        """
        Wrapper to automate the removal of exons.
        :param exons:
        :return:
        """

        for exon in exons:
            self.remove_exon(exon)

    def remove_utrs(self):
        """Method to strip a transcript from its UTRs.
        It will not execute anything if the transcript lacks a CDS or
        it has more than one ORF defined.
        """

        self.finalize()

        if self.selected_cds_length == 0:
            self.logger.debug("No CDS")
            return
        elif self.three_utr_length + self.five_utr_length == 0:
            self.logger.debug("No UTR")
            return  # No UTR to strip

        elif self.number_internal_orfs > 1:
            return
        elif re.search(r"\.orf[0-9]+$", self.id):
            self.logger.debug("Multiple ORFs already")
            return

        self.finalized = False
        cds_start, cds_end = self.combined_cds[0][0], self.combined_cds[-1][1]
        assert isinstance(cds_start, int)
        assert isinstance(cds_end, int)

        self.exons = self.combined_cds
        self.start = cds_start
        self.end = cds_end
        self.internal_orfs, self.combined_utr = [], []
        self.segments = []
        # Need to recalculate it
        self.__cdna_length = None
        self.finalize()
        assert self.combined_utr == self.three_utr == self.five_utr == [], (
            self.combined_utr, self.three_utr, self.five_utr, self.start, self.end)

    def strip_cds(self, strand_specific=True):
        """Method to completely remove CDS information from a transcript.
        Necessary for those cases where the input is malformed.

        :param strand_specific: boolean flag. If set to False and the transcript is monoexonic,
        the strand will be removed from it.
        """

        self.logger.debug("Stripping CDS from {0}".format(self.id))
        self.finalized = False
        assert len(self.exons) > 0
        if self.monoexonic is True and strand_specific is False:
            self.strand = None

        if self.feature == "mRNA":
            self.feature = "transcript"
        self.combined_cds = []
        self._combined_cds_introns = set()
        self._selected_cds_introns = set()
        self.selected_internal_orf_index = None
        self.combined_utr = []
        self.segments = []
        self.internal_orfs = []
        self.finalize()

    def copy(self):
        """
        Method to return a shallow copy of the current instance.
        :return:
        """

        return copy.copy(self)

    def deepcopy(self):
        """
        Method to return a deep copy of the current instance.
        :return:
        """

        return copy.deepcopy(self)

    def finalize(self):
        """Function to calculate the internal introns from the exons.
        In the first step, it will sort the exons by their internal coordinates.
        """

        if self.finalized is True:
            return

        finalize(self)

        return

    def unfinalize(self):

        if self.finalized is False:
            return

        self.internal_orfs = []
        self.combined_utr = []
        self.finalized = False

    def reverse_strand(self):
        """Method to reverse the strand"""
        if self.strand == "+":
            self.strand = "-"
        elif self.strand == "-":
            self.strand = "+"
        elif self.strand is None:
            pass
        return

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

        retrieval.load_information_from_db(self,
                                           json_conf,
                                           introns=introns,
                                           session=session,
                                           data_dict=data_dict)

    def load_orfs(self, candidate_orfs):

        """
        Thin layer over the load_orfs method from the retrieval module.
        :param candidate_orfs: list of candidate ORFs in BED12 format.
        :return:
        """

        retrieval.load_orfs(self, candidate_orfs)

    def find_overlapping_cds(self, candidate_orfs):

        """
        Thin wrapper for the homonym function in retrieval
        :param candidate_orfs: List of candidate ORFs
        :return:
        """

        return retrieval.find_overlapping_cds(self, candidate_orfs)

    def as_dict(self, remove_attributes=True):

        """
        Method to transform the transcript object into a JSON-friendly representation.
        :return:
        """

        state = dict()

        for key in ["chrom", "source", "start", "end", "strand", "score", "attributes"]:

            state[key] = getattr(self, key)
            if key == "attributes" and remove_attributes is True:
                for subkey in [_ for _ in state[key] if _ not in ("ID", "Parent", "Name")]:
                    del state[key][subkey]

        state["exons"] = []
        for exon in self.exons:
            state["exons"].append([exon[0], exon[1]])

        state["orfs"] = dict()
        state["selected_orf"] = self.selected_internal_orf_index
        for index, orf in enumerate(self.internal_orfs):
            state["orfs"][str(index)] = []
            for segment in orf:
                if segment[0] == "CDS":
                    to_store = [segment[0], [segment[1][0],
                                             segment[1][1]],
                                segment[2]]
                else:
                    to_store = [segment[0], [segment[1][0],
                                             segment[1][1]]]

                state["orfs"][str(index)].append(to_store)
        state["parent"] = getattr(self, "parent")
        state["id"] = getattr(self, "id")
        return state

    def load_dict(self, state):
        self.finalized = False

        for key in ["chrom", "source",
                    "start", "end", "strand", "score",
                    "parent", "id"]:
            if key not in state:
                raise CorruptIndex("Key not found for {}: {}".format(self.id, key))
            if isinstance(state[key], str):
                state[key] = intern(state[key])
            setattr(self, key, state[key])

        self.attributes = {}
        for key, val in state["attributes"].items():
            self.attributes[intern(key)] = val

        self.exons = []
        self.combined_cds = []
        for exon in state["exons"]:
            if len(exon) != 2:
                raise CorruptIndex("Invalid exonic values of {}".format(self.id))
            self.exons.append(tuple(exon))

        self.internal_orfs = []
        try:
            for orf in iter(state["orfs"][_] for _ in sorted(state["orfs"])):
                neworf = []
                for segment in orf:
                    if segment[0] == "CDS":
                        new_segment = (segment[0],
                                       tuple(segment[1]),
                                       int(segment[2]))
                        self.combined_cds.append(tuple(segment[1]))
                    else:
                        new_segment = (segment[0],
                                       tuple(segment[1]))
                    neworf.append(new_segment)

                self.internal_orfs.append(neworf)
            self.selected_internal_orf_index = state["selected_orf"]
        except (ValueError, IndexError):
            raise CorruptIndex("Invalid values for ORFs of {}".format(self.id))

        self.finalize()

    def add_derived_child(self, name):

        """Method to add a derived child (eg TAIR protein) to a transcript"""

        if not isinstance(name, (str, bytes)):
            raise ValueError("Invalid child type: {} (value: {})".format(
                type(name), name))

        self.__derived_children.add(name)

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

        Wrapper for the clique_methods functions.
        As we are interested only in the communities, not the cliques,
        this wrapper discards the cliques
        (first element of the Abstractlocus.find_communities results)
        """
        data = dict((obj, obj) for obj in objects)
        communities = find_communities(define_graph(data, inters=cls.is_intersecting))

        return communities

    @classmethod
    @functools.lru_cache(maxsize=None, typed=True)
    def get_available_metrics(cls) -> list:
        """This function retrieves all metrics available for the class."""

        metrics = [member[0] for member in inspect.getmembers(cls) if
                   "__" not in member[0] and isinstance(cls.__dict__[member[0]], Metric)]

        # metrics = list(x[0] for x in filter(
        #     lambda y: "__" not in y[0] and isinstance(cls.__dict__[y[0]], Metric),
        #     inspect.getmembers(cls)))
        # assert "tid" in metrics and "parent" in metrics and "score" in metrics
        _metrics = sorted([metric for metric in metrics])
        final_metrics = ["tid", "parent", "score"] + _metrics
        return final_metrics

    # ###################Class properties##################################

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

    @property
    def json_conf(self):
        """
        Configuration dictionary. It can be None.
        :return:
        """

        return self.__json_conf

    @json_conf.setter
    def json_conf(self, json_conf):

        """
        Setter for the configuration dictionary.
        :param json_conf: None or a dictionary
        :type json_conf: (None | dict)
        :return:
        """

        assert isinstance(json_conf, dict) or json_conf is None
        self.__json_conf = json_conf

    @logger.deleter
    def logger(self):
        """
        Destroyer for the logger. It sets the internal __logger attribute to None.
        """
        self.__logger = None

    @property
    def phases(self):
        """
        This property contains the first phase gleaned for each internal ORF from the
         GFF. Structure:

         dict[(start,end)] = phase

        where (start, end) are the coordinates of the coding segment

        :return: __phases, a list
        :rtype: dict
        """
        return self.__phases

    @phases.setter
    def phases(self, phases):
        """
        Setter for phases. The input must be a list.
        :param phases:
        :return:
        """

        assert isinstance(phases, dict)
        self.__phases = phases

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
        self.__id = intern(newid)
    # pylint: enable=invalid-name

    @property
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

    @property
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

        self.attributes["gene_id"] = self.__parent
        self.__parent = [intern(_) for _ in self.__parent]

    @property
    def gene(self):

        if "gene_id" not in self.attributes:
            self.attributes["gene_id"] = self.parent[0]

        return self.attributes["gene_id"]

    @property
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

    @property
    @functools.lru_cache(maxsize=None, typed=True)
    def available_metrics(self) -> list:
        """Return the list of available metrics, using the "get_metrics" function."""
        return self.get_available_metrics()

    @property
    def name(self):

        if "Name" not in self.attributes:
            keys = [_ for _ in self.attributes.keys() if _.lower() == "name"]
            if len(keys) > 0:
                self.name = self.attributes[keys[0]]
                for key in keys:
                    del self.attributes[key]
            else:
                self.name = None
        return self.attributes["Name"]

    @name.setter
    def name(self, name):
        if not isinstance(name, (type(None), str)):
            raise ValueError("Invalid name: {0}".format(name))

        if name is not None:
            name = intern(name)
        self.attributes["Name"] = name

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
            self.selected_internal_orf_index = None
            return tuple([])
        else:
            assert self.selected_internal_orf_index is not None
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
            return list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                        utr_segment[0] == "UTR" and utr_segment[1][0] > self.selected_cds_start)
            #
            #
            # return list(
            #     filter(lambda exon: exon[0] == "UTR" and exon[1][0] > self.selected_cds_start,
            #            self.selected_internal_orf))
        else:
            return list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                        utr_segment[0] == "UTR" and utr_segment[1][1] < self.selected_cds_start)
            #
            #
            # return list(
            #     filter(lambda exon: exon[0] == "UTR" and exon[1][1] < self.selected_cds_start,
            #            self.selected_internal_orf))

    @property
    def three_utr(self):
        """Returns the exons in the 3' UTR of the selected ORF.
        If the end codon is absent, no UTR is given."""
        if len(self.combined_cds) == 0:
            return []
        if self.strand == "-":
            return list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                        utr_segment[0] == "UTR" and utr_segment[1][1] < self.selected_cds_end)
            # filter(lambda exon: exon[0] == "UTR" and exon[1][1] < self.selected_cds_end,
            #        self.selected_internal_orf))
        else:
            return list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                        utr_segment[0] == "UTR" and utr_segment[1][0] > self.selected_cds_end)
            #
            #
            # return list(
            #     filter(lambda exon: exon[0] == "UTR" and exon[1][0] > self.selected_cds_end,
            #            self.selected_internal_orf))

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
            assert isinstance(internal_cds[0][1], tuple), internal_cds[0]
            length = sum(x[1][1] - x[1][0] + 1 for x in internal_cds if
                         x[0] == "CDS")
            lengths.append(length)
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
                segments = set([segment[1] for segment in internal_cds if
                                segment[0] == "CDS"])
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

        return self._combined_cds_introns

        # if self.number_internal_orfs < 2:
        #     return self.selected_cds_introns
        # if self.number_internal_orfs == 0 or len(self.combined_cds) < 2:
        #     return set()
        #
        # cintrons = []
        # for position in range(len(self.combined_cds) - 1):
        #     former = self.combined_cds[position]
        #     latter = self.combined_cds[position + 1]
        #     junc = intervaltree.Interval(former[1] + 1, latter[0] - 1)
        #     if junc in self.introns:
        #         cintrons.append(junc)
        # cintrons = set(cintrons)
        # return cintrons

    @property
    def selected_cds_introns(self):
        """This property returns the introns which are located between
        CDS segments in the selected ORF."""

        return self._selected_cds_introns

        # if len(self.selected_cds) < 2:
        #     return set()
        # if self.number_internal_orfs == 0 or len(self.combined_cds) < 2:
        #     return set()
        #
        # cintrons = []
        # for first, second in zip(self.selected_cds[:-1], self.selected_cds[1:]):
        #     cintrons.append(
        #         intervaltree.Interval(first[1] + 1,
        #                               second[0] - 1)
        #     )
        # cintrons = set(cintrons)
        # assert len(cintrons) > 0
        # return cintrons

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

        if ((not isinstance(combined, list)) or
                any(self.__wrong_combined_entry(comb) for comb in combined)):
            raise TypeError("Invalid value for combined CDS: {0}".format(combined))

        # if len(combined) > 0:
        #     if isinstance(combined[0], tuple):
        #         try:
        #             combined = [intervaltree.Interval(_[0], _[1]) for _ in combined]
        #         except IndexError:
        #             raise IndexError(combined)
        #     else:
        #         assert isinstance(combined[0], intervaltree.Interval)

        self.__combined_cds = combined

    @staticmethod
    def __wrong_combined_entry(to_test):
        """
        Private method to test the correctness of entries for "combined"
        data
        :param to_test:
        :return:
        """
        if not isinstance(to_test, tuple) or len(to_test) != 2:
            return True
        elif to_test[1] < to_test[0]:
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
            return None
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
            self.__selected_cds = [segment[1] for segment in self.selected_internal_orf if
                                   segment[0] == "CDS"]
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
            try:
                return self.selected_cds[-1][1]
            except IndexError as exc:
                self.logger.exception(
                    "{0}, selected CDS: {1}, combined CDS: {2}, \
index {3}, internal ORFs: {4}".format(
                        exc,
                        self.selected_cds,
                        self.combined_cds,
                        self.selected_internal_orf_index,
                        self.internal_orfs))
                raise

    @property
    def monoexonic(self):
        """
        Property. True if the transcript has only one exon, False otherwise.
        :rtype bool
        """
        if len(self.exons) == 1:
            return True
        return False

    @property
    def is_coding(self):
        """
        Simple property to investigate whether a transcript is coding or not
        :return: boolean value
        :rtype: bool
        """

        return len(self.combined_cds) > 0

    @property
    def cds_tree(self):
        """
        This property returns an interval tree of the CDS segments.
        Used to calculate the non-coding parts of the CDS.
        :rtype: intervaltree.Intervaltree
        """
        return self.__cds_tree

    @cds_tree.setter
    def cds_tree(self, segments):
        """
        Setter for CDS tree. It checks that the calculated tree is actually valid.
        :param segments: the interval tree to be set.
        :type segments: intervaltree.Intervaltree
        :return:
        """

        if segments is None:
            pass

        elif isinstance(segments, intervaltree.IntervalTree):
            assert len(segments) == len(self.combined_cds)
        else:
            raise TypeError("Invalid cds segments: %s, type %s",
                            segments, type(segments))

        self.__cds_tree = segments

    @property
    def derived_children(self):
        """
        This property stores the names of children features (eg TAIR proteins)
        that this transcript gives origin to.
        :return:
        """

        return self.__derived_children

    # ################### Class metrics ##################################

    # Disable normal checks on names and hidden methods, as
    # checkers get confused by the Metric method
    # pylint: disable=method-hidden,invalid-name

    @Metric
    def combined_cds_length(self):
        """This property return the length of the CDS part of the transcript."""
        c_length = sum([c[1] - c[0] + 1 for c in self.combined_cds])
        if len(self.combined_cds) > 0:
            assert c_length > 0
        return c_length

    combined_cds_length.category = "CDS"

    @Metric
    def combined_cds_num(self):
        """This property returns the number of non-overlapping CDS segments
        in the transcript."""
        return len(self.combined_cds)

    combined_cds_num.category = "CDS"

    @Metric
    def combined_cds_num_fraction(self):
        """This property returns the fraction of non-overlapping CDS segments
        in the transcript
        vs. the total number of exons"""
        return len(self.combined_cds) / len(self.exons)

    combined_cds_num_fraction.category = "CDS"

    @Metric
    def combined_cds_fraction(self):
        """This property return the percentage of the CDS part of the transcript
        vs. the cDNA length"""
        return self.combined_cds_length / self.cdna_length

    combined_cds_fraction.category = "CDS"

    @Metric
    def combined_utr_length(self):
        """This property return the length of the UTR part of the transcript."""
        return sum([e[1] - e[0] + 1 for e in self.combined_utr])

    combined_utr_length.category = "UTR"

    @Metric
    def combined_utr_fraction(self):
        """This property returns the fraction of the cDNA which is not coding according
        to any ORF. Complement of combined_cds_fraction"""
        return 1 - self.combined_cds_fraction

    combined_utr_fraction.__category = "UTR"

    @Metric
    def cdna_length(self):
        """This property returns the length of the transcript."""
        try:
            self.__cdna_length = sum([e[1] - e[0] + 1 for e in self.exons])
        except TypeError:
            raise TypeError(self.exons)
        return self.__cdna_length

    cdna_length.category = "cDNA"

    @Metric
    def number_internal_orfs(self):
        """This property returns the number of ORFs inside a transcript."""
        return len(self.internal_orfs)

    number_internal_orfs.category = "CDS"

    @Metric
    def selected_cds_length(self):
        """This property calculates the length of the CDS selected as best inside
        the cDNA."""
        self.finalize()
        if len(self.combined_cds) == 0:
            self.__max_internal_orf_length = 0
        else:
            self.__max_internal_orf_length = sum(
                _[1][1] - _[1][0] + 1 for _ in self.selected_internal_orf if _[0] == "CDS")

        return self.__max_internal_orf_length

    selected_cds_length.category = "CDS"

    @Metric
    def selected_cds_num(self):
        """This property calculates the number of CDS exons for the selected ORF"""
        return sum(1 for exon in self.selected_internal_orf if exon[0] == "CDS")

    selected_cds_num.category = "CDS"

    @Metric
    def selected_cds_fraction(self):
        """This property calculates the fraction of the selected CDS vs. the cDNA length."""
        return self.selected_cds_length / self.cdna_length

    selected_cds_fraction.category = "CDS"

    @Metric
    def highest_cds_exons_num(self):
        """Returns the number of CDS segments in the selected ORF
        (irrespective of the number of exons involved)"""
        return sum(1 for _ in self.selected_internal_orf if _[0] == "CDS")
        # return len(list(filter(lambda x: x[0] == "CDS", self.selected_internal_orf)))

    highest_cds_exons_num.category = "CDS"

    @Metric
    def selected_cds_exons_fraction(self):
        """Returns the fraction of CDS segments in the selected ORF
        (irrespective of the number of exons involved)"""
        return self.highest_cds_exon_number / len(self.exons)

        # return len(list(filter(lambda x: x[0] == "CDS",
        #                        self.selected_internal_orf))) / len(self.exons)

    selected_cds_exons_fraction.category = "CDS"

    @Metric
    def highest_cds_exon_number(self):
        """This property returns the maximum number of CDS segments
        among the ORFs; this number can refer to an ORF *DIFFERENT*
        from the maximal ORF."""
        if len(self.internal_orfs) == 0:
            return 0

        cds_numbers = []

        for cds in self.internal_orfs:
            cds_numbers.append(sum(1 for segment in cds if segment[0] == "CDS"))
            # len(list(filter(lambda x: x[0] == "CDS", cds))))
        return max(cds_numbers)

    highest_cds_exon_number.category = "CDS"

    @Metric
    def selected_cds_number_fraction(self):
        """This property returns the proportion of best possible CDS segments
        vs. the number of exons. See selected_cds_number."""
        return self.selected_cds_num / self.exon_num

    selected_cds_number_fraction.category = "CDS"

    @Metric
    def cds_not_maximal(self):
        """This property returns the length of the CDS excluded from the selected ORF."""
        if len(self.internal_orfs) < 2:
            return 0
        return self.combined_cds_length - self.selected_cds_length

    cds_not_maximal.category = "CDS"

    @Metric
    def cds_not_maximal_fraction(self):
        """This property returns the fraction of bases not in the selected ORF compared to
        the total number of CDS bases in the cDNA."""
        if self.combined_cds_length == 0:
            return 0
        else:
            return self.cds_not_maximal / self.combined_cds_length

    cds_not_maximal_fraction.category = "CDS"

    @Metric
    def five_utr_length(self):
        """Returns the length of the 5' UTR of the selected ORF."""
        if len(self.combined_cds) == 0:
            return 0
        return sum(utr[1] - utr[0] + 1 for utr in self.five_utr)

    five_utr_length.category = "UTR"

    @Metric
    def five_utr_num(self):
        """This property returns the number of 5' UTR segments for the selected ORF."""
        return len(self.five_utr)

    five_utr_num.category = "UTR"

    @Metric
    def five_utr_num_complete(self):
        """This property returns the number of 5' UTR segments for the selected ORF,
        considering only those which are complete exons."""
        return sum(1 for utr in self.five_utr if utr in self.exons)

    five_utr_num_complete.category = "UTR"

    @Metric
    def three_utr_length(self):
        """Returns the length of the 5' UTR of the selected ORF."""
        if len(self.combined_cds) == 0:
            return 0
        return sum(x[1] - x[0] + 1 for x in self.three_utr)

    three_utr_length.__category = "UTR"

    @Metric
    def three_utr_num(self):
        """This property returns the number of 3' UTR segments
        (referred to the selected ORF)."""
        return len(self.three_utr)

    three_utr_num.category = "UTR"

    @Metric
    def three_utr_num_complete(self):
        """This property returns the number of 3' UTR segments for the selected ORF,
        considering only those which are complete exons."""
        return sum(1 for utr in self.three_utr if utr in self.exons)

    three_utr_num_complete.category = "UTR"

    @Metric
    def utr_num(self):
        """Returns the number of UTR segments (referred to the selected ORF)."""
        return len(self.three_utr + self.five_utr)

    utr_num.category = "UTR"

    @Metric
    def utr_num_complete(self):
        """Returns the number of UTR segments which are
        complete exons (referred to the selected ORF)."""
        return self.three_utr_num_complete + self.five_utr_num_complete

    utr_num_complete.category = "UTR"

    @Metric
    def utr_fraction(self):
        """This property calculates the length of the UTR
        of the selected ORF vs. the cDNA length."""
        return 1 - self.selected_cds_fraction

    utr_fraction.category = "UTR"

    @Metric
    def utr_length(self):
        """Returns the sum of the 5'+3' UTR lengths"""
        return self.three_utr_length + self.five_utr_length

    utr_length.category = "UTR"

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

    has_start_codon.category = "CDS"

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

    has_stop_codon.category = "CDS"

    @Metric
    def is_complete(self):
        """Boolean. True if the selected ORF has both start and end."""
        return (self.__has_start_codon is True) and (self.__has_stop_codon is True)

    is_complete.category = "CDS"

    @Metric
    def exon_num(self):
        """This property returns the number of exons of the transcript."""
        return len(self.exons)

    exon_num.category = "cDNA"

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

    exon_fraction.category = "Locus"

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

    intron_fraction.category = "Locus"

    @Metric
    def max_intron_length(self):
        """This property returns the greatest intron length for the transcript."""
        if len(self.introns) == 0:
            return 0
        return max(intron[1] + 1 - intron[0] for intron in self.introns)

    max_intron_length.category = "Intron"

    @Metric
    def min_intron_length(self):
        if len(self.introns) == 0:
            return 0
        return min(intron[1] + 1 - intron[0] for intron in self.introns)

    min_intron_length.category = "Intron"

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

    start_distance_from_tss.category = "CDS"

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

    selected_start_distance_from_tss.category = "CDS"

    @Metric
    def source_score(self):

        """This metric returns a score that is assigned to the transcript
        in virtue of its origin."""

        if (self.json_conf and "pick" in self.json_conf and
                self.source in self.json_conf["pick"]["source_score"]):
            return self.json_conf["pick"]["source_score"][self.source]
        else:
            return 0

    source_score.category = "External"

    @Metric
    def selected_end_distance_from_tes(self):
        """This property returns the distance of the end of the best CDS
        from the transcript end site.
        If no CDS is defined, it defaults to 0."""
        if self.is_coding is False:
            return 0

        distance = 0
        if self.strand == "+":
            # Case 1: the stop is after the latest junction
            for exon in sorted([_ for _ in self.exons
                                if _[1] > self.selected_cds_end]):
                if exon[0] <= self.selected_cds_end <= exon[1]:
                    distance += exon[1] - self.selected_cds_end
                else:
                    distance += exon[1] - exon[0] + 1
        elif self.strand == "-":
            for exon in sorted([_ for _ in self.exons
                                if _[0] < self.selected_cds_end], reverse=True):
                if exon[0] <= self.selected_cds_end <= exon[1]:
                    distance += self.selected_cds_end - exon[0]  # Exclude end
                else:
                    distance += exon[1] - exon[0] + 1
        return distance

    selected_end_distance_from_tes.category = "CDS"

    @Metric
    def selected_end_distance_from_junction(self):
        """This metric returns the distance between the stop codon and the
        last junction of the transcript. In many eukaryotes, this distance
        cannot exceed 50-55 bps, otherwise the transcript becomes a target of NMD.
        If the transcript is not coding or there is no junction downstream of
        the stop codon, the metric returns 0."""

        if self.monoexonic is True or self.is_coding is False:
            return 0

        self.finalize()
        distance = 0
        if self.strand == "+":
            # Case 1: the stop is after the latest junction
            if self.selected_cds_end > max(self.splices):
                pass
            else:
                for exon in sorted([_ for _ in self.exons
                                    if _[1] > self.selected_cds_end])[:-1]:
                    if exon[0] <= self.selected_cds_end <= exon[1]:
                        distance += exon[1] - self.selected_cds_end  # Exclude end
                    else:
                        distance += exon[1] - exon[0] + 1
        elif self.strand == "-":
            if self.selected_cds_end < min(self.splices):
                pass
            else:
                for exon in sorted([_ for _ in self.exons
                                    if _[0] < self.selected_cds_end], reverse=True)[:-1]:
                    if exon[0] <= self.selected_cds_end <= exon[1]:
                        distance += self.selected_cds_end - exon[0]  # Exclude end
                    else:
                        distance += exon[1] - exon[0] + 1
        return distance

    selected_end_distance_from_junction.category = "CDS"

    @Metric
    def end_distance_from_junction(self):
        """This metric returns the cDNA distance between the stop codon and
        the last junction in the transcript.
        In many eukaryotes, this distance cannot exceed 50-55 bps
        otherwise the transcript becomes a target of NMD.
        If the transcript is not coding or there is no junction downstream
        of the stop codon, the metric returns 0.
        This metric considers the combined CDS end."""

        if self.monoexonic is True or self.is_coding is False:
            return 0

        distance = 0
        self.finalize()
        if self.strand == "+":
            # Case 1: the stop is after the latest junction
            if self.combined_cds_end > max(self.splices):
                pass
            else:
                for exon in sorted([_ for _ in self.exons
                                    if _[1] > self.combined_cds_end])[:-1]:
                    if exon[0] <= self.combined_cds_end <= exon[1]:
                        distance += exon[1] - self.combined_cds_end
                    else:
                        distance += exon[1] - exon[0] + 1
        elif self.strand == "-":
            if self.combined_cds_end < min(self.splices):
                pass
            else:
                for exon in sorted([_ for _ in self.exons
                                    if _[0] < self.combined_cds_end], reverse=True)[:-1]:
                    if exon[0] <= self.combined_cds_end <= exon[1]:
                        distance += self.combined_cds_end - exon[0]  # Exclude end
                    else:
                        distance += exon[1] - exon[0] + 1
        return distance

    end_distance_from_junction.category = "CDS"

    @Metric
    def end_distance_from_tes(self):
        """This property returns the distance of the end of the combined CDS
        from the transcript end site.
        If no CDS is defined, it defaults to 0."""

        if self.is_coding is False:
            return 0

        distance = 0
        if self.strand == "+":
            # Case 1: the stop is after the latest junction
            for exon in sorted([_ for _ in self.exons
                                if _[1] > self.combined_cds_end]):
                if exon[0] <= self.combined_cds_end <= exon[1]:
                    distance += exon[1] - self.combined_cds_end
                else:
                    distance += exon[1] - exon[0] + 1
        elif self.strand == "-":
            for exon in sorted([_ for _ in self.exons
                                if _[0] < self.combined_cds_end], reverse=True):
                if exon[0] <= self.combined_cds_end <= exon[1]:
                    distance += self.combined_cds_end - exon[0]  # Exclude end
                else:
                    distance += exon[1] - exon[0] + 1
        return distance

    end_distance_from_tes.category = "CDS"

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

    combined_cds_intron_fraction.category = "Locus"

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

    selected_cds_intron_fraction.category = "CDS"

    @Metric
    def retained_intron_num(self):
        """This property records the number of introns in the transcripts
        which are marked as being retained.
        See the corresponding method in the sublocus class."""
        return len(self.retained_introns)

    retained_intron_num.category = "Locus"

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

    retained_fraction.category = "Locus"

    @Metric
    def proportion_verified_introns(self):
        """This metric returns, as a fraction, how many of the transcript introns
        are validated by external data. Monoexonic transcripts are set to 1."""
        if self.monoexonic is True:
            return 1
        else:
            return len(self.verified_introns) / len(self.introns)

    proportion_verified_introns.category = "External"

    @Metric
    def non_verified_introns_num(self):
        """
        This metric returns the number of introns of the transcript which are not validated
        by external data.
        :rtype : int
        """
        num = len(set.difference(self.introns, self.verified_introns))
        if num < 0:
            # This is a clear error
            self.logger.error("Erroneous number of verified introns for %s; total %d, verified %d, subtraction %d",
                              self.id, len(self.introns), len(self.verified_introns), num)
            self.logger.error("Introns for %s: %s; verified: %s. Resetting verified to 0.",
                              self.id, self.introns, self.verified_introns)
            num = len(self.introns)

        return num

    non_verified_introns_num.category = "External"

    @Metric
    def verified_introns_num(self):
        """
        This metric returns the number of introns of the transcript which are validated
        by external data.
        :rtype : int
        """
        return len(self.verified_introns)

    verified_introns_num.category = "External"

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

    proportion_verified_introns_inlocus.category = "Locus"

    @Metric
    def num_introns_greater_than_max(self):
        """
        This metric returns the number of introns greater
        than the maximum acceptable intron size
        indicated in the constructor.
        :rtype : int
        """

        return sum(1 for intron in self.introns if
                   intron[1] - intron[0] + 1 > self.intron_range[1])
        #
        # return len(list(filter(lambda x: x[1]-x[0]+1 > self.intron_range[1],
        #                        self.introns)))

    num_introns_greater_than_max.category = "Intron"

    @Metric
    def num_introns_smaller_than_min(self):
        """
        This metric returns the number of introns smaller
        than the mininum acceptable intron size
        indicated in the constructor.
        :rtype : int
        """

        return sum(1 for intron in self.introns if
                   intron[1] - intron[0] + 1 < self.intron_range[0])

        #
        # return len(list(filter(lambda x: x[1]-x[0]+1 < self.intron_range[0],
        #                        self.introns)))

    num_introns_smaller_than_min.category = "Intron"

    @Metric
    def snowy_blast_score(self):

        """
        Metric that indicates how good a hit is compared to the competition, in terms of BLAST
        similarities.
        As in SnowyOwl, the score for each hit is calculated by taking the percentage of positive
        matches and dividing it by (2 * len(self.blast_hits)).
        IMPORTANT: when splitting transcripts by ORF, a blast hit is added to the new transcript
        only if it is contained within the new transcript.
        This WILL screw up a bit the homology score.
        :return
        """

        if len(self.blast_hits) == 0:
            self.__blast_score = 0
        elif self.__blast_score == 0 and len(self.blast_hits) > 0:
            score = 0
            for hit in self.blast_hits:
                score += hit["global_positives"]/(2 * len(self.blast_hits))
            self.__blast_score = score

        return self.__blast_score

    snowy_blast_score.category = "External"

    @Metric
    def best_bits(self):
        """Metric that returns the best BitS associated with the transcript."""

        return max([0] + [hit["bits"] for hit in self.blast_hits])

    best_bits.category = "External"

    @Metric
    def blast_score(self):
        """
        Interchangeable alias for testing different blast-related scores.
        Current: best bit score.
        :return:
        """
        # return self.snowy_blast_score
        return self.best_bits

    blast_score.category = "External"

    @Metric
    def canonical_intron_proportion(self):

        """
        This metric returns the proportion of canonical introns
         of the transcript on its total number of introns.
        :return:
        """

        return float(self.attributes.get("canonical_proportion", 0))

    canonical_intron_proportion.category = "Intron"