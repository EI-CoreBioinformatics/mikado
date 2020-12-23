# coding: utf-8

"""
This module defines the RNA objects. It also defines Metric, a property alias.
"""

# pylint: disable=too-many-lines

import copy
import functools
import inspect
import logging
import re
from ast import literal_eval
from sys import intern, maxsize
import operator
import intervaltree
from sqlalchemy import and_
from sqlalchemy import bindparam
from sqlalchemy.ext import baked
from sqlalchemy.sql.expression import desc, asc  # SQLAlchemy imports
from ..exceptions import ModificationError, InvalidTranscript, CorruptIndex
from ..parsers.GFF import GffLine
from ..parsers.GTF import GtfLine
from ..parsers.bed12 import BED12
from ..serializers.blast_serializer import Query, Hit
from ..serializers.external import External, ExternalSource
from ..serializers.orf import Orf
from ..transcripts.clique_methods import find_communities, define_graph
from ..utilities.log_utils import create_null_logger
from .transcript_methods import splitting, retrieval
from .transcript_methods.finalizing import finalize
from .transcript_methods.printing import create_lines_cds
from .transcript_methods.printing import create_lines_no_cds, create_lines_bed, as_bed12
from ..utilities.intervaltree import Interval, IntervalTree
from ..utilities.namespace import Namespace
from collections.abc import Hashable
import numpy as np
import pysam


class Metric(property):
    """
    Simple aliasing of property. All transcript metrics
    should use this alias, not "property", as a decorator.
    Any metric has two properties of its own:
    - "category": this is just a string which acts as a general label for the metric
    (e.g., "CDS", "Locus", etc.)
    - "usable_raw": a boolean flag that indicates whether the metric can be used as
    a raw score in itself. *Only metrics whose possible values are between 0 and 1 can be used
    this way*.
    """

    __category__ = None
    __name__ = "Metric"

    def category_getter(self):
        return self.__category__

    def category_setter(self, category):
        if category is not None and not isinstance(category, (str, bytes)):
            raise TypeError("Invalid category specified: {}".format(category))
        self.__category__ = category

    category = property(fget=category_getter,
                        fset=category_setter,
                        fdel=None)

    # The "raw_valid"
    __usable_raw__ = False

    def raw_getter(self):

        return self.__usable_raw__

    def raw_setter(self, boolean):
        if not isinstance(boolean, bool):
            raise ValueError("The \"usable_raw\" property must be boolean!")
        self.__usable_raw__ = boolean

    usable_raw = property(fget=raw_getter, fset=raw_setter, fdel=None)

    __rtype__ = None

    def rtype_getter(self):

        return self.__rtype__

    __valid_types = [None, "int", "float", "str", "bytes",
                     "complex", "dict", "bytearray", "set", "frozenset",
                     "list", "bool", "object", "slice", "tuple"]

    def rtype_setter(self, string):
        if not isinstance(string, (str, bytes)) and string is not None:
            raise ValueError("The \"rtype\" property must be a string!")
        elif isinstance(string, bytes):
            string = string.decode()
        if string not in self.__valid_types:
            raise ValueError("{} is an invalid type for a Metric".format(string))

        self.__rtype__ = string

    rtype = property(fget=rtype_getter,
                     fset=rtype_setter,
                     fdel=None)

    @property
    def settable(self):
        return self.fset is not None

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

    # External scores
    external_baked = bakery(lambda session: session.query(External))
    external_baked += lambda q: q.filter()

    external_sources = bakery(lambda session: session.query(ExternalSource))

    # ######## Class special methods ####################
    def __init__(self, *args,
                 source=None,
                 logger=None,
                 intron_range=(0, maxsize),
                 trust_orf=False,
                 is_reference=None,
                 accept_undefined_multi=False):

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

        :param accept_undefined_multi: optional boolean flag. If set to True, Mikado will not mark the transcript as
        invalid if it is multiexonic but without a defined strand. Normally, we would discard such transcripts.
        :type accept_undefined_multi: bool
        """

        # Mock setting of base hidden variables
        self.__id = ""
        self.__finalized = False  # Flag. We do not want to repeat the finalising more than once.
        self.__source = None
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
        self.__exon_num = None

        self.__combined_cds_length = 0
        self.__selected_cds = []
        self.__cdna_length = None
        self._combined_cds_introns = set()
        self._selected_cds_introns = set()
        self.__selected_cds_locus_fraction = 0
        self.__combined_cds_locus_fraction = 0
        self._trust_orf = trust_orf  # This is used inside finalizing.__check_internal_orf
        self.__combined_utr = []
        # pylint: enable=invalid-name
        self._selected_internal_orf_cds = []
        # This is used to set the phase if the CDS is loaded from the GFF
        self.__phases = dict()  # will contain (start, phase) for each CDS exon
        self.__blast_score = 0  # Homology score
        self.__derived_children = set()
        self._original_source = None
        self.__external_scores = Namespace(default=(0, False))
        self.__internal_orf_transcripts = []
        self.__cds_not_maximal = None

        # Starting settings for everything else
        self.__chrom = None
        self.is_reference = is_reference
        self.feature = "transcript"
        self.__start, self.__end = None, None
        self.attributes = dict()
        self.exons, self.combined_cds, self.combined_utr = [], [], []
        self.stop_codon = []
        self.logger = logger
        self.introns = set()
        self.splices = set()
        self.selected_internal_orf_index = None
        self.non_overlapping_cds = None
        self.__verified_introns = set()
        self.segments = []
        self.intron_range = intron_range
        self.internal_orfs = []
        self.blast_hits = []

        # Relative properties
        self.retained_introns = ()
        self.retained_fraction = 0
        self.__cds_disrupted_by_ri = False
        self.exon_fraction = self.intron_fraction = 1
        self.cds_intron_fraction = self.selected_cds_intron_fraction = 1

        # Json configuration
        self.__json_conf = None

        # Things that will be populated by querying the database
        self.loaded_bed12 = []
        self.engine, self.session, self.sessionmaker = None, None, None
        # Initialisation of the CDS segments used for finding retained introns
        self.__cds_tree = IntervalTree()
        self.__expandable = False
        self.__segmenttree = IntervalTree()
        self.__cds_introntree = IntervalTree()
        self._possibly_without_exons = False
        self._accept_undefined_multi = accept_undefined_multi

        if len(args) == 0:
            pass
        else:
            self.__initialize_with_line(args[0])
            self._original_source = self.source
            self.feature = intern(self.feature)

        if source is not None:
            self.source = self._original_source = source

    def __initialize_with_line(self, transcript_row):
        """
        Private method to copy the necessary attributes from
        an external GTF/GFF3 row.
        :param transcript_row:
        :return:
        """

        if isinstance(transcript_row, (str, bytes)):
            if isinstance(transcript_row, bytes):
                transcript_row = transcript_row.decode()
            _ = GffLine(transcript_row)
            if _.header is False and _.is_transcript is True and _.id is not None:
                transcript_row = _
            else:
                _ = GtfLine(transcript_row)
                if _.header is False and _.is_transcript is True and _.id is not None:
                    transcript_row = _
                else:
                    _ = BED12(transcript_row)
                    if _.header is False and _.name is not None:
                        transcript_row = _

        if isinstance(transcript_row, (GffLine, GtfLine)):
            self.__initialize_with_gf(transcript_row)
        elif isinstance(transcript_row, BED12):
            self.__initialize_with_bed12(transcript_row)
        elif isinstance(transcript_row, pysam.AlignedSegment):
            self.__initialize_with_bam(transcript_row)
        else:
            raise TypeError("Invalid data type: {0}".format(type(transcript_row)))

    def __initialize_with_bed12(self, transcript_row: BED12):

        """
        :param transcript_row:
        :type transcript_row: pysam.AlignedSegment
        :return:
        """

        if transcript_row.header is True:
            raise InvalidTranscript("I cannot initialise a valid transcript with a header (ie empty) BED line.")

        self.chrom = transcript_row.chrom
        self.name = self.id = transcript_row.name
        self.start, self.end = transcript_row.start, transcript_row.end
        self.score = transcript_row.score
        self.strand = transcript_row.strand
        exon_starts = np.array([_ + self.start for _ in transcript_row.block_starts])
        exon_ends = exon_starts + np.array(transcript_row.block_sizes) - 1
        isizes = exon_starts[1:] - exon_ends[:-1]  # This functions also for monoexonic
        if np.where(isizes < 0)[0].size > 0:
            raise InvalidTranscript("Overlapping exons found for {}!".format(self.id))
        if np.where(isizes <= 1)[0].size > 0:
            self.logger.debug("Merging touching exons")
            exon_starts = np.concatenate((np.array([exon_starts[0]]), exon_starts[np.where(isizes > 1)[0] + 1]))
            exon_ends = np.concatenate((exon_ends[np.where(isizes > 1)[0]], np.array([exon_ends[-1]])))

        self.add_exons(list(zip(list(exon_starts), list(exon_ends))))
        self.parent = getattr(transcript_row, "parent", transcript_row.id)
        self.source = getattr(transcript_row, "source", None)
        # Now we have to calculate the CDS
        cds = []
        if transcript_row.coding is True:
            for exon in self.exons:
                if exon[1] >= transcript_row.thick_start and exon[0] <= transcript_row.thick_end:
                    cds.append((int(max(exon[0], transcript_row.thick_start)),
                                int(min(exon[1], transcript_row.thick_end))))
            self.add_exons(cds, features="CDS")
        self.finalize()

    def __initialize_with_bam(self, transcript_row: pysam.AlignedSegment):

        if transcript_row.is_unmapped is True:
            raise InvalidTranscript("I cannot transform BAM unmapped reads.")

        self.chrom = str(transcript_row.reference_name)
        transcript_row.cigar = [(key, val) if key not in (7, 8) else (0, val) for key, val in transcript_row.cigar]

        try:
            start, end = transcript_row.reference_start, transcript_row.get_blocks()[-1][1]
        except IndexError:
            exc = "Invalid start and end for BAM row: {}".format(transcript_row)
            raise InvalidTranscript(exc)

        current = start
        current_exon = None
        self.score = transcript_row.mapq

        r_length = transcript_row.inferred_length  # Read length
        matches = 0
        alen = 0

        # Get the exons
        for cigar, length in transcript_row.cigartuples:
            if cigar in (4, 5):  # Insertion, clipping (soft or hard)
                continue
            elif cigar == 1:
                alen += length
            elif cigar == 3:  # Intron
                if current_exon is None:
                    continue  # Read positioned at the end/beginning of scaffold
                # assert current_exon is not None
                self.add_exon(current_exon)
                current_exon = None
                current += length
            elif cigar in (0, 2):  # Match or deletion
                if current_exon is None:
                    current_exon = (current + 1, current + length)
                else:
                    current_exon = (current_exon[0], current_exon[1] + length)
                current += length
                if cigar == 0:
                    matches += length
                    alen += length
        self.add_exon(current_exon)
        coverage = round(100 * alen / r_length, 2)
        # Get the tags
        tags = dict(transcript_row.tags)
        if "MD" in tags:
            snps = sum(len(_) for _ in re.split("[0-9]*", tags["MD"]) if not _.startswith("^"))
            identity = round(100 * (matches - snps) / r_length, 2)
            self.attributes["identity"] = identity
            del tags["MD"]

        # Set the strand
        if "XS" in tags:
            self.strand = tags["XS"]
            del tags["XS"]
        elif "ts" in tags:
            self.strand = tags["ts"]  # Minimap2
            del tags["ts"]
        else:
            if transcript_row.is_reverse:
                self.strand = "-"
            else:
                self.strand = "+"

        self.attributes.update(tags)
        self.attributes["coverage"] = coverage
        self.attributes["cigar"] = transcript_row.cigarstring
        self.id = transcript_row.query_name
        self.parent = self.attributes["gene_id"] = "{0}.gene".format(transcript_row.query_name)
        self.finalize()

    def __initialize_with_gf(self, transcript_row: (GffLine, GtfLine)):

        self.chrom = transcript_row.chrom

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

        booleans = {"True": True, "False": False, "None": None}

        for key, val in transcript_row.attributes.items():
            if not isinstance(val, Hashable):
                pass
            elif val in booleans:
                val = booleans[val]
            elif isinstance(val, bool):
                pass
            else:
                try:
                    val = int(val)
                except (ValueError, OverflowError):
                    try:
                        val = float(val)
                    except ValueError:
                        pass
                except TypeError:
                    pass
            self.attributes[intern(key)] = val

        self.blast_hits = []
        self.json_conf = None

        if transcript_row.is_transcript is False:
            if transcript_row.is_exon is False and transcript_row.feature not in ("match", "tss", "tts"):
                raise TypeError("Invalid GF line")
            self.__expandable = True

            if "cDNA_match" in transcript_row.feature and isinstance(transcript_row, GffLine):
                self.parent = transcript_row.id
                self.id = transcript_row.id
                self.__expandable = True
            elif transcript_row.feature in ("match",):
                self.parent = transcript_row.parent or transcript_row.id
                if transcript_row.id is not None:
                    self.id = transcript_row.id
                else:
                    self.id = transcript_row.transcript
                self._possibly_without_exons = True
                # if transcript_row.feature == "tts":  # Termination site, especially in Augustus
                #     if self.strand == "-":
                #         self.end = transcript_row.start
                #     else:
                #         self.end = transcript_row.end
                # elif transcript_row.feature == "tss":  # Transcription start site, especially in Augustus
                #     if self.strand == "-":
                #         self.start = transcript_row.end
                #     else:
                #         self.start = transcript_row.start
            elif transcript_row.is_exon is True and isinstance(transcript_row, GffLine) and transcript_row.feature not in ("cDNA_match", "match"):
                self.id = transcript_row.parent[0]
            else:
                self.parent = transcript_row.gene
                self.id = transcript_row.transcript
            if transcript_row.is_exon is True:
                self.add_exon(transcript_row)
        else:
            self.parent = transcript_row.parent
            if transcript_row.id is not None:
                self.id = transcript_row.id
            else:
                raise ValueError(transcript_row)
            self.feature = intern(transcript_row.feature)

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

        return \
            (
                self.start, self.end, self.combined_cds, self.internal_orfs, self.chrom, self.strand, self.exons
            ) == (
                other.start, other.end, other.combined_cds, other.internal_orfs, other.chrom, other.strand, other.exons
            )

    def __len__(self) -> int:
        """Returns the length occupied by the unspliced transcript on the genome."""
        if self.end is not None and self.start is not None:
            return self.end - self.start + 1
        else:
            return 0

    def __lt__(self, other) -> bool:
        """A transcript is lesser than another if it is on a lexicographic inferior chromosome,
        or if it begins before the other, or (in the case where they begin at the same location)
        it ends earlier than the other.
        """

        if not hasattr(other, "chrom") or not hasattr(other, "start") or not hasattr(other, "end"):
            raise TypeError("I cannot compare two different types: {} vs {}".format(
                type(self), type(other)
            ))

        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        if self == other:
            return False
        elif self.start < other.start:
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

        state = dict()
        for key, item in self.__dict__.items():
            if key.endswith("json_conf") or key in ("_Transcript__segmenttree", "_Transcript__cds_tree"):
                continue
            try:
                state[key] = copy.deepcopy(item)
            except TypeError:
                raise TypeError(key, item)

        self.logger = logger

        if hasattr(self, "json_conf") and self.json_conf is not None:
            if "json_conf" not in state:
                state["json_conf"] = dict()
                for key in self.json_conf:
                    if not isinstance(self.json_conf[key], dict):
                        if key == "compiled":
                            continue
                        if isinstance(self.json_conf[key], pysam.FastaFile):
                            state["json_conf"][key] = self.json_conf[key].filename
                        else:
                            state["json_conf"][key] = copy.deepcopy(self.json_conf[key])
                    else:
                        state["json_conf"][key] = dict()
                        for subkey in self.json_conf[key]:
                            if subkey == "compiled":
                                continue
                            else:
                                if isinstance(self.json_conf[key][subkey], pysam.FastaFile):
                                    state["json_conf"][key][subkey] = self.json_conf[key][subkey].filename
                                else:
                                    try:
                                        state["json_conf"][key][subkey] = copy.deepcopy(self.json_conf[key][subkey])
                                    except TypeError:
                                        raise TypeError((key, subkey, self.json_conf[key][subkey]))

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
        self.__json_conf = None
        self.__json_conf = state.pop("json_conf", None)
        self.__dict__.update(state)
        self._calculate_cds_tree()
        self._calculate_segment_tree()
        # Set the logger to NullHandler
        self.logger = None

    # ######## Class instance methods ####################

    def add_exon(self, gffline, feature=None, phase=None):
        """This function will append an exon/CDS feature to the object.
        :param gffline: an annotation line
        :type gffline: (Mikado.parsers.GFF.GffLine | Mikado.parsers.GTF.GtfLine | tuple | list)
        :type feature: flag to indicate what kind of feature we are adding
        """

        if self.finalized is True:
            raise ModificationError("You cannot add exons to a finalized transcript!")

        _gtype = type(gffline)

        if _gtype in (tuple, list):
            assert len(gffline) == 2
            try:
                start, end = sorted(gffline)
            except TypeError:
                raise TypeError((gffline, _gtype))
            if feature is None:
                feature = "exon"
        elif _gtype is intervaltree.Interval:
            start, end = gffline[0], gffline[1]
            if feature is None:
                feature = "exon"
        elif _gtype is Interval:
            start, end = gffline.start, gffline.end
            if feature is None:
                feature = "exon"
        elif _gtype not in (GtfLine, GffLine):
            raise InvalidTranscript("Unkwown feature type! %s", _gtype)
        else:
            start, end = sorted([gffline.start, gffline.end])
            if feature is None:
                feature = gffline.feature
            if _gtype is GffLine and "cdna_match" in gffline.feature.lower():
                gffline.parent = gffline.id

            if self.id not in gffline.parent:
                raise InvalidTranscript(
                    """Mismatch between transcript and exon:
                    {0}
                    {1}
                    {2}""".format(self.id, gffline.parent, gffline))
            assert gffline.is_exon is True, str(gffline)
            assert gffline.chrom == self.chrom, (gffline.chrom, self.chrom)
            phase = gffline.phase

        assert isinstance(start,
                          (int, np.int, np.int32, np.int64)) and isinstance(end,
                                                                            (int, np.int, np.int32, np.int64))
        if feature.upper().endswith("CDS"):
            store = self.combined_cds
            if phase is not None:
                self.phases[(start, end)] = phase

        elif "combined_utr" in feature or "UTR" in feature.upper():
            store = self.combined_utr
        elif feature.endswith("exon") or "match" in feature:
            store = self.exons
        elif feature == "start_codon":
            self.has_start_codon = True
            return
        elif feature == "stop_codon":
            # Here we have to store it, as the stop codon has to be included in the CDS for Mikado.
            self.has_stop_codon = True
            store = self.stop_codon
            if phase is not None:
                self.phases[(start, end)] = phase
            # return
        elif feature == "intron":
            store = self.introns
        else:
            raise InvalidTranscript("Unknown feature: {0}".format(gffline.feature))

        segment = tuple([int(start), int(end)])
        if segment in store:
            return
        if self.__expandable is True:
            self.start = min([self.start, start])
            self.end = max([self.end, end])
        store.append(segment)
        return

    def add_exons(self, exons, features=None, phases=None):

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
        if phases is None:
            phases = [None] * len(exons)
        elif len(phases) != len(exons):
            raise InvalidTranscript("Mismatch between exons and features! %s,\t%s",
                                    exons,
                                    features)

        for exon, feature, phase in zip(exons, features, phases):
            self.add_exon(exon, feature=feature, phase=phase)
        return

    def format(self, format_name,
               with_introns=False,
               with_cds=True,
               all_orfs=False,
               transcriptomic=False):

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
            lines = [create_lines_bed(self, transcriptomic=transcriptomic)]
        else:
            to_gtf = (format_name == "gtf")
            if with_cds is True:
                lines = create_lines_cds(self,
                                         to_gtf=to_gtf,
                                         with_introns=with_introns,
                                         all_orfs=all_orfs,
                                         transcriptomic=transcriptomic
                                         )
            else:
                lines = create_lines_no_cds(self,
                                            to_gtf=to_gtf,
                                            with_introns=with_introns,
                                            transcriptomic=transcriptomic)

        return "\n".join(lines)

    def get_internal_orf_beds(self):
        """This method will return all internal ORFs as BED12 objects"""

        if hasattr(self, "cdna"):
            seq = self.cdna
        else:
            seq = None

        row = BED12(transcriptomic=True, coding=True, start_adjustment=False, max_regression=0,
                    table=self.codon_table)
        row.header = False
        row.chrom = self.id
        row.strand = "+"
        row.start = 1
        row.end = self.cdna_length

        if not self.internal_orfs:
            row.coding = False
            row.block_count = 0
            row.thick_start = 1  # Necessary b/c I am subtracting one
            row.thick_end = 0
            row.name = self.tid
            row.block_count = 0
            row.block_starts = [0]
            row.block_sizes = [0]
            row = BED12(row, seq,
                        coding=False, transcriptomic=True, max_regression=0, start_adjustment=False,
                        table=self.codon_table)
            assert row.invalid is False, ("\n".join([str(row), row.invalid_reason]))
            yield row

        else:
            for index, iorf in enumerate(self.internal_orfs):
                new_row = row.copy()
                cds_len = 0
                cds_start = float("inf")
                cds_end = float("-inf")
                if self.strand == "-":
                    cds_start, cds_end = cds_end, cds_start
                phase = None
                cdses = sorted([(_[1][0], _[1][1], _[2]) for _ in iorf if _[0] == "CDS"],
                               key=operator.itemgetter(0, 1))
                cds_len = sum(_[1] - _[0] + 1 for _ in cdses)
                if self.strand == "-":
                    cds_start = cdses[-1][1]
                    cds_end = cdses[0][0]
                    phase = cdses[-1][2]
                else:
                    cds_start = cdses[0][0]
                    cds_end = cdses[-1][1]
                    phase = cdses[0][2]
                self.logger.debug("Phase: %s. CDSes: %s", phase, cdses)

                # Now convert to transcriptomic coordinates
                if self.strand == "-":
                    utr = sum(_[1][1] - _[1][0] + 1 for _ in iorf if _[0] == "UTR" and _[1][0] > cds_start)
                else:
                    utr = sum(_[1][1] - _[1][0] + 1 for _ in iorf if _[0] == "UTR" and _[1][1] < cds_start)
                new_row.thick_start = utr + 1
                new_row.thick_end = new_row.thick_start + cds_len - 1
                new_row.name = "{}_orf{}".format(self.tid, index)
                new_row.block_starts = [row.thick_start]
                new_row.block_sizes = [cds_len]
                new_row.phase = phase
                # self.logger.debug(new_row)
                new_row = BED12(new_row,
                                sequence=seq,
                                phase=phase,
                                coding=True, transcriptomic=True, max_regression=0, start_adjustment=False,
                                table=self.codon_table)
                if (cds_len - phase) % 3 != 0 and cds_end not in (self.start, self.end):
                    raise AssertionError("Invalid CDS length for {}:\n{}\n{}".format(self.id,
                                                                                     iorf,
                                                                                     new_row))

                if new_row.invalid is True:
                    exc = "Invalid ORF for {id} ({new_row.thick_start}, {new_row.thick_end}, phase {new_row.phase})"
                    exc += ", reason: {reason}"
                    exc = exc.format(id=self.id, reason=new_row.invalid_reason, **locals())
                    raise AssertionError(exc)

                yield new_row

    @Metric
    def is_reference(self):
        """Checks whether the transcript has been marked as reference by Mikado prepare"""

        if self.__is_reference is None and self.json_conf is None:
            return False
        elif self.__is_reference is None:
            self.__is_reference = (self.original_source in self.json_conf["prepare"]["files"]["reference"])

        return self.__is_reference

    @is_reference.setter
    def is_reference(self, value):
        if isinstance(value, str):
            if value.lower() == "true":
                value = True
            elif value.lower() == "none":
                value = None
            elif value.lower() == "false":
                value = False
            else:
                pass
        if not isinstance(value, bool) and value is not None:
            if value == 1:
                value = True
            elif value == 0:
                value = False
            else:
                raise ValueError("Invalid value: {} (type: {})".format(value, type(value)))
        self.__is_reference = value

    @is_reference.deleter
    def is_reference(self):
        self.__is_reference = None

    is_reference.category = "External"
    is_reference.rtype = "bool"

    @property
    def frames(self):
        """This property will return a dictionary with three keys - the three possible frames, 0, 1 and 2 - and within
        each, a set of the positions that are in that frame. If the transcript does not have """
        self.finalize()
        frames = {0: set(), 1: set(), 2: set()}
        for orf in self.internal_orfs:
            if self.strand == "-":
                exons = sorted([(_[1][0], _[1][1], _[2]) for _ in orf
                                if _[0] == "CDS"], key=operator.itemgetter(0, 1), reverse=True)
                for start, end, phase in exons:
                    frame = ((3 - phase) % 3 - 1) % 3
                    for pos in range(end, start - 1, -1):
                        frame = abs((frame + 1) % 3)
                        frames[frame].add(pos)
            else:
                exons = sorted([(_[1][0], _[1][1], _[2]) for _ in orf
                      if _[0] == "CDS"], key=operator.itemgetter(0, 1), reverse=False)
                for start, end, phase in exons:
                    frame = ((3 - phase) % 3 - 1) % 3  # Previous frame before beginning of the feature
                    for pos in range(start, end + 1):
                        frame = abs((frame + 1) % 3)
                        frames[frame].add(pos)
        return frames

    @property
    def framed_codons(self):

        codons = list(zip(*[sorted(self.frames[0]), sorted(self.frames[1]), sorted(self.frames[2])]))
        if self.strand == "-":
            codons = list(reversed(codons))

        return codons

    @property
    def _selected_orf_transcript(self):

        """This method will return the selected internal ORF as a transcript object."""

        self.finalize()
        if not self.is_coding:
            return []
        return self._internal_orfs_transcripts[self.selected_internal_orf_index]

    @property
    def _internal_orfs_transcripts(self):
        """This method will return all internal ORFs as transcript objects.
        Note: this will exclude the UTR part, even when the transcript only has one ORF."""

        self.finalize()
        if not self.is_coding:
            return []
        elif len(self.__internal_orf_transcripts) == len(self.internal_orfs):
            return self.__internal_orf_transcripts
        else:
            for num, orf in enumerate(self.internal_orfs, start=1):
                torf = Transcript()
                torf.chrom, torf.strand = self.chrom, self.strand
                torf.derives_from = self.id
                torf.id = "{}.orf{}".format(self.id, num)
                __exons, __phases = [], []
                for segment in [_ for _ in orf if _[0] == "CDS"]:
                    __exons.append(segment[1])
                    __phases.append(segment[2])
                torf.add_exons(__exons, features="exon", phases=None)
                torf.add_exons(__exons, features="CDS", phases=__phases)
                torf.finalize()
                self.__internal_orf_transcripts.append(torf)

        return self.__internal_orf_transcripts

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
            if ("exon", exon) in self.segments:
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

    _orf_pattern = re.compile(r"\.orf[0-9]+$", re.IGNORECASE)

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
        elif self._orf_pattern.search(self.id):
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
        self.stop_codon = []
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
        self.__internal_orf_transcripts = []
        self.combined_utr = []
        self.__cdna_length = None
        self.finalized = False

    def reverse_strand(self):
        """Method to reverse the strand"""
        if self.strand == "+":
            self.strand = "-"
        elif self.strand == "-":
            self.strand = "+"
        elif self.strand is None:
            pass
        self.strip_cds()
        self.logger.warning("Transcript %s has been assigned to the wrong strand, reversing it.",
                            self.id)
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

    def as_dict(self, remove_attributes=False):

        """
        Method to transform the transcript object into a JSON-friendly representation.
        :return:
        """

        state = dict()
        self.finalize()
        mmetrics = self.get_modifiable_metrics()
        if not mmetrics:
            raise ValueError((self.id, mmetrics))

        state["external"] = dict((key, value) for key, value in self.external_scores.items())

        for key in ["chrom", "source", "start", "end", "strand", "score", "attributes"]:

            state[key] = getattr(self, key)
            if key == "attributes" and remove_attributes is True:
                for subkey in [_ for _ in state[key] if _ not in ("ID", "Parent", "Name")]:
                    del state[key][subkey]

        state["exons"] = []
        for exon in self.exons:
            state["exons"].append([exon[0], exon[1]])
        state["introns"] = list((intron[0], intron[1]) for intron in self.introns)

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
        state["finalized"] = getattr(self, "finalized")
        state["splices"] = list(self.splices.copy())
        state["combined_utr"] = list(self.combined_utr)
        state["combined_cds"] = list(self.combined_cds)
        # Now let'add the things that we calculate with the basic lengths
        state["calculated"] = self.__get_calculated_stats()

        for metric in mmetrics:
            state[metric] = getattr(self, metric)

        return state

    def load_dict(self, state, trust_orf=False, accept_undefined_multi=False):

        for key in ["chrom", "source",
                    "start", "end", "strand", "score",
                    "parent", "id"]:
            if key not in state:
                raise CorruptIndex("Key not found for {}: {} ({})".format(self.id, key, type(state)))
            if isinstance(state[key], str):
                state[key] = intern(state[key])
            setattr(self, key, state[key])

        self.external_scores.update(state.get("external", dict()))
        self._original_source = self.source
        self.attributes = state["attributes"].copy()

        # for subkey in self.attributes:
        #     if self.attributes[subkey] == maxsize:
        #         self.attributes[subkey] = float("inf")
        #     elif self.attributes[subkey] == maxsize:
        #         self.attributes[subkey] = float("-inf")

        self.exons = []
        self.combined_cds = []
        for exon in state["exons"]:
            if len(exon) != 2:
                raise CorruptIndex("Invalid exonic values of {}".format(self.id))
            self.exons.append(tuple(exon))
        for intron in state["introns"]:
            if len(intron) != 2:
                raise CorruptIndex("Invalid intronic values of {}".format(self.id))
            self.introns.add(tuple(intron))
        self.splices = set(state["splices"])

        self._trust_orf = trust_orf
        self._accept_undefined_multi = accept_undefined_multi
        self.internal_orfs = []

        self.logger.debug("Starting to load the ORFs for %s", self.id)
        try:
            indices = dict((int(_), _) for _ in state["orfs"])
            __phases = dict()
            for index in sorted(indices.keys()):
                orf = state["orfs"][indices[index]]
                neworf = []

                for segment in orf:
                    if len(segment) == 3:
                        assert segment[0] == "CDS"

                        new_segment = (segment[0],
                                       tuple(segment[1]),
                                       int(segment[2]))
                        self.combined_cds.append(new_segment[1])
                        if index == 0:
                            __phases[new_segment[1]] = new_segment[2]
                    else:
                        assert segment[0] != "CDS"
                        new_segment = (segment[0],
                                       tuple(segment[1]))
                    neworf.append(new_segment)

                self.internal_orfs.append(neworf)
                self.phases = __phases
            self.logger.debug("ORFs: %s", " ".join(str(_) for _ in self.internal_orfs))
            self.selected_internal_orf_index = state["selected_orf"]
        except (ValueError, IndexError):
            raise CorruptIndex("Invalid values for ORFs of {}".format(self.id))

        calculated = state["calculated"]
        del state["calculated"]
        self.__load_calculated_stats(calculated)
        for metric in self.get_modifiable_metrics():
            try:
                setattr(self, metric, state[metric])
            except KeyError:
                raise KeyError((self.id, metric))
        self.finalized = state["finalized"]
        if self.finalized:
            self.combined_utr = sorted([tuple(combi) for combi in state["combined_utr"]])
        else:
            self.combined_utr = []
        self.finalize()

    def __get_calculated_stats(self):
        d = dict()
        self.finalize()
        d["exon_num"] = self.exon_num
        d["three_utr"] = self.three_utr
        d["five_utr"] = self.five_utr
        d["combined_cds_length"] = self.combined_cds_length
        d["combined_cds_fraction"] = self.combined_cds_fraction
        d["selected_cds_num"] = self.selected_cds_num
        d["selected_cds_number_fraction"] = self.selected_cds_number_fraction
        d["three_utr_length"] =self.three_utr_length
        d["five_utr_length"] = self.five_utr_length
        d["three_utr_length"] = self.three_utr_length
        d["five_utr_length"] = self.five_utr_length
        d["utr_length"] = self.utr_length
        d["utr_fraction"] = self.utr_fraction
        d["five_utr_num"] = self.five_utr_num
        d["five_utr_num_complete"] = self.five_utr_num_complete
        d["three_utr_num"] = self.three_utr_num
        d["three_utr_num_complete"] = self.three_utr_num_complete
        d["utr_num_complete"] = self.utr_num_complete
        d["utr_num"] = self.utr_num
        d["highest_cds_exons_num"] = self.highest_cds_exons_num
        d["min_intron_length"] = self.min_intron_length
        d["max_intron_length"] = self.max_intron_length
        d["highest_cds_exon_number"] = self.highest_cds_exon_number
        d["selected_cds_exons_fraction"] = self.selected_cds_exons_fraction
        d["cdna_length"] = self.cdna_length
        d["cds_not_maximal"] = self.cds_not_maximal
        d["cds_not_maximal_fraction"] = self.cds_not_maximal_fraction
        d["end_distance_from_tes"] = self.end_distance_from_tes
        d["end_distance_from_junction"] = self.end_distance_from_junction
        d["selected_end_distance_from_junction"] = self.selected_end_distance_from_junction
        d["selected_end_distance_from_tes"] = self.selected_end_distance_from_tes
        d["selected_start_distance_from_tss"] = self.selected_start_distance_from_tss
        d["start_distance_from_tss"] = self.start_distance_from_tss

        return d

    def __load_calculated_stats(self, state):
        self.__exon_num = state["exon_num"]
        self.__three_utr = state["three_utr"]
        self.__five_utr = state["five_utr"]
        self.__combined_cds_length = state["combined_cds_length"]
        self.__combined_cds_fraction = state["combined_cds_fraction"]
        self.__selected_cds_num = state["selected_cds_num"]
        self.__selected_cds_number_fraction = state["selected_cds_number_fraction"]
        self.__three_utr_length = state["three_utr_length"]
        self.__five_utr_length = state["five_utr_length"]
        self.__three_utr_length = state["three_utr_length"]
        self.__five_utr_length = state["five_utr_length"]
        self.__utr_length = state["utr_length"]
        self.__utr_fraction = state["utr_fraction"]
        self.__five_utr_num = state["five_utr_num"]
        self.__five_utr_num_complete = state["five_utr_num_complete"]
        self.__three_utr_num = state["three_utr_num"]
        self.__three_utr_num_complete = state["three_utr_num_complete"]
        self.__utr_num_complete = state["utr_num_complete"]
        self.__utr_num = state["utr_num"]
        self.__highest_cds_exons_num = state["highest_cds_exons_num"]
        self.__min_intron_length = state["min_intron_length"]
        self.__max_intron_length = state["max_intron_length"]
        self.__highest_cds_exon_number = state["highest_cds_exon_number"]
        self.__selected_cds_exons_fraction = state["selected_cds_exons_fraction"]
        self.__cdna_length = state["cdna_length"]
        self.__cds_not_maximal = state["cds_not_maximal"]
        self.__cds_not_maximal_fraction = state["cds_not_maximal_fraction"]
        self.__end_distance_from_tes = state["end_distance_from_tes"]
        self.__end_distance_from_junction = state["end_distance_from_junction"]
        self.__selected_end_distance_from_junction = state["selected_end_distance_from_junction"]
        self.__selected_end_distance_from_tes = state["selected_end_distance_from_tes"]
        self.__selected_start_distance_from_tss = state["selected_start_distance_from_tss"]
        self.__start_distance_from_tss = state["start_distance_from_tss"]

    def add_derived_child(self, name):

        """Method to add a derived child (eg TAIR protein) to a transcript"""

        if not isinstance(name, (str, bytes)):
            raise ValueError("Invalid child type: {} (value: {})".format(
                type(name), name))

        self.__derived_children.add(name)

    def find(self, *args, **kwargs):

        """Alias for search"""

        self.search(*args, **kwargs)

    def search(self, start, end, value=None, max_distance=0):

        """Method to use the segment tree to find all segments intersecting the given segment."""
        return self.segmenttree.search(start - 1, end, max_distance=max_distance, value=value)

    def find_upstream(self, start, end, offset=0):

        """Method to use the segment tree to find all segments in the transcript upstream of a given interval."""
        return self.segmenttree.upstream_of_interval(Interval(start - offset, end),
                                                     num_intervals=self.exon_num + len(self.introns),
                                                     max_dist=10**9)

    def find_downstream(self, start, end, offset=0):

        """Method to use the segment tree to find all segments in the transcript downstream of a given interval."""

        return self.segmenttree.downstream_of_interval(Interval(start - offset, end),
                                                       num_intervals=self.exon_num + len(self.introns),
                                                       max_dist=10**9)


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

        _metrics = sorted([metric for metric in metrics])
        final_metrics = ["tid", "alias", "parent", "original_source", "score"] + _metrics
        return final_metrics

    @classmethod
    @functools.lru_cache(maxsize=None, typed=True)
    def get_modifiable_metrics(cls) -> list:

        metrics = []
        for member in inspect.getmembers(cls):
            not_private = (not member[0].startswith("_" + cls.__name__ + "__") and not member[0].startswith("__"))
            in_dict = (member[0] in cls.__dict__)
            if in_dict:
                is_metric = isinstance(cls.__dict__[member[0]], Metric)
                has_fset = (getattr(cls.__dict__[member[0]], "fset", None) is not None)
            else:
                is_metric = None
                has_fset = None
            if all([not_private, in_dict, is_metric, has_fset]):
                metrics.append(member[0])
        return metrics

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
                logger = create_null_logger()
                self.__logger = logger
            else:
                pass
        else:
            assert isinstance(logger, logging.Logger), type(logger)
            self.__logger = logger
        # self.__logger.propagate = False

    @logger.deleter
    def logger(self):
        self.__logger = create_null_logger()

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
    def start(self):
        return self.__start
    
    @start.setter
    def start(self, start: [int,None]):
        if self.finalized is True:
            raise InvalidTranscript("I cannot modify the start of a finalised transcript!")
        if start is not None and not isinstance(start, (int)):
            raise TypeError("Start values can only be integers or None!")
        self.__start = start
    
    @property
    def end(self):
        return self.__end

    @end.setter
    def end(self, end: [int, None]):
        if self.finalized is True:
            raise InvalidTranscript("I cannot modify the end of a finalised transcript!")
        if end is not None and not isinstance(end, (int)):
            raise TypeError("End values can only be integers or None!")
        self.__end = end

    @property
    def chrom(self):
        return self.__chrom

    @chrom.setter
    def chrom(self, chrom: [None, str]):
        if self.finalized is True:
            raise InvalidTranscript("I cannot modify the chromosome of a finalised transcript!")
        if isinstance(chrom, bytes):
            chrom = intern(chrom.decode())
        elif chrom is not None:
            chrom = intern(str(chrom))
        self.__chrom = chrom

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

        assert isinstance(phases, dict), (phases, type(phases))
        self.__phases = phases

    #
    # pylint: disable=invalid-name
    @property
    def id(self):
        """ID of the transcript - cannot be an undefined value."""
        return self.__id

    @id.setter
    def id(self, newid):

        if not isinstance(newid, str):
            raise ValueError("Invalid value for id: {0}, type {1}".format(newid, type(newid)))
        self.__id = intern(newid)
    # pylint: enable=invalid-name

    @property
    def tid(self):
        """ID of the transcript - cannot be an undefined value. Alias of id."""
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

        _ptype = type(parent)

        if _ptype in (list, type(None)):
            self.__parent = parent
        elif _ptype is str:
            if "," in parent:
                self.__parent = parent.split(",")
            else:
                self.__parent = [parent]
        else:
            raise ValueError("Invalid value for parent: {0}, type {1}".format(
                parent, type(parent)))

        self.attributes["gene_id"] = self.__parent
        if self.__parent:
            self.__parent = [intern(str(_)) for _ in self.__parent]

    @property
    def source(self):
        """Source of the transcript. Equivalent to the second field in the GFF/GTF files."""
        return self.__source

    @source.setter
    def source(self, source):
        if source is not None and not isinstance(source, str):
            raise TypeError("Source values must be strings or None!")
        self.__source = source

    @property
    def original_source(self):
        """This property returns the original source assigned to the transcript (before Mikado assigns its own final
        source value)."""
        return self._original_source

    @property
    def gene(self):

        if "gene_id" not in self.attributes:
            self.attributes["gene_id"] = self.parent[0]

        return self.attributes["gene_id"]

    @property
    def location(self):
        """Web-apollo compatible string for the location of the transcript."""

        return "{}:{}..{}".format(self.chrom,
                                  self.start,
                                  self.end)

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
        if isinstance(name, (int, float)):
            name = str(name)
        elif isinstance(name, (bytes,)):
            name = name.decode()

        if not isinstance(name, (type(None), str)):
            raise ValueError("Invalid name: {0}".format(name))

        if name is not None:
            name = intern(name)
        self.attributes["Name"] = name

    @property
    def alias(self):
        """This property returns the alias of the transcript, if present, else its ID"""
        if "alias" in self.attributes.keys():
            if "," not in self.attributes["alias"]:
                return self.attributes["alias"]
            else:
                return self.attributes["alias"].split(",")
        elif "Alias" in self.attributes.keys():
            if "," not in self.attributes["Alias"]:
                return self.attributes["Alias"]
            else:
                return self.attributes["Alias"].split(",")
        else:
            return self.id

    @alias.setter
    def alias(self, name):
        if not isinstance(name, (type(None), str)):
            raise ValueError("Invalid name: {0}".format(name))
        self.attributes["alias"] = name

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
        index = self.selected_internal_orf_index
        if index is not None:
            return self.internal_orfs[index]
        else:
            return tuple([])

        # if self.combined_cds_length == 0:
        #     # self.__max_internal_orf_length = 0
        #     # self.selected_internal_orf_index = None
        #     return tuple([])
        # else:
        #     assert self.selected_internal_orf_index is not None
        #     return self.internal_orfs[self.selected_internal_orf_index]

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
        return self.__five_utr

    @property
    def three_utr(self):
        """Returns the exons in the 3' UTR of the selected ORF.
        If the end codon is absent, no UTR is given."""
        return self.__three_utr

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
            self.__max_internal_orf_length = 0
            self.__selected_cds = []
            return
        if not isinstance(index, int):
            raise TypeError()
        if index < 0 or index >= len(self.internal_orfs):
            raise IndexError("No ORF corresponding to this index: {0}".format(index))
        self.__max_internal_orf_index = index
        selected_cds = [segment[1] for segment in self.selected_internal_orf if
                               segment[0] == "CDS"]
        self.__selected_cds = selected_cds
        ar = np.array(list(zip(*[segment[1] for segment in self.selected_internal_orf if
                               segment[0] == "CDS"])))
        self.__max_internal_orf_length = int(np.subtract(ar[1], ar[0] - 1).sum())

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
    def finalized(self):
        return self.__finalized

    @finalized.setter
    def finalized(self, finalized):
        if finalized not in (True, False):
            raise ValueError("This is a boolean flag")
        self.__finalized = finalized

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

        if self.finalized is True:
            raise NotImplementedError("I cannot reset the exons in a finalised transcript.")

        if not isinstance(args[0], (set, list)):
            raise TypeError(type(args[0]))

        for exon in args[0]:
            if (not isinstance(exon, tuple) or len(exon) != 2 or
                    not (isinstance(exon[0], int) and isinstance(exon[1], int))):
                raise TypeError("Exons must be tuples of two integers!")

        self.__exons = list(sorted(args[0]))

    def _set_exons(self, exons):
        """Private method that bypasses the checks within the direct setter, for speed purposes in finalising."""
        if self.finalized is True:
            raise NotImplementedError("I cannot reset the exons in a finalised transcript.")

        self.__exons = exons

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

        if len(combined) > 0:
            self.__combined_cds_length = sum([_[1] - _[0] + 1 for _ in combined])
        else:
            self.__combined_cds_length = 0

        self.__combined_cds = combined

    @property
    def coding_exons(self):

        """This property returns a list of exons which have at least a coding part."""

        @functools.lru_cache(maxsize=None, typed=True)
        def _wrapped(combined_cds, combined_cds_length, exons):
            if combined_cds_length == 0:
                return []
            coding = []
            for exon in exons:
                for ccd in combined_cds:
                    if self.overlap(exon, ccd) > 0:
                        coding.append(exon)
                        break
                    continue
            return coding

        return _wrapped(tuple(self.combined_cds), self.combined_cds_length, tuple(self.exons))

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

        # if len(self.combined_cds) == 0:
        #     self.__selected_cds = []
        # else:
        #     self.__selected_cds = [segment[1] for segment in self.selected_internal_orf if
        #                            segment[0] == "CDS"]

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
        """
        if len(self.__segmenttree) != len(self.combined_cds) + len(self.combined_cds_introns):
            self._calculate_segment_tree()

        return self.__cds_tree

    def _calculate_cds_tree(self):

        """
        :rtype: intervaltree.IntervalTree
        """

        self.__cds_tree = IntervalTree()

        for exon in self.combined_cds:
            self.__cds_tree.add(exon[0], exon[1], value=Interval(exon[0], exon[1], value="CDS"))

        for intron in self.combined_cds_introns:
            self.__cds_tree.add(intron[0], intron[1], value=Interval(intron[0], intron[1], value="intron"))

        return

    @property
    def codon_table(self):
        """This property returns the codon table for the project. Default: 0 (Standard, but only ATG is considered
        a valid start codon)."""

        if self.json_conf is None:
            return 0
        return self.json_conf.get("serialise", {}).get("codon_table", 0)

    @property
    def segmenttree(self):

        """
        :rtype: intervaltree.IntervalTree
        """

        if len(self.__segmenttree) != self.exon_num + len(self.introns):
            self._calculate_segment_tree()

        return self.__segmenttree

    def _calculate_segment_tree(self):

        self.__segmenttree = IntervalTree()
        for exon in self.exons:
            try:
                self.__segmenttree.add(exon[0], exon[1], value="exon")
            except AssertionError as exc:
                raise AssertionError(f"Exon for {self.id} invalid: {exon}\n{exc}")

        for intron in self.introns:
            try:
                self.__segmenttree.add(intron[0], intron[1], value="intron")
            except AssertionError as exc:
                raise AssertionError(f"Intron for {self.id} invalid: {intron}\n{exc}")

    @property
    def derived_children(self):
        """
        This property stores the names of children features (eg TAIR proteins)
        that this transcript gives origin to.
        :return:
        """

        return self.__derived_children

    @property
    def external_scores(self):

        """
        **SPECIAL** this Namespace contains all the information regarding external scores
        for the transcript. If an absent property is not defined in the Namespace,
        Mikado will set a default value of 0 into the Namespace and return it.
        """

        return self.__external_scores

    @property
    def external(self):

        """
        **SPECIAL** this Namespace contains all the information regarding external scores
        for the transcript. If an absent property is not defined in the Namespace,
        Mikado will set a default value of 0 into the Namespace and return it.
        """
        return self.__external_scores

    # ################### Class metrics ##################################

    # Disable normal checks on names and hidden methods, as
    # checkers get confused by the Metric method
    # pylint: disable=method-hidden,invalid-name

    @Metric
    def source_score(self):

        """This metric returns a score that is assigned to the transcript
        in virtue of its origin."""

        if self.json_conf is not None:
            return self.json_conf.get("prepare", {}).get("files", {}).get(
                "source_score", {}).get(self.original_source, 0)
        else:
            return 0

    source_score.category = "External"
    source_score.rtype = "float"

    @Metric
    def combined_cds_length(self):
        """This property return the length of the CDS part of the transcript."""

        # if self.__combined_cds_length == 0 and len(self.combined_cds) > 0:
        #     ar = np.array(list(*zip(self.combined_cds)))
        #     c_length = int(np.subtract(ar[1], ar[0] - 1).sum())
        #     assert c_length > 0
        return self.__combined_cds_length

    combined_cds_length.category = "CDS"
    combined_cds_length.rtype = "int"

    @Metric
    def combined_cds_num(self):
        """This property returns the number of non-overlapping CDS segments
        in the transcript."""
        return len(self.combined_cds)

    combined_cds_num.category = "CDS"
    combined_cds_num.rtype = "int"

    @Metric
    def combined_cds_num_fraction(self):
        """This property returns the fraction of non-overlapping CDS segments
        in the transcript
        vs. the total number of exons"""
        return len(self.combined_cds) / len(self.exons)

    combined_cds_num_fraction.category = "CDS"
    combined_cds_num_fraction.rtype = "float"
    combined_cds_num_fraction.usable_raw = True

    @Metric
    def combined_cds_fraction(self):
        """This property return the percentage of the CDS part of the transcript
        vs. the cDNA length"""
        return self.__combined_cds_fraction

    combined_cds_fraction.category = "CDS"
    combined_cds_fraction.rtype = "float"
    combined_cds_fraction.usable_raw = True

    @Metric
    def combined_utr_length(self):
        """This property return the length of the UTR part of the transcript."""
        return sum([e[1] - e[0] + 1 for e in self.combined_utr])

    combined_utr_length.category = "UTR"
    combined_utr_length.rtype = "int"

    @Metric
    def combined_utr_fraction(self):
        """This property returns the fraction of the cDNA which is not coding according
        to any ORF. Complement of combined_cds_fraction"""
        return 1 - self.combined_cds_fraction

    combined_utr_fraction.category = "UTR"
    combined_utr_fraction.usable_raw = True
    combined_utr_fraction.rtype = "float"

    def __calculate_cdna_length(self):
        self.__cdna_length = int(sum([_[1] - _[0] + 1 for _ in self.exons]))

    @Metric
    def cdna_length(self):
        """This property returns the length of the transcript."""
        if self.__cdna_length is None and self.finalized is True:
            raise AssertionError
        if self.finalized is False or self.__cdna_length is None:
            self.__calculate_cdna_length()

        return self.__cdna_length

    cdna_length.category = "cDNA"
    cdna_length.rtype = "int"

    @Metric
    def number_internal_orfs(self):
        """This property returns the number of ORFs inside a transcript."""
        return len(self.internal_orfs)

    number_internal_orfs.category = "CDS"
    number_internal_orfs.rtype = "int"

    @Metric
    def selected_cds_length(self):
        """This property calculates the length of the CDS selected as best inside
        the cDNA."""

        return self.__max_internal_orf_length

    selected_cds_length.category = "CDS"
    selected_cds_length.rtype = "int"

    @Metric
    def selected_cds_num(self):
        """This property calculates the number of CDS exons for the selected ORF"""
        return self.__selected_cds_num

    selected_cds_num.category = "CDS"
    selected_cds_num.rtype = "int"

    @Metric
    def selected_cds_fraction(self):
        """This property calculates the fraction of the selected CDS vs. the cDNA length."""
        return self.selected_cds_length / self.cdna_length

    selected_cds_fraction.category = "CDS"
    selected_cds_fraction.usable_raw = True
    selected_cds_fraction.rtype = "float"

    @Metric
    def highest_cds_exons_num(self):
        """Returns the number of CDS segments in the selected ORF
        (irrespective of the number of exons involved)"""
        return self.__highest_cds_exons_num
        # return len(list(filter(lambda x: x[0] == "CDS", self.selected_internal_orf)))

    highest_cds_exons_num.category = "CDS"
    highest_cds_exons_num.rtype = "int"

    @Metric
    def selected_cds_exons_fraction(self):
        """Returns the fraction of CDS segments in the selected ORF
        (irrespective of the number of exons involved)"""
        return self.__selected_cds_exons_fraction

    selected_cds_exons_fraction.category = "CDS"
    selected_cds_exons_fraction.usable_raw = True
    selected_cds_exons_fraction.rtype = "float"

    def __calculate_selected_cds_exons_fraction(self):
        self.__selected_cds_exons_fraction = self.highest_cds_exon_number / self.exon_num

    @Metric
    def highest_cds_exon_number(self):
        """This property returns the maximum number of CDS segments
        among the ORFs; this number can refer to an ORF *DIFFERENT*
        from the maximal ORF."""
        return self.__highest_cds_exon_number

    highest_cds_exon_number.category = "CDS"
    highest_cds_exon_number.rtype = "int"

    def __calculate_highest_cds_exon_number(self):
        if len(self.internal_orfs) == 0:
            self.__highest_cds_exon_number = 0
        else:
            cds_numbers = []
            for cds in self.internal_orfs:
                cds_numbers.append(sum(1 for segment in cds if segment[0] == "CDS"))
                # len(list(filter(lambda x: x[0] == "CDS", cds))))
            self.__highest_cds_exon_number = max(cds_numbers)

    @Metric
    def selected_cds_number_fraction(self):
        """This property returns the proportion of best possible CDS segments
        vs. the number of exons. See selected_cds_number."""
        return self.__selected_cds_number_fraction

    selected_cds_number_fraction.category = "CDS"
    selected_cds_number_fraction.rtype = "float"

    @Metric
    def cds_not_maximal(self):
        """This property returns the length of the CDS excluded from the selected ORF."""
        return self.__cds_not_maximal

    def _set_cds_not_maximal(self):
        if len(self.internal_orfs) < 2:
            self.__cds_not_maximal = 0
        else:
            self.__cds_not_maximal = self.combined_cds_length - self.selected_cds_length

    cds_not_maximal.category = "CDS"
    cds_not_maximal.rtype = "int"

    @Metric
    def cds_not_maximal_fraction(self):
        """This property returns the fraction of bases not in the selected ORF compared to
        the total number of CDS bases in the cDNA."""
        return self.__cds_not_maximal_fraction

    def _set_cds_not_maximal_fraction(self):
        self._set_cds_not_maximal()
        if self.combined_cds_length == 0:
            self.__cds_not_maximal_fraction = 0
        else:
            self.__cds_not_maximal_fraction = self.cds_not_maximal / self.combined_cds_length

    cds_not_maximal_fraction.category = "CDS"
    cds_not_maximal_fraction.usable_raw = True
    cds_not_maximal_fraction.rtype = "float"

    @Metric
    def five_utr_length(self):
        """Returns the length of the 5' UTR of the selected ORF."""
        return self.__five_utr_length

    five_utr_length.category = "UTR"

    @Metric
    def five_utr_num(self):
        """This property returns the number of 5' UTR segments for the selected ORF."""
        return self.__five_utr_num

    five_utr_num.category = "UTR"
    five_utr_num.rtype = "int"

    @Metric
    def five_utr_num_complete(self):
        """This property returns the number of 5' UTR segments for the selected ORF,
        considering only those which are complete exons."""
        return self.__five_utr_num_complete

    five_utr_num_complete.category = "UTR"
    five_utr_num_complete.rtype = "int"

    @Metric
    def three_utr_length(self):
        """Returns the length of the 5' UTR of the selected ORF."""
        return self.__three_utr_length

    three_utr_length.__category = "UTR"
    three_utr_length.rtype = "int"

    @Metric
    def three_utr_num(self):
        """This property returns the number of 3' UTR segments
        (referred to the selected ORF)."""
        return self.__three_utr_num

    three_utr_num.category = "UTR"
    three_utr_num.rtype = "int"

    @Metric
    def three_utr_num_complete(self):
        """This property returns the number of 3' UTR segments for the selected ORF,
        considering only those which are complete exons."""
        return self.__three_utr_num_complete

    three_utr_num_complete.category = "UTR"
    three_utr_num_complete.rtype = "int"

    @Metric
    def utr_num(self):
        """Returns the number of UTR segments (referred to the selected ORF)."""
        return self.__utr_num

    utr_num.category = "UTR"
    utr_num.rtype = "int"

    @Metric
    def utr_num_complete(self):
        """Returns the number of UTR segments which are
        complete exons (referred to the selected ORF)."""
        return self.__utr_num_complete

    utr_num_complete.category = "UTR"
    utr_num_complete.rtype = "int"

    @Metric
    def utr_fraction(self):
        """This property calculates the length of the UTR
        of the selected ORF vs. the cDNA length."""
        return self.__utr_fraction

    utr_fraction.category = "UTR"
    utr_fraction.usable_raw = True
    utr_fraction.rtype = "float"

    @Metric
    def utr_length(self):
        """Returns the sum of the 5'+3' UTR lengths"""
        return self.__utr_length

    utr_length.category = "UTR"
    utr_length.rtype = "int"

    def _set_basic_lengths(self):

        self.__exon_num = len(self.exons)
        if len(self.combined_cds) == 0:
            self.__three_utr = []
            self.__five_utr = []
        if self.strand == "-":
            self.__three_utr = list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                                    utr_segment[0] == "UTR" and utr_segment[1][1] < self.selected_cds_end)
            self.__five_utr = list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                                   utr_segment[0] == "UTR" and utr_segment[1][0] > self.selected_cds_start)
        else:
            self.__three_utr = list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                                    utr_segment[0] == "UTR" and utr_segment[1][0] > self.selected_cds_end)
            self.__five_utr = list(utr_segment[1] for utr_segment in self.selected_internal_orf if
                                   utr_segment[0] == "UTR" and utr_segment[1][1] < self.selected_cds_start)

        self.__combined_cds_length = sum([e[1] - e[0] + 1 for e in self.combined_cds])
        self.__combined_cds_fraction = self.combined_cds_length / self.cdna_length
        self.__selected_cds_num = sum(1 for exon in self.selected_internal_orf if exon[0] == "CDS")
        self.__selected_cds_number_fraction = self.selected_cds_num / self.exon_num
        self.__three_utr_length = sum(x[1] - x[0] + 1 for x in self.three_utr)
        self.__five_utr_length = sum(x[1] - x[0] + 1 for x in self.five_utr)
        self.__utr_length = self.__three_utr_length + self.__five_utr_length
        self.__utr_fraction = 1 - self.selected_cds_fraction
        self.__five_utr_num = len(self.five_utr)
        self.__five_utr_num_complete = sum(1 for utr in self.five_utr if utr in self.exons)
        self.__three_utr_num = len(self.three_utr)
        self.__three_utr_num_complete = sum(1 for utr in self.three_utr if utr in self.exons)
        self.__utr_num_complete = self.three_utr_num_complete + self.five_utr_num_complete
        self.__utr_num = self.three_utr_num + self.five_utr_num
        self.__highest_cds_exons_num = sum(1 for _ in self.selected_internal_orf if _[0] == "CDS")
        self.__calculate_min_intron_length()
        self.__calculate_max_intron_length()
        self.__calculate_highest_cds_exon_number()
        self.__calculate_selected_cds_exons_fraction()
        self._set_cds_not_maximal()
        self._set_cds_not_maximal_fraction()

    @Metric
    def has_start_codon(self):
        """Boolean. True if the selected ORF has a start codon."""
        return self.__has_start_codon and (self.combined_cds_length > 0)

    @has_start_codon.setter
    def has_start_codon(self, value):
        """Setter. Checks that the argument is boolean."""

        if not isinstance(value, bool) and value is not None:
            if value == 1:
                value = True
            elif value == 0:
                value = False
            else:
                raise TypeError(
                    "Invalid value for has_start_codon: {0}".format(type(value)))
        self.__has_start_codon = value

    has_start_codon.category = "CDS"
    has_start_codon.rtype = "bool"

    @Metric
    def has_stop_codon(self):
        """Boolean. True if the selected ORF has a stop codon."""
        return self.__has_stop_codon and (self.combined_cds_length > 0)

    @has_stop_codon.setter
    def has_stop_codon(self, value):
        """Setter. Checks that the argument is boolean."""

        if not isinstance(value, bool) and value is not None:
            if value == 1:
                value = True
            elif value == 0:
                value = False
            else:
                raise TypeError(
                    "Invalid value for has_stop_codon: {0}".format(type(value)))

        self.__has_stop_codon = value

    has_stop_codon.category = "CDS"
    has_stop_codon.rtype = "bool"

    @Metric
    def is_complete(self):
        """Boolean. True if the selected ORF has both start and end."""
        return (self.has_start_codon is True) and (self.has_stop_codon is True)

    is_complete.category = "CDS"
    is_complete.rtype = "bool"

    @Metric
    def exon_num(self):
        """This property returns the number of exons of the transcript."""

        return self.__exon_num or len(self.exons)

    exon_num.category = "cDNA"
    exon_num.rtype = "int"

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

        if not isinstance(args[0], (float, int)):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        elif args[0] <= 0 or args[0] > 1:
            raise ValueError("Invalid value for the fraction: {0}".format(args[0]))
        self.__exon_fraction = args[0]

    exon_fraction.category = "Locus"
    exon_fraction.usable_raw = True
    exon_fraction.rtype = "float"

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
            raise ValueError("""It is impossible that the intron fraction is null \
when the transcript has at least one intron!""")
        self.__intron_fraction = args[0]

    intron_fraction.category = "Locus"
    intron_fraction.usable_raw = True
    intron_fraction.rtype = "float"

    @Metric
    def combined_cds_locus_fraction(self):
        """This metric returns the fraction of CDS bases of the transcript
        vs. the total of CDS bases in the locus."""
        return self.__combined_cds_locus_fraction

    @combined_cds_locus_fraction.setter
    def combined_cds_locus_fraction(self, value):
        if not isinstance(value, (int, float)) and 0 <= value <= 1:
            raise TypeError("The fraction should be a number between 0 and 1")
        elif self.combined_cds_length == 0 and value > 0:
            raise ValueError("{} has no CDS, its CDS fraction cannot be greater than 0!".format(self.id))
        self.__combined_cds_locus_fraction = value

    combined_cds_locus_fraction.category = "Locus"
    combined_cds_locus_fraction.usable_raw = True
    combined_cds_locus_fraction.rtype = "float"

    @Metric
    def selected_cds_locus_fraction(self):
        """This metric returns the fraction of CDS bases of the transcript
        vs. the total of CDS bases in the locus."""
        return self.__selected_cds_locus_fraction

    @selected_cds_locus_fraction.setter
    def selected_cds_locus_fraction(self, value):
        if not isinstance(value, (int, float)) and 0 <= value <= 1:
            raise TypeError("The fraction should be a number between 0 and 1")
        elif self.selected_cds_length == 0 and value > 0:
            raise ValueError("{} has no CDS, its CDS fraction cannot be greater than 0!")
        self.__selected_cds_locus_fraction = value

    selected_cds_locus_fraction.category = "Locus"
    selected_cds_locus_fraction.usable_raw = True
    selected_cds_locus_fraction.rtype = "float"

    @Metric
    def max_intron_length(self):
        """This property returns the greatest intron length for the transcript."""
        return self.__max_intron_length

    max_intron_length.category = "Intron"
    max_intron_length.rtype = "int"

    def __calculate_max_intron_length(self):
        if len(self.introns) == 0:
            self.__max_intron_length = 0
        else:
            self.__max_intron_length = max(intron[1] + 1 - intron[0] for intron in self.introns)

    @Metric
    def min_intron_length(self):
        """This property returns the smallest intron length for the transcript."""
        return self.__min_intron_length

    def __calculate_min_intron_length(self):

        if len(self.introns) == 0:
            self.__min_intron_length = 0
        else:
            self.__min_intron_length = min(intron[1] + 1 - intron[0] for intron in self.introns)

    min_intron_length.category = "Intron"
    min_intron_length.rtype = "int"

    @Metric
    def start_distance_from_tss(self):
        """This property returns the distance of the start of the combined CDS
        from the transcript start site.
        If no CDS is defined, it defaults to 0."""

        return self.__start_distance_from_tss

    start_distance_from_tss.category = "CDS"
    start_distance_from_tss.rtype = "int"

    def __calculate_start_distance_from_tss(self):

        if len(self.internal_orfs) < 2:
            self.__start_distance_from_tss = self.selected_start_distance_from_tss
        else:
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
            self.__start_distance_from_tss = distance

    # pylint: disable=invalid-name
    @Metric
    def selected_start_distance_from_tss(self):
        """This property returns the distance of the start of the best CDS
        from the transcript start site.
        If no CDS is defined, it defaults to 0."""
        return self.__selected_start_distance_from_tss

    selected_start_distance_from_tss.category = "CDS"
    selected_start_distance_from_tss.rtype = "int"

    def __calculate_selected_start_distance_from_tss(self):
        if len(self.combined_cds) == 0:
            self.__selected_start_distance_from_tss = 0
            return
        else:
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
            self.__selected_start_distance_from_tss = distance

    @Metric
    def selected_end_distance_from_tes(self):
        """This property returns the distance of the end of the best CDS
        from the transcript end site.
        If no CDS is defined, it defaults to 0."""
        return self.__selected_end_distance_from_tes

    selected_end_distance_from_tes.category = "CDS"
    selected_end_distance_from_tes.rtype = "int"

    def __calculate_selected_end_distance_from_tes(self):
        if self.is_coding is False:
            self.__selected_end_distance_from_tes = 0
        else:
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
            self.__selected_end_distance_from_tes = distance

    @Metric
    def selected_end_distance_from_junction(self):
        """This metric returns the distance between the stop codon and the
        last junction of the transcript. In many eukaryotes, this distance
        cannot exceed 50-55 bps, otherwise the transcript becomes a target of NMD.
        If the transcript is not coding or there is no junction downstream of
        the stop codon, the metric returns 0."""

        return self.__selected_end_distance_from_junction

    selected_end_distance_from_junction.category = "CDS"
    selected_end_distance_from_junction.rtype = "int"

    def __calculate_selected_end_distance_from_junction(self):

        if self.monoexonic is True or self.is_coding is False:
            self.__selected_end_distance_from_junction = 0
        else:
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
            self.__selected_end_distance_from_junction = distance

    @Metric
    def end_distance_from_junction(self):
        """This metric returns the cDNA distance between the stop codon and
        the last junction in the transcript.
        In many eukaryotes, this distance cannot exceed 50-55 bps
        otherwise the transcript becomes a target of NMD.
        If the transcript is not coding or there is no junction downstream
        of the stop codon, the metric returns 0.
        This metric considers the combined CDS end."""

        return self.__end_distance_from_junction

    end_distance_from_junction.category = "CDS"
    end_distance_from_junction.rtype = "int"

    def __calculate_end_distance_from_junction(self):

        if self.monoexonic is True or self.is_coding is False:
            self.__end_distance_from_junction = 0
        else:
            distance = 0
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
            self.__end_distance_from_junction = distance

    @Metric
    def end_distance_from_tes(self):
        """This property returns the distance of the end of the combined CDS
        from the transcript end site.
        If no CDS is defined, it defaults to 0."""

        return self.__end_distance_from_tes

    end_distance_from_tes.category = "CDS"
    end_distance_from_tes.rtype = "int"

    def __calculate_end_distance_from_tes(self):
        if self.is_coding is False:
            self.__end_distance_from_tes = 0
        else:
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
            self.__end_distance_from_tes = distance

    def _set_distances(self):
        self.__calculate_end_distance_from_tes()
        self.__calculate_end_distance_from_junction()
        self.__calculate_selected_end_distance_from_junction()
        self.__calculate_selected_end_distance_from_tes()
        self.__calculate_selected_start_distance_from_tss()
        self.__calculate_start_distance_from_tss()

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
    combined_cds_intron_fraction.usable_raw = True

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
    selected_cds_intron_fraction.usable_raw = True
    selected_cds_intron_fraction.rtype = "float"

    @Metric
    def retained_intron_num(self):
        """This property records the number of introns in the transcripts
        which are marked as being retained.
        See the corresponding method in the sublocus class."""
        return len(self.retained_introns)

    retained_intron_num.category = "Locus"
    retained_intron_num.rtype = "int"

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
    retained_fraction.usable_raw = True
    retained_fraction.rtype = "float"

    @Metric
    def cds_disrupted_by_ri(self):
        """This property describes whether the CDS is interrupted within a retained intron."""

        return self.__cds_disrupted_by_ri

    @cds_disrupted_by_ri.setter
    def cds_disrupted_by_ri(self, boolean):
        if boolean not in (True, False):
            raise TypeError("Invalid value for cds_disrupted_by_ri, it should be boolean")
        self.__cds_disrupted_by_ri = boolean

    cds_disrupted_by_ri.category = "Locus"
    cds_disrupted_by_ri.usable_raw = True
    cds_disrupted_by_ri.rtype = "bool"

    @Metric
    def proportion_verified_introns(self):
        """This metric returns, as a fraction, how many of the transcript introns
        are validated by external data. Monoexonic transcripts are set to 1."""
        if self.monoexonic is True:
            return 0
        else:
            return len(self.verified_introns) / len(self.introns)

    proportion_verified_introns.category = "External"
    proportion_verified_introns.usable_raw = True
    proportion_verified_introns.rtype = "float"

    @Metric
    def non_verified_introns_num(self):
        """
        This metric returns the number of introns of the transcript which are not validated
        by external data."""

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
    non_verified_introns_num.rtype = "int"

    @property
    def verified_introns(self):
        """This property holds the verified introns in a set. It also verifies that the introns are contained
        within the transcript."""

        if set.difference(self.__verified_introns, self.introns):
            self.logger.debug("Invalid verified junctions found for %s, removing them", self.id)
            self.__verified_introns = set.intersection(self.introns, self.__verified_introns)
        return self.__verified_introns

    @verified_introns.setter
    def verified_introns(self, item):

        if not isinstance(item, set):
            raise TypeError("Invalid type for verified junctions, expected set, got {}".format(
                type(item)))
        self.__verified_introns = set.intersection(item, self.introns)

    @Metric
    def verified_introns_num(self):
        """
        This metric returns the number of introns of the transcript which are validated
        by external data."""

        assert len(self.verified_introns) <= len(self.introns)
        return len(self.verified_introns)

    verified_introns_num.category = "External"
    verified_introns_num.rtype = "int"

    @Metric
    def proportion_verified_introns_inlocus(self):
        """This metric returns, as a fraction, how many of the
        verified introns inside the Locus are contained inside the transcript."""
        return self.__proportion_verified_introns_inlocus

    @proportion_verified_introns_inlocus.setter
    def proportion_verified_introns_inlocus(self, *args):
        """Setter for retained_intron_fraction."""

        if not isinstance(args[0], (float, int)) or (args[0] < 0 or args[0] > 1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))

        value = args[0]
        if value == 0:
            assert len(self.verified_introns) == 0
        assert 0 <= value <= 1
        self.__proportion_verified_introns_inlocus = value

    proportion_verified_introns_inlocus.category = "Locus"
    proportion_verified_introns_inlocus.rtype = "float"
    proportion_verified_introns_inlocus.usable_raw = True

    @Metric
    def num_introns_greater_than_max(self):
        """
        This metric returns the number of introns greater
        than the maximum acceptable intron size
        indicated in the constructor."""

        return sum(1 for intron in self.introns if
                   intron[1] - intron[0] + 1 > self.intron_range[1])
        #
        # return len(list(filter(lambda x: x[1]-x[0]+1 > self.intron_range[1],
        #                        self.introns)))

    num_introns_greater_than_max.category = "Intron"
    num_introns_greater_than_max.rtype = "int"

    @Metric
    def num_introns_smaller_than_min(self):
        """
        This metric returns the number of introns smaller
        than the mininum acceptable intron size
        indicated in the constructor."""

        return sum(1 for intron in self.introns if
                   intron[1] - intron[0] + 1 < self.intron_range[0])

        #
        # return len(list(filter(lambda x: x[1]-x[0]+1 < self.intron_range[0],
        #                        self.introns)))

    num_introns_smaller_than_min.category = "Intron"
    num_introns_smaller_than_min.rtype = "int"

    @Metric
    def snowy_blast_score(self):

        """
        Metric that indicates how good a hit is compared to the competition, in terms of BLAST
        similarities.
        As in SnowyOwl, the score for each hit is calculated by taking the coverage of the target
        and dividing it by (2 * len(self.blast_hits)).
        IMPORTANT: when splitting transcripts by ORF, a blast hit is added to the new transcript
        only if it is contained within the new transcript.
        This WILL screw up a bit the homology score.
        """

        if len(self.blast_hits) == 0:
            self.__blast_score = 0
        elif self.__blast_score == 0 and len(self.blast_hits) > 0:

            score = 0
            for hit in self.blast_hits:
                score += (hit["target_cov"]) / (2 * len(self.blast_hits))

            self.__blast_score = score

        return self.__blast_score

    snowy_blast_score.category = "External"
    snowy_blast_score.rtype = "float"

    @Metric
    def best_bits(self):
        """Metric that returns the best BitS associated with the transcript."""

        return max([0] + [hit["bits"] for hit in self.blast_hits])

    best_bits.category = "External"
    best_bits.rtype = "float"

    @Metric
    def blast_score(self):
        """
        Interchangeable alias for testing different blast-related scores.
        Current: best bit score."""

        # return self.snowy_blast_score
        return self.best_bits

    blast_score.category = "External"
    blast_score.rtype = "float"

    @Metric
    def blast_query_coverage(self):

        """
        This metric will return the **query** coverage for the best BLAST hit according to the evalue.
        If no BLAST hits are available for the sequence, it will return 0.
        :return:
        """

        if len(self.blast_hits) == 0:
            return 0
        else:
            hits = sorted(self.blast_hits,
                          key=operator.itemgetter("evalue"),
                          reverse=False)
            return hits[0]["query_cov"]

    blast_query_coverage.category = "External"
    blast_query_coverage.rtype = "float"
    blast_query_coverage.usable_raw = True

    @Metric
    def blast_identity(self):

        """
        This metric will return the alignment identity for the best BLAST hit according to the evalue.
        If no BLAST hits are available for the sequence, it will return 0.
        :return:
        :return:
        """

        if len(self.blast_hits) == 0:
            return 0
        else:
            hits = sorted(self.blast_hits,
                          key=operator.itemgetter("evalue"),
                          reverse=False)
            assert hits[0]["global_identity"] is not None, hits[0]

            return hits[0]["global_identity"] / 100

    blast_identity.category = "External"
    blast_identity.rtype = "float"
    blast_identity.usable_raw = True

    @Metric
    def blast_target_coverage(self):

        """
        This metric will return the **target** coverage for the best BLAST hit according to the evalue.
        If no BLAST hits are available for the sequence, it will return 0.
        :return:
        :return:
        """

        if len(self.blast_hits) == 0:
            return 0
        else:
            hits = sorted(self.blast_hits,
                          key=operator.itemgetter("evalue"),
                          reverse=False)
            return hits[0]["target_cov"]

    blast_target_coverage.category = "External"
    blast_target_coverage.rtype = "float"
    blast_target_coverage.usable_raw = True

    @Metric
    def canonical_intron_proportion(self):

        """
        This metric returns the proportion of canonical introns
         of the transcript on its total number of introns."""

        return float(self.attributes.get("canonical_proportion", 0))

    canonical_intron_proportion.category = "Intron"
    canonical_intron_proportion.usable_raw = True
    canonical_intron_proportion.rtype = "float"

    @Metric
    def suspicious_splicing(self):

        """This metric will return True if the transcript either has canonical introns
        on both strands (probably a chimeric artifact between two neighbouring loci,
        or if it has no canonical splicing event but it would if it were assigned to the opposite strand
        (probably a strand misassignment on the part of the assembler/predictor)."""

        mixed = bool(self.attributes.get("mixed_splices", False))
        canonical_on_reverse = self.attributes.get("canonical_on_reverse_strand", False)
        if not isinstance(canonical_on_reverse, bool):
            canonical_on_reverse = literal_eval(canonical_on_reverse)
            self.attributes["canonical_on_reverse_strand"] = canonical_on_reverse

        return self.monoexonic is False and (canonical_on_reverse or mixed)

    suspicious_splicing.category = "Intron"
    suspicious_splicing.rtype = "bool"

    @Metric
    def only_non_canonical_splicing(self):

        """This metric will return True if the canonical_number is 0"""

        return self.monoexonic is False and (int(self.attributes.get("canonical_number", 1)) == 0)

    only_non_canonical_splicing.category = "Intron"
    only_non_canonical_splicing.rtype = "bool"

    @Metric
    def max_exon_length(self):
        """This metric will return the length of the biggest exon in the transcript."""

        if len(self.exons) == 0:
            return 0
        else:
            return max([_[1] - _[0] + 1 for _ in self.exons])

    max_exon_length.category = "cDNA"
    max_exon_length.rtype = "int"

    @Metric
    def min_exon_length(self):
        """This metric will return the length of the biggest exon in the transcript."""

        if len(self.exons) == 0:
            return 0
        else:
            return min([_[1] - _[0] + 1 for _ in self.exons])

    max_exon_length.category = "cDNA"
    max_exon_length.rtype = "int"

