#!/usr/bin/env python3
# coding: utf-8

"""
Superlocus module. The class here defined is the uppermost container for transcripts
and is used to define all the possible children (subloci, monoloci, loci, etc.)
"""

# Core imports
import collections
from sys import version_info
import networkx
from sqlalchemy import bindparam
from sqlalchemy.engine import Engine
from sqlalchemy.ext import baked
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.orm import session as sasession
from sqlalchemy.sql.expression import and_
from .abstractlocus import Abstractlocus
from .monosublocusholder import MonosublocusHolder
from .sublocus import Sublocus
from ..exceptions import NotInLocusError
from ..parsers.GFF import GffLine
from ..serializers.blast_serializer import Hit, Query, Target
from ..serializers.external import External
from ..serializers.junction import Junction, Chrom
from ..serializers.orf import Orf
from ..utilities import dbutils, grouper
from ..scales.assignment.assigner import Assigner
import bisect
from sys import maxsize
import functools

if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict
from .locus import Locus
from .excluded import Excluded
from typing import Union
from ..utilities.intervaltree import Interval, IntervalTree
from itertools import combinations
import random

# The number of attributes is something I need
# pylint: disable=too-many-instance-attributes


class Superlocus(Abstractlocus):
    """The superlocus class is used to define overlapping regions
    on the genome, and it receives as input transcript class instances.
    """

    __name__ = "superlocus"

    bakery = baked.bakery()
    db_baked = bakery(lambda session: session.query(Chrom))
    db_baked += lambda q: q.filter(Chrom.name == bindparam("chrom_name"))

    junction_baked = bakery(lambda session: session.query(Junction))
    junction_baked += lambda q: q.filter(and_(
        Junction.chrom == bindparam("chrom"),
        Junction.junction_start >= bindparam("junctionStart"),
        Junction.junction_end <= bindparam("junctionEnd"),
        Junction.strand == bindparam("strand")
    ))

    hit_baked = bakery(lambda session: session.query(Hit))
    hit_baked += lambda q: q.filter(and_(
        Hit.query_id == bindparam("query_id"),
        Hit.evalue <= bindparam("evalue"),
        Hit.hit_number <= bindparam("hit_number")
    ))

    _complex_limit = (5000, 5000)

    # Junction.strand == bindparam("strand")))

    # ###### Special methods ############

    def __init__(self,
                 transcript_instance,
                 verified_introns=None,
                 stranded=True,
                 json_conf=None,
                 source="",
                 logger=None,
                 **kwargs):

        """

        :param transcript_instance: an instance of the Transcript class
        :type transcript_instance: [Transcript|None]
        :param stranded: boolean flag that indicates whether
        the Locus should use or ignore strand information
        :type stranded: bool
        :param json_conf: a configuration dictionary derived from JSON/YAML config files
        :type json_conf: dict
        :param source: optional source for the locus
        :type source: str
        :param logger: the logger for the class
        :type logger: logging.Logger

        The superlocus class is instantiated from a transcript_instance class,
        which it copies in its entirety.

        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - introns - a *set* which contains the positions of each
        *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus

        The constructor method takes the following keyword arguments:
        - stranded    True if all transcripts inside the superlocus are required
        to be on the same strand
        - json_conf    Required. A dictionary with the coniguration necessary
        for scoring transcripts.
        - purge        Flag. If True, all loci holding only transcripts with a 0
        score will be deleted
        from further consideration.
        """

        super().__init__(source=source,
                         transcript_instance=None,
                         verified_introns=verified_introns,
                         json_conf=json_conf,
                         logger=logger,
                         **kwargs)
        self.approximation_level = 0
        self.feature = self.__name__
        # if json_conf is None or not isinstance(json_conf, dict):
        #     raise NoJsonConfigError("I am missing the configuration for prioritizing transcripts!")
        self.__regressor = None
        self.stranded = stranded
        self.splices = set(self.splices)
        self.introns = set(self.introns)

        # Flags
        self.subloci_defined = False
        self.monosubloci_defined = False
        self.loci_defined = False
        self.monosubloci_metrics_calculated = False

        # Objects used during computation
        self.subloci = []
        self.loci = SortedDict()
        self.monosubloci = dict()
        self.monoholders = []

        # Connection objects
        self.engine = self.sessionmaker = self.session = None
        # Excluded object
        self.excluded = Excluded(json_conf=self.json_conf)
        self.__retained_sources = set()
        self.__data_loaded = False
        self.__lost = dict()
        if transcript_instance is not None:
            super().add_transcript_to_locus(transcript_instance)
            assert transcript_instance.monoexonic is True or len(self.introns) > 0
            if self.stranded is True:
                self.strand = transcript_instance.strand
        else:
            self.logger.warning("Creating an empty superlocus as no transcript was provided!")

    def __create_locus_lines(self, superlocus_line, new_id, print_cds=True):

        """
        Private method to prepare the lines for printing out loci
        into GFF/GTF files.
        """

        lines = []
        self.define_loci()
        if len(self.loci) > 0:
            source = "{0}_loci".format(self.source)
            superlocus_line.source = source
            lines.append(str(superlocus_line))
            found = dict()

            for _, locus_instance in self.loci.items():
                locus_instance.source = source
                locus_instance.parent = new_id
                if locus_instance.id in found:
                    found[locus_instance.id] += 1
                    locus_instance.counter = found[locus_instance.id]
                else:
                    found[locus_instance.id] = 0
                lines.append(locus_instance.__str__(print_cds=print_cds).rstrip())
        return lines

    def __create_monolocus_holder_lines(self, superlocus_line, new_id, print_cds=True):

        """
        Private method to prepare the lines for printing out monosubloci
        into GFF/GTF files.
        """

        lines = []
        self.define_monosubloci()
        if len(self.monoholders) > 0:
            source = "{0}_monosubloci".format(self.source)
            superlocus_line.source = source
            lines.append(str(superlocus_line))
            found = dict()
            for monosublocus_instance in self.monoholders:
                monosublocus_instance.source = source
                monosublocus_instance.parent = new_id
                if monosublocus_instance.id in found:
                    found[monosublocus_instance.id] += 1
                    monosublocus_instance.counter = found[monosublocus_instance.id]
                else:
                    found[monosublocus_instance.id] = 0

                lines.append(monosublocus_instance.__str__(print_cds=print_cds).rstrip())

        return lines

    def __create_monolocus_lines(self, superlocus_line, new_id, print_cds=True):

        """
        Private method to prepare the lines for printing out monosubloci
        into GFF/GTF files.
        """

        lines = []
        self.define_monosubloci()
        if len(self.monosubloci) > 0:
            source = "{0}_monosubloci".format(self.source)
            superlocus_line.source = source
            lines.append(str(superlocus_line))
            found = dict()
            for monosublocus_instance in self.monosubloci:
                monosublocus_instance.source = source
                monosublocus_instance.parent = new_id
                if monosublocus_instance.id in found:
                    found[monosublocus_instance.id] += 1
                    monosublocus_instance.counter = found[monosublocus_instance.id]
                else:
                    found[monosublocus_instance.id] = 0

                lines.append(monosublocus_instance.__str__(print_cds=print_cds).rstrip())

        return lines

    def __create_sublocus_lines(self, superlocus_line, new_id, print_cds=True):
        """
        Private method to prepare the lines for printing out subloci
        into GFF/GTF files.
        """

        source = "{0}_subloci".format(self.source)
        superlocus_line.source = source
        lines = [str(superlocus_line)]
        self.define_subloci()
        found = dict()
        for sublocus_instance in self.subloci:
            try:
                sublocus_instance.source = source
            except AttributeError:
                raise AttributeError(sublocus_instance)
            sublocus_instance.parent = new_id
            if sublocus_instance.id in found:
                found[sublocus_instance.id] += 1
                sublocus_instance.counter = found[sublocus_instance.id]
            else:
                found[sublocus_instance.id] = 0
            lines.append(sublocus_instance.__str__(print_cds=print_cds).rstrip())
        return lines

    # This discrepancy with the base class is necessary
    # pylint: disable=arguments-differ

    def format(self, print_cds=True, level=None):

        """
        Alias for __str__.
        :param print_cds: Boolean. It indicates whether to print the CDS features or not.
        :param level: level which we wish to print for. Can be "loci", "subloci", "monosubloci"
        :return: formatted GFF strings
        """
        return self.__str__(print_cds=print_cds,
                            level=level)

    def __str__(self, level=None, print_cds=True):

        """
        :param level: level which we wish to print for. Can be "loci", "subloci", "monosubloci"
        :type level: str
        :param print_cds: flag. If set to False, only the exonic information will be printed.
        :type print_cds: bool

        This function will return the desired level of children loci.
        The keyword level accepts the following four values:
        - "None" - print whatever is available.
        - "loci" - print the final loci
        - "monosubloci" - print the monosubloci
        - "subloci" - print the subloci.

        The function will then return the desired location in GFF3-compliant format.
        """

        if abs(self.start) == float("inf") or abs(self.start) == maxsize:
            return ''

        superlocus_line = GffLine('')
        superlocus_line.chrom = self.chrom
        superlocus_line.feature = self.__name__
        superlocus_line.start, \
            superlocus_line.end, \
            superlocus_line.score = self.start, self.end, "."
        superlocus_line.strand = self.strand
        superlocus_line.phase, superlocus_line.score = None, None
        new_id = "{0}_{1}".format(self.source, self.id)
        superlocus_line.id, superlocus_line.name = new_id, self.name
        if self.approximation_level > 0:
            superlocus_line.attributes["approximation_level"] = self.approximation_level
        if len(self.__retained_sources) > 0:
            superlocus_line.attributes["retained_sources"] = ",".join(
                sorted(list(self.__retained_sources))
            )

        lines = []
        if level not in (None, "loci", "subloci", "monosubloci"):
            raise ValueError("Unrecognized level: {0}".format(level))

        elif level == "loci" or (level is None and self.loci_defined is True):
            lines = self.__create_locus_lines(
                superlocus_line,
                new_id,
                print_cds=print_cds
            )
        elif level == "monosubloci" or (level is None and self.monosubloci_defined is True):
            lines = self.__create_monolocus_holder_lines(superlocus_line,
                                                         new_id,
                                                         print_cds=print_cds)
            # lines = self.__create_monolocus_lines(superlocus_line,
            #                                       new_id,
            #                                       print_cds=print_cds)
        elif level == "subloci" or (level is None and self.monosubloci_defined is False):
            lines = self.__create_sublocus_lines(superlocus_line,
                                                 new_id,
                                                 print_cds=print_cds)
        if len(lines) > 0:
            lines.append("###")
        return "\n".join([line for line in lines if line is not None and line != ''])
    # pylint: enable=arguments-differ

    def add_locus(self, locus: Locus):

        for transcript in locus.transcripts.values():
            self.add_transcript_to_locus(transcript, check_in_locus=False)

        # self.source = locus.source
        self.json_conf = locus.json_conf
        self.loci[locus.id] = locus
        self.loci_defined = True

    def as_dict(self, with_subloci=True, with_monoholders=True):

        state = super().as_dict()
        state["start"], state["end"] = self.start, self.end
        if with_subloci is True:
            state["subloci"] = [sublocus.as_dict() for sublocus in self.subloci]
        else:
            state["subloci"] = []
        state["loci"] = dict((lid, locus.as_dict()) for lid, locus in self.loci.items())
        if with_monoholders is True:
            state["monoholders"] = [mono.as_dict() for mono in self.monoholders]
        else:
            state["monoholders"] = []
        state["excluded"] = self.excluded.as_dict()
        return state

    def load_dict(self, state, print_subloci=True, print_monoloci=True, load_transcripts=True):
        super().load_dict(state, load_transcripts=load_transcripts)
        # self.json_conf = state["json_conf"]
        # self.chrom, self.start, self.end, self.strand = (state["chrom"], state["start"], state["end"], state["strand"])
        self.loci = dict()
        for lid, stat in state["loci"].items():
            locus = Locus()
            locus.load_dict(stat)
            self.loci[lid] = locus

        if print_subloci is True or print_monoloci is True:

            if print_subloci is True:
                self.excluded = Excluded()
                self.excluded.load_dict(state["excluded"])
                self.subloci = []
                for stat in state["subloci"]:
                    sub = Sublocus(json_conf=self.json_conf)
                    sub.load_dict(stat)
                    assert isinstance(sub, Sublocus)
                    self.subloci.append(sub)
            if print_monoloci is True:
                self.monoholders = []
                for stat in state["monoholders"]:
                    sub = MonosublocusHolder()
                    sub.load_dict(stat)
                    self.monoholders.append(sub)
        self.chrom, self.strand, self.start, self.end = state["chrom"], state["strand"], state["start"], state["end"]
        if len(self.loci) > 0 or len(self.transcripts) > 0:
            assert self.start != maxsize
            assert not self.id.endswith(str(maxsize))

    # ########### Class instance methods ############

    def split_strands(self):
        """This method will divide the superlocus on the basis of the strand.
        The rationale is to parse a GFF file without regard for the
        strand, in order to find all intersecting loci;
        and subsequently break the superlocus into the different components.
        Notice that each strand might generate more than one superlocus,
        if genes on a different strand link what are
        two different superloci.
        """

        self.logger.debug("Splitting by strand for {0}".format(self.id))
        if self.stranded is True:
            self.logger.warning("Trying to split by strand a stranded Locus, {0}!".format(self.id))
            yield self

        else:
            plus, minus, nones = [], [], []
            for cdna_id in self.transcripts:
                cdna = self.transcripts[cdna_id]
                self.logger.debug("{0}: strand {1}".format(cdna_id, cdna.strand))
                if cdna.strand == "+":
                    plus.append(cdna)
                elif cdna.strand == "-":
                    minus.append(cdna)
                elif cdna.strand is None:
                    nones.append(cdna)

            new_loci = []
            for strand in plus, minus, nones:
                if len(strand) > 0:
                    strand = sorted(strand)
                    new_locus = Superlocus(strand[0],
                                           stranded=True,
                                           json_conf=self.json_conf,
                                           source=self.source,
                                           logger=self.logger
                                           )
                    assert len(new_locus.introns) > 0 or new_locus.monoexonic is True
                    for cdna in strand[1:]:
                        if new_locus.in_locus(new_locus, cdna):
                            new_locus.add_transcript_to_locus(cdna)
                        else:
                            assert len(new_locus.introns) > 0 or new_locus.monoexonic is True
                            new_loci.append(new_locus)
                            new_locus = Superlocus(cdna,
                                                   stranded=True,
                                                   json_conf=self.json_conf,
                                                   source=self.source,
                                                   logger=self.logger
                                                   )
                    assert len(new_locus.introns) > 0 or new_locus.monoexonic is True
                    new_loci.append(new_locus)

            self.logger.debug(
                "Defined %d superloci by splitting by strand at %s.",
                len(new_loci), self.id)
            for new_locus in iter(sorted(new_loci)):
                if self.regressor is not None:
                    new_locus.regressor = self.regressor
                yield new_locus
        # raise StopIteration

    # @profile
    def connect_to_db(self, engine):

        """
        :param engine: the connection pool
        :type engine: Engine

        This method will connect to the database using the information
        contained in the JSON configuration.
        """

        if engine is None:
            self.engine = dbutils.connect(self.json_conf)
        else:
            self.engine = engine

        self.sessionmaker = sessionmaker()
        self.sessionmaker.configure(bind=self.engine)
        self.session = self.sessionmaker()

    # @asyncio.coroutine
    def load_transcript_data(self, tid, data_dict):
        """
        :param tid: the name of the transcript to retrieve data for.
        :type tid: str

        :param data_dict: the dictionary to use for data retrieval, if specified.
        If None, a DB connection will be established to retrieve the necessary data.
        :type data_dict: (None | dict)

        This routine is used to load data for a single transcript."""

        self.logger.debug("Retrieving data for {0}".format(tid))
        self.transcripts[tid].logger = self.logger
        self.transcripts[tid].json_conf = self.json_conf
        self.transcripts[tid].load_information_from_db(self.json_conf,
                                                       introns=self.locus_verified_introns,
                                                       session=self.session,
                                                       data_dict=data_dict)
        to_remove, to_add = False, dict()

        if self.json_conf["pick"]["chimera_split"]["execute"] is True:
            if self.transcripts[tid].number_internal_orfs > 1:
                new_tr = list(self.transcripts[tid].split_by_cds())
                if len(new_tr) > 1:
                    for new in new_tr:
                        assert new.id not in to_add
                        to_add[new.id] = new
                    to_remove = True
                    self.logger.info("%s has been split into %d different transcripts.",
                                     tid, len(new_tr))

        del data_dict
        return to_remove, to_add
        # @profile

    def _load_introns(self, data_dict):

        """Private method to load the intron data into the locus.
        :param data_dict: Dictionary containing the preloaded data, if available.
        :return:
        """

        if len(self.introns) == 0:
            if self.monoexonic is False:
                raise ValueError("%s is multiexonic but has no introns defined!",
                                 self.id)
            self.logger.debug("No introns for %s", self.id)
            return

        self.logger.debug("Querying the DB for introns, %d total", len(self.introns))
        if data_dict is None:
            if self.json_conf["db_settings"]["db"] is None or self.json_conf["db_settings"]["db"] == "":
                return  # No data to load
            # dbquery = self.db_baked(self.session).params(chrom_name=self.chrom).all()

            ver_introns = self.engine.execute(" ".join([
                "select junction_start, junction_end, strand from junctions where",
                "chrom_id = (select chrom_id from chrom where name = \"{chrom}\")",
                "and junction_start > {start} and junction_end < {end}"]).format(
                    chrom=self.chrom, start=self.start, end=self.end
            ))
            ver_introns = dict(((junc.junction_start, junc.junction_end), junc.strand)
                               for junc in ver_introns)

            self.logger.debug("Found %d verifiable introns for %s",
                              len(ver_introns), self.id)

            for intron in self.introns:
                self.logger.debug("Checking %s%s:%d-%d",
                                  self.chrom, self.strand, intron[0], intron[1])
                if (intron[0], intron[1]) in ver_introns:
                    self.logger.debug("Verified intron %s:%d-%d",
                                      self.chrom, intron[0], intron[1])
                    self.locus_verified_introns.add((intron[0],
                                                     intron[1],
                                                     ver_introns[(intron[0], intron[1])]))
        else:
            for intron in self.introns:
                self.logger.debug("Checking %s%s:%d-%d",
                                  self.chrom, self.strand, intron[0], intron[1])
                key = (self.chrom, intron[0], intron[1])
                # Ignore the strand 'cause we are loading BEFORE splitting by strand
                if key in data_dict["junctions"]:
                    self.logger.debug("Verified intron %s%s:%d-%d",
                                      self.chrom,
                                      data_dict["junctions"][key],
                                      intron[0], intron[1])
                    # Start, Stop, Strand
                    self.locus_verified_introns.add((intron[0],
                                                     intron[1],
                                                     data_dict["junctions"][key]))

    def _create_data_dict(self, engine, tid_keys):

        assert engine is not None

        data_dict = dict()
        self.logger.debug("Starting to load hits and orfs for %d transcripts",
                          len(tid_keys))
        data_dict["hits"] = collections.defaultdict(list)
        data_dict["orfs"] = collections.defaultdict(list)
        data_dict["external"] = collections.defaultdict(dict)

        for tid_group in grouper(tid_keys, 100):
            query_ids = dict((query.query_id, query) for query in
                             self.session.query(Query).filter(
                                 Query.query_name.in_(tid_group)))
            # Retrieve the external scores

            if query_ids:
                external = self.session.query(External).filter(External.query_id.in_(query_ids.keys()))
            else:
                external = []

            for ext in external:
                if ext.rtype == "int":
                    score = int(ext.score)
                elif ext.rtype == "float":
                    score = float(ext.score)
                elif ext.rtype == "bool":
                    score = bool(int(ext.score))
                else:
                    raise ValueError("Invalid rtype: {}".format(ext.rtype))

                data_dict["external"][ext.query][ext.source] = (score, ext.valid_raw)

            # Load the ORFs from the table
            if query_ids:
                orfs = self.session.query(Orf).filter(Orf.query_id.in_(query_ids.keys()))
            else:
                orfs = []

            for orf in orfs:
                data_dict["orfs"][orf.query].append(orf.as_bed12())

            # Now retrieve the HSPs from the BLAST HSP table
            hsp_command = " ".join([
                "select * from hsp where",
                "hsp_evalue <= {0} and query_id in {1} order by query_id;"]).format(
                self.json_conf["pick"]["chimera_split"]["blast_params"]["hsp_evalue"],
                "({0})".format(", ".join([str(_) for _ in query_ids.keys()]))
            )

            hsps = dict()
            targets = set()

            for hsp in engine.execute(hsp_command):
                if hsp.query_id not in hsps:
                    hsps[hsp.query_id] = collections.defaultdict(list)
                hsps[hsp.query_id][hsp.target_id].append(hsp)
                targets.add(hsp.target_id)

            # Now that we have the HSPs, load the corresponding HITs
            hit_command = " ".join([
                "select * from hit where evalue <= {0}",
                "and hit_number <= {1} and query_id in {2}",
                "order by query_id, evalue asc;"
            ]).format(
                self.json_conf["pick"]["chimera_split"]["blast_params"]["evalue"],
                self.json_conf["pick"]["chimera_split"]["blast_params"]["max_target_seqs"],
                "({0})".format(", ".join([str(_) for _ in query_ids.keys()])))

            if len(targets) > 0:
                target_ids = dict()
                for target_group in grouper(targets, 100):
                    target_ids.update(dict((target.target_id, target) for target in
                                      self.session.query(Target).filter(
                                          Target.target_id.in_(target_group))))
            else:
                target_ids = dict()
            current_hit = None
            for hit in engine.execute(hit_command):
                if current_hit != hit.query_id:
                    current_hit = hit.query_id
                current_counter = 0

                current_counter += 1

                my_query = query_ids[hit.query_id]
                my_target = target_ids[hit.target_id]

                data_dict["hits"][my_query.query_name].append(
                    Hit.as_full_dict_static(
                        hit,
                        hsps[hit.query_id][hit.target_id],
                        my_query,
                        my_target
                    )
                )

        return data_dict

    def load_all_transcript_data(self, engine=None, data_dict=None):

        """
        This method will load data into the transcripts instances,
        and perform the split_by_cds if required
        by the configuration.
        Asyncio coroutines are used to decrease runtime.

        :param engine: a connection engine
        :type engine: Engine

        :param data_dict: the dictionary to use for data retrieval, if specified.
        If None, a DB connection will be established to retrieve the necessary data.
        :type data_dict: (None | dict)

        """

        if self.__data_loaded is True:
            return

        # Before starting to load: let's reduce the loci.
        if len(self.transcripts) >= self._complex_limit[0]:
            self.approximation_level = 1
            self.reduce_method_one(None)

        if data_dict is None and engine is None:
            raise ValueError("Both engine and data_dict are void.")
        elif data_dict is None:
            self.connect_to_db(engine)

        self.logger.debug("Type of data dict: %s",
                          type(data_dict))
        if isinstance(data_dict, dict):
            self.logger.debug("Length of data dict: %s", len(data_dict))

        tid_keys = list(self.transcripts.keys())
        to_remove, to_add = set(), dict()
        # This will function even if data_dict is None
        self._load_introns(data_dict)

        if data_dict is None:
            data_dict = self._create_data_dict(engine, tid_keys)

            self.logger.debug("Verified %d introns for %s",
                              len(self.locus_verified_introns),
                              self.id)

            self.logger.debug("Finished retrieving data for %d transcripts",
                              len(tid_keys))
            self.session.close()
            sasession.close_all_sessions()

        for tid in tid_keys:
            remove_flag, new_transcripts = self.load_transcript_data(tid, data_dict)
            if remove_flag is True:
                to_remove.add(tid)
                to_add.update(new_transcripts)

        assert len(to_remove) < len(to_add) + len(self.transcripts)

        if len(to_remove) > 0:
            self.logger.debug("Rebuilding the superlocus as %d transcripts have been excluded or split.",
                              len(to_remove))
            for tid in to_remove:
                self.remove_transcript_from_locus(tid)
            for tid, transcript in to_add.items():
                self.logger.debug("Adding %s to %s", tid, self.id)
                self.add_transcript_to_locus(transcript, check_in_locus=False)

        elif len(to_remove) == len(self.transcripts):
            self.logger.warning("No transcripts left for %s", self.name)
            [self.excluded.add_transcript_to_locus(_) for _ in to_remove]
            self._remove_all()

        del data_dict

        num_coding = 0
        for tid in self.transcripts:
            if self.transcripts[tid].combined_cds_length > 0:
                num_coding += 1
            else:
                self.transcripts[tid].feature = "ncRNA"

        self.logger.debug(
            "Found %d coding transcripts out of %d in %s",
            num_coding,
            len(self.transcripts),
            self.id)

        self.session = None
        self.__data_loaded = True
        self.sessionmaker = None
        self.stranded = False

    # ##### Sublocus-related steps ######

    def reduce_complex_loci(self, transcript_graph: networkx.Graph):

        """
        Method which checks whether a locus has too many transcripts and tries to reduce them.

        :param transcript_graph: the transcript graph to analyse for redundancies
        :return:
        """

        max_edges = max([d for n, d in transcript_graph.degree])
        self.approximation_level = 0
        if len(transcript_graph) < self._complex_limit[0] and max_edges < self._complex_limit[1]:
            return transcript_graph
        self.logger.warning("Complex superlocus with %d nodes \
        with the most connected having %d edges",
                            len(transcript_graph), max_edges)

        self.approximation_level = 1
        transcript_graph, max_edges = self.reduce_method_one(transcript_graph)

        if len(transcript_graph) < self._complex_limit[0] and max_edges < self._complex_limit[1]:
            self.logger.warning("Approximation level 1 for %s", self.id)
            self.logger.debug("Remaining transcripts: %s", ", ".join(self.transcripts.keys()))
            return transcript_graph

        self.logger.warning(
            "Still %d nodes with the most connected with %d edges after approximation 1",
            len(transcript_graph), max_edges)

        self.approximation_level = 2
        transcript_graph, max_edges = self.reduce_method_two(transcript_graph)
        # if len(transcript_graph) < self._complex_limit[0] and max_edges < self._complex_limit[1]:
        self.logger.warning("Approximation level 2 for %s", self.id)
        return transcript_graph

    def reduce_method_one(self, transcript_graph: Union[None, networkx.Graph]) -> [networkx.Graph, int]:

        """Approximation level one: we are going to group together all transcripts that have identical intron chains,
        and remove any that is completely contained within the longest ones. Reference transcripts get an automatic
        pass."""

        ichains = collections.defaultdict(list)
        cds_only = self.json_conf["pick"]["clustering"]["cds_only"]
        for transcript in self.transcripts.values():
            # Ordered by coordinates
            if cds_only:
                bisect.insort(ichains[tuple(sorted(transcript.selected_cds_introns))], transcript)
            else:
                bisect.insort(ichains[tuple(sorted(transcript.introns))], transcript)

        to_remove = set()
        for ichain, transcripts in ichains.items():
            current_coords = [transcripts[0].start, transcripts[0].end]
            current_id = transcripts[0].id
            for transcript in transcripts[1:]:  # we know that they are sorted left to right
                self.logger.debug("Comparing %s to %s", current_id, transcript.id)
                if transcript_graph and transcript.id not in transcript_graph.neighbors(current_id):
                    self.logger.debug("%s and %s are not neighbours.", current_id, transcript.id)
                    continue
                if transcript.end <= current_coords[1] and transcript.start > current_coords[0]:
                    if transcript.is_reference:
                        continue
                    else:
                        self.logger.debug("Removing %s as it is contained", transcript.id)
                        to_remove.add(transcript.id)
                elif transcript.end == current_coords[1] and transcript.start == current_coords[0]:
                    if transcript.is_reference:
                        current_id = transcript.id
                        if not self.transcripts[current_id].is_reference:
                            self.logger.debug("Removing %s as it is identical to a reference", current_id)
                            to_remove.add(current_id)
                    elif self.transcripts[current_id].is_reference:
                        self.logger.debug("Removing %s as it is identical to a reference", transcript.id)
                        to_remove.add(transcript.id)
                    else:
                        # chosen = np.random.choice(sorted([current_id, transcript.id]))
                        chosen = random.choice(sorted([current_id, transcript.id]))
                        if current_id == chosen:
                            to_remove.add(transcript.id)
                        else:
                            to_remove.add(chosen)
                            current_id = transcript.id
                elif transcript.end > current_coords[1]:
                    self.logger.debug("Removing %s as it is contained in %s", current_id, transcript.id)
                    if transcript.start == current_coords[0] and not self.transcripts[current_id].is_reference:
                        to_remove.add(current_id)
                    current_coords, current_id = (transcript.start, transcript.end), transcript.id

        if to_remove:
            if transcript_graph:
                transcript_graph.remove_nodes_from(to_remove)
            [self.excluded.add_transcript_to_locus(self.transcripts[tid]) for tid in to_remove]
            self.logger.debug("Removing the following transcripts from %s: %s",
                              self.id, ", ".join(to_remove))
            for tid in to_remove:
                self.remove_transcript_from_locus(tid)

        max_edges = max([d for n, d in transcript_graph.degree]) if transcript_graph else 0
        return transcript_graph, max_edges

    def reduce_method_two(self, transcript_graph: networkx.Graph) -> [networkx.Graph, int]:

        """Approximation level two: we are going to remove transcripts that are contained within others. So e.g.
        a transcript with two exons might be seen as redundant with, and therefore removed, with a transcript with
        three or more exons."""

        to_remove = set()
        done = set()
        cds_only = self.json_conf["pick"]["clustering"]["cds_only"]
        if cds_only:
            order = sorted(list(transcript_graph.nodes), key=lambda node: self.transcripts[node].selected_cds_start)
        else:
            order = sorted(list(transcript_graph.nodes), key=lambda node: self.transcripts[node].start)

        for tid in order:
            current = self.transcripts[tid]
            for neighbour in transcript_graph.neighbors(tid):
                couple = tuple(sorted((neighbour, tid)))
                if couple in done:
                    continue
                done.add(couple)
                self.logger.debug("Comparing %s to %s", tid, neighbour)
                if neighbour in to_remove:
                    continue
                neighbour = self.transcripts[neighbour]
                inters = set.intersection(current.introns, neighbour.introns)
                if inters == current.introns:
                    self.logger.debug("Evaluating %s (template %s) for removal", current.id, neighbour.id)
                    if current.is_reference:
                        self.logger.debug("Ignoring %s for removal as it is a reference transcript", current.id)
                        continue
                    if cds_only:
                        comparison = Assigner.compare(current.copy().remove_utrs(),
                                                      neighbour.copy().remove_utrs())[0]
                    else:
                        comparison = Assigner.compare(current, neighbour)[0]
                    self.logger.debug(
                        "Evaluating %s (template %s) for removal (%s)",
                        current.id, neighbour.id, (comparison.n_prec[0], comparison.n_recall[0]))
                    if comparison.n_prec[0] == 100:
                        self.logger.debug("Removing %s as it is contained in %s", current.id, neighbour.id)
                        to_remove.add(current.id)
                        break
                    self.logger.debug("%s is not perfectly contained in %s", current.id, neighbour.id)
                if inters == neighbour.introns:
                    if neighbour.is_reference:
                        self.logger.debug("Ignoring %s for removal as it is a reference transcript", neighbour.id)
                        continue
                    if cds_only:
                        comparison = Assigner.compare(neighbour.copy().remove_utrs(),
                                                      current.copy().remove_utrs())[0]
                    else:
                        comparison = Assigner.compare(neighbour, current)[0]
                    self.logger.debug(
                        "Evaluating %s (template %s) for removal (%s)",
                        current.id, neighbour.id, (comparison.n_prec[0], comparison.n_recall[0]))
                    if comparison.n_prec[0] == 100:
                        self.logger.debug("Removing %s as it is contained in %s", neighbour.id, current.id)
                        to_remove.add(neighbour.id)
                        break
                    self.logger.debug("%s is not perfectly contained in %s", neighbour.id, current.id)
                else:
                    self.logger.debug("%s and %s have differing introns, ignoring the comparison",
                                      neighbour.id, current.id)
                    continue

        if to_remove:
            transcript_graph.remove_nodes_from(to_remove)
            [self.excluded.add_transcript_to_locus(self.transcripts[tid]) for tid in to_remove]
            self.logger.debug("Removing the following transcripts from %s: %s", self.id, ", ".join(to_remove))
            for tid in to_remove:
                self.logger.debug("Removing %s", tid)
                self.remove_transcript_from_locus(tid)

        max_edges = max([d for n, d in transcript_graph.degree])
        return transcript_graph, max_edges

    def define_subloci(self):
        """This method will define all subloci inside the superlocus.
        Steps:
            - Call the BronKerbosch algorithm to define cliques
            - Call the "merge_cliques" algorithm the merge the cliques.
            - Create "sublocus" objects from the merged cliques
            and store them inside the instance store "subloci"
        """

        if self.subloci_defined is True:
            return

        self.compile_requirements()
        self.subloci = []

        # First, check whether we need to remove CDS from anything.
        for tid in super()._check_not_passing(section="cds_requirements"):
            self.transcripts[tid].strip_cds(strand_specific=True)
            self.metrics_calculated = False

        # Check whether there is something to remove
        self._check_requirements()
        excluded_tids = list(self._excluded_transcripts.keys())
        for excluded_tid in excluded_tids:
            to_remove = self._excluded_transcripts[excluded_tid]
            self.excluded.add_transcript_to_locus(to_remove, check_in_locus=False)
            self.remove_transcript_from_locus(to_remove.id)
            if excluded_tid in self._excluded_transcripts:
                del self._excluded_transcripts[excluded_tid]

        if len(self.transcripts) == 0:
            # we have removed all transcripts from the Locus. Set the flag to True and exit.
            self.subloci_defined = True
            return

        # cds_only = self.json_conf["pick"]["clustering"]["cds_only"]
        self.logger.debug("Calculating the transcript graph for %d transcripts", len(self.transcripts))
        transcript_graph = self.define_graph()

        transcript_graph = self.reduce_complex_loci(transcript_graph)
        if len(self.transcripts) > len(transcript_graph) and self.reference_update is False:
            self.logger.warning("Discarded %d transcripts from %s due to approximation level %d",
                                len(self.transcripts) - len(transcript_graph),
                                self.id,
                                self.approximation_level)
            for tid in set.difference(set(self.transcripts.keys()), set(transcript_graph.nodes())):
                self.logger.debug("Discarding %s from %s", tid, self.id)
                del self.transcripts[tid]

        if len(self.transcripts) == 0:
            # we have removed all transcripts from the Locus. Set the flag to True and exit.
            self.logger.warning("Discarded all transcripts from %s", self.id)
            self.subloci_defined = True
            return

        self.logger.debug("Calculated the transcript graph for %d transcripts: %s",
                          len(self.transcripts),
                          str(transcript_graph))
        self.logger.debug("Calculating the transcript communities")
        subloci = self.find_communities(transcript_graph)
        self.logger.debug("Calculated the transcript communities")

        # Now we should define each sublocus and store it in a permanent structure of the class
        for subl in subloci:
            if len(subl) == 0:
                continue
            subl = [self.transcripts[x] for x in subl]
            subl = sorted(subl)
            new_sublocus = Sublocus(subl[0],
                                    json_conf=self.json_conf,
                                    logger=self.logger,
                                    use_transcript_scores=self._use_transcript_scores
                                    )
            new_sublocus.logger = self.logger
            if self.regressor is not None:
                new_sublocus.regressor = self.regressor
            for ttt in subl[1:]:
                try:
                    new_sublocus.add_transcript_to_locus(ttt)
                except NotInLocusError as orig_exc:
                    exc_text = """Sublocus: {0}
                    Offending transcript:{1}
                    In locus manual check: {2}
                    Original exception: {3}""".format(
                        "{0} {1}:{2}-{3} {4}".format(
                            subl[0].id, subl[0].chrom, subl[0].start,
                            subl[0].end, subl[0].exons),
                        "{0} {1}:{2}-{3} {4}".format(ttt.id, ttt.chrom,
                                                     ttt.start, ttt.end, ttt.exons),
                        "Chrom {0} Strand {1} overlap {2}".format(
                            new_sublocus.chrom == ttt.chrom,
                            "{0}/{1}/{2}".format(
                                new_sublocus.strand,
                                ttt.strand,
                                new_sublocus.strand == ttt.strand
                            ),
                            self.overlap((subl[0].start, subl[1].end),
                                         (ttt.start, ttt.end)) > 0
                        ),
                        orig_exc
                    )
                    raise NotInLocusError(exc_text)

            new_sublocus.parent = self.id
            new_sublocus.metrics_calculated = False
            new_sublocus.get_metrics()
            self.subloci.append(new_sublocus)

        self.subloci = sorted(self.subloci)
        self.subloci_defined = True

    def define_monosubloci(self):

        """This is a wrapper method that defines the monosubloci for each sublocus.
        """
        if self.monosubloci_defined is True:
            return

        self.logger.debug("Calculating subloci for %s, %d transcripts",
                          self.id, len(self.transcripts))
        self.define_subloci()
        self.logger.debug("Calculated subloci for %s, %d transcripts",
                          self.id, len(self.transcripts))
        self.monosubloci = dict()
        # Extract the relevant transcripts
        for sublocus_instance in sorted(self.subloci):
            sublocus_instance.define_monosubloci(purge=self.purge)
            for transcript in sublocus_instance.excluded.transcripts.values():
                self.excluded.add_transcript_to_locus(transcript)
            for tid in sublocus_instance.transcripts:
                # Update the score
                self.transcripts[tid].score = sublocus_instance.transcripts[tid].score
            for monosubl in sublocus_instance.monosubloci:
                monosubl.parent = self.id
                # self.monosubloci.append(monosubl)
                self.monosubloci[monosubl.tid] = monosubl

        # self.monosubloci = sorted(self.monosubloci)
        if self.logger.level == 10:  # DEBUG
            self.logger.debug("Monosubloci for %s:\n\t\t%s",
                              self.id,
                              "\n\t\t".join(
                                  ["{}, transcript: {}".format(
                                      self.monosubloci[_].id, _) for _ in self.monosubloci]))

        self.monosubloci_defined = True

    def print_subloci_metrics(self):
        """Wrapper method to create a csv.DictWriter instance and call
        the sublocus.print_metrics method
        on it for each sublocus."""

        # self.get_sublocus_metrics()

        for slocus in self.subloci:
            for row in slocus.print_metrics():
                yield row
        if self.excluded is not None:
            yield from self.excluded.print_metrics()

    def print_subloci_scores(self):
        """Wrapper method to create a csv.DictWriter instance and call the
        sublocus.print_metrics method
        on it for each sublocus."""

        # self.get_sublocus_metrics()

        for slocus in self.subloci:
            for row in slocus.print_scores():
                yield row

    def print_monoholder_metrics(self):

        """Wrapper method to create a csv.DictWriter instance and call the
        MonosublocusHolder.print_metrics method
        on it."""

        self.define_monosubloci()

        # self.available_monolocus_metrics = set(self.monoholder.available_metrics)
        if len(self.monoholders) == 0:
            return
        for monoholder in self.monoholders:
            for row in monoholder.print_metrics():
                yield row

    def print_monoholder_scores(self):

        """Wrapper method to create a csv.DictWriter instance and call
        the MonosublocusHolder.print_scores method on it."""

        self.define_monosubloci()

        # self.available_monolocus_metrics = set(self.monoholder.available_metrics)
        if len(self.monoholders) == 0:
            return
        for monoholder in self.monoholders:
            for row in monoholder.print_scores():
                yield row

    def print_loci_metrics(self):

        self.define_loci()

        if len(self.loci) == 0:
            return []
        for locus in self.loci:
            for row in self.loci[locus].print_metrics():
                yield row

    def print_loci_scores(self):

        """Wrapper method to create a csv.DictWriter instance and call
        the Locus.print_scores method on it."""

        self.define_loci()

        # self.available_monolocus_metrics = set(self.monoholder.available_metrics)
        if len(self.loci) == 0:
            return
        for locus in self.loci:
            for row in self.loci[locus].print_scores():
                yield row

    def define_loci(self):
        """This is the final method in the pipeline. It creates a container
        for all the monosubloci (an instance of the class MonosublocusHolder)
        and retrieves the loci it calculates internally."""

        if self.loci_defined is True:
            return

        self.logger.debug("Calculating monosubloci for %s, %d transcripts",
                          self.id, len(self.transcripts))
        self.define_monosubloci()
        self.logger.debug("Calculated monosubloci for %s, %d transcripts",
                          self.id, len(self.transcripts))
        self.calculate_mono_metrics()

        self.loci = SortedDict()
        if len(self.monoholders) == 0:
            if len(self.transcripts) == 0:
                # We have already purged
                self.logger.debug("No locus retained for %s", self.id)
            else:
                self.logger.warning("No locus retained for %s", self.id)
            self.loci_defined = True
            return

        loci = []
        for monoholder in self.monoholders:
            monoholder.define_loci(purge=self.purge)
            for locus_instance in monoholder.loci:
                monoholder.loci[locus_instance].parent = self.id
                loci.append(monoholder.loci[locus_instance])

        for locus in sorted(loci):
            self.loci[locus.id] = locus
            self.loci[locus.id].logger = self.logger
            self.loci[locus.id].set_json_conf(self.json_conf)

        self.loci_defined = True

        self.logger.debug("Looking for AS events in %s: %s",
                          self.id,
                          self.json_conf["pick"]["alternative_splicing"]["report"])
        if self.json_conf["pick"]["alternative_splicing"]["report"] is True:
            self.define_alternative_splicing()

        self.__find_lost_transcripts()
        while len(self.lost_transcripts) > 0:
            new_locus = None
            for transcript in self.lost_transcripts.values():
                if new_locus is None:
                    new_locus = Superlocus(transcript,
                                           json_conf=self.json_conf,
                                           use_transcript_scores=self._use_transcript_scores,
                                           stranded=self.stranded,
                                           verified_introns=self.locus_verified_introns,
                                           logger=self.logger,
                                           source=self.source
                                           )
                else:
                    new_locus.add_transcript_to_locus(transcript,
                                                      check_in_locus=False)
            new_locus.define_loci()
            self.loci.update(new_locus.loci)
            self.__lost = new_locus.lost_transcripts

        if self.json_conf["pick"]["run_options"]["only_reference_update"] is True:
            lids = list(self.loci.keys())[:]
            for lid in lids:
                if self.loci[lid].has_reference_transcript is False:
                    self.logger.debug("Removing %s (primary: %s) as it has no reference transcripts.",
                                      lid, self.loci[lid].primary_transcript_id)
                    del self.loci[lid]
            self.logger.debug("Remaining loci in %s: %s", self.id, ",".join(list(self.loci.keys())))

        return

    def __find_lost_transcripts(self):

        self.__lost = dict()
        cds_only = self.json_conf["pick"]["clustering"]["cds_only"]
        simple_overlap = self.json_conf["pick"]["clustering"]["simple_overlap_for_monoexonic"]
        cdna_overlap = self.json_conf["pick"]["clustering"]["min_cdna_overlap"]
        cds_overlap = self.json_conf["pick"]["clustering"]["min_cds_overlap"]

        loci_transcripts = set()
        for locus in self.loci.values():
            loci_transcripts.update(set([_ for _ in locus.transcripts.keys()]))

        not_loci_transcripts = set.difference({_ for _ in self.transcripts.keys()
                                               if _ not in self._excluded_transcripts},
                                              loci_transcripts)

        for locus in self.loci.values():
            to_remove = set()
            for tid in not_loci_transcripts:
                found = False
                for ltid in locus.transcripts:
                    found = found or MonosublocusHolder.is_intersecting(
                        self.transcripts[tid],
                        locus[ltid],
                        cds_only=cds_only,
                        logger=self.logger,
                        min_cdna_overlap=cdna_overlap,
                        min_cds_overlap=cds_overlap,
                        simple_overlap_for_monoexonic=simple_overlap
                    )
                    if found:
                        self.logger.debug("%s intersects %s in %s, not lost.",
                                          tid, ltid, locus.id)
                        break
                if found:
                    to_remove.add(tid)
            not_loci_transcripts = set.difference(not_loci_transcripts, to_remove)

        self.__lost.update({tid: self.transcripts[tid] for tid in not_loci_transcripts})

        if len(self.__lost):
            self.logger.debug("Lost %s transcripts from %s; starting the recovery process: TIDs: %s",
                              len(self.lost_transcripts), self.id, ", ".join(self.__lost))

    def define_alternative_splicing(self):

        """
         This method will consider all possible candidates for alternative splicing
         for each of the final loci, after excluding transcripts which potentially map
         to more than one Locus (in order to remove chimeras).
         It will then call the add_transcript_to_locus method to try to add
         the transcript to the relevant Locus container.
        """

        # First off, define genes

        self.define_loci()

        candidates = collections.defaultdict(set)
        primary_transcripts = set(locus.primary_transcript_id for locus in self.loci.values())
        self.logger.debug("Primary transcripts: %s", primary_transcripts)

        cds_only = self.json_conf["pick"]["clustering"]["cds_only"]
        cds_overlap = self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"]
        cdna_overlap = self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"]

        self.logger.debug("Defining the transcript graph")
        t_graph = self.define_as_graph(inters=MonosublocusHolder.is_intersecting,
                                       cds_only=cds_only,
                                       min_cdna_overlap=cdna_overlap,
                                       min_cds_overlap=cds_overlap,
                                       simple_overlap_for_monoexonic=False)

        self.logger.debug("Defined the transcript graph")

        loci_cliques = dict()
        for lid, locus_instance in self.loci.items():
            try:
                neighbors = set(t_graph.neighbors(locus_instance.primary_transcript_id))
            except networkx.exception.NetworkXError:
                raise networkx.exception.NetworkXError(
                    "{} {}".format(
                    # locus_instance.primary_transcript.attributes["Alias"],
                    locus_instance.primary_transcript_id,
                    list(t_graph.nodes)
                ))
            loci_cliques[lid] = neighbors
            self.logger.debug("Neighbours for %s: %s", lid, neighbors)

        for tid in iter(tid for tid in self.transcripts if tid not in primary_transcripts):
            loci_in = list(llid for llid in loci_cliques if
                           tid in loci_cliques[llid])
            if len(loci_in) == 1:
                candidates[loci_in[0]].add(tid)
            elif len(loci_in) > 1:
                # These are transcripts that match more than one locus
                continue
            else:
                # These transcripts have been lost.
                self.__lost.update({tid: self[tid]})

        for lid in candidates:
            for tid in sorted(candidates[lid],
                              key=lambda ttid: self.transcripts[ttid].score,
                              reverse=True):
                self.logger.debug("Adding %s to %s", tid, lid)
                self.loci[lid].add_transcript_to_locus(self.transcripts[tid])
            self.loci[lid].finalize_alternative_splicing()

        # Now we have to recheck that no AS event is linking more than one locus.
        to_remove = collections.defaultdict(list)
        for lid in self.loci:
            for tid, transcript in [_ for _ in self.loci[lid].transcripts.items() if
                                    _[0] != self.loci[lid].primary_transcript_id]:
                for olid in [_ for _ in self.loci if _ != lid]:
                    is_compatible = MonosublocusHolder.in_locus(self.loci[olid],
                                                                transcript)
                    if is_compatible is True:
                        self.logger.debug("%s is compatible with more than one locus. Removing it.", tid)
                        to_remove[lid].append(tid)

        for lid in to_remove:
            for tid in to_remove[lid]:
                self.loci[lid].remove_transcript_from_locus(tid)

            self.loci[lid].finalize_alternative_splicing()

        # This is *necessary* because the names of the loci *MIGHT HAVE CHANGED*!
        new = dict()
        for lid in self.loci:
            new[lid] = self.loci[lid]

        self.loci = new
        return

    def calculate_mono_metrics(self):
        """Wrapper to calculate the metrics for the monosubloci."""
        self.monoholders = []
        self.define_monosubloci()

        mono_graph = super().define_graph(
            self.monosubloci,
            inters=MonosublocusHolder.in_locus,
            logger=self.logger,
            cds_only=self.json_conf["pick"]["clustering"]["cds_only"],
            min_cdna_overlap=self.json_conf["pick"]["clustering"]["min_cdna_overlap"],
            min_cds_overlap=self.json_conf["pick"]["clustering"]["min_cds_overlap"],
            simple_overlap_for_monoexonic=self.json_conf["pick"]["clustering"]["simple_overlap_for_monoexonic"])

        assert len(mono_graph.nodes()) == len(self.monosubloci)

        communities = self.find_communities(mono_graph)

        for community in communities:
            community = set(community)
            monosub = self.monosubloci[community.pop()]
            holder = MonosublocusHolder(monosub,
                                        json_conf=self.json_conf,
                                        logger=self.logger,
                                        use_transcript_scores=self._use_transcript_scores)
            while len(community) > 0:
                holder.add_monosublocus(self.monosubloci[community.pop()],
                                        check_in_locus=False)
            self.monoholders.append(holder)

        for monoholder in self.monoholders:
            # monoholder.scores_calculated = False
            if self.regressor is not None:
                monoholder.regressor = self.regressor
            monoholder.filter_and_calculate_scores()

    def compile_requirements(self):
        """Quick function to evaluate the filtering expression, if it is present."""

        for section in ["requirements", "as_requirements", "cds_requirements"]:
            if section in self.json_conf:
                if "compiled" in self.json_conf[section]:
                    return
                else:
                    self.json_conf[section]["compiled"] = compile(
                        self.json_conf[section]["expression"],
                        "<json>", "eval")
                    return
            else:
                continue

    # ############ Class methods ###########

    # The discrepancy is by design
    # pylint: disable=arguments-differ

    def define_graph(self) -> networkx.Graph:

        """Overloading of the original method to avoid being a O(n**2) algorithm."""

        graph = networkx.Graph()

        # As we are using intern for transcripts, this should prevent
        # memory usage to increase too much
        keys = set.difference(set(self.transcripts.keys()), set(self._excluded_transcripts))

        # objects = dict((tid, self.transcripts[tid]) for tid in keys)
        graph.add_nodes_from(keys)

        monos = []
        intronic = collections.defaultdict(set)

        for tid in keys:
            transcript = self.transcripts[tid]
            if self._cds_only is True and transcript.is_coding:
                if transcript.selected_cds_introns:
                    for intron in transcript.selected_cds_introns:
                        intronic[intron].add(tid)
                else:
                    start, end = sorted([transcript.selected_cds_start, transcript.selected_cds_end])
                    interval = Interval(start, end, transcript.id)
                    monos.append(interval)
                    monos.append(interval)
            else:
                if transcript.introns:
                    for intron in transcript.introns:
                        intronic[intron].add(tid)
                else:
                    interval = Interval(transcript.start, transcript.end, transcript.id)
                    monos.append(interval)

        edges = set()
        for intron in intronic:
            edges.update(set(combinations(intronic[intron], 2)))

        # Now the monoexonic
        monoexonic = IntervalTree()
        [monoexonic.add_interval(interval) for interval in monos]
        monos = sorted(monos)
        for mono in monos:
            edges.update(set((mono.value, omono.value) for omono in
                             (other for other in monoexonic.find(mono[0], mono[1], strict=False)
                              if mono.value != other.value)))
        graph.add_edges_from(edges)

        return graph

    def define_as_graph(self,
                        inters=MonosublocusHolder.is_intersecting,
                        cds_only=False,
                        min_cdna_overlap=0.2,
                        min_cds_overlap=0.2,
                        simple_overlap_for_monoexonic=True):

        """This method will try to build the AS graph using a O(nlogn) rather than O(n^2) algorithm."""

        method = functools.partial(inters, cds_only=cds_only, min_cdna_overlap=min_cdna_overlap,
                                   min_cds_overlap=min_cds_overlap,
                                   simple_overlap_for_monoexonic=simple_overlap_for_monoexonic,
                                   logger=self.logger)

        graph = networkx.Graph()
        graph.add_nodes_from(self.transcripts.keys())

        itree = IntervalTree()
        primaries = set()

        for lid in self.loci:
            transcript = self.loci[lid].primary_transcript
            primaries.add(transcript.id)
            if transcript.is_coding and cds_only:
                start, end = sorted([transcript.selected_cds_start, transcript.selected_cds_end])
            else:
                start, end = transcript.start, transcript.end
            itree.add_interval(Interval(start, end, transcript.id))

        for tid in self.transcripts:
            if tid in primaries:
                continue
            transcript = self[tid]
            if transcript.is_coding and cds_only:
                start, end = sorted([transcript.selected_cds_start, transcript.selected_cds_end])
            else:
                start, end = transcript.start, transcript.end
            for found in itree.find(start, end, strict=False):
                locus_transcript = found.value
                if method(self[locus_transcript], transcript):
                    graph.add_edge(tid, locus_transcript)

        return graph

    @classmethod
    def is_intersecting(cls, transcript, other, cds_only=False):
        """
        :rtype : bool
        :param transcript: a transcript for which we wish to verify
        whether it is intersecting with another transcript or not.
        :type transcript: Mikado.loci_objects.transcript.Transcript
        :param other: the transcript which will be used for the comparison.
        :type other: Mikado.loci_objects.transcript.Transcript

        :param cds_only: boolean flag. If enabled, only CDS exons/intron
        will be considered when deciding whether two transcripts are part
        of the same Locus or not.
        :type cds_only: bool


        When comparing two transcripts, for the definition of subloci inside
        superloci we follow these rules:

        If both are multiexonic, the function verifies whether there is at
        least one intron in common.
        If both are monoexonic, the function verifies whether there is some overlap between them.
        If one is monoexonic and the other is not, the function will return False by definition.
        """

        transcript.finalize()
        other.finalize()
        if transcript == other:
            return False  # We do not want intersection with oneself

        if transcript.monoexonic is False and other.monoexonic is False:
            if cds_only is False or transcript.is_coding is False or other.is_coding is False:
                intersection = set.intersection(transcript.introns, other.introns)
            else:
                intersection = set.intersection(transcript.selected_cds_introns,
                                                other.selected_cds_introns)
            intersecting = (len(intersection) > 0)

        elif transcript.monoexonic is True and other.monoexonic is True:

            if cds_only is False or transcript.is_coding is False or other.is_coding is False:
                intersecting = (cls.overlap(
                    (transcript.start, transcript.end),
                    (other.start, other.end), positive=False) > 0)
            else:
                # intersecting = any([cls.overlap(cds_comb[0],
                #                                 cds_comb[1],
                #                                 positive=False) > 0] for cds_comb in itertools.product(
                #     transcript.internal_orf_boundaries,
                #     other.internal_orf_boundaries))
                intersecting = (cls.overlap(
                    (transcript.selected_cds_start, transcript.selected_cds_end),
                    (other.selected_cds_start, other.selected_cds_end), positive=False) > 0)
        else:
            intersecting = False

        return intersecting
    # pylint: enable=arguments-differ

    # ############## Properties ############
    @property
    def id(self) -> str:
        """
        This is a generic string generator for all inherited children.
        :rtype : str
        """
        if self.stranded is True:
            strand = self.strand
        else:
            strand = "mixed"
        return "{0}:{1}{2}:{3}-{4}".format(
            self.__name__,
            self.chrom,
            strand,
            self.start,
            self.end)

    @property
    def lost_transcripts(self):
        return self.__lost.copy()
