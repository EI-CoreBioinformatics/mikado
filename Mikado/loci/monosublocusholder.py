# coding: utf-8

"""
This module defines a holder of Monosublocus instances. It is the last step
before the definition of real loci.
"""

import itertools
from sys import version_info
from ..transcripts.transcript import Transcript
from .abstractlocus import Abstractlocus
from .locus import Locus
from .monosublocus import Monosublocus
from .sublocus import Sublocus
from ..parsers.GFF import GffLine
from ..utilities import overlap
if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict

# Resolution order is important here!
# pylint: disable=too-many-instance-attributes
class MonosublocusHolder(Sublocus, Abstractlocus):
    """This is a container that groups together the transcripts
    surviving the selection for the Monosublocus.
    The class inherits from both sublocus and Abstractlocus
    (the main abstract class) in order to be able to reuse
    some of the code present in the former.
    Internally, the most important method is define_loci -
    which will select the best transcript(s) and remove all the overlapping ones.
    The intersection function for this object is quite laxer than in previous stages,
    and so are the requirements for the inclusion.
    """

    __name__ = "monosubloci_holder"

    # pylint: disable=super-init-not-called
    def __init__(self, transcript_instance=None,
                 json_conf=None, logger=None,
                 verified_introns=None, **kwargs):

        # I know what I am doing by NOT calling the Sublocus super but rather
        # Abstractlocus
        Abstractlocus.__init__(self,
                               transcript_instance=None,
                               verified_introns=verified_introns,
                               json_conf=json_conf,
                               logger=logger,
                               **kwargs)
        self._not_passing = set()
        self.splitted = False
        self.metrics_calculated = False
        self.excluded = None
        self.feature = "MonosublocusHolder"
        # Add the transcript to the Locus
        self.locus_verified_introns = set()
        if transcript_instance is not None:
            self.add_monosublocus(transcript_instance, check_in_locus=False)
        self.loci = SortedDict()
        self.attributes = dict()

    # Overriding is correct here
    # pylint: disable=arguments-differ
    def add_transcript_to_locus(self, transcript, check_in_locus=True):
        """Override of the sublocus method, and reversal to the original
        method in the Abstractlocus class.
        The check_in_locus boolean flag is used to decide
        whether to check if the transcript is in the Locus or not.
        This should be set to False for the first transcript, and True afterwards.

        :param transcript: a Transcript instance
        :type transcript: Transcript

        :param check_in_locus: optional flag to pass to the basic method.
        If set to False it disables checks on whether the transcript
        is really in the locus
        :type check_in_locus: bool
        """

        #
        # monosublocus: Abstractlocus,
        # transcript: Transcript,
        # flank = 0,
        # logger = None,
        # cds_only = False,
        # min_cdna_overlap = 0.2,
        # min_cds_overlap = 0.2,
        # classic_method = False

        if check_in_locus is True and self.in_locus(
                self,
                transcript,
                flank=self.json_conf["pick"]["clustering"]["flank"],
                logger=self.logger,
                cds_only=self.json_conf["pick"]["clustering"]["cds_only"],
                min_cdna_overlap=self.json_conf["pick"]["clustering"]["min_cdna_overlap"],
                min_cds_overlap=self.json_conf["pick"]["clustering"]["min_cds_overlap"],
                simple_overlap_for_monoexonic=self.json_conf["pick"]["clustering"]["simple_overlap_for_monoexonic"]
        ) is False:

                self.logger.debug("%s is not a valid intersection for %s", transcript.id, self.id)
                return False

        Abstractlocus.add_transcript_to_locus(self, transcript, check_in_locus=False)
        self.locus_verified_introns = set.union(self.locus_verified_introns,
                                                transcript.verified_introns)

    # pylint: enable=arguments-differ
    @classmethod
    def is_intersecting(cls,
                        transcript,
                        other,
                        cds_only=False,
                        logger=None,
                        min_cdna_overlap=0.2,
                        min_cds_overlap=0.2,
                        simple_overlap_for_monoexonic=True) -> bool:
        return Sublocus.is_intersecting(transcript, other, cds_only=cds_only,
                                        logger=logger, min_cdna_overlap=min_cdna_overlap,
                                        min_cds_overlap=min_cds_overlap,
                                        simple_overlap_for_monoexonic=simple_overlap_for_monoexonic)

    def add_monosublocus(self, monosublocus_instance: Monosublocus,
                         check_in_locus=True):
        """Wrapper to extract the transcript from the monosubloci and pass it to the constructor.

        :param monosublocus_instance
        :type monosublocus_instance: Monosublocus
        """
        assert len(monosublocus_instance.transcripts) == 1
        # if len(self.transcripts) == 0:
        #     check_in_locus = False
        # else:
        #     check_in_locus = True
        for tid in monosublocus_instance.transcripts:
            self.add_transcript_to_locus(monosublocus_instance.transcripts[tid],
                                         check_in_locus=check_in_locus)

    def __str__(self, print_cds=False, source_in_name=True):
        """This special method is explicitly *not* implemented;
        this Locus object is not meant for printing, only for computation!

        :param print_cds: flag. Ignored.
        """

        lines = []

        self_line = GffLine('')
        for attr in ["chrom", 'feature', 'source', 'start', 'end', 'strand']:
            setattr(self_line, attr, getattr(self, attr))
        self.filter_and_calculate_scores()

        self_line.phase, self_line.score = None, self.score
        if source_in_name is True:
            self_line.id = "{0}_{1}".format(self.source, self.id)
        else:
            self_line.id = self.id
        self_line.name = self.name
        self_line.parent = self.parent
        self_line.attributes.update(self.attributes)
        self_line.attributes["multiexonic"] = (not self.monoexonic)
        lines.append(str(self_line))

        for tid in self.transcripts:
            transcript_instance = self.transcripts[tid]
            transcript_instance.source = self.source
            transcript_instance.parent = self_line.id
            self.logger.debug(self.attributes)
            for attribute in self.attributes:
                if attribute not in transcript_instance.attributes:
                    if attribute == "is_fragment" and self.attributes[attribute] is False:
                        continue
                    transcript_instance.attributes[attribute] = self.attributes[attribute]

            lines.append(transcript_instance.format(
                "gff",
                all_orfs=self.json_conf["pick"]["output_format"]["report_all_orfs"],
                with_cds=print_cds).rstrip())

        return "\n".join(lines)

    def define_monosubloci(self, **kwargs):
        """Overriden and set to NotImplemented to avoid cross-calling it when inappropriate."""
        raise NotImplementedError("Monosubloci are the input of this object, not the output.")

    def define_loci(self, purge=False, excluded=None):
        """This is the main function of the class. It is analogous
        to the define_subloci class defined for sublocus objects,
        but it returns "Locus" objects (not "Monosublocus").

        Optional parameters:

        :param purge: flag. If set to True, all loci whose transcripts
        have scores of 0 will be thrown out
        into an Excluded holder.
        :type purge: bool

        :param excluded
        :type excluded: Excluded

        """
        if self.splitted is True:
            return

        self.excluded = excluded

        self.filter_and_calculate_scores()

        graph = self.define_graph(
            self.transcripts,
            inters=self.is_intersecting,
            logger=self.logger,
            cds_only=self.json_conf["pick"]["clustering"]["cds_only"],
            min_cdna_overlap=self.json_conf["pick"]["clustering"]["min_cdna_overlap"],
            min_cds_overlap=self.json_conf["pick"]["clustering"]["min_cds_overlap"],
            simple_overlap_for_monoexonic=self.json_conf["pick"]["clustering"]["simple_overlap_for_monoexonic"]
        )

        loci = []
        while len(graph) > 0:
            # cliques = self.find_cliques(graph)
            communities = self.find_communities(graph)
            to_remove = set()
            for locus_comm in communities:
                locus_comm = dict((x, self.transcripts[x]) for x in locus_comm)
                selected_tid = self.choose_best(locus_comm)
                selected_transcript = self.transcripts[selected_tid]
                to_remove.add(selected_tid)
                to_remove.update(set(graph.neighbors(selected_tid)))
                if purge is False or selected_transcript.score > 0:
                    new_locus = Locus(selected_transcript, logger=self.logger, json_conf=self.json_conf,
                                      use_transcript_scores=self._use_transcript_scores)
                    loci.append(new_locus)
            self.logger.debug("Removing {0} transcripts from {1}".format(len(to_remove), self.id))
            graph.remove_nodes_from(to_remove)  # Remove nodes from graph, iterate

        for locus in sorted(loci):
            self.loci[locus.id] = locus
            if self.regressor is not None:
                self.loci[locus.id].regressor = self.regressor
        self.splitted = True
        return

    @classmethod
    def in_locus(cls, monosublocus: Abstractlocus,
                 transcript: Transcript,
                 flank=0,
                 logger=None,
                 cds_only=False,
                 min_cdna_overlap=0.2,
                 min_cds_overlap=0.2,
                 simple_overlap_for_monoexonic=False) -> bool:

        """This method checks whether a transcript / monosbulocus
        falls inside the Locus coordinates.

        :rtype: bool

        :param monosublocus: Monosublocus instance
        :type monosublocus: Monosublocus

        :param transcript: the transcript to be compared
        :type transcript: Transcript

        :param flank: optional flank argument
        :type flank: int
        """

        if hasattr(transcript, "transcripts"):
            assert len(transcript.transcripts) == 1
            transcript = transcript.transcripts[list(transcript.transcripts.keys())[0]]
            assert hasattr(transcript, "finalize")
        if Abstractlocus.in_locus(monosublocus, transcript, flank=flank) is True:
            is_in_locus = False
            for tran in monosublocus.transcripts:
                tran = monosublocus.transcripts[tran]
                is_in_locus = cls.is_intersecting(tran,
                                                  transcript,
                                                  logger=logger,
                                                  cds_only=cds_only,
                                                  min_cds_overlap=min_cds_overlap,
                                                  min_cdna_overlap=min_cdna_overlap,
                                                  simple_overlap_for_monoexonic=simple_overlap_for_monoexonic
                                                  )
                if is_in_locus is True:
                    break
            return is_in_locus
        else:
            return False

    @property
    def id(self):
        """
        Wrapper for the id method of abstractlocus. Necessary to redefine the name.
        """

        return Abstractlocus.id.fget(self)  # @UndefinedVariable

    def as_dict(self):
        state = Abstractlocus.as_dict(self)
        return state

    def load_dict(self, state):
        Abstractlocus.load_dict(self, state)
