# coding: utf-8

"""
This class is used to store all the data calculated by the Assigner class, in order to produce
the RefMap/Stats files.
"""

import argparse
import collections
import logging
import operator
from logging import handlers as log_handlers

from Mikado.transcripts.transcript import Transcript
from . import calc_f1
from .resultstorer import ResultStorer


# noinspection PyPropertyAccess,PyPropertyAccess,PyPropertyAccess
# pylint: disable=too-many-instance-attributes
class Accountant:
    """This class stores the data necessary to calculate the final statistics
     - base and exon Sn/Sp/F1 etc."""

    def __init__(self, genes: dict, args: argparse.Namespace):

        """Class constructor. It requires:
        :param genes: a dictionary
        :type genes: dict
        :param args: A namespace (like those provided by argparse)
        containing the parameters for the run.
        """

        self.args = args
        self.queue_handler = None
        self.logger = None
        self.__setup_logger()
        self.logger.debug("Started with stat printing")

        self.introns = dict()
        self.exons = dict()
        self.starts = dict()
        self.ends = dict()
        self.intron_chains = collections.Counter()
        self.monoexonic_matches = (set(), set())
        self.ref_genes = dict()
        self.pred_genes = dict()
        self.__setup_reference_data(genes)
        self.self_analysis = False
        if hasattr(args, "self") and args.self is True:
            self.self_analysis = True

    def __setup_reference_data(self, genes):

        """
        Private method that prepares the reference data into the data structure
        that will be used to compare each prediction with a reference transcript.
        """

        for gene in genes:
            self.ref_genes[gene] = dict()
            for transcr in genes[gene]:
                # 0b000, indicating:
                # First bit: stringent match (100% F1)
                # Second bit: normal match (95% F1)
                # Third bit: lenient match (80% F1)
                self.ref_genes[gene][transcr.id] = 0b1000
                #                 self.ref_transcript_num+=1
                if transcr.chrom not in self.introns:
                    self.exons[transcr.chrom] = dict([("+", dict()), ("-", dict())])
                    self.starts[transcr.chrom] = dict([("+", dict()), ("-", dict())])
                    self.ends[transcr.chrom] = dict([("+", dict()), ("-", dict())])

                    self.introns[transcr.chrom] = dict([("+", dict()), ("-", dict())])
                    self.intron_chains[transcr.chrom] = dict([("+", dict()), ("-", dict())])
                if transcr.strand is None:
                    strand = "+"
                else:
                    strand = transcr.strand

                # 0b000001: in reference
                # 0b000010: in prediction
                # 0b000100: single
                # 0b001000: internal
                # 0b010000: border
                # 0b100000: single match lenient

                if transcr.exon_num == 1:
                    exon = tuple([transcr.exons[0][0], transcr.exons[0][1]])
                    self.exons[transcr.chrom][strand][exon] = 0b00000
                    self.exons[transcr.chrom][strand][exon] |= 0b1
                    self.exons[transcr.chrom][strand][exon] |= 0b100
                else:
                    self.__store_multiexonic_reference(transcr, strand)
        return

    def __store_multiexonic_reference(self, transcr, strand):

        """
        Private method to store a multiexonic reference transcript.
        :param transcr:
        :param strand:
        :return:
        """

        intron_chain = []
        for intron in sorted(transcr.introns):
            intron = tuple([intron[0], intron[1]])
            intron_chain.append(intron)
            self.introns[transcr.chrom][strand][intron] = 0b01
        intron_chain = tuple(intron_chain)
        if intron_chain not in self.intron_chains[transcr.chrom][strand]:
            self.intron_chains[transcr.chrom][
                strand][intron_chain] = [set(), set()]
        self.intron_chains[transcr.chrom][
            strand][intron_chain][0].add(transcr.id)

        for index, exon in enumerate(transcr.exons):
            exon = tuple([exon[0], exon[1]])
            if exon not in self.exons[transcr.chrom][strand]:
                self.exons[transcr.chrom][strand][exon] = 0b00000
            self.exons[transcr.chrom][strand][exon] |= 0b01
            if index == 0:
                self.exons[transcr.chrom][strand][exon] |= 0b010000
                self.starts[transcr.chrom][strand][exon[1]] = 0b01
            elif index == transcr.exon_num - 1:
                self.exons[transcr.chrom][strand][exon] |= 0b010000
                self.ends[transcr.chrom][strand][exon[0]] = 0b01
            else:
                self.exons[transcr.chrom][strand][exon] |= 0b01000

    def __setup_logger(self):

        """
        Private method to set up the logger using indications in the
        args namespace.
        """

        if hasattr(self.args, "log_queue"):
            # noinspection PyUnresolvedReferences
            self.queue_handler = log_handlers.QueueHandler(self.args.log_queue)
        else:
            self.queue_handler = logging.NullHandler
        self.logger = logging.getLogger("stat_logger")
        self.logger.addHandler(self.queue_handler)
        # noinspection PyUnresolvedReferences
        if self.args.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        self.logger.propagate = False
        return

    # pylint: disable=too-many-locals
    def __retrieve_gene_stats(self, store_name="ref"):

        """
        This private method recovers from the pred/ref store the statistics regarding the
        number of genes and transcripts found and missed,
        both stringently and leniently.
        The keyword store_name, which must evaluate to either
        "ref"
        or
        "pred",
        indicates whether we are intersted in prediction of reference data.
        :param store_name: str
        :return: { "stringent": (transcript_stringent, gene_stringent),
                   "lenient": (transcript_lenient, gene_lenient),
                   "total": total_transcripts,
                   "private": (private_transcripts, private_genes)
                   }
        """

        found_transcripts_stringent = 0  # 100% match
        found_transcripts = 0  # 95% match
        found_transcripts_lenient = 0  # 80% match
        found_genes_stringent = 0  # 100% match
        found_genes = 0  # 95% match
        found_genes_lenient = 0  # 80% match
        private_genes = 0
        private_transcripts = 0
        total_transcripts = 0

        if store_name == "ref":
            store = self.ref_genes
        elif store_name == "pred":
            store = self.pred_genes
        else:
            raise ValueError("Invalid store selected: {0}".format(store_name))

        # 0b000, indicating:
        # First bit: stringent match (100% F1)
        # Second bit: normal match (95% F1)
        # Third bit: lenient match (80% F1)
        # Fourth bit: private

        for gene in store:
            gene_match = 0b000
            gene_not_found = 0b1
            for _, val in store[gene].items():
                total_transcripts += 1
                found_transcripts_lenient += (0b100 & val) >> 2
                found_transcripts += (0b10 & val) >> 1
                found_transcripts_stringent += 0b1 & val
                # Will evaluate to 1 if the transcript has at least one match
                private_transcripts += (0b1000 & val) >> 3
                gene_not_found &= (0b1000 & val) >> 3
                gene_match |= val

            found_genes_stringent += gene_match & 0b1
            # Shift back by 1
            found_genes += (gene_match & 0b10) >> 1
            # Shift back by 2
            found_genes_lenient += (gene_match & 0b100) >> 2
            private_genes += gene_not_found

        self.logger.debug("""Found %s transcripts:
        \tstringent\t%s
        \tlenient\t%s
        \ttotal\t%s""",
                          store_name,
                          found_transcripts_stringent,
                          found_transcripts_lenient,
                          total_transcripts)

        result = dict()
        result["stringent"] = (found_transcripts_stringent, found_genes_stringent)
        result["standard"] = (found_transcripts, found_genes)
        result["lenient"] = (found_transcripts_lenient, found_genes_lenient)
        result["total"] = total_transcripts
        result["private"] = (private_transcripts, private_genes)

        return result
    # pylint: enable=too-many-locals

    def __calculate_gene_stats(self):

        """
         Private method that calculates the number of found transcripts
         in the reference and prediction sets. Wrapper around a double
         call to __retrieve_gene_stats
         :return: double the result of __retrieve_gene_stats
         """

        ref_results = self.__retrieve_gene_stats(store_name="ref")

        pred_results = self.__retrieve_gene_stats(store_name="pred")

        result = {"ref": ref_results,
                  "pred": pred_results}

        return result

    def __calculate_intron_stats(self):

        """
        This private method calculates the raw numbers regarding intron and intron chain statistics.
        :return: result_dictionary
        :rtype: dict
        """

        intron_stats = [0, 0, 0]  # Common, prediction, reference

        for chrom in self.introns:
            for strand in self.introns[chrom]:
                for intron in self.introns[chrom][strand]:
                    intron_stats[0] += (0b01 & self.introns[chrom][strand][intron]) & \
                                     ((0b10 & self.introns[chrom][strand][intron]) >> 1)
                    intron_stats[1] += (0b10 & self.introns[chrom][strand][intron]) >> 1
                    intron_stats[2] += 0b01 & self.introns[chrom][strand][intron]

        # Common, prediction, reference
        intron_chains_nonred = [0, 0, 0]
        intron_chains_common_ref = set()
        intron_chains_common_pred = set()
        intron_chains_pred = 0
        intron_chains_ref = 0

        for chrom in self.intron_chains:
            for strand in self.intron_chains[chrom]:
                for _, intron_val in self.intron_chains[chrom][strand].items():
                    if len(intron_val[0]) > 0 and len(intron_val[1]) > 0:
                        intron_chains_nonred[0] += 1
                        intron_chains_common_ref.update(intron_val[0])
                        intron_chains_common_pred.update(intron_val[1])
                    # In reference
                    if len(intron_val[0]) > 0:
                        intron_chains_nonred[2] += 1
                        # intron_chains_ref_nonred += 1
                        intron_chains_ref += len(intron_val[0])
                    # In prediction
                    if len(intron_val[1]) > 0:
                        intron_chains_nonred[1] += 1
                        # intron_chains_pred_nonred += 1
                        intron_chains_pred += len(intron_val[1])

        intron_chains_common_ref = len(intron_chains_common_ref)
        intron_chains_common_pred = len(intron_chains_common_pred)

        self.logger.debug("""Intron chains:\
        \treference %d\t%d
        \tprediction\t%d
        \t%d
        \tcommon\t%d %d %d""",
                          intron_chains_ref,
                          intron_chains_nonred[2],
                          intron_chains_pred,
                          intron_chains_nonred[1],
                          intron_chains_nonred[0],
                          intron_chains_common_ref,
                          intron_chains_common_pred)

        result_dictionary = dict()
        result_dictionary["introns"] = intron_stats
        result_dictionary["intron_chains"] = dict()
        result_dictionary["intron_chains"]["non_redundant"] = intron_chains_nonred

        result_dictionary["intron_chains"]["redundant"] = [intron_chains_common_pred,
                                                           intron_chains_common_ref]
        return result_dictionary

    def __extract_terminal_stats(self, result_dictionary):
        """
        Private method to
        :param result_dictionary: the results dictionary
        :type result_dictionary: dict
        :return: updated result dictionary
        :rtype: dict
        """

        # 0b000001: in reference
        # 0b000010: in prediction
        # 0b000100: single
        # 0b001000: internal
        # 0b010000: border
        # 0b100000: single match lenient

        exon_common_lenient = 0
        exon_ref_lenient = 0
        exon_pred_lenient = 0

        for chrom in self.starts:
            for strand in self.starts[chrom]:
                for start in self.starts[chrom][strand]:
                    exon_common_lenient += (0b1 & self.starts[chrom][strand][start]) & \
                                           ((0b10 & self.starts[chrom][strand][start]) >> 1)
                    exon_ref_lenient += 0b01 & self.starts[chrom][strand][start]
                    exon_pred_lenient += (0b10 & self.starts[chrom][strand][start]) >> 1
        starts_common = exon_common_lenient
        starts_ref = exon_ref_lenient
        starts_pred = exon_pred_lenient

        self.logger.debug("Starts %s", [starts_common, starts_ref, starts_pred])

        for chrom in self.ends:
            for strand in self.ends[chrom]:
                for end in self.ends[chrom][strand]:
                    exon_common_lenient += (0b01 & self.ends[chrom][strand][end]) & \
                                           ((0b10 & self.ends[chrom][strand][end]) >> 1)
                    exon_ref_lenient += 0b01 & self.ends[chrom][strand][end]
                    exon_pred_lenient += (0b10 & self.ends[chrom][strand][end]) >> 1

        ends_common = exon_common_lenient - starts_common
        ends_ref = exon_ref_lenient - starts_ref
        ends_pred = exon_pred_lenient - starts_pred
        self.logger.debug("Ends %s", [ends_common, ends_ref, ends_pred])
        result_dictionary["ends"] = [ends_common, ends_pred, ends_ref]
        result_dictionary["starts"] = [starts_common, starts_pred, starts_ref]
        result_dictionary["exons"] = dict()
        result_dictionary["exons"]["lenient"] = [exon_common_lenient,
                                                 exon_pred_lenient,
                                                 exon_ref_lenient]

        return result_dictionary

    def __extract_internal_stats(self, result_dictionary):
        """
        Private method that extracts the internal stats from the stored
        data and updates the result_dictionary dictionary.
        :param result_dictionary: the dictionary that stores the results
        :return:
        """

        # Re-extract the previous stats
        # Common, prediction, reference
        exon_lenient = result_dictionary["exons"]["lenient"]

        # Common, prediction, reference
        bases = [0, 0, 0]

        # Common, prediction, reference
        exon_stringent = [0, 0, 0]

        for chrom in self.exons:
            for strand in self.exons[chrom]:
                # 0b000001: in reference
                # 0b000010: in prediction
                # 0b000100: single
                # 0b001000: internal
                # 0b010000: border
                # 0b100000: single match lenient
                curr_pred_bases = set()
                curr_ref_bases = set()
                curr_exon = (-1, -1)

                for exon in sorted(self.exons[chrom][strand], key=operator.itemgetter(0)):
                    # Out of previous overlap
                    if exon[0] > curr_exon[1]:
                        bases[0] += len(set.intersection(curr_pred_bases, curr_ref_bases))
                        bases[1] += len(curr_pred_bases)
                        bases[2] += len(curr_ref_bases)
                        curr_exon = exon
                        curr_pred_bases = set()
                        curr_ref_bases = set()
                    else:
                        curr_exon = (curr_exon[0], exon[1])

                    # Exon is in reference
                    if (0b01 & self.exons[chrom][strand][exon]) == 0b1:
                        curr_ref_bases.update(set(range(exon[0], exon[1])))
                        exon_stringent[2] += 1
                        # Internal (first condition)
                        # OR
                        # Single exon
                        if (0b001000 & self.exons[chrom][strand][exon]) == 0b001000:
                            exon_lenient[2] += 1
                        elif 0b100000 & self.exons[chrom][strand][exon]:
                            exon_lenient[2] += 1

                    # Exon is in prediction
                    if (0b10 & self.exons[chrom][strand][exon]) == 0b10:
                        curr_pred_bases.update(set(range(exon[0], exon[1])))
                        exon_stringent[1] += 1
                        # Either internal, or a single exon which
                        # has not a lenient single match with the reference annotation
                        if (0b001000 & self.exons[chrom][strand][exon]) == 0b001000:
                            exon_lenient[1] += 1
                        elif 0b100000 & self.exons[chrom][strand][exon]:
                            exon_lenient[1] += 1

                    # In both
                    if (0b11 & self.exons[chrom][strand][exon]) == 0b11:
                        exon_stringent[0] += 1
                        if (0b001000 & self.exons[chrom][strand][exon]) == 0b001000:
                            exon_lenient[0] += 1
                        elif 0b100000 & self.exons[chrom][strand][exon]:
                            exon_lenient[0] += 1

                bases[0] += len(set.intersection(curr_pred_bases, curr_ref_bases))
                bases[1] += len(curr_pred_bases)
                bases[2] += len(curr_ref_bases)

        # result_dictionary["bases"] = dict()

        assert bases[0] <= min(bases[1], bases[2]), bases
        result_dictionary["bases"] = bases
        result_dictionary["exons"]["stringent"] = exon_stringent
        result_dictionary["exons"]["lenient"] = exon_lenient
        return result_dictionary

    def __calculate_exon_stats(self):

        """
        This private method calculates the raw numbers regarding base and exon statistics.
        :return:
        """

        result = dict()
        result = self.__extract_terminal_stats(result)
        result = self.__extract_internal_stats(result)
        return result

    def __store_multiexonic_result(self, transcr, strand, result):

        """
        Private procedure to store the results regarding a multiexonic
        transcript compared to the reference.
        :param transcr: a transcript instance
        :param strand: the strand to be used for storing
        :param result: a ResultStorer instance
        :return:
        """

        # assert transcr.monoexonic is False
        ic_key = tuple([tuple([intron[0], intron[1]]) for intron in sorted(transcr.introns)])
        # if result.ccode == ("=",):
        #     assert ic_key in self.intron_chains[transcr.chrom][strand]
        #     assert result.ref_id[0] in self.intron_chains[transcr.chrom][strand][ic_key][0]

        if ic_key not in self.intron_chains[transcr.chrom][strand]:
            self.intron_chains[transcr.chrom][strand][ic_key] = [set(), set()]
        self.intron_chains[transcr.chrom][strand][ic_key][1].add(transcr.id)

        for intron in transcr.introns:
            intron = tuple([intron[0], intron[1]])
            if intron not in self.introns[transcr.chrom][strand]:
                self.introns[transcr.chrom][strand][intron] = 0b0
            self.introns[transcr.chrom][strand][intron] |= 0b10

        for index, exon in enumerate(transcr.exons):
            exon = tuple([exon[0], exon[1]])
            if exon not in self.exons[transcr.chrom][strand]:
                self.exons[transcr.chrom][strand][exon] = 0b0
            self.exons[transcr.chrom][strand][exon] |= 0b10  # set it as "in prediction"
            if index == 0:
                self.exons[transcr.chrom][strand][exon] |= 0b010000
                if exon[1] not in self.starts[transcr.chrom][strand]:
                    self.starts[transcr.chrom][strand][exon[1]] = 0b0
                self.starts[transcr.chrom][strand][exon[1]] |= 0b10
            elif index == transcr.exon_num - 1:
                self.exons[transcr.chrom][strand][exon] |= 0b010000
                if exon[0] not in self.ends[transcr.chrom][strand]:
                    self.ends[transcr.chrom][strand][exon[0]] = 0b00
                self.ends[transcr.chrom][strand][exon[0]] |= 0b10
            else:
                self.exons[transcr.chrom][strand][exon] |= 0b01000

    def __store_monoexonic_result(self, transcr, strand, result: ResultStorer, other_exon=None):
        """
        Private procedure to store the comparison of a monoexonic
        transcript with the reference.
        :param transcr: the monoexonic transcript
        :param strand: the strand to be used for storing
        :param other_exon: the reference exon this transcript is assigned to
        :return:
        """
        assert transcr.monoexonic is True
        exon = tuple([transcr.exons[0][0], transcr.exons[0][1]])
        if exon not in self.exons[transcr.chrom][strand]:
            self.exons[transcr.chrom][strand][exon] = 0b0
        self.exons[transcr.chrom][strand][exon] |= 0b10
        self.exons[transcr.chrom][strand][exon] |= 0b100
        if result.ccode == ("_",):
            self.monoexonic_matches[0].add(result.ref_id)
            self.monoexonic_matches[1].add(transcr.id)

        if other_exon is not None:
            # other_exon = tuple([other_exon[0], other_exon[1]])
            assert isinstance(other_exon, tuple)
            assert other_exon in self.exons[transcr.chrom][strand], (transcr.id,
                                                                     transcr.exons,
                                                                     other_exon)
            if not 0b100 & self.exons[transcr.chrom][strand][other_exon]:
                # Check the other exon is marked as single
                self.exons[transcr.chrom][strand][other_exon] |= 0b100

            # self.exons[transcr.chrom][strand][exon] |= 0b100000
            self.exons[transcr.chrom][strand][other_exon] |= 0b100000

    def store(self, transcr: Transcript, result: ResultStorer, other_exon):

        """Add exons introns intron chains etc. to the storage.
        :param transcr: a "Transcript" instance
        :param result: Instance of "result_storer"
        :param other_exon: either None or an integer duplex
        :type other_exon: None | (int,int)

        A transcript is considered a perfect match if it has:
            junction_f1==100 and nucleotide_f1==100;
        a lenient match if it has junction_f1==100 and nucleotide_f1>95,
        i.e. min(nucleotide_precision, nucleotide_recall)>90.4
        """

        for parent in transcr.parent:
            if parent not in self.pred_genes:
                self.pred_genes[parent] = dict()
            if transcr.id not in self.pred_genes[parent]:
                self.pred_genes[parent][transcr.id] = 0b1000

        if transcr.strand is None:
            strand = "+"
        else:
            strand = transcr.strand

        if result.ccode != ("u",):
            found_one = False

            for refid, refgene, nf1 in zip(result.ref_id, result.ref_gene, result.n_f1):
                # self.ref_genes[refgene][refid] |= 0b1000
                if nf1 > 0:
                    self.ref_genes[refgene][refid] &= ~(1 << 3)  # Clear the fourth bit
                    found_one = True
            if found_one:
                for parent in transcr.parent:
                    self.pred_genes[parent][transcr.id] &= ~(1 << 3)  # Clear the fourth bit

            if len(result.j_f1) == 1:
                zipper = zip(result.ref_id, result.ref_gene, result.j_f1, result.n_f1)
                for refid, refgene, junc_f1, nucl_f1 in zipper:
                    if junc_f1 == 100 and nucl_f1 >= 80:
                        for parent in transcr.parent:
                            self.pred_genes[parent][transcr.id] |= 0b100
                            if nucl_f1 >= 95:
                                self.pred_genes[parent][transcr.id] |= 0b10
                            if nucl_f1 == 100:
                                self.pred_genes[parent][transcr.id] |= 0b01
                        # Unset the "private" mark

                        self.ref_genes[refgene][refid] |= 0b100
                        if nucl_f1 >= 95:
                            self.ref_genes[refgene][refid] |= 0b10
                        if nucl_f1 == 100:
                            self.ref_genes[refgene][refid] |= 0b01
        else:
            for parent in transcr.parent:
                self.pred_genes[parent][transcr.id] &= 0b1000

        if transcr.chrom not in self.exons:
            self.exons[transcr.chrom] = dict([("+", dict()), ("-", dict())])
            self.starts[transcr.chrom] = dict([("+", dict()), ("-", dict())])
            self.ends[transcr.chrom] = dict([("+", dict()), ("-", dict())])
            self.introns[transcr.chrom] = dict([("+", dict()), ("-", dict())])
            self.intron_chains[transcr.chrom] = dict([("+", dict()), ("-", dict())])

        if transcr.exon_num > 1:
            self.__store_multiexonic_result(transcr, strand, result)
        else:
            self.__store_monoexonic_result(transcr, strand, result, other_exon=other_exon)

    @staticmethod
    def __calculate_statistics(common, pred, ref):
        """
        Function to calculate precision, recall and F1 given
        the numbers in common, of the prediction, and of the reference
        :param common:
        :param pred:
        :param ref:
        :return: precision, recall, f1stat
        """

        if ref > 0:
            recall = common / ref
        else:
            recall = 0
        if pred > 0:
            precision = common / pred
        else:
            precision = 0
        f1stat = calc_f1(recall, precision)

        return precision, recall, f1stat

    @staticmethod
    def __redundant_stats(common_pred, common_ref, pred, ref):

        """
        Function to calculate precision, recall and F1 given the numbers of:
        - predictions which are matches
        - reference which are matches
        - total predictions
        - total references

        :param common_pred: number of predictions which are a match
        :param common_ref: number of references which are a match
        :param pred: total number of predictions
        :param ref: total number of references
        :return:
        """

        if ref > 0:
            recall = common_ref / ref
        else:
            recall = 0

        if pred > 0:
            precision = common_pred / pred
        else:
            precision = 0

        f1stat = calc_f1(recall, precision)

        return precision, recall, f1stat

    @staticmethod
    def __format_rowname(stru):
        """
        Private method to format the name of the rows in the stats table
        :param stru:
        :return:
        """

        total_length = 35
        return " " * (total_length-len(stru)-1) + "{}:".format(stru.rstrip(":"))

    @classmethod
    def __format_comparison_line(cls, name, private, total):
        """
        Method to create the string that will be printed out to screen for the stats file.
        :return:
        """

        if total > 0:
            perc = 100 * private / total
        else:
            perc = 0

        string = "{0} {1}/{2}  ({3:.2f}%)".format(cls.__format_rowname(name),
                                                  private,
                                                  total,
                                                  perc)
        return string

    @classmethod
    def __format_comparison_segments(cls, name, private, common):
        """
        Method to create the string that will be printed out to screen for the stats file.
        :return:
        """

        if common > 0:
            perc = 100 * (private - common) / private
        else:
            perc = 0

        string = "{0} {1}/{2}  ({3:.2f}%)".format(cls.__format_rowname(name),
                                                  private-common,
                                                  private,
                                                  perc)
        return string

    def print_stats(self):

        """
        This method calculates and prints the final stats file.
        It is called when the prediction file
        parsing has been terminated correctly.
        """

        # num_ref_transcripts = sum(len(self.ref_genes[gid]) for gid in self.ref_genes)
        # num_pred_transcripts = sum(len(self.pred_genes[gid]) for gid in self.pred_genes)

        # Dictionary with the necessary data
        bases_exon_result = self.__calculate_exon_stats()
        intron_results = self.__calculate_intron_stats()
        gene_transcript_results = self.__calculate_gene_stats()

        # self.logger.debug("Missed exons lenient: %s", bases_exon_result["missed"])
        self.logger.debug("Exon stringent %s", bases_exon_result["exons"]["stringent"])
        self.logger.debug("Exon lenient %s", bases_exon_result["exons"]["lenient"])

        self.logger.debug(
            """Bases:
            \tcommon\t%d
            \tprediction\t%d
            \treference\t%d""",
            *bases_exon_result["bases"])

        bases_prec, bases_recall, bases_f1 = self.__calculate_statistics(
            *bases_exon_result["bases"])

        (exon_stringent_precision, exon_stringent_recall, exon_stringent_f1) = (
            self.__calculate_statistics(*bases_exon_result["exons"]["stringent"]))

        (exon_lenient_precision, exon_lenient_recall, exon_lenient_f1) = (
            self.__calculate_statistics(*bases_exon_result["exons"]["lenient"]))

        intron_precision, intron_recall, intron_f1 = self.__calculate_statistics(
            *intron_results["introns"]
        )

        intron_stats = self.__calculate_statistics(
            *intron_results["intron_chains"]["non_redundant"])

        intron_chains_precision, intron_chains_recall, intron_chains_f1 = intron_stats

        # Transcript level
        # noinspection PyTypeChecker
        (tr_precision_stringent, tr_recall_stringent, tr_f1_stringent) = (
            self.__redundant_stats(
                gene_transcript_results["pred"]["stringent"][0],
                gene_transcript_results["ref"]["stringent"][0],
                gene_transcript_results["pred"]["total"],
                gene_transcript_results["ref"]["total"]))

        (tr_precision_standard, tr_recall_standard, tr_f1_standard) = (
            self.__redundant_stats(
                gene_transcript_results["pred"]["standard"][0],
                gene_transcript_results["ref"]["standard"][0],
                gene_transcript_results["pred"]["total"],
                gene_transcript_results["ref"]["total"]))

        # noinspection PyTypeChecker
        (tr_precision_lenient, tr_recall_lenient, tr_f1_lenient) = (
            self.__redundant_stats(
                gene_transcript_results["pred"]["lenient"][0],
                gene_transcript_results["ref"]["lenient"][0],
                gene_transcript_results["pred"]["total"],
                gene_transcript_results["ref"]["total"]))

        # Gene level
        ref_genes = len(self.ref_genes)
        pred_genes = len(self.pred_genes)
        # noinspection PyTypeChecker
        (gene_precision_stringent, gene_recall_stringent, gene_f1_stringent) = (
            self.__redundant_stats(
                gene_transcript_results["pred"]["stringent"][1],
                gene_transcript_results["ref"]["stringent"][1],
                pred_genes,
                ref_genes))

        # noinspection PyTypeChecker
        (gene_precision_standard, gene_recall_standard, gene_f1_standard) = (
            self.__redundant_stats(
                gene_transcript_results["pred"]["standard"][1],
                gene_transcript_results["ref"]["standard"][1],
                pred_genes,
                ref_genes))

        # noinspection PyTypeChecker
        (gene_precision_lenient, gene_recall_lenient, gene_f1_lenient) = (
            self.__redundant_stats(
                gene_transcript_results["pred"]["lenient"][1],
                gene_transcript_results["ref"]["lenient"][1],
                pred_genes,
                ref_genes))

        # noinspection PyUnresolvedReferences
        with open("{0}.stats".format(self.args.out), 'wt') as out:

            # noinspection PyUnresolvedReferences
            print("Command line:\n{0:>10}".format(self.args.commandline), file=out)
            print(gene_transcript_results["ref"]["total"], "reference RNAs in",
                  len(self.ref_genes), "genes", file=out)
            print(gene_transcript_results["pred"]["total"], "predicted RNAs in ",
                  len(self.pred_genes), "genes", file=out)

            print("-" * 33, "|   Sn |   Pr |   F1 |", file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Base level"),
                bases_recall * 100,
                bases_prec * 100,
                bases_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Exon level (stringent)"),
                exon_stringent_recall * 100,
                exon_stringent_precision * 100,
                exon_stringent_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Exon level (lenient)"),
                exon_lenient_recall * 100,
                exon_lenient_precision * 100,
                exon_lenient_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Intron level"),
                intron_recall * 100,
                intron_precision * 100, intron_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Intron chain level"),
                intron_chains_recall * 100,
                intron_chains_precision * 100,
                intron_chains_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Transcript level (stringent)"),
                tr_recall_stringent * 100,
                tr_precision_stringent * 100,
                tr_f1_stringent * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Transcript level (>=95% base F1)"),
                tr_recall_standard * 100,
                tr_precision_standard * 100,
                tr_f1_standard * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Transcript level (>=80% base F1)"),
                tr_recall_lenient * 100,
                tr_precision_lenient * 100,
                tr_f1_lenient * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Gene level (100% base F1)"),
                gene_recall_stringent * 100,
                gene_precision_stringent * 100,
                gene_f1_stringent * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Gene level (>=95% base F1)"),
                gene_recall_standard * 100,
                gene_precision_standard * 100,
                gene_f1_standard * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Gene level (>=80% base F1)"),
                gene_recall_lenient * 100,
                gene_precision_lenient * 100,
                gene_f1_lenient * 100), file=out)
            print(file=out)
            print("#   Matching: in prediction; matched: in reference.", file=out)
            print("", file=out)
            print("{0} {1}".format(
                self.__format_rowname("Matching intron chains"),
                intron_results["intron_chains"]["redundant"][0]), file=out)
            print("{0} {1}".format(
                self.__format_rowname("Matched intron chains"),
                intron_results["intron_chains"]["redundant"][1]), file=out)
            print("{0} {1}".format(
                self.__format_rowname("Matching monoexonic transcripts"),
                len(self.monoexonic_matches[1])), file=out)
            print("{0} {1}".format(
                self.__format_rowname("Matched monoexonic transcripts"),
                len(self.monoexonic_matches[0])), file=out)
            print("{0} {1}".format(
                self.__format_rowname("Total matching transcripts"),
                intron_results["intron_chains"]["redundant"][0] + len(self.monoexonic_matches[1])),
                file=out)
            print("{0} {1}".format(
                self.__format_rowname("Total matched transcripts"),
                intron_results["intron_chains"]["redundant"][1] + len(self.monoexonic_matches[0])),
                file=out)

            print("", file=out)

            print(self.__format_comparison_segments("Missed exons (stringent)",
                                                    bases_exon_result["exons"]["stringent"][2],
                                                    bases_exon_result["exons"]["stringent"][0]),
                  file=out)
            print(self.__format_comparison_segments("Novel exons (stringent)",
                                                    bases_exon_result["exons"]["stringent"][1],
                                                    bases_exon_result["exons"]["stringent"][0]),
                  file=out)

            print(self.__format_comparison_segments("Missed exons (lenient)",
                                                    bases_exon_result["exons"]["lenient"][2],
                                                    bases_exon_result["exons"]["lenient"][0]),
                  file=out)
            print(self.__format_comparison_segments("Novel exons (lenient)",
                                                    bases_exon_result["exons"]["lenient"][1],
                                                    bases_exon_result["exons"]["lenient"][0]),
                  file=out)

            print(self.__format_comparison_segments("Missed introns",
                                                    intron_results["introns"][2],
                                                    intron_results["introns"][0]),
                  file=out)
            print(self.__format_comparison_segments("Novel introns",
                                                    intron_results["introns"][1],
                                                    intron_results["introns"][0]),
                  file=out)
            print("", file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Missed transcripts",
                                                gene_transcript_results["ref"]["private"][0],
                                                gene_transcript_results["ref"]["total"]),
                  file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Novel transcripts",
                                                gene_transcript_results["pred"]["private"][0],
                                                gene_transcript_results["pred"]["total"]),
                  file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Missed genes",
                                                gene_transcript_results["ref"]["private"][1],
                                                len(self.ref_genes)), file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Novel genes",
                                                gene_transcript_results["pred"]["private"][1],
                                                len(self.pred_genes)), file=out)

        self.logger.removeHandler(self.queue_handler)
        # self.queue_handler.close()
        return
