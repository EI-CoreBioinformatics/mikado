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
from ..transcripts import Transcript, Namespace
from ..utilities.f1 import calc_f1
from .resultstorer import ResultStorer
from ..utilities.intervaltree import Interval, IntervalTree
import networkx as nx
import numpy as np


# noinspection PyPropertyAccess,PyPropertyAccess,PyPropertyAccess
# pylint: disable=too-many-instance-attributes
class Accountant:
    """This class stores the data necessary to calculate the final statistics
     - base and exon Sn/Sp/F1 etc."""

    def __init__(self, genes, args: argparse.Namespace, counter=None, fuzzymatch=0, load_ref=False):

        """Class constructor. It requires:
        :param genes: a dictionary
        :type genes: [dict|GeneDict]
        :param args: A namespace (like those provided by argparse)
        containing the parameters for the run.
        """

        self.args = args
        self.queue_handler = None
        self.logger = None
        self._counter = counter
        self.__setup_logger()
        self.logger.debug("Started with stat printing")

        self.introns = dict()
        self.exons = dict()
        self.starts = dict()
        self.ends = dict()
        self.intron_chains = dict()
        self.monoexonic_matches = (set(), set())
        self.ref_genes = dict()
        self.pred_genes = dict()
        if load_ref is True:
            self.logger.info("Starting loading the reference for the accountant %s", self._counter)
            self.__setup_reference_data(genes)
        self.self_analysis = False
        self.__fuzzymatch = fuzzymatch
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
            if intron not in self.introns[transcr.chrom][strand]:
                self.introns[transcr.chrom][strand][intron] = [0, 0]
            self.introns[transcr.chrom][strand][intron][0] += 1
        intron_chain = tuple(intron_chain)
        if intron_chain not in self.intron_chains[transcr.chrom][strand]:
            self.intron_chains[transcr.chrom][strand][intron_chain] = [set(), set()]
        self.intron_chains[transcr.chrom][strand][intron_chain][0].add(transcr.id)

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
        if self._counter is None:
            self.logger = logging.getLogger("stat_logger")
        else:
            self.logger = logging.getLogger("stat_logger-{}".format(self._counter))
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

    def __calculate_chains_stats(self, chains):

        """
        This takes up a dictionary of the form (intron_chain: set(reference ids), set(prediction ids))

        :param chains:
        :return:
        """

        graph = nx.Graph()
        positions = IntervalTree()

        evaluator = dict()
        intron_chains_common_ref = set()
        intron_chains_common_pred = set()
        intron_chains_pred = 0
        intron_chains_ref = 0

        intron_chains_nonred = np.array([0, 0, 0])

        for chain in chains:
            ichain = IntervalTree.from_tuples(chain)
            evaluator[ichain] = chain
            graph.add_node(ichain)
            positions.insert(ichain.start, ichain.end, value=ichain)

        for ichain, chain in evaluator.items():
            found = positions.find(ichain.start, ichain.end)
            for ochain in found:
                if ichain.fuzzy_equal(ochain.value, self.__fuzzymatch):  # Define which chains are broadly similar
                    graph.add_edge(ichain, ochain.value)

        communities = nx.algorithms.connected_components(graph)
        for community in communities:  # Now we have to go back to the original sets ...
            ref_vals = set()
            pred_vals = set()
            for chain in community:
                ref_vals.update(chains[evaluator[chain]][0])
                pred_vals.update(chains[evaluator[chain]][1])
            if ref_vals and pred_vals:
                intron_chains_common_ref.update(ref_vals)
                intron_chains_common_pred.update(pred_vals)
                intron_chains_nonred[0] += 1

            if ref_vals:
                intron_chains_nonred[2] += 1
                intron_chains_ref += len(ref_vals)
            if pred_vals:
                intron_chains_nonred[1] += 1
                intron_chains_pred += len(pred_vals)

        intron_chains_common_ref = len(intron_chains_common_ref)
        intron_chains_common_pred = len(intron_chains_common_pred)
        intron_chains_red = np.array([intron_chains_common_pred, intron_chains_common_ref,
                                       intron_chains_pred, intron_chains_ref])

        return intron_chains_red, intron_chains_nonred

    def __calculate_intron_stats(self):

        """
        This private method calculates the raw numbers regarding intron and intron chain statistics.
        :return: result_dictionary
        :rtype: dict
        """

        intron_stats = [0, 0, 0, 0]  # Common - pred side, common - ref side, prediction, reference
        intron_nr_stats = [0, 0, 0]  # Common, prediction, reference

        splice_stats = [0, 0, 0]  # Common, prediction, reference

        # Common, prediction, reference
        intron_chains_red = np.array([0, 0, 0, 0])  # Common - pred side, common - ref side, prediction, reference
        intron_chains_nonred = np.array([0, 0, 0])

        for chrom in self.intron_chains:
            for strand in self.intron_chains[chrom]:
                splice_starts = [set(), set()]  # reference, prediction
                splice_ends = [set(), set()]
                intron_nrs = [collections.Counter(), collections.Counter()]
                for chain, intron_val in self.intron_chains[chrom][strand].items():  # Two sets, refs and preds
                    __ref_ic_val = len(intron_val[0])
                    __pred_ic_val = len(intron_val[1])
                    for intron in chain:
                        if __ref_ic_val:  # reference
                            splice_starts[0].add(intron[0])
                            splice_ends[0].add(intron[1])
                        if __pred_ic_val:  # prediction
                            splice_starts[1].add(intron[0])
                            splice_ends[1].add(intron[1])
                        intron_nrs[0][intron] += __ref_ic_val
                        intron_nrs[1][intron] += __pred_ic_val

                __red, __nonred = self.__calculate_chains_stats(self.intron_chains[chrom][strand])
                intron_chains_red += __red
                intron_chains_nonred += __nonred

                splice_stats[0] += len(set.intersection(*splice_starts)) + len(set.intersection(*splice_ends))
                splice_stats[1] += len(splice_starts[1]) + len(splice_ends[1])  # prediction
                splice_stats[2] += len(splice_starts[0]) + len(splice_ends[0])  # reference
                for key in set.union(set(intron_nrs[0].keys()), set(intron_nrs[1].keys())):
                    __ref_val, __pred_val = intron_nrs[0][key], intron_nrs[1][key]
                    if min(__ref_val, __pred_val) > 0:
                        intron_stats[1] += __ref_val
                        intron_stats[0] += __pred_val
                        intron_nr_stats[0] += 1

                    intron_stats[2] += __pred_val
                    intron_stats[3] += __ref_val
                    intron_nr_stats[1] += (__pred_val > 0)
                    intron_nr_stats[2] += (__ref_val > 0)

        self.logger.debug("Intron stats: {}".format(", ".join([str(_) for _ in intron_stats])))
        # intron_nr_stats = [len(set.intersection(*intron_nrs)), len(intron_nrs[1]), len(intron_nrs[0])]
        self.logger.debug("Intron NR stats: {}".format(", ".join([str(_) for _ in intron_nr_stats])))


        self.logger.debug("""Intron chains:\
        \treference %d\t%d
        \tprediction\t%d
        \t%d
        \tcommon\t%d %d %d""",
                          intron_chains_red[3],
                          intron_chains_nonred[2],
                          intron_chains_red[2],
                          intron_chains_nonred[1],
                          intron_chains_nonred[0],
                          intron_chains_red[1],
                          intron_chains_red[0])

        result_dictionary = dict()
        result_dictionary["introns"] = dict()
        result_dictionary["introns"]["non_redundant"] = intron_nr_stats
        result_dictionary["introns"]["redundant"] = intron_stats
        result_dictionary["splices"] = splice_stats
        result_dictionary["intron_chains"] = dict()
        result_dictionary["intron_chains"]["non_redundant"] = intron_chains_nonred

        result_dictionary["intron_chains"]["redundant"] = intron_chains_red
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

    def serialize(self):

        simplified = Namespace()
        simplified.attributes = ("pred_genes", "ref_genes", "intron_chains", "monoexonic_matches",
                                 "introns", "exons", "starts", "ends")

        for attr in simplified.attributes:
            setattr(simplified, attr, getattr(self, attr))

        return simplified

    def merge_into(self, accountant):

        for parent in accountant.pred_genes:
            if parent not in self.pred_genes:
                self.pred_genes[parent] = accountant.pred_genes[parent].copy()
            else:
                for tid in accountant.pred_genes[parent]:
                    if tid in self.pred_genes[parent]:
                        self.pred_genes[parent][tid] |= accountant.pred_genes[parent][tid]
                    else:
                        self.pred_genes[parent][tid] = accountant.pred_genes[parent][tid]

        def clearBit(int_type, offset):
            mask = ~(1 << offset)
            return (int_type & mask)

        def getBit(int_type, offset):
            mask = 1 << offset
            if (int_type & mask):
                return 1
            else:
                return 0

        def setBit(int_type, offset):
            mask = 1 << offset
            return (int_type | mask)

        for ref_gene in accountant.ref_genes:
            for refid in accountant.ref_genes[ref_gene]:
                val = self.ref_genes[ref_gene][refid]
                oval = accountant.ref_genes[ref_gene][refid]
                for pos in range(3):
                    bit = getBit(val, pos) | getBit(oval, pos)
                    if bit:
                        val = setBit(val, pos)
                bit = getBit(val, 3) & getBit(oval, 3)
                if not bit:
                    val = clearBit(val, 3)
                self.ref_genes[ref_gene][refid] = val

        for chrom in accountant.intron_chains:
            if chrom not in self.intron_chains:
                self.intron_chains[chrom] = dict()
            for strand in accountant.intron_chains[chrom]:
                if strand not in self.intron_chains[chrom]:
                    self.intron_chains[chrom][strand] = dict()
                for ic_key in accountant.intron_chains[chrom][strand]:
                    if ic_key not in self.intron_chains[chrom][strand]:
                        self.intron_chains[chrom][strand][ic_key] = [set(), set()]
                    self.intron_chains[chrom][strand][ic_key][0].update(
                        accountant.intron_chains[chrom][strand][ic_key][0]
                    )
                    self.intron_chains[chrom][strand][ic_key][1].update(
                        accountant.intron_chains[chrom][strand][ic_key][1]
                    )

        self.monoexonic_matches[0].update(accountant.monoexonic_matches[0])
        self.monoexonic_matches[1].update(accountant.monoexonic_matches[1])

        for feature in ("exons", "starts", "ends"):
            store = getattr(accountant, feature)
            my_store = getattr(self, feature)
            for chrom in store:
                if chrom not in my_store:
                    my_store[chrom] = dict()
                for strand in store[chrom]:
                    if strand not in my_store[chrom]:
                        my_store[chrom][strand] = dict()
                    for feat in store[chrom][strand]:
                        if feat not in my_store[chrom][strand]:
                            my_store[chrom][strand][feat] = 0b0
                        my_store[chrom][strand][feat] |= store[chrom][strand][feat]

    def __store_multiexonic_result(self, transcr, strand):

        """
        Private procedure to store the results regarding a multiexonic
        transcript compared to the reference.
        :param transcr: a transcript instance
        :param strand: the strand to be used for storing
        :return:
        """

        # assert transcr.monoexonic is False
        ic_key = tuple([tuple([intron[0], intron[1]]) for intron in sorted(transcr.introns)])
        # We register the variables here because multiple access slows the program down!
        chrom = transcr.chrom[:]
        exon_num = transcr.exon_num
        if ic_key not in self.intron_chains[chrom][strand]:
            self.intron_chains[chrom][strand][ic_key] = [set(), set()]
        self.intron_chains[chrom][strand][ic_key][1].add(transcr.id)

        for index, exon in enumerate(transcr.exons):
            exon = tuple([exon[0], exon[1]])
            if exon not in self.exons[chrom][strand]:
                self.exons[chrom][strand][exon] = 0b0
            self.exons[chrom][strand][exon] |= 0b10  # set it as "in prediction"
            if index == 0:
                self.exons[chrom][strand][exon] |= 0b010000
                if exon[1] not in self.starts[transcr.chrom][strand]:
                    self.starts[chrom][strand][exon[1]] = 0b0
                self.starts[chrom][strand][exon[1]] |= 0b10
            elif index == exon_num - 1:
                self.exons[chrom][strand][exon] |= 0b010000
                if exon[0] not in self.ends[chrom][strand]:
                    self.ends[chrom][strand][exon[0]] = 0b00
                self.ends[chrom][strand][exon[0]] |= 0b10
            else:
                self.exons[chrom][strand][exon] |= 0b01000

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
            if other_exon not in self.exons[transcr.chrom][strand]:
                self.exons[transcr.chrom][strand][other_exon] = 0b1  # Necessary as we are trying to avoid preloading

            assert other_exon in self.exons[transcr.chrom][strand]

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
                if refgene not in self.ref_genes:
                    self.ref_genes[refgene] = dict()
                if refid not in self.ref_genes[refgene]:
                    self.ref_genes[refgene][refid] = 0b1000

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
            self.__store_multiexonic_result(transcr, strand)
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
    def __calculate_redundant_statistics(common_pred, common_ref, pred, ref):
        """
        Function to calculate precision, recall and F1 for non-redundant statistics given
        the numbers in common for the prediction, for the reference, and then the numbers of the reference and
        of the prediction.
        :param common:
        :param pred:
        :param ref:
        :return: precision, recall, f1stat
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
        if self.self_analysis is True:
            self.logger.info("No general statistics file will be printed for a self-analysis run.")
            self.logger.removeHandler(self.queue_handler)
            return

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

        (splice_precision, splice_recall, splice_f1) = self.__calculate_statistics(
            *intron_results["splices"]
        )

        intron_precision, intron_recall, intron_f1 = self.__calculate_redundant_statistics(
            *intron_results["introns"]["redundant"]
        )

        intron_nr_precision, intron_nr_recall, intron_nr_f1 = self.__calculate_statistics(
            *intron_results["introns"]["non_redundant"]
        )

        # TODO: revise
        intron_chain_stats = self.__calculate_redundant_statistics(
            *intron_results["intron_chains"]["redundant"])

        intron_chains_precision, intron_chains_recall, intron_chains_f1 = intron_chain_stats

        intron_chain_nr_stats = self.__calculate_statistics(
            *intron_results["intron_chains"]["non_redundant"])

        intron_chains_nr_precision, intron_chains_nr_recall, intron_chains_nr_f1 = intron_chain_nr_stats

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
                self.__format_rowname("Splice site level"),
                splice_recall * 100,
                splice_precision * 100,
                splice_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Intron level"),
                intron_recall * 100,
                intron_precision * 100, intron_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Intron level (NR)"),
                intron_nr_recall * 100,
                intron_nr_precision * 100, intron_nr_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Intron chain level"),
                intron_chains_recall * 100,
                intron_chains_precision * 100,
                intron_chains_f1 * 100), file=out)
            print("{0} {1:.2f}  {2:.2f}  {3:.2f}".format(
                self.__format_rowname("Intron chain level (NR)"),
                intron_chains_nr_recall * 100,
                intron_chains_nr_precision * 100,
                intron_chains_nr_f1 * 100), file=out)
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
                                                    intron_results["introns"]["non_redundant"][2],
                                                    intron_results["introns"]["non_redundant"][0]),
                  file=out)
            print(self.__format_comparison_segments("Novel introns",
                                                    intron_results["introns"]["non_redundant"][1],
                                                    intron_results["introns"]["non_redundant"][0]),
                  file=out)
            print("", file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Missed transcripts (0% nF1)",
                                                gene_transcript_results["ref"]["private"][0],
                                                gene_transcript_results["ref"]["total"]),
                  file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Novel transcripts (0% nF1)",
                                                gene_transcript_results["pred"]["private"][0],
                                                gene_transcript_results["pred"]["total"]),
                  file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Missed genes (0% nF1)",
                                                gene_transcript_results["ref"]["private"][1],
                                                len(self.ref_genes)), file=out)

            # noinspection PyTypeChecker
            print(self.__format_comparison_line("Novel genes (0% nF1)",
                                                gene_transcript_results["pred"]["private"][1],
                                                len(self.pred_genes)), file=out)

        self.logger.removeHandler(self.queue_handler)
        # self.queue_handler.close()
        return
