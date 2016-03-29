# coding: utf-8

"""
This class is the main workhorse of the compare.py utility.
"""

import sys
import csv
from intervaltree import IntervalTree
from logging import handlers as log_handlers
import queue
import logging
import collections
import argparse
import operator
from collections import namedtuple
from .resultstorer import ResultStorer
from . import calc_f1
from ..loci.transcript import Transcript
from ..exceptions import InvalidTranscript, InvalidCDS
from .accountant import Accountant
from ..utilities import overlap
from .contrast import compare as c_compare


# noinspection PyPropertyAccess,PyPropertyAccess
class Assigner:

    """
    This class has the purpose of assiging each prediction transcript to its best match
    among the reference transcripts.
    """

    def __init__(self,
                 genes: dict,
                 positions: collections.defaultdict,
                 args: argparse.Namespace,
                 stat_calculator: Accountant):

        """

        :param genes: a dictionary which contains
        the gene containers for the reference transcript objects.
        :type genes: dict

        :param positions: a defaultdict which is used for fast lookup of genomic positions
        :type positions: collections.defaultdict

        :param args: the parameters passed through the command line
        :type args: None | argparse.Namespace

        :param stat_calculator: an instance of Accountant, to which the results will be sent to.
        :type stat_calculator: Accountant
        """

        if args is None:
            self.args = argparse.Namespace()
            self.args.out = sys.stdout
            self.args.distance = 2000
            self.args.protein_coding = False
            self.args.loq_queue = queue.Queue()
            self.args.exclude_utr = False
            self.args.verbose = False
        else:
            self.args = args

        # noinspection PyUnresolvedReferences
        # pylint: disable=no-member
        self.queue_handler = log_handlers.QueueHandler(self.args.log_queue)
        # pylint: enable=no-member
        self.logger = logging.getLogger("Assigner")
        self.logger.addHandler(self.queue_handler)
        # noinspection PyUnresolvedReferences
        if args.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        self.logger.propagate = False

        self.genes = genes
        self.positions = positions
        self.gene_matches = collections.defaultdict(dict)
        for gid in self.genes:
            for tid in self.genes[gid].transcripts:
                self.gene_matches[gid][tid] = []

        #
        # for chrom in positions:
        #     for key in positions[chrom]:
        #         for gid in positions[chrom][key]:
        #             self.gene_matches[gid] = dict()
        #             for tid in self.genes[gid]:
        #                 self.gene_matches[gid][tid.id] = []
        #

        self.indexer = collections.defaultdict(list).fromkeys(self.positions)

        for chrom in positions:
            self.indexer[chrom] = IntervalTree.from_tuples(self.positions[chrom].keys())

        # noinspection PyUnresolvedReferences
        self.tmap_out = open("{0}.tmap".format(args.out), 'wt')
        self.tmap_rower = csv.DictWriter(self.tmap_out, ResultStorer.__slots__, delimiter="\t")
        self.tmap_rower.writeheader()
        self.done = 0
        self.stat_calculator = stat_calculator

    def add_to_refmap(self, result: ResultStorer) -> None:
        """
        :param result: the result of the compare function

        This method verifies whether the comparison in input is the best available for the
        candidate reference gene, and loads it into the refmap dictionary if such is the case.

        """

        if result is not None and result.ccode != ("u",):
            self.gene_matches[result.ref_gene[0]][result.ref_id[0]].append(result)
        return

    @classmethod
    def find_neighbours(cls, keys, position, distance=2000):
        """
        This class method is used to find the possible matches of a given prediction key.
        :param keys: the start
        :type keys: intervaltree.IntervalTree

        :param position: the position of my prediction in the genome
        :type position: (int, int)

        :param distance: optional maximum distance of a prediction
        from reference before being called an unknown
        :type distance: int

        :return: a list of distances, with the following format:
                 (start, end), distance
        """

        start, end = position

        # This should happen only if we are analysing a prediction from a scaffold
        # with no reference annotation on it
        if len(keys) == 0:
            return []

        found = keys.search(start-distance, end+distance)

        distances = []
        for key in found:
            # Append the key (to be used later for retrieval) and the distance
            distances.append(
                ((key.begin, key.end),
                 max(0,
                     max(start, key.begin) - min(end, key.end))))

        # Sort by distance
        distances = sorted(distances, key=operator.itemgetter(1))

        return distances

    @staticmethod
    def dubious_getter(dubious_result):
        """
        Function used to perform the sorting of the matches.
        :param dubious_result: a result
        :type dubious_result: ResultStorer
        """
        getter = operator.attrgetter("j_f1", "n_f1")
        return getter(dubious_result[0])

    @staticmethod
    def get_f1(curr_result):
        """
        Simple getter for F1 statistics (N_F1 and J_F1)
        :param curr_result: a result storer
        :type curr_result: ResultStorer

        :return: (J_F1, N_F1)
        :rtype (float, float)
        """
        return curr_result.j_f1[0], curr_result.n_f1[0]

    def __prepare_transcript(self, prediction: Transcript):
        """
        Private method that checks that a prediction transcript is OK
        before starting to analyse its concordance with the reference.
        :param prediction:
        :return:
        """

        # Prepare the prediction to be analysed
        prediction.logger = self.logger
        try:
            prediction.finalize()
            # noinspection PyUnresolvedReferences
            if self.args.exclude_utr is True:
                prediction.remove_utrs()
        except InvalidCDS:
            try:
                prediction.strip_cds()
            except InvalidTranscript as err:
                self.logger.warn("Invalid transcript (due to CDS): %s",
                                 prediction.id)
                self.logger.warn("Error message: %s", err)
                self.done += 1
                self.print_tmap(None)
                return None
        except InvalidTranscript as err:
            #         args.queue.put_nowait("mock")
            self.logger.warn("Invalid transcript: %s", prediction.id)
            self.logger.warn("Error message: %s", err)
            self.done += 1
            self.print_tmap(None)
            return None
        return prediction

    def __check_for_fusions(self, prediction, matches):

        """
        This private method checks whether a transcript with
        multiple matches at distance 0 is indeed a fusion, or not.
        :param prediction:
        :param matches:
        :return:
        """

        new_matches = collections.defaultdict(list)

        # Matches format: [(start,end), distance]
        for match in matches:
            best = None
            gene_matches = [self.genes[_] for _ in
                            self.positions[prediction.chrom][match[0]]]

            self.logger.debug("Match for %s: %s",
                              match,
                              [_.id for _ in gene_matches])
            for gene_match in gene_matches:
                __res = sorted([self.calc_and_store_compare(prediction, tra) for tra in gene_match],
                               reverse=True, key=self.get_f1)
                new_matches[(match[0], gene_match.strand)].append(__res)

            #     best_res = (gene_match.strand, __res)
            #     if best is None:
            #         best = best_res
            #     else:
            #         best = sorted([best, best_res],
            #                       reverse=True,
            #                       key=lambda res: self.get_f1(res[1][0]))[0]
            # if best is not None:
            #     new_matches.append(best)

        strands = set(_[1] for _ in new_matches)
        # strands = set(_[0] for _ in matches)
        if len(strands) > 1 and prediction.strand in strands:
            matches = dict((_[0][0],
                            sorted(_[1], reverse=True,
                                   key=lambda res: self.get_f1(res[1][0]))[0])
                           for _ in new_matches.items() if
                           _[0][1] == prediction.strand)
            # matches = list(mmatch for mmatch in matches
            #                if mmatch[0] == prediction.strand)
            if len(matches) == 0:
                raise ValueError("I filtered out all matches. This is wrong!")
        else:
            matches = dict((_[0][0],
                            sorted(_[1], reverse=True,
                                   key=lambda res: self.get_f1(res[1][0]))[0])
                           for _ in new_matches.items())

        best_fusion_results = []
        results = []  # Final results
        dubious = []  # Necessary for a double check.

        for match in matches:



            m_res = match[1]
            # A fusion is called only if I have one of the following conditions:
            # the transcript gets one of the junctions of the other transcript
            # the exonic overlap is >=10% (n_recall)_

            if m_res[0].j_f1[0] == 0 and m_res[0].n_recall[0] < 10:
                dubious.append(m_res)
                continue
            # List of ResultStorer instances
            results.extend(m_res)
            best_fusion_results.append(m_res[0])

        if len(results) == 0:
            self.logger.debug("Filtered out all results for %s, using the dubious ones",
                              prediction.id)
            # I have filtered out all the results,
            # because I only match partially the reference genes
            dubious = sorted(dubious, key=self.dubious_getter)
            results = dubious[0]
            best_fusion_results = [results[0]]

        fused_group = tuple(sorted([_.ref_gene[0] for _ in best_fusion_results]))

        values = []
        if len(fused_group) > 1:
            for result in results:
                result.ccode = ("f", result.ccode[0])
            for key in ResultStorer.__slots__:
                if key in ["gid", "tid", "distance"]:
                    values.append(getattr(best_fusion_results[0], key))
                elif key == "ccode":
                    values.append(
                        tuple(["f"] + [getattr(_, "ccode")[1] for _ in best_fusion_results])
                    )
                else:
                    val = tuple([getattr(x, key)[0] for x in best_fusion_results])
                    values.append(val)
            best_result = ResultStorer(*values)
        else:
            best_result = sorted(best_fusion_results, reverse=True, key=self.get_f1)[0]
        # Finished creating the fusion result

        return results, best_result

    def __prepare_result(self, prediction, distances):

        """
        This private method prepares the matching result for cases where the
        minimum distance from a reference gene/transcript is 0.
        :param prediction:
        :param distances:
        :return:
        """

        matches = list(distance for distance in distances if distance[1] == 0)
        self.logger.debug("Matches for %s: %s",
                          prediction.id,
                          matches)

        if len(matches) > 1 and prediction.strand is not None:
                correct = list()
                for match in matches:
                    for gid in self.positions[prediction.chrom][match[0]]:
                        if any([self.genes[gid].strand in (None, prediction.strand)]):
                            correct.append(match)
                            break
                if len(correct) > 0:
                    matches = correct[:]
                del correct

        if len(matches) > 1:
            self.logger.debug("More than one match for %s: %s",
                              prediction.id,
                              matches)
            results, best_result = self.__check_for_fusions(prediction, matches)

        else:
            matches = [self.genes[_] for _
                       in self.positions[prediction.chrom][matches[0][0]]]
            results = []
            for match in matches:
                self.logger.debug("%s: type %s", repr(match), type(match))
                results.extend([self.calc_and_store_compare(prediction, tra) for tra in match])

            results = sorted(results, reverse=True,
                             key=operator.attrgetter("j_f1", "n_f1"))

            best_result = results[0]

        return results, best_result

    def get_best(self, prediction: Transcript):

        """
        :param prediction: the candidate transcript to be analysed
        :type prediction: Transcript

        This function will get the best possible assignment for each transcript.
        Fusion genes are called when the following conditions are verified:
        - the prediction intersects (at least) two transcripts in (at least)
        two different loci
        - the suspected fusion transcript lies on the same strand of all
        candidate fused genes
        - each candidate transcript has at least one fusion or
        10% of its nucleotides covered by the fusion transcript.
        The 10% threshold is hard-coded in the function.
        """

        self.logger.debug("Started with %s", prediction.id)

        # Quickly check that the transcript is OK
        prediction = self.__prepare_transcript(prediction)
        if prediction is None:
            return None

        # Ignore non-coding RNAs if we are interested in protein-coding transcripts only
        # noinspection PyUnresolvedReferences
        if self.args.protein_coding is True and prediction.combined_cds_length == 0:
            #         args.queue.put_nowait("mock")
            self.logger.debug("No CDS for %s. Ignoring.", prediction.id)
            self.done += 1
            self.print_tmap(None)
            return None

        if prediction.chrom in self.indexer:
            keys = self.indexer[prediction.chrom]
        else:
            keys = IntervalTree()

        # noinspection PyUnresolvedReferences
        distances = self.find_neighbours(keys,
                                         (prediction.start, prediction.end),
                                         distance=self.args.distance)
        self.logger.debug("Distances for %s: %s",
                          prediction.id,
                          distances)

        # Unknown transcript
        # noinspection PyUnresolvedReferences
        if len(distances) == 0 or distances[0][1] > self.args.distance:
            ccode = "u"
            # noinspection PyTypeChecker,PyUnresolvedReferences
            best_result = ResultStorer("-", "-",
                                       ccode,
                                       prediction.id,
                                       ",".join(prediction.parent), *[0] * 9 + ["-"])
            self.stat_calculator.store(prediction, best_result, None)
            results = [best_result]
        elif distances[0][1] > 0:
            # Polymerase run-on
            match = self.genes[self.positions[prediction.chrom][distances[0][0]][0]]
            results = [self.calc_and_store_compare(prediction, reference) for reference in match]
            best_result = sorted(results,
                                 key=operator.attrgetter("distance"))[0]
        else:
            # All other cases
            results, best_result = self.__prepare_result(prediction, distances)

        for result in results:
            self.add_to_refmap(result)

        self.logger.debug("Finished with %s", prediction.id)
        self.print_tmap(best_result)
        self.done += 1
        return best_result

    def finish(self):
        """
        This function is called upon completion of the input file parsing.
        It will print out the reference file assignments into the refmap
        file, and the final statistics into the stats file.
        """
        self.logger.info("Finished parsing, total: %d transcript%s.",
                         self.done, "s" if self.done > 1 else "")
        self.refmap_printer()
        self.stat_calculator.print_stats()
        self.tmap_out.close()

    def calc_and_store_compare(self, prediction: Transcript, reference: Transcript) -> ResultStorer:
        """Thin layer around the calc_and_store_compare class method.

        :param prediction: a Transcript instance.

        :param reference: a Transcript instance to which the prediction is compared to.

        :rtype ResultStorer

        """

        result, reference_exon = c_compare(prediction, reference)

        assert reference_exon is None or reference_exon in reference.exons
        self.stat_calculator.store(prediction, result, reference_exon)

        return result

    @staticmethod
    def compare(prediction: Transcript, reference: Transcript) -> (ResultStorer, tuple):

        """Function to compare two transcripts and determine a ccode. Thin wrapper
        around the compare cython code.

        :param prediction: the transcript query
        :type prediction: Transcript

        :param reference: the reference transcript against which we desire to
        calculate the ccode and other stats.
        :type reference: Transcript

        :rtype (ResultStorer, (int,int)) | (ResultStorer, None)

        Available ccodes (from Cufflinks documentation):

        - =    Complete intron chain match
        - c    Contained (perfect junction recall and precision, imperfect recall)
        - j    Potentially novel isoform (fragment): at least one splice junction is shared
        with a reference transcript
        - e    Single exon transfrag overlapping a reference exon and at least
        10 bp of a reference intron, indicating a possible pre-mRNA fragment.
        - i    A *monoexonic* transfrag falling entirely within a reference intron
        - o    Generic exonic overlap with a reference transcript
        - p    Possible polymerase run-on fragment (within 2Kbases of a reference transcript)
        - u    Unknown, intergenic transcript
        - x    Exonic overlap with reference on the opposite strand (class codes e, o, m, c, _)
        - X    Overlap on the opposite strand, with some junctions in common (probably a serious mistake,
               unless non-canonical splicing junctions are involved).

        Please note that the description for i is changed from Cufflinks.

        We also provide the following additional classifications:

        - f    gene fusion - in this case, this ccode will be followed by the
        ccodes of the matches for each gene, separated by comma
        - _    Complete match, for monoexonic transcripts
        (nucleotide F1>=80% - i.e. min(precision,recall)>=66.7%
        - m    Exon overlap between two monoexonic transcripts
        - n    Potential extension of the reference - we have added new splice junctions
        *outside* the boundaries of the transcript itself
        - C    Contained transcript with overextensions on either side
        (perfect junction recall, imperfect nucleotide specificity)
        - J    Potentially novel isoform, where all the known junctions
        have been confirmed and we have added others as well *externally*
        - I    *multiexonic* transcript falling completely inside a known transcript
        - h    AS event in which at least a couple of introns overlaps but without any
               junction in common.
        - O    Reverse generic overlap - the reference is monoexonic while the prediction isn't
        - P    Possible polymerase run-on fragment
        - mo   Monoexonic overlap - the prediction is monoexonic and the reference is multiexonic
        (within 2K bases of a reference transcript), on the opposite strand

        This is a class method, and can therefore be used outside of a class instance.
        """

        return c_compare(prediction, reference)

    def print_tmap(self, res: ResultStorer):
        """
        This method will print a ResultStorer instance onto the TMAP file.
        :param res: result from compare
        :type res: ResultStorer | None
        """
        if self.done % 10000 == 0 and self.done > 0:
            self.logger.info("Done %d transcripts", self.done)
        elif self.done % 1000 == 0 and self.done > 0:
            self.logger.debug("Done %d transcripts", self.done)
        if res is not None:
            if not isinstance(res, ResultStorer):
                self.logger.exception("Wrong type for res: %s", str(type(res)))
                self.logger.exception(repr(res))
                raise ValueError
            else:
                self.tmap_rower.writerow(res.as_dict())

    @staticmethod
    def __result_sorter(result):

        """
        Method to sort the results for the refmap. Order:
        - CCode does not contain "x", "P", "p" (i.e. fragments on opposite strand or
        polymerase run-on fragments)
        - Exonic F1 (e_f1)
        - Junction F1 (j_f1)
        - "f" in ccode (i.e. transcript is a fusion)
        - Nucleotide F1 (n_f1)

        :param result: a resultStorer object
        :type result: ResultStorer
        :return: (int, float, float, float)
        """

        bad_ccodes = ["x", "P", "p"]
        bad_ccodes = set(bad_ccodes)

        orderer = (len(set.intersection(bad_ccodes, set(result.ccode))) == 0,
                   result.j_f1, result.e_f1,
                   result.n_f1,
                   "f" in result.ccode)

        return orderer

    @classmethod
    def result_sorter(cls, result):

        """
        Public interface of __result_sorter
        :param result: a result
        :return:
        """

        return cls.__result_sorter(result)

    def refmap_printer(self) -> None:

        """Function to print out the best match for each gene."""
        self.logger.info("Starting printing RefMap")
        # noinspection PyUnresolvedReferences
        with open("{0}.refmap".format(self.args.out), 'wt') as out:
            fields = ["ref_id", "ccode", "tid", "gid",
                      "ref_gene", "best_ccode", "best_tid", "best_gid"]
            out_tuple = namedtuple("refmap", fields)

            rower = csv.DictWriter(out, fields, delimiter="\t")
            rower.writeheader()

            for gid in sorted(self.gene_matches.keys()):

                rows = []
                best_picks = []
                assert len(self.gene_matches[gid].keys()) > 0
                for tid in sorted(self.gene_matches[gid].keys()):
                    if len(self.gene_matches[gid][tid]) == 0:
                        row = tuple([tid, gid, "NA", "NA", "NA"])
                    else:
                        if any(True if (x.j_f1[0] > 0 or x.n_f1[0] > 0) else False
                               for x in self.gene_matches[gid][tid]):
                            best = sorted(self.gene_matches[gid][tid],
                                          key=self.__result_sorter, reverse=True)[0]
                        else:
                            best = sorted(self.gene_matches[gid][tid],
                                          key=operator.attrgetter("distance"),
                                          reverse=False)[0]
                        best_picks.append(best)
                        row = tuple([tid, gid, ",".join(best.ccode),
                                     best.tid, best.gid])

                    rows.append(row)

                if len(best_picks) > 0:
                    best_pick = sorted(best_picks,
                                       key=self.__result_sorter,
                                       reverse=True)[0]
                else:
                    best_pick = None

                for row in rows:
                    if best_pick is not None:
                        assert row[2] != "NA", row
                        row = out_tuple(row[0], row[2],
                                        row[3], row[4],
                                        row[1], ",".join(best_pick.ccode),
                                        best_pick.tid, best_pick.gid)
                    else:
                        row = out_tuple(row[0], "NA", "NA", "NA", row[1], "NA", "NA", "NA")
                    # noinspection PyProtectedMember,PyProtectedMember
                    rower.writerow(row._asdict())
        self.logger.info("Finished printing RefMap")
        return None
