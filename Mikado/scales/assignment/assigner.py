# coding: utf-8

"""
This class is the main workhorse of the compare.py utility.
"""

import argparse
import collections
import csv
import gzip
import logging
import operator
import queue
import sys
from collections import namedtuple
from functools import partial
from logging import handlers as log_handlers
from ...transcripts.transcript import Transcript, Namespace
from ..accountant import Accountant
from ..contrast import compare as c_compare
from ..resultstorer import ResultStorer
from ...exceptions import InvalidTranscript, InvalidCDS
from ...utilities.intervaltree import IntervalTree
import msgpack
import tempfile
from ..reference_preparation.gene_dict import GeneDict


# noinspection PyPropertyAccess,PyPropertyAccess

def msgpack_default(o):
    """Function to convert the types for msgpack"""
    if isinstance(o, tuple):
        return {'__type__': 'tuple', 'value': list(o)} 
    elif isinstance(o, set):
        return {'__type__': 'set', 'value': list(o)}
    elif isinstance(o, ResultStorer):
        return {'__type__': 'rstor', 'value': o.as_dict()}
    else:
        return o


def msgpack_convert(o):
    """Function to re-extract the data from msgpack"""
    if isinstance(o, dict) and o.get('__type__', None):
        if o['__type__'] == "set":
            return set(o['value'])
        elif o['__type__'] == b'list':
            return list(o[b'value'])
        elif o['__type__'] == 'tuple':
            return tuple(o['value'])
        elif o['__type__'] == 'rstor':
            return ResultStorer(state=o['value'])
        else:
            return o['value']
    else:
        return o


class Assigner:

    """
    This class has the purpose of assiging each prediction transcript to its best match
    among the reference transcripts.
    """

    def __init__(self,
                 index: str,
                 args: argparse.Namespace,
                 results=(),
                 printout_tmap=True,
                 fuzzymatch=0,
                 counter=None):

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

        :param printout_tmap: boolean value. If set to True, the object will print each row.
        :type printout_tmap: bool
        """

        if args is None:
            self.args = argparse.Namespace()
            self.args.out = sys.stdout
            self.args.distance = 2000
            self.args.protein_coding = False
            self.args.loq_queue = queue.Queue()
            self.args.exclude_utr = False
            self.args.verbose = False
            self.lenient = False
            self.use_prediction_alias = False
            self.__report_fusions = True
        else:
            self.args = args
            if not hasattr(args, "out"):
                self.args.out = sys.stdout
            if not hasattr(args, "distance"):
                self.args.distance = 2000
            if not hasattr(args, "protein_coding"):
                self.args.protein_coding = False
            if not hasattr(args, "loq_queue"):
                self.args.loq_queue = queue.Queue()
            if not hasattr(args, "exclude_utr"):
                self.args.exclude_utr = False
            if not hasattr(args, "verbose"):
                self.args.verbose = False
            self.lenient = getattr(args, "lenient", False)
            self.use_prediction_alias = getattr(args, "use_prediction_alias", False)
            self.__report_fusions = getattr(args, "report_fusions", True)

        # noinspection PyUnresolvedReferences
        # pylint: disable=no-member
        self.queue_handler = log_handlers.QueueHandler(self.args.log_queue)

        # pylint: enable=no-member
        if counter is None:
            self.logger = logging.getLogger("Assigner")
        else:
            self.logger = logging.getLogger("Assigner-{}".format(counter))
        self.logger.addHandler(self.queue_handler)
        # noinspection PyUnresolvedReferences
        if args.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        if hasattr(args, "fuzzymatch") and args.fuzzymatch not in (0, None):
            if fuzzymatch in (0, None):
                fuzzymatch = args.fuzzymatch

        if not isinstance(fuzzymatch, int) and not (fuzzymatch is None):
            raise TypeError(fuzzymatch)

        self.__fuzzymatch = fuzzymatch

        self.logger.propagate = True
        self.dbname = index
        self.genes = GeneDict(self.dbname, logger=self.logger,
                              exclude_utr=self.args.exclude_utr, protein_coding=self.args.protein_coding)
        self.positions = collections.defaultdict(dict)
        self.indexer = collections.defaultdict(list)
        self._load_positions()
        self.printout_tmap = printout_tmap
        self.__done = 0
        # noinspection PyUnresolvedReferences
        if self.printout_tmap:
            if self.args.gzip is False:
                self.tmap_out = open("{0}.tmap".format(args.out), 'wt')
            else:
                self.tmap_out = gzip.open("{0}.tmap.gz".format(args.out), 'wt')
            self.tmap_rower = csv.DictWriter(self.tmap_out, ResultStorer.__slots__, delimiter="\t")
            self.tmap_rower.writeheader()

        self.gene_matches = collections.defaultdict(dict)
        self.done = 0

        if self.printout_tmap is True:
            for gid in self.genes:
                for tid in self.genes[gid].transcripts:
                    self.gene_matches[gid][tid] = []

        self.stat_calculator = Accountant(self.genes, args=args,
                                          fuzzymatch=self.__fuzzymatch, counter=counter,
                                          load_ref=self.printout_tmap)
        self.self_analysis = self.stat_calculator.self_analysis
        self.__merged = False

        if results:
            self.load_from_results(results)
            self.__merged = True
            return

    def _load_positions(self):

        for row in self.genes.positions:
            chrom, start, end, gid = row
            key = (start, end)
            if key not in self.positions[chrom]:
                self.positions[chrom][key] = []
            self.positions[chrom][key].append(gid)

        for chrom in self.positions:
            self.indexer[chrom] = IntervalTree.from_tuples(self.positions[chrom].keys())

    def load_result(self, refmap, stats):

        gene_matches = msgpack.loads(refmap,
                                     raw=False, use_list=False,
                                     strict_map_key=False, object_hook=msgpack_convert)

        for gid, gene_match in gene_matches.items():
            for tid in gene_match:
                for match in gene_match[tid]:
                    self.gene_matches[gid][tid].append(ResultStorer(state=match))

        temp_stats = Namespace()
        for attribute, stat in stats:
            stat = msgpack.loads(stat,
                                 raw=False,
                                 use_list=False,
                                 strict_map_key=False,
                                 object_hook=msgpack_convert)
            setattr(temp_stats, attribute, stat)
        self.stat_calculator.merge_into(temp_stats)

    def dump(self):

        """Method to dump all results into the database"""

        refmap = msgpack.dumps(self.gene_matches, default=msgpack_default, strict_types=False)
        simplified = self.stat_calculator.serialize()
        stats = []
        for attribute in simplified.attributes:
            stat = msgpack.dumps(getattr(simplified, attribute),
                          default=msgpack_default, strict_types=True)
            stats.append((attribute, stat))

        return refmap, stats

    def add_to_refmap(self, result: ResultStorer) -> None:
        """
        :param result: the result of the compare function

        This method verifies whether the comparison in input is the best available for the
        candidate reference gene, and loads it into the refmap dictionary if such is the case.

        """

        if result is not None and result.ccode != ("u",):
            gid, tid = result.ref_gene[0], result.ref_id[0]
            if gid not in self.gene_matches:
                self.gene_matches[gid] = dict()
            if tid not in self.gene_matches[gid]:
                self.gene_matches[gid][tid] = list()
            self.gene_matches[gid][tid].append(result)
        return

    @classmethod
    def find_neighbours(cls, keys, position, distance=2000):
        """
        This class method is used to find the possible matches of a given prediction key.
        :param keys: the start
        :type keys: Mikado.utilities.intervaltree.IntervalTree

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

        found = keys.search(start, end, max_distance=distance)

        distances = []
        for key in found:
            # Append the key (to be used later for retrieval) and the distance
            distances.append(
                ((key.start, key.end),
                 max(0,
                     max(start, key.start) - min(end, key.end))))

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
        Simple getter for F1 statistics (N_F1, J_F1 and distance).
        :param curr_result: a result storer
        :type curr_result: ResultStorer

        :return: (J_F1, N_F1, distance)
        :rtype (float, float, int)
        """
        return (curr_result.j_f1[0],
                curr_result.n_f1[0],
                -curr_result.distance[0],
                (curr_result.ccode[0] not in ("x", "X")))

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
                self.logger.warning("Invalid transcript (due to CDS): %s", prediction.id)
                self.logger.warning("Error message: %s", err)
                # self.done += 1
                self.print_tmap(None)
                return None
        except InvalidTranscript as err:
            self.logger.warning("Invalid transcript: %s", prediction.id)
            self.logger.warning("Error message: %s", err)
            # self.done += 1
            self.print_tmap(None)
            return None
        return prediction

    def __check_for_fusions(self, prediction, matches, fuzzymatch=None):

        """
        This private method checks whether a transcript with
        multiple matches at distance 0 is indeed a fusion, or not.
        :param prediction:
        :param matches:
        :type matches: list
        :return:
        """

        # Input: dict[(start, end), distance]: [gene1, gene2 ...]
        # What we need to achieve in this section of the program
        # Detect which genes in range *are* fused
        # Assign a "fused" tag to each transcript in the fused gene that passes the criteria

        results = []
        new_matches = dict()
        match_to_gene = dict()

        strands = collections.defaultdict(set)

        # Get all the results for the single gene
        # We *want* to do the calculation for all hits
        for match in matches:
            match_to_gene[match[0]] = self.positions[prediction.chrom][match[0]]
            for gene, gene_match in iter((gene, self.genes[gene]) for gene in
                                         match_to_gene[match[0]]):
                strands[gene_match.strand].add(gene)
                new_matches[gene] = sorted(
                    [self.calc_and_store_compare(prediction, tra, fuzzymatch=fuzzymatch) for tra in gene_match],
                    key=self.get_f1, reverse=True)

        # If we have candidates for the fusion which are on its same strand
        # Keep only those. Otherwise in compact genomes we might call as "fusions"
        # all transcripts which are assigned to genes overlapping on the same strand!
        if prediction.strand is not None and len(strands[prediction.strand]) > 0:
            for strand in iter(key for key in strands.keys() if key != prediction.strand):
                for gene in strands[strand]:
                    # Add to the final results the data for the gene
                    results.extend(new_matches[gene])
                    del new_matches[gene]

        if len(new_matches) == 0:
            error = AssertionError("Filtered all results for %s. This is wrong!",
                                   prediction.id)
            self.logger.error(error)
            raise error

        fused = collections.defaultdict(list)
        dubious = set()
        best = set()

        for match, genes in match_to_gene.items():
            local_best = []
            for gene in genes:
                if gene not in new_matches:
                    continue
                elif new_matches[gene] == []:
                    continue
                # The best match is the first
                if new_matches[gene][0].j_f1[0] == 0 and new_matches[gene][0].n_recall[0] < 10:
                    dubious.add(gene)
                    continue
                else:
                    local_best.append(new_matches[gene][0])
                    fused[match].append(gene)
            if len(local_best) > 0:
                best.add((match, sorted(local_best, key=self.get_f1, reverse=True)[0]))

        if len(best) == 0:
            self.logger.debug("Filtered out all results for %s, selecting the best dubious one",
                              prediction.id)
            # I have filtered out all the results,
            # because I only match partially the reference genes
            best_result = sorted([new_matches[gene][0] for gene in dubious],
                                 key=self.get_f1, reverse=True)[0]

        else:
            best = [_[1] for _ in sorted(best, key=lambda res: (res[0][0], res[0][1]))]
            if len(best) == 1:
                best_result = best.pop()
            elif len(best) > 1 and self.__report_fusions is False:
                best_result = sorted(best, key=operator.attrgetter("j_f1", "n_f1")).pop()
            else:  # We have a fusion
                # Now retrieve the results according to their order on the genome
                # Keep only the result, not their position

                chrom = prediction.chrom
                start = min([prediction.start] + [self.genes[_.ref_gene[0]][_.ref_id[0]].start for _ in best])
                end = max([prediction.end] + [self.genes[_.ref_gene[0]][_.ref_id[0]].end for _ in best])
                location = "{}:{}..{}".format(chrom, start, end)

                best_result = []
                for result in best:
                    new_result = []
                    for key in ResultStorer.__slots__:
                        if key == "location":
                            new_result.append(location)
                        elif key == "ccode":
                            tup = tuple(["f"] + [getattr(result, key)[0]])
                            new_result.append(tup)
                        else:
                            new_result.append(getattr(result, key))
                    new_result = ResultStorer(*new_result)
                    best_result.append(new_result)

                for match, genes in fused.items():
                    for gene in iter(_ for _ in genes if _ not in dubious):
                        for position, result in enumerate(new_matches[gene]):
                            if result.j_f1[0] > 0 or result.n_recall[0] > 10:
                                result.ccode = ("f", result.ccode[0])
                                new_matches[gene][position] = result

        # Finally add the results
        for gene in new_matches:
            results.extend(new_matches[gene])
        self.logger.debug("\n".join([str(result) for result in results]))

        return results, best_result

    def __prepare_result(self, prediction, distances, fuzzymatch=0):

        """
        This private method prepares the matching result for cases where the
        minimum distance from a reference gene/transcript is 0.
        :param prediction:
        :param distances:
        :return:
        """

        matches = list(distance for distance in distances if distance[1] == 0)
        self.logger.debug("Matches for %s: %s", prediction.id,
                          ", ".join([_.id for _ in self.genes[self.positions[prediction.chrom][distances[0][0]][0]]]))

        same_strand = False
        if len(matches) > 1 and prediction.strand is not None:
            correct = list()
            for match in matches:
                for gid in self.positions[prediction.chrom][match[0]]:
                    if any([self.genes[gid].strand in (None, prediction.strand)]):
                        correct.append(match)
                        break
            if len(correct) > 0:
                matches = correct[:]
                same_strand = True
            del correct
        elif len(matches) > 1 and prediction.strand is None:
            same_strand = True

        if len(matches) > 1 and same_strand is True:
            self.logger.debug("More than one match for %s: %s",
                              prediction.id,
                              matches)
            results, best_result = self.__check_for_fusions(prediction, matches, fuzzymatch=fuzzymatch)
        else:
            matches = [self.genes[_] for _
                       in self.positions[prediction.chrom][matches[0][0]]]
            results = []
            for match in matches:
                self.logger.debug("%s: type %s", repr(match), type(match))
                for tra in match:
                    try:
                        results.append(self.calc_and_store_compare(prediction, tra, fuzzymatch=fuzzymatch))
                    except TypeError:
                        failed = []
                        for intron in prediction.introns:
                            if not (isinstance(intron[0], int) and isinstance(intron[1], int)):
                                failed.append(intron)

                        raise TypeError((failed, type(prediction), prediction.introns))

            results = sorted(results, reverse=True,
                             key=self.get_f1)

            best_result = results[0]

        return results, best_result

    def self_analyse_prediction(self, prediction: Transcript, distances, fuzzymatch=0):

        """This method will be invoked during a self analysis run."""

        assert len(distances) >= 1 and distances[0][1] == 0

        genes = []
        __gene_removed = False
        __gene_found = False
        for distance in distances:
            for posis in self.positions[prediction.chrom][distance[0]]:
                gene = self.genes[posis]
                if prediction.id in gene.transcripts.keys():
                    __gene_found = True
                    if len(gene.transcripts) == 1:
                        __gene_removed = True
                        continue
                genes.append((self.genes[posis], distance[1]))

        if __gene_found is False:
            raise AssertionError("Parent for {} not found in the index!".format(prediction.id))

        genes = sorted(genes, key=operator.itemgetter(1))

        # noinspection PyUnresolvedReferences
        if len(genes) == 0 or genes[0][1] > self.args.distance:
            # noinspection PyTypeChecker,PyUnresolvedReferences
            best_result = ResultStorer("-", "-",
                                       "u",
                                       prediction.id,
                                       prediction.parent[0],
                                       prediction.exon_num,
                                       "-",
                                       *[0] * 9 + ["-"] + [prediction.location])
            self.stat_calculator.store(prediction, best_result, None)
            results = [best_result]
        elif genes[0][1] > 0:  # Up or downstream
            results = [self.calc_and_store_compare(prediction, reference, fuzzymatch=fuzzymatch)
                       for reference in genes[0][0]]
            best_result = sorted(results,
                                 key=operator.attrgetter("distance"))[0]
        else:
            # All other cases
            genes = [_[0] for _ in genes if _[1] == 0]
            same_strand = False
            if prediction.strand is not None:
                same_strand_genes = [_ for _ in genes if _.strand == prediction.strand]
                if len(same_strand_genes) >= 1:
                    genes = same_strand_genes
                    same_strand = True

            if len(genes) == 1:
                results = [self.calc_and_store_compare(prediction, reference, fuzzymatch=fuzzymatch)
                           for reference in genes[0] if reference.id != prediction.id]
                assert len(results) > 0, (genes[0].transcripts.keys(), prediction.id)
                best_result = sorted(results,
                                     key=operator.attrgetter("distance"))[0]
            else:
                result_dict = dict()
                results = []

                for gene in genes:
                    if gene.id in prediction.parent:
                        if prediction.id not in gene.transcripts:
                            raise AssertionError("{} not found in {}!".format(prediction.id, gene.id))
                        elif len(gene.transcripts) == 1:
                            raise AssertionError(
                                "We should have ignored {}! Gene removed: {}".format(
                                    gene.id, __gene_removed))

                    result_dict[gene.id] = sorted(
                        [self.calc_and_store_compare(prediction, reference, fuzzymatch=fuzzymatch)
                         for reference in gene if reference.id != prediction.id],
                        key=self.get_f1,
                        reverse=True)
                    if len(result_dict[gene.id]) == 0:
                        raise AssertionError(
                            "Nothing found for {} vs. {} (transcripts: {})".format(
                                gene.id, prediction.id, ", ".join(list(gene.transcripts.keys()))
                            ))

                if same_strand is True:
                    # This is a fusion, period
                    results = []
                    best_results = []
                    for gene in result_dict:
                        if result_dict[gene][0].n_recall[0] > 10 or result_dict[gene][0].j_f1[0] > 0:
                            results.extend(result_dict[gene])
                            best_results.append(result_dict[gene][0])
                    if len(best_results) == 1:
                        best_result = best_results[0]
                    elif len(best_results) > 1 and self.__report_fusions is True:
                        values = []
                        for key in ResultStorer.__slots__:
                            if key in ["gid", "tid", "distance", "tid_num_exons"]:
                                values.append(getattr(best_results[0], key))
                            elif key == "ccode":
                                values.append(tuple(["f"] + [_.ccode[0] for _ in best_results]))
                            else:
                                val = tuple([getattr(result, key)[0] for result in best_results])
                                values.append(val)
                        best_result = ResultStorer(*values)
                        for result in results:
                            if result.j_f1[0] > 0 or result.n_recall[0] > 10:
                                result.ccode = ("f", result.ccode[0])
                    elif len(best_results) > 1 and self.__report_fusions is False:
                        best_result = [_ for _ in sorted(best_results, key=operator.attrgetter("j_f1", "n_f1",))].pop()

                # Check how many
                if not results:
                    # This is not a fusion. Let us just select the best match.
                    best = None
                    for gene in result_dict:
                        if best is None:
                            best = [gene, result_dict[gene][0]]
                        else:
                            if sorted([best[1], result_dict[gene][0]],
                                      key=self.get_f1, reverse=True)[0] == result_dict[gene][0]:
                                best = [gene, result_dict[gene][0]]
                    results = result_dict[best[0]]
                    best_result = best[1]

        return results, best_result

    def get_best(self, prediction: Transcript, fuzzymatch=None):

        """
        :param prediction: the candidate transcript to be analysed
        :type prediction: Transcript

        :param fuzzymatch: leniency in determining whether two introns are related
        :type fuzzymatch: (int|None)

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

        if fuzzymatch is None:
            fuzzymatch = self.__fuzzymatch

        # Quickly check that the transcript is OK
        prediction = self.__prepare_transcript(prediction)
        if prediction is None:
            return None
        elif prediction.finalized is False:
            self.logger.error("%s failed to be finalised. Ignoring it.", prediction.id)
            return None

        # Ignore non-coding RNAs if we are interested in protein-coding transcripts only
        # noinspection PyUnresolvedReferences
        if self.args.protein_coding is True and prediction.combined_cds_length == 0:
            self.logger.debug("No CDS for %s. Ignoring.", prediction.id)
            # self.done += 1
            self.print_tmap(None)
            return None

        if prediction.chrom in self.indexer:
            keys = self.indexer[prediction.chrom]
        else:
            keys = IntervalTree()

        # noinspection PyUnresolvedReferences
        # Here I am getting all the hits within striking distance for my hit
        distances = self.find_neighbours(keys,
                                         (prediction.start, prediction.end),
                                         distance=self.args.distance)
        if self.self_analysis is True:
            results, best_result = self.self_analyse_prediction(prediction, distances,
                                                                fuzzymatch=fuzzymatch)
        else:
            self.logger.debug("Distances for %s: %s", prediction.id, distances)
            # Unknown transcript
            # noinspection PyUnresolvedReferences
            if len(distances) == 0 or distances[0][1] > self.args.distance:
                ccode = "u"
                # noinspection PyTypeChecker,PyUnresolvedReferences
                best_result = ResultStorer("-", "-",
                                           ccode,
                                           prediction.id,
                                           prediction.parent[0],
                                           prediction.exon_num,
                                           "-",
                                           *[0] * 9 + ["-"] + [prediction.location])
                self.stat_calculator.store(prediction, best_result, None)
                results = [best_result]
            elif distances[0][1] > 0:
                # Polymerase run-on
                match = self.genes[self.positions[prediction.chrom][distances[0][0]][0]]
                results = [self.calc_and_store_compare(prediction, reference) for reference in match]
                best_result = sorted(results, key=operator.attrgetter("distance"))[0]
            else:
                # All other cases
                try:
                    results, best_result = self.__prepare_result(prediction, distances, fuzzymatch=fuzzymatch)
                except (InvalidTranscript, ZeroDivisionError, AssertionError) as exc:
                    self.logger.error("Something went wrong with %s. Ignoring it. Error: %s",
                                      prediction.id, exc)
                    self.logger.debug(exc)
                    return None

        if self.use_prediction_alias is True:
            if isinstance(best_result, list):
                for idx in range(len(best_result)):
                    best_result[idx].tid = prediction.alias
            else:
                best_result.tid = prediction.alias
            for idx in range(len(results)):
                results[idx].tid = prediction.alias

        for result in results:
            self.add_to_refmap(result)

        self.logger.debug("Finished with %s", prediction.id)
        if isinstance(best_result, list):
            for result in best_result:
                self.print_tmap(result)
        else:
            self.print_tmap(best_result)

        return best_result

    def finish(self):
        """
        This function is called upon completion of the input file parsing.
        It will print out the reference file assignments into the refmap
        file, and the final statistics into the stats file.
        """

        self.logger.debug("Finished parsing, total: %d transcript%s.", self.done, "s" if self.done > 1 else "")

        self.print_refmap()
        self.stat_calculator.print_stats()
        self.logger.info("Finished printing final stats")
        self.tmap_out.close()
        self.logger.info("Closed output files")

    def calc_and_store_compare(self, prediction: Transcript, reference: Transcript, fuzzymatch=0) -> ResultStorer:
        """Thin layer around the calc_and_store_compare class method.

        :param prediction: a Transcript instance.

        :param reference: a Transcript instance to which the prediction is compared to.

        :rtype ResultStorer

        """

        result, reference_exon = c_compare(prediction, reference,
                                           lenient=self.lenient, fuzzymatch=fuzzymatch)

        assert reference_exon is None or reference_exon in reference.exons
        self.stat_calculator.store(prediction, result, reference_exon)

        return result

    @staticmethod
    def compare(prediction: Transcript,
                reference: Transcript,
                lenient=False,
                strict_strandedness=False,
                fuzzy_match=0) -> (ResultStorer, tuple):

        """Function to compare two transcripts and determine a ccode. Thin wrapper
        around the compare cython code.

        :param prediction: the transcript query
        :type prediction: Transcript

        :param reference: the reference transcript against which we desire to
        calculate the ccode and other stats.
        :type reference: Transcript

        :param lenient: boolean flag. If switched to True, the comparison will be performed by considering
        only the *internal* boundaries of the transcripts.
        :type lenient: bool

        :param strict_strandedness: boolean flag. If switched to True, transcripts without a strand
        will be considered as stranded when calculating the type of overlap with the transcript.
        Eg. if the option is set to True and the transcript has a class code of e, it will become an x.
        :type strict_strandedness: bool

        :rtype (ResultStorer, (int,int)) | (ResultStorer, None)

        Please see the class_codes subpackage for details on the available class codes.

        This is a class method, and can therefore be used outside of a class instance.
        """

        return c_compare(prediction,
                         reference,
                         lenient=lenient,
                         strict_strandedness=strict_strandedness,
                         fuzzymatch=fuzzy_match)

    def print_tmap(self, res):
        """
        This method will print a ResultStorer instance onto the TMAP file.
        :param res: result from compare
        :type res: (ResultStorer | None)
        """

        if res is not None:
            if not isinstance(res, ResultStorer):
                self.logger.exception("Wrong type for res: %s", str(type(res)))
                self.logger.exception(repr(res))
                raise ValueError
            else:
                if self.printout_tmap is False:
                    pass
                else:
                    self.__done += 1
                    self.tmap_rower.writerow(res.as_dict())

    @staticmethod
    def result_sorter(result):

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

        bad_ccodes = ["x", "X", "P", "p"]
        bad_ccodes = set(bad_ccodes)

        orderer = (len(set.intersection(bad_ccodes, set(result.ccode))) == 0,
                   result.j_f1, result.e_f1,
                   result.n_f1,
                   "f" in result.ccode)

        return orderer

    def print_refmap(self) -> None:

        """Function to print out the best match for each gene."""
        self.logger.info("Starting printing RefMap")
        # noinspection PyUnresolvedReferences
        if self.args.gzip is False:
            opening_function = partial(open, "{0}.refmap".format(self.args.out), 'wt')
        else:
            opening_function = partial(gzip.open, "{0}.refmap.gz".format(self.args.out), 'wt')

        with opening_function() as out:
            if self.args.extended_refmap is True:
                fields = ["ref_id", "ccode", "tid", "gid",
                          "nRecall", "nPrecision", "nF1",
                          "jRecall", "jPrecision", "jF1",
                          "eRecall", "ePrecision", "eF1",
                          "ref_gene",
                          "best_ccode", "best_tid", "best_gid",
                          "best_nRecall", "best_nPrecision", "best_nF1",
                          "best_jRecall", "best_jPrecision", "best_jF1",
                          "best_eRecall", "best_ePrecision", "best_eF1",
                          "location"]
            else:
                fields = ["ref_id", "ccode", "tid", "gid",
                          "nF1", "jF1", "eF1",
                          "ref_gene",
                          "best_ccode", "best_tid", "best_gid",
                          "best_nF1", "best_jF1", "best_eF1",
                          "location"]
            out_tuple = namedtuple("refmap", fields)

            rower = csv.DictWriter(out, fields, delimiter="\t")
            rower.writeheader()

            for gid in sorted(self.gene_matches.keys()):

                rows = []
                best_picks = []
                assert len(self.gene_matches[gid].keys()) > 0
                for tid in sorted(self.gene_matches[gid].keys()):
                    if len(self.gene_matches[gid][tid]) == 0:
                        # First part of the tuple
                        row = tuple([tid, gid] + ["NA"] * 6)
                    else:
                        # Choose the best hit for the transcript
                        if any((x.j_f1[0] > 0 or x.n_f1[0] > 0) for x in self.gene_matches[gid][tid]):
                                best = sorted(self.gene_matches[gid][tid],
                                              key=self.result_sorter, reverse=True)[0]
                        else:
                            best = sorted(self.gene_matches[gid][tid],
                                          key=operator.attrgetter("distance"),
                                          reverse=False)[0]
                        best_picks.append(best)
                        # Store the result for the transcript
                        if self.args.extended_refmap is True:
                            row = tuple([tid, gid, ",".join(best.ccode),
                                         best.tid, best.gid,
                                         best.n_recall[0], best.n_prec[0], best.n_f1[0],
                                         best.j_recall[0], best.j_prec[0], best.j_f1[0],
                                         best.e_recall[0], best.e_prec[0], best.e_f1[0],
                                         best.location])
                        else:
                            row = tuple([tid, gid, ",".join(best.ccode),
                                         best.tid, best.gid,
                                         best.n_f1[0], best.j_f1[0], best.e_f1[0],
                                         best.location])

                    rows.append(row)

                if len(best_picks) > 0:
                    best_pick = sorted(best_picks,
                                       key=self.result_sorter,
                                       reverse=True)[0]
                else:
                    best_pick = None

                for row in rows:
                    if best_pick is not None:
                        if self.self_analysis is False:
                            assert row[2] != "NA", row
                        if self.args.extended_refmap is True:
                            row = out_tuple(row[0],  # Ref TID
                                            row[2],  # class code
                                            row[3],  # Pred TID
                                            row[4],  # Pred GID
                                            row[5], row[6], row[7],  # N
                                            row[8], row[9], row[10], #J
                                            row[11], row[12], row[13], #E
                                            row[1],
                                            ",".join(best_pick.ccode),
                                            best_pick.tid,
                                            best_pick.gid,
                                            best_pick.n_recall[0], best_pick.n_prec[0], best_pick.n_f1[0],
                                            best_pick.j_recall[0], best_pick.j_prec[0], best_pick.j_f1[0],
                                            best_pick.e_recall[0], best_pick.e_prec[0], best_pick.e_f1[0],
                                            row[14]  #Location
                                            )
                        else:
                            # fields = ["ref_id", "ccode", "tid", "gid",
                            #                           "nF1", "jF1", "eF1",
                            #                           "ref_gene",
                            #                           "best_ccode", "best_tid", "best_gid",
                            #                           "best_nF1", "best_jF1", "best_eF1",
                            #                           "location"]
                            try:
                                row = out_tuple(ref_id=row[0],  # Ref TID
                                                ccode=row[2],  # class code
                                                tid=row[3],  # Pred TID
                                                gid=row[4],  # Pred GID
                                                nF1=row[5], jF1=row[6], eF1=row[7],  # Pred F1
                                                ref_gene=row[1],
                                                best_ccode=",".join(best_pick.ccode),
                                                best_tid=best_pick.tid,
                                                best_gid=best_pick.gid,
                                                best_nF1=best_pick.n_f1[0],
                                                best_jF1=best_pick.j_f1[0],
                                                best_eF1=best_pick.e_f1[0],
                                                location=row[8]  # Location
                                                )
                            except IndexError:
                                self.logger.critical("Error in creating the refmap output")
                                self.logger.critical(row)
                                raise
                    else:
                        if self.args.extended_refmap is True:
                            row = out_tuple(*[row[0]] + ["NA"] * 12 + [row[1]] + ["NA"] * 12 + [
                                self.genes[gid][row[0]].location])
                        else:
                            row = out_tuple(*[row[0]] + ["NA"] * 6 + [row[1]] + ["NA"] * 6 + [
                                self.genes[gid][row[0]].location])
                    # noinspection PyProtectedMember,PyProtectedMember
                    rower.writerow(row._asdict())
        self.logger.info("Finished printing RefMap")
        return None
