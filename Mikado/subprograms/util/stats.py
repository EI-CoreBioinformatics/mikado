#!/usr/bin/env python3
# coding: utf-8

""" Script to calculate statistics about an annotation file.
It can take both GTF and GFF files as input."""

import sys
import argparse
import csv
from ...exceptions import InvalidCDS
from .. import to_gff
from ...loci import Transcript, Gene
from ...parsers import GFF
import numpy
from collections import namedtuple, Counter
from ...utilities.log_utils import create_default_logger
from collections import defaultdict

__author__ = "Luca Venturini"

# pylint: disable=E1101
numpy.seterr(all="ignore")  # Suppress warnings
numpy.warnings.filterwarnings("ignore")
# pylint: enable=E1101


def itemize(counter):
    """
    Private static method to convert a counter into a 2-dimensional numpy array.
    :param counter: the counter to transform
    :type counter: Counter
    :return: Numpy array of the counter
    :rtype: array
    """

    if len(counter) == 0:
        return numpy.array([[], []])  # Empty 2dimensional array
    else:
        return numpy.array(list(zip(*counter.items())))


def weighted_percentile(a, percentile=numpy.array([75, 25]), weights=None):
    """
    O(nlgn) implementation for weighted_percentile.
    From: http://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
    Kudos to SO user Nayyary http://stackoverflow.com/users/2004093/nayyarv

    :param a: array or Counter
    :type a: (Counter|list|numpy.array|set|tuple)

    :param percentile: the percentiles to calculate.
    :type percentile: (numpy.array|list|tuple)

    """

    percentile = numpy.array(percentile)/100.0

    if isinstance(a, Counter):
        a, weigths = itemize(a)
    else:
        assert isinstance(a, (list, set, tuple, numpy.ndarray)), (a, type(a))
        if not isinstance(a, type(numpy.array)):
            a = numpy.array(a)
        if weights is None:
            weights = numpy.ones(len(a))

    a_indsort = numpy.argsort(a)
    a_sort = a[a_indsort]
    weights_sort = weights[a_indsort]
    ecdf = numpy.cumsum(weights_sort)

    percentile_index_positions = percentile * (weights.sum() - 1) + 1
    # need the 1 offset at the end due to ecdf not starting at 0
    locations = numpy.searchsorted(ecdf, percentile_index_positions)

    out_percentiles = numpy.zeros(len(percentile_index_positions))

    for i, empiricalLocation in enumerate(locations):
        # iterate across the requested percentiles
        if ecdf[empiricalLocation-1] == numpy.floor(percentile_index_positions[i]):
            # i.e. is the percentile in between 2 separate values
            uppWeight = percentile_index_positions[i] - ecdf[empiricalLocation-1]
            lowWeight = 1 - uppWeight

            out_percentiles[i] = a_sort[empiricalLocation-1] * lowWeight + \
                                 a_sort[empiricalLocation] * uppWeight
        else:
            # i.e. the percentile is entirely in one bin
            out_percentiles[i] = a_sort[empiricalLocation]

    return out_percentiles


class TranscriptComputer(Transcript):
    """
    Class that is used to calculate and store basic statistics about a transcript object.
    """

    data_fields = ["parent", 'chrom',
                   'start', 'end',
                   'introns', 'exons',
                   'exon_lengths', 'intron_lengths',
                   'cdna_length', 'selected_cds_length',
                   'cds_intron_lengths', 'cds_exon_lengths',
                   "five_utr_length", "three_utr_length",
                   "five_utr_num", "three_utr_num",
                   "selected_end_distance_from_junction"]
    data_tuple = namedtuple("transcript_data", data_fields, verbose=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.exon_lengths = []
        self.cds_exon_lengths = []
        self.utr_exon_lengths = []

        self.intron_lengths = []
        self.cds_intron_lengths = []
        self.utr_intron_lengths = []

    def finalize(self):
        """
        Method to be called when all exons/features have been
        added to the transcript. It will call the parent's finalize method,
        followed by calculation of the necessary statistics.
        """
        try:
            super().finalize()
        except InvalidCDS:
            super().strip_cds()

        self.exon_lengths = [e[1] - e[0] + 1 for e in self.exons]
        self.cds_exon_lengths = [c[1] - c[0] + 1 for c in self.selected_cds]
        self.utr_exon_lengths = [u[1] - u[0] + 1 for u in self.three_utr + self.five_utr]

        self.intron_lengths = [i[1] - i[0] + 1 for i in self.introns]
        self.cds_intron_lengths = [i[1] - i[0] for i in self.selected_cds_introns]
        self.utr_intron_lengths = [i[1] - i[0] for i in self.introns if
                                   i not in self.selected_cds_introns]

    def as_tuple(self):
        """Method to build a namedtuple containing only the basic information for stat building.

        We want to analyze the following:
        - cDNA length
        - CDS length
        - Exons (number and length)
        - CDS Exons (number and length)
        - Introns (number and length)
        - CDS Introns (number and length)
        """

        self.finalize()
        constructor = dict()
        for field in self.data_fields:
            constructor[field] = getattr(self, field)

        return self.data_tuple(**constructor)


def finalize_wrapper(finalizable):

    """
    Wrapper around the finalize method for the object
    :param finalizable: an object with the finalize method
    :return:
    """

    finalizable.finalize()
    return finalizable


class Calculator:

    """
    This class has the purpose of parsing a reference file,
    calculating the statistics, and printing them out.
    """

    def __init__(self, parsed_args):

        """Constructor function"""

        self.gff = parsed_args.gff
        if isinstance(self.gff, GFF.GFF3):
            self.is_gff = True
        else:
            self.is_gff = False

        self.__logger = create_default_logger("calculator")
        if parsed_args.verbose is True:
            self.__logger.setLevel("DEBUG")
        self.only_coding = parsed_args.only_coding
        self.tab_handle = parsed_args.tab_stats
        if self.tab_handle is not None:
            # TID, Parent, Exon number, cDNA length, CDS length, #coding exons, cDNA/CDS ratio
            # 5'UTR, 5'UTR exons, 3'UTR, 3'UTR exons
            __fieldnames = ["TID", "GID",
                            "Exon number",
                            "cDNA length",
                            "CDS length",
                            "# coding exons",
                            "cDNA/CDS ratio",
                            "5'UTR",
                            "5'UTR exons",
                            "3'UTR",
                            "3'UTR exons"]
            self.tab_writer = csv.DictWriter(self.tab_handle, __fieldnames, delimiter="\t")
            self.tab_writer.writeheader()

        self.out = parsed_args.out
        self.genes = dict()
        self.coding_genes = []
        self.__distances = numpy.array([])
        self.__positions = defaultdict(list)
        self.__coding_positions = defaultdict(list)
        self.__coding_distances = numpy.array([])
        self.__fieldnames = ['Stat', 'Total', 'Average', 'Mode', 'Min',
                             '1%', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', '99%', 'Max']
        self.__rower = csv.DictWriter(self.out,
                                      self.__fieldnames,
                                      delimiter="\t")
        self.__arrays = dict()
        self.__stores = dict()
        self.__prepare_stores()

    def parse_input(self):
        """
        Method to parse the input GTF/GFF file.
        """
        transcript2gene = dict()

        derived_features = set()

        current_gene = None

        for record in self.gff:
            if record.is_gene is True:
                self.__store_gene(current_gene)
                current_gene = Gene(record,
                                    only_coding=self.only_coding,
                                    logger=self.__logger)
            elif record.is_transcript is True or (self.is_gff is True and record.feature == "match"):
                if record.parent is None and record.feature != "match":
                    raise TypeError("No parent found for:\n{0}".format(str(record)))
                elif record.feature == "match":
                    gid = record.id
                else:
                    gid = record.parent[0]
                    # assert current_gene is not None, record
                if current_gene is None or gid != current_gene.id:
                    # Create a gene record
                    self.__store_gene(current_gene)
                    new_record = record.copy()
                    # if record.feature != "match":
                    #     new_record.feature = "gene"

                    current_gene = Gene(
                        new_record,
                        gid=gid,
                        only_coding=self.only_coding,
                        logger=self.__logger)
                transcript2gene[record.id] = gid
                current_gene.transcripts[record.id] = TranscriptComputer(record,
                                                                         logger=self.__logger)
            elif record.is_derived is True and record.is_exon is False:
                derived_features.add(record.id)
            elif record.is_exon is True:
                if self.is_gff is False:
                    if current_gene is None or record.gene != current_gene.id:
                        self.__store_gene(current_gene)
                        new_record = record.copy()
                        new_record.feature = "gene"
                        new_record.id = new_record.gene
                        current_gene = Gene(
                            new_record,
                            only_coding=self.only_coding,
                            logger=self.__logger)
                        record.id = record.transcript
                        transcript2gene[record.transcript] = record.gene
                        current_gene.transcripts[record.transcript] = TranscriptComputer(record,
                                                                                         logger=self.__logger)
                    elif record.transcript not in current_gene:
                        assert record.transcript not in transcript2gene, record.transcript
                        transcript2gene[record.transcript] = record.gene
                        current_gene.transcripts[record.transcript] = TranscriptComputer(record,
                                                                                         logger=self.__logger)
                    else:
                        current_gene.transcripts[record.transcript].add_exon(record)
                elif self.is_gff is True and "cDNA_match" in record.feature:
                    # Here we assume that we only have "cDNA_match" lines, with no parents
                    record.parent = record.id
                    if current_gene is None or record.id != current_gene.id:
                        self.__store_gene(current_gene)
                        new_record = record.copy()
                        current_gene = Gene(
                            new_record,
                            gid=new_record.id,
                            only_coding=self.only_coding,
                            logger=self.__logger)
                        transcript2gene[record.id] = record.id
                        current_gene.transcripts[record.id] = TranscriptComputer(record,
                                                                                 logger=self.__logger,
                                                                                 )
                    elif record.id not in current_gene:
                        raise ValueError(
                            "cDNA_match instances should not have more than one transcript per \"gene\"!")
                    else:
                        current_gene.transcripts[record.id].add_exon(record)

                elif self.is_gff is True:
                    current_gene.add_exon(record)

            elif record.header is True:
                continue
            else:
                continue

        self.__store_gene(current_gene)

    def __call__(self):

        self.parse_input()
        self.gff.close()
        distances = []
        for chromosome in self.__positions:
            __ordered = sorted(self.__positions[chromosome])
            for index, position in enumerate(__ordered[:-1]):
                distances.append(__ordered[index + 1][0] - position[1])

        self.__distances = numpy.array(distances)

        distances = []
        for chromosome in self.__positions:
            __ordered = sorted(self.__coding_positions[chromosome])
            for index, position in enumerate(__ordered[:-1]):
                distances.append(__ordered[index + 1][0] - position[1])

        self.__coding_distances = numpy.array(distances)
        self.writer()
        self.out.close()

    @staticmethod
    def get_stats(row: dict, array: numpy.array) -> dict:
        """
        Method to calculate the necessary statistic from a row of values.
        :param row: the output dictionary row.
        :type row: dict

        :param array: an array of values.

        :rtype : dict
        """

        # Decimal to second digit precision

        if array is None or len(array) == 0:
            return row

        array, weights = array

        if len(weights) > 0 and sum(weights) != 0:
            row["Average"] = "{0:,.2f}".format(round(
                sum(_[0] * _[1] for _ in zip(array, weights)) / sum(weights), 2))
            sorter = numpy.argsort(weights)
            try:
                moder = array[sorter][weights[sorter].searchsorted(weights.max()):]
            except TypeError as exc:
                raise TypeError((exc, array, weights, sorter))
            row["Mode"] = ";".join(str(x) for x in sorted(moder))
        else:
            row["Average"] = "NA"
            row["Mode"] = "NA"

        keys = ['Min', '1%', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', '99%', 'Max']
        if len(array) == 0:
            quantiles = ["NA"]*len(keys)
        else:
            quantiles = weighted_percentile(array,
                                            [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100],
                                            weights)
        for key, val in zip(keys, quantiles):
            try:
                row[key] = "{0:,.0f}".format(val)  # No decimal
            except KeyError:
                row[key] = val
            except ValueError:
                row[key] = val
            except TypeError:
                row[key] = val
            except Exception:
                raise
        return row

    def __prepare_stores(self):

        self.__stores["genes"] = set()
        self.__stores["coding_genes"] = set()
        self.__stores["monoexonic_genes"] = set()
        self.__stores["transcripts_per_gene"] = Counter()
        self.__stores["coding_transcripts_per_gene"] = Counter()
        self.__stores["exons"], self.__stores["exons_coding"] = Counter(), Counter()
        self.__stores["exon_num"], self.__stores["exon_num_coding"] = Counter(), Counter()
        self.__stores["introns"], self.__stores["introns_coding"] = Counter(), Counter()
        self.__stores["cds_introns"] = Counter()
        self.__stores["cds_exons"] = Counter()
        self.__stores["cds_exon_num"] = Counter()
        self.__stores["cds_exon_num_coding"] = Counter()
        self.__stores["cdna_lengths"] = Counter()  # Done
        self.__stores["cdna_lengths_coding"] = Counter()
        self.__stores["cds_lengths"] = Counter()  # Done
        self.__stores["cds_lengths_coding"] = Counter()  # Done
        self.__stores["cds_ratio"] = Counter()
        self.__stores["monoexonic_lengths"] = Counter()
        self.__stores["multiexonic_lengths"] = Counter()
        self.__stores["monocds_lengths"] = Counter()

        self.__stores["five_utr_lengths"] = Counter()
        self.__stores["five_utr_nums"] = Counter()
        self.__stores["three_utr_lengths"] = Counter()
        self.__stores["three_utr_nums"] = Counter()
        self.__stores["end_distance_from_junction"] = Counter()

    def __store_gene(self, gene):

        """

        :param gene:
        :type gene: (None|Mikado.loci.Gene)
        :return:
        """

        if gene is None:
            return

        gene.finalize()
        if self.only_coding is True and gene.is_coding is False:
            return

        self.__positions[gene.chrom].append((gene.start, gene.end, str(gene.strand)))
        self.__stores["genes"].add(gene.id)
        if gene.is_coding is True:
            self.__coding_positions[gene.chrom].append((gene.start, gene.end, gene.strand))
            self.__stores["coding_genes"].add(gene.id)

        self.__stores["transcripts_per_gene"].update([gene.num_transcripts])
        self.__stores["coding_transcripts_per_gene"].update([gene.num_coding_transcripts])

        if gene.monoexonic is True:
            self.__stores["monoexonic_genes"].add(gene.id)

        for tid in gene.transcripts:

            if self.tab_handle is not None:
                # __fieldnames = ["TID", "GID",
                #                 "Exon number",
                #                 "cDNA length",
                #                 "CDS length",
                #                 "# coding exons",
                #                 "cDNA/CDS ratio",
                #                 "5'UTR",
                #                 "5'UTR exons",
                #                 "3'UTR",
                #                 "3'UTR exons"]
                row = dict()
                row["TID"] = tid
                row["GID"] = gene.id
                row["Exon number"] = len(gene.transcripts[tid].exon_lengths)
                row["cDNA length"] = gene.transcripts[tid].cdna_length
                row["CDS length"] = gene.transcripts[tid].selected_cds_length
                row["# coding exons"] = len(gene.transcripts[tid].cds_exon_lengths)
                row["cDNA/CDS ratio"] = "{:.2f}".format(100 * gene.transcripts[tid].selected_cds_length / gene.transcripts[tid].cdna_length)
                row["5'UTR"] = gene.transcripts[tid].five_utr_length
                row["5'UTR exons"] = gene.transcripts[tid].five_utr_num
                row["3'UTR"] = gene.transcripts[tid].three_utr_length
                row["3'UTR exons"] = gene.transcripts[tid].three_utr_num
                self.tab_writer.writerow(row)

            self.__stores["exons"].update(
                gene.transcripts[tid].exon_lengths)
            exon_number = len(gene.transcripts[tid].exon_lengths)
            self.__stores["exon_num"].update([exon_number])
            if exon_number == 1:
                self.__stores["monoexonic_lengths"].update(
                    [gene.transcripts[tid].cdna_length])
            else:
                self.__stores["multiexonic_lengths"].update(
                    [gene.transcripts[tid].cdna_length])
            self.__stores["introns"].update(gene.transcripts[tid].intron_lengths)
            self.__stores["cds_introns"].update(gene.transcripts[tid].cds_intron_lengths)
            self.__stores["cds_exons"].update(gene.transcripts[tid].cds_exon_lengths)
            cds_num = len(gene.transcripts[tid].cds_exon_lengths)
            if cds_num == 1:
                self.__stores["monocds_lengths"].update(
                    [gene.transcripts[tid].selected_cds_length])
            elif exon_number == 1:
                if gene.transcripts[tid].selected_cds_length > 0:
                    self.__stores["monocds_lengths"].update(
                        [gene.transcripts[tid].selected_cds_length])

            self.__stores["cds_exon_num"].update([cds_num])
            self.__stores["cdna_lengths"].update(
                [gene.transcripts[tid].cdna_length])
            self.__stores["cds_lengths"].update(
                [gene.transcripts[tid].selected_cds_length])
            if gene.transcripts[tid].selected_cds_length > 0:
                __cds_length = gene.transcripts[tid].selected_cds_length
                __cdna_length = gene.transcripts[tid].cdna_length
                assert __cds_length > 0
                self.__stores["cds_lengths_coding"].update(
                    [gene.transcripts[tid].selected_cds_length])
                self.__stores["five_utr_lengths"].update(
                    [gene.transcripts[tid].five_utr_length])
                self.__stores["three_utr_lengths"].update(
                    [gene.transcripts[tid].three_utr_length])
                self.__stores["five_utr_nums"].update(
                    [gene.transcripts[tid].five_utr_num])
                self.__stores["three_utr_nums"].update(
                    [gene.transcripts[tid].three_utr_num])
                self.__stores["end_distance_from_junction"].update(
                    [gene.transcripts[tid].selected_end_distance_from_junction])
                self.__stores["cds_ratio"].update([100 * __cds_length / __cdna_length])

            if self.only_coding is False:
                if gene.transcripts[tid].selected_cds_length > 0:
                    self.__stores["cdna_lengths_coding"].update(
                        [gene.transcripts[tid].cdna_length])
                    self.__stores["exons_coding"].update(gene.transcripts[tid].exon_lengths)
                    self.__stores["exon_num_coding"].update(
                        [len(gene.transcripts[tid].exon_lengths)])
                    self.__stores["cds_exon_num_coding"].update(
                        [len(gene.transcripts[tid].cds_exon_lengths)])
                    self.__stores["introns_coding"].update(
                        gene.transcripts[tid].intron_lengths)
        return

    def __finalize_arrays(self):

        self.__arrays["Transcripts per gene"] = itemize(self.__stores["transcripts_per_gene"])

        self.__arrays["Coding transcripts per gene"] = itemize(self.__stores["coding_transcripts_per_gene"])
        self.__arrays["Intergenic distances"] = itemize(Counter(self.__distances))
        self.__arrays["Intergenic distances (coding)"] = itemize(Counter(self.__coding_distances))

        self.__arrays['CDNA lengths'] = itemize(self.__stores["cdna_lengths"])
        self.__arrays["CDNA lengths (mRNAs)"] = itemize(self.__stores["cdna_lengths_coding"])
        self.__arrays['CDS lengths'] = itemize(self.__stores["cds_lengths"])
        if self.only_coding is False:
            self.__arrays["CDS lengths (mRNAs)"] = itemize(self.__stores["cds_lengths_coding"])
            self.__arrays['Exons per transcript (mRNAs)'] = itemize(self.__stores["exon_num_coding"])
            self.__arrays['Exon lengths (mRNAs)'] = itemize(self.__stores["exons_coding"])
            self.__arrays["CDS exons per transcript (mRNAs)"] = itemize(
                self.__stores["cds_exon_num_coding"])

        self.__arrays['Monoexonic transcripts'] = itemize(self.__stores["monoexonic_lengths"])
        self.__arrays['MonoCDS transcripts'] = itemize(self.__stores["monocds_lengths"])
        self.__arrays['Exons per transcript'] = itemize(self.__stores["exon_num"])
        self.__arrays['Exon lengths'] = itemize(self.__stores["exons"])
        self.__arrays["Intron lengths"] = itemize(self.__stores["introns"])
        self.__arrays["Intron lengths (mRNAs)"] = itemize(self.__stores["introns_coding"])
        self.__arrays["CDS exons per transcript"] = itemize(self.__stores["cds_exon_num"])
        self.__arrays["CDS exon lengths"] = itemize(self.__stores["cds_exons"])
        self.__arrays["CDS Intron lengths"] = itemize(self.__stores["cds_introns"])
        self.__arrays["5'UTR exon number"] = itemize(self.__stores["five_utr_nums"])
        self.__arrays["3'UTR exon number"] = itemize(self.__stores["three_utr_nums"])
        self.__arrays["5'UTR length"] = itemize(self.__stores["five_utr_lengths"])
        self.__arrays["3'UTR length"] = itemize(self.__stores["three_utr_lengths"])
        self.__arrays["Stop distance from junction"] = itemize(
            self.__stores["end_distance_from_junction"])
        self.__arrays["CDS/cDNA ratio"] = itemize(self.__stores["cds_ratio"])
    # pylint: enable=too-many-locals,too-many-statements

    def writer(self):
        """Method which creates the final output"""

        self.__finalize_arrays()
        self.__rower.writeheader()
        self.__write_statrow('Number of genes',
                             len(self.__stores["genes"]))
        self.__write_statrow("Number of genes (coding)",
                             len(self.__stores["coding_genes"]))
        self.__write_statrow("Number of monoexonic genes",
                             len(self.__stores["monoexonic_genes"])
                             )
        self.__write_statrow('Transcripts per gene',
                             total=numpy.dot)
        self.__write_statrow("Coding transcripts per gene", total=numpy.dot)

        self.__write_statrow('CDNA lengths', total=numpy.dot)
        self.__write_statrow("CDNA lengths (mRNAs)", total=numpy.dot)
        self.__write_statrow('CDS lengths', total=numpy.dot)
        if self.only_coding is False:
            self.__write_statrow("CDS lengths (mRNAs)", total=False)

        self.__write_statrow("CDS/cDNA ratio", total=False)

        self.__write_statrow('Monoexonic transcripts', total=sum)
        self.__write_statrow('MonoCDS transcripts', total=sum)
        self.__write_statrow('Exons per transcript', total=numpy.dot)

        if self.only_coding is False:
            self.__write_statrow('Exons per transcript (mRNAs)',
                                 total="Exon lengths (mRNAs)")
        self.__write_statrow('Exon lengths', total=False)
        if self.only_coding is False:
            self.__write_statrow('Exon lengths (mRNAs)', total=False)
        self.__write_statrow("Intron lengths", total=False)
        if self.only_coding is False:
            self.__write_statrow("Intron lengths (mRNAs)", total=False)

        self.__write_statrow("CDS exons per transcript",
                             total="CDS exon lengths")

        if self.only_coding is False:
            self.__write_statrow("CDS exons per transcript (mRNAs)",
                                 total="CDS exon lengths")

        self.__write_statrow("CDS exon lengths", total=numpy.dot)
        self.__write_statrow("CDS Intron lengths", total=numpy.dot)
        self.__write_statrow("5'UTR exon number", total=sum)
        self.__write_statrow("3'UTR exon number", total=sum)
        self.__write_statrow("5'UTR length", total=numpy.dot)
        self.__write_statrow("3'UTR length", total=numpy.dot)
        self.__write_statrow("Stop distance from junction",
                             total=False)
        self.__write_statrow("Intergenic distances", total=False)
        self.__write_statrow("Intergenic distances (coding)", total=False)

    def __write_statrow(self, stat, total=True):
        """
        Static method to write out a statistic to the
        output file.
        :param stat: the name of the row
        :type stat: str
        :param total: value to display in the "Total" column
        :type total: str | int | sum | bool
        """
        row = dict()
        for key in self.__fieldnames:
            row[key] = "NA"
        row["Stat"] = stat
        if total is False:
            total = "NA"
        elif total is True:
            if len(self.__arrays[stat]) == 2 and isinstance(
                        self.__arrays[stat], numpy.ndarray):
                total = len(self.__arrays[stat][0])
            else:
                total = len(self.__arrays[stat])
        else:
            if total is sum:
                try:
                    _, weights = self.__arrays[stat]
                    total = sum(weights)
                except ValueError:
                    total = "NA"
            elif total is numpy.dot:
                total = numpy.dot(self.__arrays[stat][0],
                                  self.__arrays[stat][1])

            elif not isinstance(total, int):
                assert total in self.__arrays
                if len(self.__arrays[total]) == 2 and isinstance(
                        self.__arrays[total], numpy.ndarray):

                    total = len(self.__arrays[total][0])
                else:
                    total = len(self.__arrays[total])
            else:
                pass  # Just keep the total that was passed from the external

        row["Total"] = total
        if stat in self.__arrays:
            current_array = self.__arrays[stat]
            # assert isinstance(current_array, Counter), type(current_array)
        else:
            current_array = None
        row = self.get_stats(row, current_array)
        self.__rower.writerow(row)
# pylint: enable=too-many-instance-attributes


def launch(args):

    """
    Very simple launcher function, calls Calculator from this module.

    :param args: the argparse Namespace.
    """

    calculator = Calculator(args)
    calculator()


def stats_parser():

    """
    Argument parser.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--only-coding', dest="only_coding", action="store_true", default=False)
    parser.add_argument('--tab-stats',
                        dest="tab_stats",
                        default=None,
                        type=argparse.FileType('w'),
                        help="Optional tabular file to write statistics for each transcript.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument('gff', type=to_gff, help="GFF file to parse.")
    parser.add_argument('out', type=argparse.FileType('w'), default=sys.stdout, nargs='?')
    parser.set_defaults(func=launch)
    return parser
