#!/usr/bin/env python3
# coding: utf-8

""" Script to calculate statistics about an annotation file.
It can take both GTF and GFF files as input."""

import sys
import argparse
import re
import csv
from ...exceptions import InvalidTranscript, InvalidCDS
from .. import to_gff
from ...loci_objects.transcript import Transcript
from ...parsers import GFF
# pylint: disable=E1101,no-name-in-module
from numpy import array as num_array
from numpy import mean, percentile
import numpy
# pylint: enable=E1101
# import numpy
from collections import namedtuple, Counter
from array import array as c_array

__author__ = "Luca Venturini"

# pylint: disable=E1101
numpy.seterr(all="ignore")  # Suppress warnings
numpy.warnings.filterwarnings("ignore")
# pylint: enable=E1101


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
        self.utr_exon_lengths = [u[2] - u[1] + 1 for u in self.three_utr + self.five_utr]

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


# pylint: disable=too-many-instance-attributes
class GeneObject:

    """
    Basic gene container for the annotation file parsing.
    """

    def __init__(self, record, only_coding=False):

        self.transcripts = dict()
        self.chrom = record.chrom
        self.start = record.start
        self.end = record.end
        self.strand = record.strand
        # pylint: disable=invalid-name
        self.id = record.id
        # pylint: enable=invalid-name

        self.only_coding = only_coding
        self.coding_transcripts = set()

    def add_transcript(self, tcomputer: TranscriptComputer):
        """
        Simplified addition of a transcript to the gene.
        :param tcomputer: the candidate transcript
        :type tcomputer: TranscriptComputer
        """
        self.transcripts[tcomputer.id] = tcomputer
        self.start = min(self.start, tcomputer.start)
        self.end = max(self.end, tcomputer.end)

    def finalize(self):
        """
        Wrapper for the Transcript.finalize method, plus some attribute setting for the
        instance as well.
        """
        if len(self.transcripts) == 0:
            return

        to_remove = set()
        for tid in self.transcripts:
            try:
                self.transcripts[tid] = self.transcripts[tid].as_tuple()
                if self.only_coding is True and self.transcripts[tid].selected_cds_length == 0:
                    to_remove.add(tid)
                if self.transcripts[tid].selected_cds_length > 0:
                    self.coding_transcripts.add(tid)
            except InvalidTranscript as _:
                print("Invalid transcript: {0}".format(tid), file=sys.stderr)
                to_remove.add(tid)
        if len(to_remove) == len(self.transcripts):
            self.transcripts = dict()
        else:
            for tid_to_remove in to_remove:
                del self.transcripts[tid_to_remove]

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        else:
            if self.start != other.start:
                return self.start < other.start
            else:
                return self.end < other.end

    def __eq__(self, other):
        if not isinstance(self, type(other)):
            return False
        if self.chrom == other.chrom:
            if self.strand == other.strand:
                if self.start == other.start and self.end == other.end:
                    return True
        return False

    @property
    def introns(self):
        """
        It returns the set of all introns in the container.
        :rtype : set
        """

        return set(self.transcripts[tid].introns for tid in self.transcripts)

    @property
    def exons(self):
        """
        It returns the set of all exons in the container.
        :rtype : set
        """
        return set.union(*[set(self.transcripts[tid].exons) for tid in self.transcripts])

    @property
    def has_monoexonic(self):
        """
        True if any of the transcripts is monoexonic.
        :rtype : bool
        """
        return any(len(self.transcripts[tid].introns) == 0 for tid in self.transcripts.keys())

    @property
    def num_transcripts(self):
        """
        Number of transcripts.
        :rtype : int
        """
        return len(self.transcripts)

    @property
    def is_coding(self):
        """
        Property. It evaluates to True if at least one transcript is coding, False otherwise.
        """
        return any(self.transcripts[tid].selected_cds_length > 0 for tid in self.transcripts.keys())


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
        self.only_coding = parsed_args.only_coding
        self.out = parsed_args.out
        self.genes = dict()
        self.coding_genes = []
        self.__distances = num_array([])
        self.__coding_distances = num_array([])
        self.__fieldnames = ['Stat', 'Total', 'Average', 'Mode', 'Min',
                             '5%', '10%', '25%', 'Median', '75%', '90%', '95%', 'Max']
        self.__rower = csv.DictWriter(self.out,
                                      self.__fieldnames,
                                      delimiter="\t")
        self.__arrays = dict()

    @staticmethod
    def __is_gene(record):
        """
        Private method to estimate whether a record is a gene or not
        :param record:
        :return: bool
        """

        if "locus" == record.feature or record.is_gene is True:
            return True
        elif record.is_parent is True and record.is_transcript is False:
            return True
        return False

    def parse_input(self):
        """
        Method to parse the input GTF/GFF file.
        """
        transcript2gene = dict()

        derived_features = set()

        for record in self.gff:
            if record.header is True or re.search(r"[^^]locus", record.feature):
                continue
            elif record.is_transcript is True:
                if record.parent is None:
                    raise TypeError("No parent found for:\n{0}".format(str(record)))
                transcript2gene[record.id] = record.parent[0]
                if self.is_gff is False:
                    new_record = record.copy()
                    new_record.feature = "gene"
                    if new_record.gene not in self.genes:
                        self.genes[new_record.gene] = GeneObject(
                            new_record,
                            only_coding=self.only_coding)
                self.genes[record.parent[0]].transcripts[record.id] = TranscriptComputer(record)
            elif record.is_derived is True and record.is_gene is False:
                derived_features.add(record.id)
            elif self.__is_gene(record) is True:
                self.genes[record.id] = GeneObject(
                    record, only_coding=self.only_coding)
            else:
                for parent in iter(pparent for pparent in record.parent if
                                   pparent not in derived_features):
                    try:
                        gid = transcript2gene[parent]
                    except KeyError as err:
                        raise KeyError("{0}, line: {1}".format(err, record))
                    self.genes[gid].transcripts[parent].add_exon(record)

        for gid in self.genes:
            self.genes[gid].finalize()

    def __call__(self):

        self.parse_input()
        ordered_genes = sorted(self.genes.values())
        self.coding_genes = list(g for g in ordered_genes if g.is_coding is True)
        distances = []
        for index, gene in enumerate(ordered_genes[:-1]):
            next_gene = ordered_genes[index + 1]
            if next_gene.chrom != gene.chrom:
                continue
            distances.append(next_gene.start + 1 - gene.end)

        self.__distances = num_array(distances)

        distances = []
        for index, gene in enumerate(self.coding_genes[:-1]):
            next_gene = self.coding_genes[index + 1]
            if next_gene.chrom != gene.chrom:
                continue
            distances.append(next_gene.start + 1 - gene.end)

        self.__coding_distances = num_array(distances)

        self.writer()

    @staticmethod
    def get_stats(row: dict, array: num_array) -> dict:
        """
        Method to calculate the necessary statistic from a row of values.
        :param row: the output dictionary row.
        :type row: dict

        :param array: an array of values.

        :rtype : dict
        """

        # Decimal to second digit precision
        if array is None:
            return row
        row["Average"] = "{0:,.2f}".format(round(mean(array), 2))

        counter_object = Counter(array)
        moder = [x for x in counter_object if
                 counter_object[x] == counter_object.most_common(1)[0][1]]
        row["Mode"] = ";".join(str(x) for x in moder)
        keys = ['Min', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', 'Max']
        if len(array) == 0:
            quantiles = ["NA"]*len(keys)
        else:
            quantiles = [percentile(array, x) for x in [0, 5, 10, 25, 50, 75, 90, 95, 100]]
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

    # pylint: disable=too-many-locals,too-many-statements
    def __prepare_arrays(self):
        """
        This private method prepares the arrays which will be used to calculate
        all the necessary statistics.
        Such arrays include therefore intergenic distances, exon lengths, etc.,
        casted into integer/float arrays whenever possible.
        :return:
        """

        self.__arrays["Number of transcripts"] = num_array(
            [self.genes[x].num_transcripts for x in self.genes])
        self.__arrays["Transcripts per gene"] = num_array(
            [self.genes[_].num_transcripts for _ in self.genes])
        self.__arrays["Coding transcripts per gene"] = num_array(
            [self.genes[_].num_transcripts for _ in self.genes])
        self.__arrays["Intergenic distances"] = self.__distances
        self.__arrays["Intergenic distances (coding)"] = self.__coding_distances

        exons, exons_coding = c_array('i'), c_array('i')
        exon_num, exon_num_coding = c_array('i'), c_array('i')
        introns, introns_coding = c_array('i'), c_array('i')
        cds_introns = c_array('i')
        cds_exons = c_array('i')
        cds_exon_num = c_array('i')
        cds_exon_num_coding = c_array('i')
        cdna_lengths = c_array('i')  # Done
        cdna_lengths_coding = c_array('i')
        cds_lengths = c_array('i')  # Done
        monoexonic_lengths, multiexonic_lengths = c_array('i'), c_array('i')
        monocds_lengths = c_array('i')

        five_utr_lengths, five_utr_nums = c_array('i'), c_array('i')
        three_utr_lengths, three_utr_nums = c_array('i'), c_array('i')

        end_distance_from_junction = c_array('i')

        for gene in self.genes:
            for tid in self.genes[gene].transcripts:
                exons.extend(self.genes[gene].transcripts[tid].exon_lengths)
                exon_number = len(self.genes[gene].transcripts[tid].exon_lengths)
                exon_num.append(exon_number)
                if exon_number == 1:
                    monoexonic_lengths.append(self.genes[gene].transcripts[tid].cdna_length)
                else:
                    multiexonic_lengths.append(self.genes[gene].transcripts[tid].cdna_length)
                introns.extend(self.genes[gene].transcripts[tid].intron_lengths)
                cds_introns.extend(self.genes[gene].transcripts[tid].cds_intron_lengths)
                cds_exons.extend(self.genes[gene].transcripts[tid].cds_exon_lengths)
                cds_num = len(self.genes[gene].transcripts[tid].cds_exon_lengths)
                if cds_num == 1:
                    monocds_lengths.append(
                        self.genes[gene].transcripts[tid].selected_cds_length)
                elif exon_num == 1:
                    if self.genes[gene].transcripts[tid].selected_cds_length > 0:
                        monocds_lengths.append(
                            self.genes[gene].transcripts[tid].selected_cds_length)

                cds_exon_num.append(cds_num)
                cdna_lengths.append(self.genes[gene].transcripts[tid].cdna_length)
                cds_lengths.append(self.genes[gene].transcripts[tid].selected_cds_length)
                if self.genes[gene].transcripts[tid].selected_cds_length > 0:
                    five_utr_lengths.append(
                        self.genes[gene].transcripts[tid].five_utr_length)
                    three_utr_lengths.append(
                        self.genes[gene].transcripts[tid].three_utr_length)
                    five_utr_nums.append(
                        self.genes[gene].transcripts[tid].five_utr_num)
                    three_utr_nums.append(
                        self.genes[gene].transcripts[tid].three_utr_num)
                    end_distance_from_junction.append(
                        self.genes[gene].transcripts[tid].selected_end_distance_from_junction)

                if self.only_coding is False:
                    if self.genes[gene].transcripts[tid].selected_cds_length > 0:
                        cdna_lengths_coding.append(
                            self.genes[gene].transcripts[tid].cdna_length)
                        exons_coding.extend(self.genes[gene].transcripts[tid].exon_lengths)
                        exon_num_coding.append(
                            len(self.genes[gene].transcripts[tid].exon_lengths))
                        cds_exon_num_coding.append(
                            len(self.genes[gene].transcripts[tid].cds_exon_lengths))
                        introns_coding.extend(
                            self.genes[gene].transcripts[tid].intron_lengths)

        self.__arrays['CDNA lengths'] = cdna_lengths
        self.__arrays["CDNA lengths (mRNAs)"] = cdna_lengths_coding
        self.__arrays['CDS lengths'] = cds_lengths
        if self.only_coding is False:
            self.__arrays["CDS lengths (mRNAs)"] = c_array('i',
                                                           [_ for _ in cds_lengths if _ > 0])
            self.__arrays['Exons per transcript (mRNAs)'] = exon_num_coding
            self.__arrays['Exon lengths (mRNAs)'] = exons_coding
            self.__arrays["CDS exons per transcript (mRNAs)"] = cds_exon_num_coding

        self.__arrays['Monoexonic transcripts'] = monoexonic_lengths
        self.__arrays['MonoCDS transcripts'] = monocds_lengths
        self.__arrays['Exons per transcript'] = exon_num
        self.__arrays['Exon lengths'] = exons
        self.__arrays["Intron lengths"] = introns
        self.__arrays["Intron lengths (mRNAs)"] = introns_coding
        self.__arrays["CDS exons per transcript"] = cds_exon_num
        self.__arrays["CDS exon lengths"] = cds_exons
        self.__arrays["CDS Intron lengths"] = cds_introns
        self.__arrays["5'UTR exon number"] = five_utr_nums
        self.__arrays["3'UTR exon number"] = three_utr_nums
        self.__arrays["5'UTR length"] = five_utr_lengths
        self.__arrays["3'UTR length"] = three_utr_lengths
        self.__arrays["Stop distance from junction"] = end_distance_from_junction
    # pylint: enable=too-many-locals,too-many-statements

    def writer(self):
        """Method which creates the final output"""

        self.__prepare_arrays()
        self.__rower.writeheader()
        self.__write_statrow('Number of genes',
                             len(self.genes))
        self.__write_statrow("Number of genes (coding)",
                             len(self.coding_genes))

        self.__write_statrow('Number of transcripts',
                             total=sum)
        self.__write_statrow('Transcripts per gene',
                             total=sum(self.genes[x].num_transcripts for x in self.genes))
        self.__write_statrow("Number of coding transcripts",
                             total=sum(len(x.coding_transcripts) for x in self.coding_genes))
        self.__write_statrow("Coding transcripts per gene",
                             total=sum(len(x.coding_transcripts) for x in self.coding_genes))

        self.__write_statrow('CDNA lengths', total=False)
        self.__write_statrow("CDNA lengths (mRNAs)", total=False)
        self.__write_statrow('CDS lengths', total=False)
        if self.only_coding is False:
            self.__write_statrow("CDS lengths (mRNAs)", total=False)

        self.__write_statrow('Monoexonic transcripts')
        self.__write_statrow('MonoCDS transcripts')
        self.__write_statrow('Exons per transcript', total='Exon lengths')

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

        self.__write_statrow("CDS exon lengths", total=sum)
        self.__write_statrow("CDS Intron lengths", total=sum)
        self.__write_statrow("5'UTR exon number", total=sum)
        self.__write_statrow("3'UTR exon number", total=sum)
        self.__write_statrow("5'UTR length", total=sum, )
        self.__write_statrow("3'UTR length", total=sum)
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
            total = len(self.__arrays[stat])
        else:
            if total is sum:
                total = sum(self.__arrays[stat])
            elif not isinstance(total, int):
                assert total in self.__arrays
                total = len(self.__arrays[total])

        row["Total"] = total
        if stat in self.__arrays:
            current_array = self.__arrays[stat]
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
    parser.add_argument('gff', type=to_gff, help="GFF file to parse.")
    parser.add_argument('out', type=argparse.FileType('w'), default=sys.stdout, nargs='?')
    parser.set_defaults(func=launch)
    return parser
