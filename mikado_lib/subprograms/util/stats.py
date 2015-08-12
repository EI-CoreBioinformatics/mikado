#!/usr/bin/env python3
# coding: utf-8

""" Script to calculate statistics about an annotation file. It can take both GTF and GFF files as input."""

import sys, argparse
from mikado_lib.exceptions import InvalidCDS
from mikado_lib.subprograms import to_gff
from mikado_lib.loci_objects import transcript
from mikado_lib.parsers import GFF
import numpy
import mikado_lib
from collections import namedtuple, Counter
from array import array as cArray

numpy.seterr(all="ignore")  # Suppress warnings
numpy.warnings.filterwarnings("ignore")


class TranscriptComputer(transcript.Transcript):
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
                   "selected_end_distance_from_junction"
                   ]
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
        Method to be called when all exons/features have been added to the transcript.
        It will call the parent's finalize method, followed by calculation of the necessary statistics.
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
        self.utr_intron_lengths = [i[1] - i[0] for i in
                                   filter(lambda j: j not in self.selected_cds_introns, self.introns)]

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
        self.id = record.id

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
            except mikado_lib.exceptions.InvalidTranscript as _:
                print("Invalid transcript: {0}".format(tid), file=sys.stderr)
                to_remove.add(tid)
        if len(to_remove) == len(self.transcripts):
            self.transcripts = dict()
        else:
            for x in to_remove:
                del self.transcripts[x]

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        else:
            if self.start != other.start:
                return self.start < other.start
            else:
                return self.end < other.end

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        if self.chrom == other.chrom and self.strand == other.strand and \
                self.start == other.start and self.end == other.end:
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
    This class has the purpose of parsing a reference file, calculating the statistics, and printing them out.
    """

    def __init__(self, parsed_args: argparse.Namespace):

        """Constructor function"""

        self.gff = parsed_args.gff
        if type(self.gff) is GFF.GFF3:
            self.is_gff = True
        else:
            self.is_gff = False
        self.only_coding = parsed_args.only_coding
        self.out = parsed_args.out
        self.genes = dict()

    def parse_input(self):
        """
        Method to parse the input GTF/GFF file.
        """
        transcript2gene = dict()

        derived_features = set()

        for record in self.gff:
            if record.header is True:
                continue

            if record.is_derived is True and record.is_gene is False:
                derived_features.add(record.id)
            elif record.feature == "Locus" or record.is_gene is True or (
                    record.is_parent is True and record.is_transcript is False):
                self.genes[record.id] = GeneObject(record, only_coding=self.only_coding)
            elif record.is_transcript is True:
                if record.parent is None:
                    raise TypeError("No parent found for:\n{0}".format(str(record)))
                transcript2gene[record.id] = record.parent[0]
                if self.is_gff is False:
                    new_record = record.copy()
                    new_record.feature = "gene"
                    if new_record.gene not in self.genes:
                        self.genes[new_record.gene] = GeneObject(new_record, only_coding=self.only_coding)
                self.genes[record.parent[0]].transcripts[record.id] = TranscriptComputer(record)
            else:
                for parent in filter(lambda pparent: pparent not in derived_features, record.parent):
                    gid = transcript2gene[parent]
                    self.genes[gid].transcripts[parent].add_exon(record)

        for gid in self.genes:
            self.genes[gid].finalize()

    def __call__(self):

        self.parse_input()
        ordered_genes = sorted(self.genes.values())
        self.coding_genes = list(filter(lambda g: g.is_coding is True, ordered_genes))
        distances = []
        for index, gene in enumerate(ordered_genes[:-1]):
            next_gene = ordered_genes[index + 1]
            if next_gene.chrom != gene.chrom:
                continue
            distances.append(next_gene.start + 1 - gene.end)

        self.distances = numpy.array(distances)

        distances = []
        for index, gene in enumerate(self.coding_genes[:-1]):
            next_gene = self.coding_genes[index + 1]
            if next_gene.chrom != gene.chrom:
                continue
            distances.append(next_gene.start + 1 - gene.end)

        self.coding_distances = numpy.array(distances)

        self.writer()

    @staticmethod
    def get_stats(row: dict, array: numpy.array) -> dict:
        """
        Method to calculate the necessary statistic from a row of values.
        :param row: the output dictionary row.
        :type row: dict

        :param array: an array of values.

        :rtype : dict
        """

        row["Average"] = "{0:,.2f}".format(round(numpy.mean(array), 2))  # Decimal to second digit precision

        counter_object = Counter(array)
        moder = [x for x in counter_object if counter_object[x] == counter_object.most_common(1)[0][1]]
        row["Mode"] = ";".join(str(x) for x in moder)
        quantiles = [numpy.percentile(array, x) for x in [0, 5, 10, 25, 50, 75, 90, 95, 100]]
        for key, val in zip(['Min', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', 'Max'], quantiles):
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

    def writer(self):
        """Method which creates the final output"""

        import csv
        fieldnames = ['Stat', 'Total', 'Average', 'Mode', 'Min', '5%', '10%', '25%', 'Median', '75%', '90%', '95%',
                      'Max']
        rower = csv.DictWriter(self.out,
                               fieldnames,
                               delimiter="\t")

        rower.writeheader()
        row = dict()
        for key in fieldnames:
            row[key] = "NA"
        row["Stat"] = 'Number of genes'
        row['Total'] = len(self.genes)
        rower.writerow(row)

        row["Stat"] = "Number of genes (coding)"
        row["Total"] = len(self.coding_genes)
        rower.writerow(row)

        row["Stat"] = 'Number of transcripts'
        row['Total'] = sum(self.genes[x].num_transcripts for x in self.genes)
        rower.writerow(row)

        row["Stat"] = 'Transcripts per gene'
        t_per_g = numpy.array(list(self.genes[x].num_transcripts for x in self.genes))
        row = self.get_stats(row, t_per_g)
        rower.writerow(row)

        row["Stat"] = "Number of coding transcripts"
        row["Total"] = sum(len(x.coding_transcripts) for x in self.coding_genes)
        rower.writerow(row)

        row["Stat"] = "Coding transcripts per gene"
        t_per_g = numpy.array(list(len(x.coding_transcripts) for x in self.coding_genes))
        row = self.get_stats(row, t_per_g)
        rower.writerow(row)

        exons = cArray('i')  # Done
        exons_coding = cArray('i')
        exon_num = cArray('i')  # Done
        exon_num_coding = cArray('i')
        introns = cArray('i')  # Done
        introns_coding = cArray('i')
        cds_introns = cArray('i')  # Done
        cds_exons = cArray('i')  # Done
        cds_exon_num = cArray('i')  # Done
        cds_exon_num_coding = cArray('i')
        cdna_lengths = cArray('i')  # Done
        cdna_lengths_coding = cArray('i')
        cds_lengths = cArray('i')  # Done
        monoexonic_lengths = cArray('i')
        multiexonic_lengths = cArray('i')
        monocds_lengths = cArray('i')

        five_utr_lengths = cArray('i')
        five_utr_nums = cArray('i')

        three_utr_lengths = cArray('i')
        three_utr_nums = cArray('i')

        end_distance_from_junction = cArray('i')

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
                if cds_num == 1 or (exon_num == 1 and self.genes[gene].transcripts[tid].selected_cds_length > 0):
                    monocds_lengths.append(self.genes[gene].transcripts[tid].selected_cds_length)
                cds_exon_num.append(cds_num)
                cdna_lengths.append(self.genes[gene].transcripts[tid].cdna_length)
                cds_lengths.append(self.genes[gene].transcripts[tid].selected_cds_length)
                if self.genes[gene].transcripts[tid].selected_cds_length > 0:
                    five_utr_lengths.append(self.genes[gene].transcripts[tid].five_utr_length)
                    three_utr_lengths.append(self.genes[gene].transcripts[tid].three_utr_length)
                    five_utr_nums.append(self.genes[gene].transcripts[tid].five_utr_num)
                    three_utr_nums.append(self.genes[gene].transcripts[tid].three_utr_num)
                    end_distance_from_junction.append(
                        self.genes[gene].transcripts[tid].selected_end_distance_from_junction)

                if self.only_coding is False and self.genes[gene].transcripts[tid].selected_cds_length > 0:
                    cdna_lengths_coding.append(self.genes[gene].transcripts[tid].cdna_length)
                    exons_coding.extend(self.genes[gene].transcripts[tid].exon_lengths)
                    exon_num_coding.append(len(self.genes[gene].transcripts[tid].exon_lengths))
                    cds_exon_num_coding.append(len(self.genes[gene].transcripts[tid].cds_exon_lengths))
                    introns_coding.extend(self.genes[gene].transcripts[tid].intron_lengths)

        row["Stat"] = 'CDNA lengths'
        row["Total"] = 'NA'
        ar = cdna_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)

        row["Stat"] = 'CDNA lengths (mRNAs)'
        row["Total"] = 'NA'
        ar = cdna_lengths_coding
        row = self.get_stats(row, ar)
        rower.writerow(row)

        row["Stat"] = 'CDS lengths'
        row["Total"] = 'NA'
        ar = cds_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = "CDS lengths (mRNAs)"
            row["Total"] = 'NA'
            ar = cArray('i', list(filter(lambda x: x > 0, cds_lengths)))
            row = self.get_stats(row, ar)
            rower.writerow(row)

        row["Stat"] = 'Monoexonic transcripts'
        row['Total'] = len(monoexonic_lengths)
        ar = monoexonic_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)

        row["Stat"] = 'MonoCDS transcripts'
        row['Total'] = len(monocds_lengths)
        ar = monoexonic_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)

        row["Stat"] = 'Exons per transcript'
        ar = exon_num
        row = self.get_stats(row, ar)
        row["Total"] = len(exons)
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = 'Exons per transcript (mRNAs)'
            ar = exon_num_coding
            row = self.get_stats(row, ar)
            row["Total"] = len(exons)
            rower.writerow(row)

        row["Stat"] = 'Exon lengths'
        ar = exons
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = 'Exon lengths (mRNAs)'
            ar = exons_coding
            row = self.get_stats(row, ar)
            row["Total"] = "NA"
            rower.writerow(row)

        row["Stat"] = "Intron lengths"
        ar = introns
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = "Intron lengths (mRNAs)"
            ar = introns_coding
            row = self.get_stats(row, ar)
            row["Total"] = "NA"
            rower.writerow(row)

        row["Stat"] = "CDS exons per transcript"
        ar = cds_exon_num
        row = self.get_stats(row, ar)
        row["Total"] = len(cds_exons)
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = "CDS exons per transcript (mRNAs)"
            ar = cds_exon_num_coding
            row = self.get_stats(row, ar)
            row["Total"] = len(cds_exons)
            rower.writerow(row)

        row["Stat"] = "CDS exon lengths"
        ar = cds_exons
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "CDS Intron lengths"
        ar = cds_introns
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "5'UTR exon number"
        ar = five_utr_nums
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "3'UTR exon number"
        ar = three_utr_nums
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "5'UTR length"
        ar = five_utr_lengths
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "3'UTR length"
        ar = three_utr_lengths
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "Stop distance from junction"
        ar = end_distance_from_junction
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "Intergenic distances"
        ar = self.distances
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)

        row["Stat"] = "Intergenic distances (coding)"
        ar = self.coding_distances
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)


def launch(args):

    """
    Very simple launcher function, calls Calculator from this module.
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
