#!/usr/local/bin/env python3

"""
This module describes an abstract class which underlies the parsing of
GTF/GFF files.
"""

import abc
import copy
from sys import intern

__author__ = 'Luca Venturini'

[intern(_) for _ in ["+", "-", "?"]]


# This class has exactly how many attributes I need it to have
# pylint: disable=too-many-instance-attributes
class GFAnnotation(metaclass=abc.ABCMeta):

    __negative_order = []
    __positive_order = []


    """
    This abstract class describes a generic GTF/GFF annotation line.
    The parsers for those two type of files inherit from this abstract class,
    which defines common methods and properties.
    """

    @abc.abstractmethod
    def __init__(self, line, my_line='', header=False):
        self.attributes = dict()
        self.header = True
        self.chrom, self.source, self.feature = None, None, None
        # pylint: disable=invalid-name
        self.id = None
        # pylint: enable=invalid-name
        self.parent = []
        self.start, self.end = None, None
        self.__score = None
        self.__strand = None
        self.__phase = None
        self._line = "NA"
        self.__gene = None
        self._transcript = None

        self.attribute_order = []
        if line is None:  # Empty constructor
            return
        if line == '' and my_line == '':
            return

        if line == '' and my_line != "":
            self._line = my_line
        else:
            self._line = line

        self._line = self._line.strip()

        self._fields = self._line.split('\t')
        self.header = header

        if self.header or len(self._fields) != 9 or self._line == '' or self._line[0] == "#":
            self.feature = None
            self.header = True
            return

        self.chrom, self.source, self.feature = self._fields[0:3]
        [intern(_) for _ in (self.chrom, self.source, self.feature)]
        self.start, self.end = tuple(int(i) for i in self._fields[3:5])

        self.score = self._fields[5]
        self.strand = self._fields[6]
        self.phase = self._fields[7]

        self._attr = self._fields[8]
        self._parse_attributes()

    def __str__(self):
        if not self.feature:
            return self._line.rstrip()
        score, strand, phase = self.__format_middle()
        attributes = self._format_attributes()
        if self.source is None:
            self.source = "Mikado"
        line = [self.chrom, self.source, self.feature,
                self.start, self.end, score,
                strand, phase, attributes]
        return "\t".join([str(_) for _ in line])

    def __len__(self):
        if self.end is not None:
            return self.end - self.start + 1
        else:
            return 0

    @abc.abstractmethod
    def _parse_attributes(self):
        """
        Abstract method. It is used by the children class to serialise the data
         inside the ninth field, which is one of the biggest differences between
         GFFs and GTFs.
        :return:
        """
        raise NotImplementedError("This is only an abstract method!")

    @abc.abstractmethod
    def _format_attributes(self):
        """ Abstract method. Each children class should implement its own version,
        which allow to prepare the ninth field for printing.
        :return:
        """
        raise NotImplementedError("This is only an abstract method!")

    def __format_middle(self):
        """
        Private method to format the middle fields (score, strand, phase)
        for printing.
        :return: score, strand, phase
        :rtype: str, str, str
        """

        if self.score is not None:
            score = str(int(round(self.score, 0)))
        else:
            score = "."
        if self.strand is None:
            strand = '.'
        else:
            strand = self.strand
        if self.phase is not None:
            phase = str(self.phase)
        else:
            phase = "."
        return score, strand, phase

    def copy(self):
        """
        Wrapper around the copy.copy function.
        """
        return copy.copy(self)

    # Instance properties
    @property
    def strand(self):

        """
        Strand attribute. One of None, +, -
        :rtype str | None
        """

        return self.__strand

    @strand.setter
    def strand(self, strand):
        """
        Strand setter. It verifies that strand is a valid value (None, +, -)
        :param strand: new strand value
        :type strand: str | None
        """

        if strand in ("+", "-"):
            self.__strand = strand
        elif strand in (None, ".", "?"):
            self.__strand = None
        else:
            raise ValueError("Invalid value for strand: {0}".format(strand))

    @property
    def score(self):
        """
        Score value. Either None or float.
        :rtype float | None
        """

        return self.__score

    @score.setter
    def score(self, *args):
        """
        Score setter. It verifies that the score is a valid value (None or float).
        :param args: List of arguments. Only the first is used.
        :type args: list | None | float
        """

        if isinstance(args[0], (float, int)):
            self.__score = args[0]
        elif args[0] is None or args[0] == '.':
            self.__score = None
        elif isinstance(args[0], str):
            self.__score = float(args[0])
        else:
            raise TypeError(args[0])

    @property
    def phase(self):
        """
        Property. Stores the phase of the feature.
        Valid values: None, 0, 1 ,2
        :return:
        """

        return self.__phase

    @phase.setter
    def phase(self, value):
        """
        Setter for the phase attribute.
        :param value:
        :return:
        """

        if value in (None, '.', '?'):
            self.__phase = None
        elif isinstance(value, (str, int, float)):
            value = int(value)
            if value not in range(3):
                raise ValueError("Invalid value for phase: {0}".format(value))
            self.__phase = value
        else:
            raise ValueError("Invalid phase: {0}".format(value))

    @property
    def is_gene(self):
        """
        Property. True if the feature ends with "gene" and has no parents.
        :rtype bool
        """

        if self.feature is not None:
            if self.feature.endswith("gene") and self.feature != "mRNA_TE_gene":
                # Hack for EnsEMBL GFFs
                if "transcript:" in self.id:
                    return False
                return True
            elif self.id is not None and self.id.startswith("gene:"):
                # Hack for EnsEMBL
                return True
        return False

    @property
    def is_exon(self):
        """
        Property. True if the feature is CDS/exon/UTR/start or stop codon.
        :rtype bool
        """

        if self.feature is None:
            return False
        _ = self.feature.lower()
        if ("cds" in _ or _.endswith("exon") or "utr" in _ or "codon" in _
            or _ == "cdna_match" or _ == "match_part"):
            return True
        return False

    @property
    def is_cds(self):
        """This property evaluates to True if the row describes a CDS/UTR segment,
        False otherwise."""
        if self.is_exon is False:
            return False
        _ = self.feature.lower()
        if "cds" in _:
            return True
        elif "utr" in _:
            return True
        elif _.endswith("codon"):
            return True
        return False

    def _sort_feature(self, feature):
        """
        Private method that sorts features according to the normal order in a GF file.
        :param feature:
        :return: numeric sort index
        """

        if self.strand == "-":
            order = self.__negative_order
        else:
            order = self.__positive_order
        if feature not in order:
            return float("inf")
        else:
            return order.index(feature)

    def __lt__(self, other):

        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        else:
            if self.start != other.start:
                return self.start < other.start
            elif self.end != other.end:
                return self.end < other.end
            elif self.feature != other.feature:
                return self._sort_feature(self.feature) < self._sort_feature(other.feature)
            else:
                return False

    def __eq__(self, other):

        if self.is_exon is True:
            return (self.chrom == other.chrom and
                    self.feature == other.feature and
                    self.start == other.start and
                    self.end == other.end)
        else:
            return (self.id == other.id and
                    self.feature == other.feature and
                    self.chrom == other.chrom and
                    self.start == other.start and
                    self.end == other.end)
