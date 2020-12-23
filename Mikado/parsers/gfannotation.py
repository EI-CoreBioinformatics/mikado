#!/usr/local/bin/env python3

"""
This module describes an abstract class which underlies the parsing of
GTF/GFF files.
"""

import abc
import copy
from sys import intern
import re


__author__ = 'Luca Venturini'

[intern(_) for _ in ["+", "-", "?", "true", "True", "false", "False"]]


def _attribute_definition(val):
    try:
        val = float(val)
        if val.is_integer():
            return int(val)
        return val
    except (ValueError, TypeError):
        if val.lower() in ("true", "false"):
            val = val.capitalize()
            if val == "True":
                return True
            else:
                return False
        return val


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
        self.__strand = None
        self.header = True
        self.chrom, self.source, self.feature = None, None, None
        # pylint: disable=invalid-name
        self.id = None
        # pylint: enable=invalid-name
        self.parent = []
        self.start, self.end = None, None
        self.__score = None
        self.__phase = None
        self.__frame = None
        self._line = "NA"
        self.__gene = None
        self._transcript = None
        self.__feature = None

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
            self.__feature = None
            self.header = True
            return

        self.chrom, self.source = self._fields[0:2]
        try:
            self.start, self.end = tuple(int(i) for i in self._fields[3:5])
        except (ValueError, TypeError):
            error = "Invalid start and end values: {}\n".format(" ".join(self._fields[3:5]))
            error += "Line: {}".format(self._line)
            raise ValueError(error)

        self.score = self._fields[5]
        self.strand = self._fields[6]
        self.phase = self._fields[7]

        self._attr = self._fields[8]
        self._parse_attributes()
        self.feature = self._fields[2]
        self.__is_exon, self.__is_gene, self.__is_cds = None, None, None
        [intern(_) for _ in (str(self.chrom), str(self.source), str(self.feature))]

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
        if self.header is False and all([_ is not None for _ in [self.end, self.start]]):
            return self.end - self.start + 1
        else:
            return 0

    def __getstate__(self):
        state = {
            "chrom": self.chrom,
            "source": self.source,
            "feature": self.feature,
            "start": self.start,
            "end": self.end,
            "score": self.score,
            "strand": self.strand,
            "phase": self.phase,
            "attributes": self.attributes,
            "is_exon": self.is_exon,
            "is_cds": self.is_cds,
            "is_gene": self.is_gene,
            "transcript": self._transcript,
            "gene": self.__gene,
            "header": self.header,
            "_line": self._line,
            "_fields": self._fields
        }
        return state

    def __setstate__(self, state):
        self._transcript = state.pop("transcript")
        self.__gene = state.pop("gene")
        self.__is_gene = state.pop("is_gene")
        self.__is_cds = state.pop("is_cds")
        self.__is_exon = state.pop("is_exon")
        self.__dict__.update(state)

    def as_dict(self):
        return self.__getstate__()

    def load_dict(self, state):
        self.__setstate__(state)

    @classmethod
    def string_from_dict(cls, data, **kwargs):

        line = [data["chrom"],
                data["source"] if data["source"] else "Mikado",
                data["feature"],
                data["start"],
                data["end"],
                data["score"] if data["score"] is not None else ".",
                data["strand"] if data["strand"] is not None else ".",
                data["phase"] if data["phase"] is not None else "."
                ]

        line = [str(item) for item in line]
        attrs = cls._format_attributes_dict(data["attributes"], **kwargs)
        return "\t".join(line + [attrs])

    @abc.abstractmethod
    def _parse_attributes(self):
        """
        Abstract method. It is used by the children class to serialise the data
         inside the ninth field, which is one of the biggest differences between
         GFFs and GTFs.
        :return:
        """
        raise NotImplementedError("This is only an abstract method!")

    @staticmethod
    def _attribute_definition(val):
        return _attribute_definition(val)

    @abc.abstractmethod
    def _format_attributes(self):
        """ Abstract method. Each children class should implement its own version,
        which allow to prepare the ninth field for printing.
        :return:
        """
        raise NotImplementedError("This is only an abstract method!")

    @staticmethod
    @abc.abstractmethod
    def _format_attributes_dict(attributes, **kwargs):
        """Method to define how to get the attribute string."""

    def __format_middle(self):
        """
        Private method to format the middle fields (score, strand, phase)
        for printing.
        :return: score, strand, phase
        :rtype: str, str, str
        """

        score = self.score
        strand = self.strand
        phase = self.phase

        if score is None:
            score = "."
        else:
            score = int(score)
        if strand is None:
            strand = '.'
        else:
            strand = strand
        if phase is None:
            phase = "."
        return score, strand, phase

    def copy(self):
        """
        Wrapper around the copy.copy function.
        """
        return copy.copy(self)

    @property
    def feature(self):
        try:
            return self.__feature
        except AttributeError:
            self.__feature = self.__dict__["feature"]
            return self.__feature

    @feature.setter
    def feature(self, feature):
        self.__feature = feature
        self.__is_exon, self.__is_cds, self._is_gene = None, None, None

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
        try:
            return self.__score
        except AttributeError:
            self.__score = self.__dict__["score"]
            return self.__score

    @score.setter
    def score(self, *args):
        """
        Score setter. It verifies that the score is a valid value (None or float).
        :param args: List of arguments. Only the first is used.
        :type args: list | None | float
        """

        score = args[0]
        if score == ".":
            score = None
        elif score is not None:
            try:
                score = float(args[0])
            except ValueError:
                score = None
            except SystemError:
                raise ValueError("Invalid score value: {}".format(args[0]))
        if score is not None and not isinstance(score, float):
            raise TypeError(score)
        self.__score = score

    @property
    def phase(self):
        """
        Property. Stores the phase of the feature.
        Valid values: None, 0, 1 ,2
        :return:
        """

        try:
            return self.__phase
        except AttributeError:
            self.__phase = self.__dict__["phase"]
            return self.__phase

    @phase.setter
    def phase(self, value):
        """
        Setter for the phase attribute.
        :param value:
        :return:
        """

        if value is None:
            self.__phase = None
        else:
            if value in (None, '.', '?'):
                self.__phase = None
            else:
                try:
                    phase = int(value)
                except (TypeError, ValueError):
                    raise ValueError("Invalid phase {0} (type: {1})\nLine:{2}".format(value, type(value),
                                                                                      self._line))
                if phase in (-1, 0, 1, 2):
                    phase = max(phase, 0)
                    self.__phase = phase
                else:
                    raise ValueError("Invalid phase {0} (type: {1})\nLine:{2}".format(value, type(value),
                                                                                      self._line))
        self._set_frame()

    def _set_frame(self):
        if self.phase is not None:
            self.__frame = (3 - self.__phase) % 3
        else:
            self.__frame = None

    @property
    def frame(self):
        return self.__frame

    @property
    def is_gene(self):
        """
        Property. True if the feature ends with "gene" and has no parents.
        :rtype bool
        """
        if self._is_gene is None:
            self._is_gene = self._set_is_gene()
        return self._is_gene

    def _set_is_gene(self):
        if self.feature is not None:
            if self.feature.endswith("gene") and self.feature != "mRNA_TE_gene":
                # Hack for EnsEMBL GFFs
                try:
                    if self.id is None or "transcript:" not in self.id:
                        return True
                    else:
                        return False
                except TypeError:
                    raise TypeError((self.id, type(self.id)))
            elif self.id is not None and self.id.startswith("gene:"):
                # Hack for EnsEMBL
                return True
            elif self.feature in ("sublocus", "monosublocus", "monosublocusholder"):
                return True
        return False

    @property
    def is_exon(self):
        """
        Property. True if the feature is CDS/exon/UTR/start or stop codon.
        :rtype bool
        """
        if self.__is_exon is None:
            self.__is_exon = self._set_is_exon()

        return self.__is_exon

    exon_pattern = re.compile("(cds|exon$|utr|codon|^cdna_match$|^match_part$)", flags=re.IGNORECASE)

    def _set_is_exon(self):

        feature = self.feature
        return (feature is not None) and (re.search(self.exon_pattern, feature) is not None)

    @property
    def is_cds(self):
        """This property evaluates to True if the row describes a CDS/UTR segment,
        False otherwise."""
        if self.__is_cds is None:
            self.__is_cds = self._set_is_cds()
        return self.__is_cds

    cds_pattern = re.compile("(cds|utr|codon)", flags=re.IGNORECASE)

    def _set_is_cds(self):
        return self.is_exon and (re.search(self.cds_pattern, self.feature) is not None)

    @property
    @abc.abstractmethod
    def _negative_order(self):
        raise NotImplementedError("This is implemented in the children classes")

    @property
    @abc.abstractmethod
    def _positive_order(self):
        raise NotImplementedError("This is implemented in the children classes")

    def _sort_feature(self, feature):
        """
        Private method that sorts features according to the normal order in a GF file.
        :param feature:
        :return: numeric sort index
        """

        try:
            strand = self.strand
        except AttributeError:
            strand = self.__dict__["strand"]

        if strand == "-":
            order = self._negative_order
        else:
            order = self._positive_order
        assert order != []

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

    def remove_attribute(self, attribute):
        """Method to remove attributes from the internal dictionary."""

        if attribute in self.attributes:
            del self.attributes[attribute]
        try:
            self.attribute_order.remove(attribute)
        except ValueError:
            pass
        return

    def add_attribute(self, attribute, value):

        self.attributes[attribute] = value

        if attribute not in self.attribute_order:
            self.attribute_order.append(attribute)

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
