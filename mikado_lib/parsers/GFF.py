#!/usr/bin/env python3
# coding: utf_8


"""
Module to serialize GFF files.
"""

from . import Parser
from .gfannotation import GFAnnotation


# This class has exactly how many attributes I need it to have
# pylint: disable=too-many-instance-attributes
class GffLine(GFAnnotation):
    """Object which serializes a GFF line."""

    def __init__(self, line: str, my_line='', header=False):
        """
        Constructor method.
        :param line: the GFF line to be serialised
        :type line: str,None

        :param my_line: optional string to be passed along. For exotic uses of the constructor.
        :type my_line: str

        :param header: boolean flag that indicates whether the instance will be a header or not.
        :type header: bool
        """

        GFAnnotation.__init__(self, line, my_line, header=header)
        _ = self.name  # Set the name

    def _parse_attributes(self):

        """
        Private method that parses the last field of the GFF line.
        :return:
        """

        self.attribute_order = []

        for item in iter(x for x in self._attr.rstrip().split(';') if x != ''):
            itemized = item.strip().split('=')
            try:
                if itemized[0].lower() == "parent":
                    self.parent = itemized[1].split(",")

                elif itemized[0].upper() == "ID":
                    self.id = itemized[1]
                else:
                    self.attributes[itemized[0]] = itemized[1]
                    self.attribute_order.append(itemized[0])
            except IndexError:
                pass

    def _sort_feature(self, feature):
        """
        Private method that sorts features according to the normal order in a GF file.
        :param feature:
        :return: numeric sort index
        """

        if self.strand == "-":
            order = ["three_prime_utr",
                     "exon",
                     "stop_codon",
                     "CDS",
                     "start_codon",
                     "five_prime_utr"]
        else:
            order = ["five_prime_utr",
                     "exon",
                     "start_codon",
                     "CDS",
                     "stop_codon",
                     "three_prime_utr"]
        if feature not in order:
            return float("inf")
        else:
            return order.index(feature)

    def _format_attributes(self):
        """
        Implementation of the abstract method for formatting the
        ninth field of a GFF.
        :return:
        """

        attrs = []
        if self.id is not None:
            attrs.append("ID={0}".format(self.id))
        if len(self.parent) > 0:
            attrs.append("Parent={0}".format(",".join(self.parent)))
        if "Name" in self.attributes and self.attributes["Name"] is not None:
            if "name" in self.attributes:
                del self.attributes["name"]
            if "name" in self.attribute_order:
                self.attribute_order.remove("name")
            attrs.append("Name={0}".format(self.name))
        elif "name" in self.attributes and self.attributes["name"] is not None:
            self.name = self.attributes["name"]
            del self.attributes["name"]
            if "name" in self.attribute_order:
                self.attribute_order.remove("name")
            attrs.append("Name={0}".format(self.name))

        if not self.attribute_order:
            self.attribute_order = sorted(list(key for key in self.attributes if
                                               key not in ["ID", "Parent", "Name"]))

        for att in iter(key for key in self.attribute_order if
                        key not in ["ID", "Parent", "Name"]):
            if att in ("gene_id", "transcript_id"):
                continue  # These are carryovers from GTF files
            if self.attributes[att] is not None:
                # try:
                attrs.append("{0}={1}".format(att.lower(), self.attributes[att]))
                # except KeyError:
                #     # Hack for those times when we modify the attributes at runtime
                #     continue
        attrs = ";".join(attrs)
        return attrs

    # id is and id remains
    # pylint: disable=invalid-name
    @property
    def id(self):
        """
        Returns the ID of the feature.
        :rtype str
        """
        return self.attributes["ID"]

    @id.setter
    def id(self, newid):
        """
        :param newid: new id of the feature
        :type newid: str
        Setter for instance id.
        """

        self.attributes["ID"] = newid
    # pylint: enable=invalid-name

    @property
    def parent(self):
        """This property looks up the "Parent" field in the "attributes" dictionary.
        Contrary to other attributes, this property returns a *list*, not a string.
        This is due to the fact that GFF files support
        multiple inheritance by separating the parent entries with a comma.

        :rtype list

        """
        if "Parent" not in self.attributes:
            self.parent = None
        return self.attributes["Parent"]

    @parent.setter
    def parent(self, parent):
        """
        Setter for the parent attribute. It will be converted into a list if it is not already.

        :param parent: the new parent attribute
        :type parent: str
        :type parent: list
        """

        if isinstance(parent, str):
            parent = parent.split(",")
        self.attributes["Parent"] = parent

    @property
    def name(self):
        """
        Returns the name of the feature. It defaults to the ID if missing.

        :rtype str
        """

        if "Name" not in self.attributes:
            self.name = self.id
        return self.attributes["Name"]

    @name.setter
    def name(self, name):
        """
        Setter for the name attribute.
        :param name
        :type name: str
        """

        self.attributes["Name"] = name

    @property
    def is_transcript(self):
        """
        Property. True if the feature is an RNA, false otherwise.
        :rtype bool
        """

        if self.feature is None:
            return False
        if self.is_gene is True:
            return False
        if self.feature.endswith("transcript") or "RNA" in self.feature.upper():
            return True
        return False

    @property
    def is_gene(self):
        """
        Property. True if the feature ends with "gene" and has no parents.
        :rtype bool
        """

        if self.feature is not None:
            if self.feature.endswith("gene"):
                return True
        return False

    @property
    def is_parent(self):
        """
        Property. True if the feature has no parent defined.
        :rtype bool
        """

        if len(self.parent) == 0:
            return True
        return False

    @property
    def gene(self):
        """
        Property. If the feature is a transcript,
        returns the parent; else, it returns None
        (as only transcript lines have the gene among the attributes)
        :rtype str
        """

        if self.is_transcript is True:
            return self.parent[0]
        else:
            return None

    @property
    def transcript(self):
        """
        Property. If the feature is a transcript, it returns the id;
        if it is an exon, it returns the parent; else, it returns None

        :rtype str
        :rtype list
        :rtype None
        """
        if self.is_exon is True:
            return self.parent
        elif self.is_transcript is True:
            return self.id
        else:
            return None

    @transcript.setter
    def transcript(self, string):
        """
        Setter for the transcript property. It performs a light sanity check on the new
        assigned value.
        :param string:
        :return:
        """
        if self.is_exon is True:
            self.parent = string
        elif self.is_transcript is True:
            self.id = string
        else:
            raise TypeError("Cannot set the transcript name of a non-transcript feature!")

    @property
    def is_derived(self):
        """
        Property. It checks whether there is a "Derives_from" attribute among the line attributes.
        :rtype bool
        """
        return "Derives_from" in self.attributes

    @property
    def derived_from(self):
        """
        Property, Set to True if the key "Derives_from" is present among the
        instance attributes.
        """

        if self.is_derived is False:
            return None
        else:
            return self.attributes["Derives_from"].split(",")


class GFF3(Parser):
    """
    Class that is used to parse a GFF file.
    """

    __annot_type__ = "gff3"

    def __init__(self, handle):
        """
        Constructor method.
        :param handle: the input file. It can be a file handle or a file name.
        :type handle: io.TextIOWrapper | str
        """
        super().__init__(handle)
        self.header = False

    def __next__(self):

        if self.closed:
            raise StopIteration

        line = self._handle.readline()
        if line == '':
            raise StopIteration

        if line[0] == "#":
            return GffLine(line, header=True)

        line = GffLine(line)
        return line
