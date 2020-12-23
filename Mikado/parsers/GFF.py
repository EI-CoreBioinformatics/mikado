#!/usr/bin/env python3
# coding: utf_8


"""
Module to serialize GFF files.
"""

from . import Parser
from .gfannotation import GFAnnotation, _attribute_definition
from sys import intern
import re


# This class has exactly how many attributes I need it to have
# pylint: disable=too-many-instance-attributes
class GffLine(GFAnnotation):
    """Object which serializes a GFF line."""

    # The (?:;|$) means "match, but **do not capture**, either semicolon or end of the line.
    _attribute_pattern = re.compile(r"([^;]*)=([^$=]*)(?:;|$)")

    def __init__(self, line, my_line='', header=False):
        """
        Constructor method.
        :param line: the GFF line to be serialised
        :type line: (str|None)

        :param my_line: optional string to be passed along. For exotic uses of the constructor.
        :type my_line: str

        :param header: boolean flag that indicates whether the instance will be a header or not.
        :type header: bool
        """

        GFAnnotation.__init__(self, line, my_line, header=header)
        self.__set_name()  # Set the name

    def _parse_attributes(self):

        """
        Private method that parses the last field of the GFF line.
        :return:
        """

        infolist = self._attribute_pattern.findall(self._attr.rstrip().rstrip(";"))
        attribute_order = [key for key, val in infolist if key not in ("Parent", "parent", "id", "ID", "Id")]
        attributes = dict((key, _attribute_definition(val)) for key, val in infolist)

        # Ensure that the "parent" attribute is a string
        if "Parent" in attributes:
            attributes["Parent"] = str(attributes["Parent"])
            self.parent = attributes["Parent"]
        elif "parent" in attributes:
            attributes["parent"] = str(attributes["parent"])
            self.parent = attributes["parent"]

        # Ensure that the "ID" attribute is a string
        if "ID" in attributes:
            attributes["ID"] = str(attributes["ID"])
            self.id = attributes["ID"]
        elif "id" in attributes:
            attributes["id"] = str(attributes["id"])
            self.id = attributes["id"]
        elif "Id" in attributes:
            attributes["Id"] = str(attributes["Id"])
            self.id = attributes["Id"]

        self.attributes.update(attributes)
        self.attribute_order = attribute_order

    def _format_attributes(self):
        """
        Implementation of the abstract method for formatting the
        ninth field of a GFF.
        :return:
        """

        attrs = self._format_attributes_dict(
            self.attributes,
            parent=self.parent,
            mid=self.id,
            name=self.name,
            attribute_order=self.attribute_order)

        return attrs

    @staticmethod
    def _format_attributes_dict(attributes, parent=None, mid=None, name=None, attribute_order=None):

        if not parent:
            try:
                parent = attributes["parent"]
            except KeyError:
                try:
                    parent = attributes["Parent"]
                except KeyError:
                    parent = None

        if not mid:
            try:
                mid = attributes["id"]
            except KeyError:
                try:
                    mid = attributes["ID"]
                except KeyError:
                    mid = None
        if not name:
            name = attributes.get("name", None)
        attrs = []
        if mid is not None:
            attrs.append("ID={0}".format(mid))
        if parent is not None:
            if isinstance(parent, str):
                parent = [parent]
            if len(parent) > 1:
                attrs.append("Parent={0}".format(",".join(parent)))
            elif len(parent) == 1:
                attrs.append("Parent={0}".format(parent[0]))
        if name is not None:
            attrs.append("Name={0}".format(name))

        if not attribute_order:
            attribute_order = sorted(list(key for key in attributes if
                                               key not in ["ID", "Parent", "Name",
                                                           "parent", "id", "name"]))

        attrs += ["{0}={1}".format(att.lower(), attributes[att]) for att in attribute_order
                  if (att not in ("ID", "Parent", "Name", "name", "gene_id", "transcript_id") and
                      attributes[att] is not None)]

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
        if "ID" in self.attributes:
            return str(self.attributes["ID"])
        else:
            return None

    @id.setter
    def id(self, newid):
        """
        :param newid: new id of the feature
        :type newid: str
        Setter for instance id.
        """

        self.attributes["ID"] = newid
        self._is_gene = None
    # pylint: enable=invalid-name

    @property
    def parent(self):
        """This property looks up the "Parent" field in the "attributes" dictionary.
        Contrary to other attributes, this property returns a *list*, not a string.
        This is due to the fact that GFF files support
        multiple inheritance by separating the parent entries with a comma.

        :rtype list

        """
        return self.__parent

    def __set_parent(self):

        if "Parent" not in self.attributes:
            self.parent = None
        self.parent = self.attributes["Parent"]

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
        self.__parent = parent

    @property
    def name(self):
        """
        Returns the name of the feature. It defaults to the ID if missing.

        :rtype str
        """

        return self.__name

    def __set_name(self):

        if "Name" not in self.attributes:
            self.name = None
        self.__name = self.attributes["Name"]

    @name.setter
    def name(self, name):
        """
        Setter for the name attribute.
        :param name
        :type name: str
        """

        if not isinstance(name, (type(None), str)):
            raise TypeError("Invalid value for name: {0}".format(name))
        self.attributes["Name"] = name
        self.__name = name

    @property
    def is_transcript(self):
        """
        Property. True if the feature is an RNA, false otherwise.
        :rtype bool
        """

        if self.feature is None:
            return False
        if self.is_gene is True and self.feature != "mRNA_TE_gene":
            return False
        if self.feature.endswith("transcript") or "RNA" in self.feature.upper():
            return True
        elif self.id is not None and "transcript:" in self.id and self.parent is not None:
            return True
        elif self.feature.endswith("_gene_segment"):  # Necessary for V_gene_segment, C_gene_segment, etc.
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

        if self.is_transcript is True and self.parent:
            return self.parent[0]
        elif self.is_gene:
            return self.id
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
        return any(("derives_from" in _.lower()) for _ in self.attributes.keys())

    @property
    def derived_from(self):
        """
        Property, Set to True if the key "Derives_from" is present among the
        instance attributes.
        """

        if self.is_derived is False:
            return None
        else:
            if "Derives_from" not in self.attributes:
                key = [_ for _ in self.attributes.keys() if _.lower() == "derives_from"]
                assert len(key) == 1, (str(self), key)
                self.attributes["Derives_from"] = self.attributes[key[0]]

            return self.attributes["Derives_from"].split(",")

    @property
    def _negative_order(self):
        negative_order = [intern(_) for _ in ["three_prime_utr",
                                              "exon",
                                              "stop_codon",
                                              "CDS",
                                              "start_codon",
                                              "five_prime_utr"]]
        return negative_order

    @property
    def _positive_order(self):
        positive_order = [intern(_) for _ in
                          ["five_prime_utr",
                           "exon",
                           "start_codon",
                           "CDS",
                           "stop_codon",
                           "three_prime_utr"]]
        return positive_order


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
        line = next(self._handle)

        if line[0] == "#":
            return GffLine(line, header=True)

        try:
            gff_line = GffLine(line)
        except Exception:
            error = "Invalid line for file {}, position {}:\n{}".format(
                self.name, self._handle.tell(), line)
            raise ValueError(error)
        return gff_line

    @property
    def file_format(self):
        return self.__annot_type__
