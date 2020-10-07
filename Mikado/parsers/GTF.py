# coding: utf_8

"""
Generic parser for GTF files.
"""

from . import Parser
from .gfannotation import GFAnnotation
import re


# This class has exactly how many attributes I need it to have
# pylint: disable=too-many-instance-attributes


class GtfLine(GFAnnotation):
    """This class defines a typical GTF line, with some added functionality
    to make it useful in e.g. parsing cufflinks GTF files or
    creating GTF lines from scratch.
    Fields:
    - chrom
    - source
    - feature
    - start,end
    - score
    - phase
    - strand
    - info (a dictionary containing all the annotations contained in the last field)

    - gene: the gene_id
    - transcript: the transcript_id

    For cufflinks files, also:
    - nearest_ref
    - tss_id
    - ccode"""

    # _slots=['chrom','source','feature','start',\
    # 'end','score','strand','phase','info']

    _attribute_pattern = re.compile(r"([^;\s]*) \"([^\"]*)\"(?:;|$)")

    def __init__(self, line, my_line='', header=False):

        self.__is_transcript = False
        super().__init__(line, my_line, header=header)
        self.__is_derived, self.__derived_from = None, False  # Placeholders
        self.__is_gene = self._set_is_gene()
        if not self.is_gene:
            self.__is_transcript = self._set_is_transcript()
            self.__set_parent()

    def _parse_attributes(self):
        """
        Method to retrieve the attributes from the last field
        of the GTF line.
        :return:
        """

        attributes = dict((key, self._attribute_definition(val))
                          for key, val in self._attribute_pattern.findall(self._attr.rstrip()))
        self.attributes.update(attributes)
        self.__name = self.__set_name()
        self.__set_gene()
        assert "gene_id" in self.attributes
        assert self.gene is not None, self.attributes

    def _format_attributes(self):

        """
        Private method to format the last field of the GTF line
        prior to printing.
        :return:
        """

        return self._format_attributes_dict(
            self.attributes,
            gene=self.gene,
            transcript=self.transcript)

    @staticmethod
    def _format_attributes_dict(attributes, gene=None, transcript=None):
        """Method to define how to get the attribute string."""

        if not gene:
            pass
            # assert attributes["gene_id"], attributes
        else:
            attributes['gene_id'] = gene
        if "gene_id" in attributes and isinstance(attributes["gene_id"], list):
            attributes["gene_id"] = ",".join(attributes["gene_id"])

        if not transcript:
            attributes["transcript_id"] = attributes.pop("transcript_id", None)
        else:
            attributes["transcript_id"] = transcript

        # assert attributes["transcript_id"]

        order = ['gene_id', 'transcript_id', 'exon_number', 'gene_name', 'transcript_name']

        info_list = []

        for tag in order:
            if tag in attributes:
                if isinstance(attributes[tag], list):
                    val = ",".join(attributes[tag])
                else:
                    val = attributes[tag]
                info_list.append("{0} \"{1}\"".format(tag, val))

        done = set()
        for info in iter(key for key in attributes if
                         attributes[key] not in (None, "", []) and key not in order):
            if info in ("Parent", "ID", "parent", "id"):
                continue

            if isinstance(attributes[info], list):
                val = ",".join(attributes[info])
            else:
                val = attributes[info]
            if info == "name":
                info = "Name"
            if info in done:
                continue
            done.add(info)
            info_list.append("{0} \"{1}\"".format(info, val))
        attrs = "; ".join(info_list) + ";"
        return attrs

    @property
    def name(self):
        """
        Returns the name of the feature. It defaults to the ID if missing.
        :rtype str
        """
        return self.__name

    def __set_name(self):

        if "Name" not in self.attributes:
            return self.id
        return self.attributes["Name"]

    @name.setter
    def name(self, *args):
        """
        Setter for name. The argument must be a string.
        :param args:
        :type args: list[(str)] | str
        """

        if not isinstance(args[0], (type(None), str)):
            raise TypeError("Invalid value for name: {0}".format(args[0]))
        self.attributes["Name"] = args[0]
        self.__name = args[0]

    @property
    def is_transcript(self):
        """
        Flag. True if feature is "transcript" or contains "RNA", False in all other cases.
        :rtype : bool
        """
        try:
            return self.__is_transcript
        except AttributeError:
            self.__is_transcript = self._set_is_transcript()
            return self.__is_transcript

    transcript_pattern = re.compile("(^transcript$|rna$)", re.IGNORECASE)

    def _set_is_transcript(self):

        # if self.feature is None:
        #     return False
        # if "transcript" == self.feature or "RNA" in self.feature:
        #     return True
        # return False

        return (self.feature is not None) and (self.transcript_pattern.search(
            self.feature) is not None)

    def __set_transcript(self):
        if self.header is True:
            self._transcript = None
        else:
            try:
                self._transcript = str(self.attributes["transcript_id"])
            except KeyError:
                raise KeyError(self.attributes)

    @property
    def parent(self):
        """This property looks up the "Parent" field in the "attributes" dictionary.
        If the line is a transcript line, it returns the gene field.
        Otherwise, it returns the transcript field.
        In order to maintain interface consistency with
        the GFF objects and contrary to other attributes,
        this property returns a *list*, not a string. This is due
        to the fact that GFF files support
        multiple inheritance by separating the parent entries with a comma.

        :rtype : list

        """

        return self.__parent

    def __set_parent(self):

        if self.is_gene:
            self.__parent = None
        self.__set_transcript()
        if self.is_transcript is True and self.gene is not None:
            self.__parent = [self.gene]
        elif self.is_transcript is True:
            raise ValueError("No gene")
        else:
            if self.transcript is not None:
                self.__parent = [self.transcript]

    @parent.setter
    def parent(self, parent):
        """
        Setter for parent. Acceptable values
        :param parent: the new parent value
        :type parent: list | str
        """
        if isinstance(parent, str):
            parent = parent.split(",")
        elif isinstance(parent, (int, float)):
            parent = [str(parent)]
        elif isinstance(parent, (bytes,)):
            parent = [parent.decode()]

        self.attributes["Parent"] = parent
        self.__parent = parent
        if self.is_transcript is True:
            self.gene = parent

    # pylint: disable=invalid-name
    @property
    def id(self):
        """
        ID of the line. "transcript_id" for transcript features, None in all other cases.
        :rtype : str | None
        """
        if self.is_transcript is True:
            return self.transcript
        else:
            return None

    @id.setter
    def id(self, newid):
        """
        Setter for id. Only transcript features can have their ID set.
        :param newid: the new ID
        """
        if self.is_transcript is True:
            self.transcript = newid
        else:
            pass
    # pylint: enable=invalid-name

    @property
    def gene(self):
        """
        Return the "gene_id" field.
        :rtype : str | None
        """

        # if "gene_id" not in self.attributes and self.is_transcript is True:
        #     self.attributes["gene_id"] = self.parent[0]
        if self.__gene is None:
            self.__set_gene()
        return self.__gene

    def __set_gene(self):
        try:
            self.__gene = str(self.attributes["gene_id"])
            if self.is_transcript:
                self.__parent = [self.__gene]
        except KeyError:
            pass

    @gene.setter
    def gene(self, gene):
        """
        Setter for gene.
        :param gene:
        :rtype : str

        """

        self.attributes["gene_id"] = self.__gene = gene
        if self.is_transcript:
            self.__parent = [gene]

    @property
    def transcript(self):
        """
        This property returns the "transcript_id" field of the GTF line.
        :rtype : str
        """
        return self._transcript

    @transcript.setter
    def transcript(self, transcript):
        """
        Setter for the transcript attribute. It also modifies the "transcript_id" field.
        :param transcript:
        :type transcript: str
        """

        self.attributes["transcript_id"] = self._transcript = str(transcript)
        if self.is_transcript is True:
            self.attributes["ID"] = transcript
        else:
            self.parent = [transcript]
        self._transcript = transcript

    @property
    def is_parent(self):
        """
        True if we are looking at a transcript, False otherwise.
        :rtype : bool
        """
        if any([self.is_transcript, self.is_gene]):
            return True
        return False

    @property
    def is_derived(self):
        """
        Property. It checks whether there is a "Derives_from" attribute among the line attributes.
        :rtype bool
        """
        if self.__is_derived is None:
            self.__is_derived = self.__set_is_derived()

        return self.__is_derived

    derived_pattern = re.compile("^derives_from$", flags=re.IGNORECASE)

    def __set_is_derived(self):
        return any(self.derived_pattern.search(key) is not None for key in self.attributes)

    @property
    def derived_from(self):
        """
        Boolean property. True if the GTF line has a "derives_from" tag,
        False otherwise.
        """
        if self.__derived_from is False:
            return self.__derived_from

    def __set_derived_from(self):

        if self.is_derived is False:
            return None
        else:
            for key in self.attributes:
                if re.search(self.derived_pattern, key) is not None:
                    return self.attributes[key].split(",")

    gene_pattern = re.compile("gene", re.IGNORECASE)

    @property
    def _negative_order(self):
        return ["3UTR",
                "exon",
                "stop_codon",
                "CDS",
                "start_codon",
                "5UTR"]

    @property
    def _positive_order(self):
        return ["5UTR",
                "exon",
                "start_codon",
                "CDS",
                "stop_codon",
                "3UTR"]


class GTF(Parser):
    """The parsing class."""

    __annot_type__ = "gtf"

    def __init__(self, handle):
        """
        Constructor for the parser.
        :param handle: either the filename or the handle for the file to parse.
        :return:
        """

        super().__init__(handle)

    def __next__(self):
        line = next(self._handle)
        try:
            return GtfLine(line)
        except Exception:
            error = "Invalid line for file {}, position {}:\n{}".format(
                self.name, self._handle.tell(), line)
            raise ValueError(error)

    @property
    def file_format(self):
        return self.__annot_type__
