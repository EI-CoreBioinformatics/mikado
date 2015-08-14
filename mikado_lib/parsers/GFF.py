#!/usr/bin/env python3
# coding: utf_8


"""
Module to serialize GFF files.
"""

from mikado_lib.parsers import Parser


class GffLine(object):
    """Object which serializes a GFF line.
    Parameters:
    :param _line: the original line
    :type _line: str

    :param _fields: the splitted line
    :type _fields: list

    :param chrom: the chromosome
    :type chrom: str

    :param source: source field, where the data originated from.
    :type source: str

    :param feature: mRNA, gene, exon, start/stop_codon, etc.
    :type feature: str

    :param start: start of the feature
    :type start: int

    :param end: stop of the feature
    :type end: int

    :param strand: strand of the feature
    :type strand: str
    :type strand: None

    :param score: Score assigned to the feature
    :type score: None
    :type score: float

    :param phase: Codon phase for protein-coding CDS/exon features
    :type phase: int
    :type phase: None

    :param attributes: a dictionary which contains the extra information.
    :type attributes: dict

    """

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

        self.attributes = dict()
        self.id = None
        self.parent = []
        self.__score = None
        self.__strand = None

        self.attributeOrder = []
        if line is None:  # Empty constructor
            return
        if line == '' and my_line == '':
            return

        if line == '' and my_line != "":
            self._line = my_line
        else:
            self._line = line

        self._fields = line.rstrip().split('\t')
        self.header = header

        if self.header or len(self._fields) != 9 or self._line == '':
            self.feature = None
            self.header = True
            return

        if len(self._fields) != 9:
            self.header = True
            return
        self.chrom, self.source, self.feature = self._fields[0:3]
        self.start, self.end = tuple(int(i) for i in self._fields[3:5])

        if self._fields[5] == '.':
            self.score = None
        else:
            self.score = float(self._fields[5])

        self.strand = self._fields[6]

        if self._fields[7] == '.':
            self.phase = None
        else:
            try:
                self.phase = int(self._fields[7])
                assert self.phase in (0, 1, 2)
            except:
                raise

        self._Attr = self._fields[8]

        self.attributeOrder = []

        for item in [x for x in self._Attr.rstrip().split(';') if x != '']:
            itemized = item.split('=')
            try:
                if itemized[0].lower() == "parent":
                    self.parent = itemized[1]

                elif itemized[0].upper() == "ID":
                    self.id = itemized[1]
                else:
                    self.attributes[itemized[0]] = itemized[1]
                    self.attributeOrder.append(itemized[0])
            except IndexError:
                pass

        _ = self.name  # Set the name

    def __str__(self):
        if not self.feature:
            return self._line.rstrip()

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
        attrs = []
        if self.id is not None:
            attrs.append("ID={0}".format(self.id))
        if len(self.parent) > 0:
            assert type(self.parent) is list, "{0}\n{1}\n{2}".format(self.parent, self.attributes, self._line)
            attrs.append("Parent={0}".format(",".join(self.parent)))
        if not self.attributeOrder:
            self.attributeOrder = sorted(list(filter(lambda x: x not in ["ID", "Parent"], self.attributes.keys())))
        for att in filter(lambda x: x not in ["ID", "Parent"], self.attributeOrder):
            if self.attributes[att] is not None:
                try:
                    attrs.append("{0}={1}".format(att, self.attributes[att]))
                except KeyError:
                    continue  # Hack for those times when we modify the attributes at runtime

        line = '\t'.join(
            [self.chrom, self.source,
             self.feature, str(self.start), str(self.end),
             str(score), strand, phase,
             ";".join(attrs)]
        )
        return line

    def __len__(self):
        if "end" in self.__dict__:
            return self.end - self.start + 1
        else:
            return 0

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

    @property
    def parent(self):
        """This property looks up the "Parent" field in the "attributes" dictionary. Contrary to other attributes,
        this property returns a *list*, not a string. This is due to the fact that GFF files support
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

        if parent is None:
            self.attributes["Parent"] = None
        elif type(parent) is str:
            new_parent = parent.split(",")
            self.attributes["Parent"] = new_parent
        elif type(parent) is list:
            self.attributes["Parent"] = parent
        else:
            raise TypeError(parent, type(parent))
        assert type(self.parent) is list or self.parent is None

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
    def strand(self):

        """
        Strand attribute. One of None, +, -
        :rtype str
        :rtype None
        """

        return self.__strand

    @strand.setter
    def strand(self, strand):
        """
        Strand setter. It verifies that strand is a valid value (None, +, -)
        :param strand: new strand value
        :type strand: str
        :type strand: None
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
        :rtype float
        :rtype None
        """

        return self.__score

    @score.setter
    def score(self, *args):
        """
        Score setter. It verifies that the score is a valid value (None or float).
        :param args: List of arguments. Only the first is used.
        :type args: list
        :type args: None
        :type args: float
        """

        if type(args[0]) in (float, int):
            self.__score = args[0]
        elif args[0] is None or args[0] == '.':
            self.__score = None
        elif type(args[0]) is str:
            self.__score = float(args[0])
        else:
            raise TypeError(args[0])

    @property
    def is_exon(self):
        """
        Property. True if the feature is CDS/exon/UTR/start or stop codon.
        :rtype bool
        """

        if self.feature is None:
            return False
        f = self.feature.lower()
        if f.endswith("cds") or f.endswith("exon") or "utr" in f or "codon" in f:
            return True
        return False

    @property
    def is_transcript(self):
        """
        Property. True if the feature is an RNA, false otherwise.
        :rtype bool
        """

        if self.feature is None:
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

        if self.feature is not None and self.feature.endswith("gene") and self.is_parent is True:
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
        returns the parent; else, it returns None (as only transcript lines have the gene among the attributes)
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

    @property
    def is_derived(self):
        """
        Property. It checks whether there is a "Derives_from" attribute among the line attributes.
        :rtype bool
        """
        return "Derives_from" in self.attributes

    @property
    def derived_from(self):
        if self.is_derived is False:
            return None
        else:
            return self.attributes["Derives_from"].split(",")


class GFF3(Parser):
    """
    Class that is used to parse a GFF file.
    """

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
