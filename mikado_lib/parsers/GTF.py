# coding: utf_8

"""
Generic parser for GTF files.
"""

from mikado_lib.parsers import Parser
from mikado_lib.parsers.gfannotation import GFAnnotation


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

    def __init__(self, line, my_line='', header=False):

        GFAnnotation.__init__(self, line, my_line, header=header)

    def _parse_attributes(self):
        """
        Method to retrieve the attributes from the last field
        of the GTF line.
        :return:
        """

        info_list = []
        for info in [x for x in self._attr.rstrip().split(';') if x != '']:
            info_list.append(info)
            info = info.lstrip().split(' ')
            try:
                self.attributes[info[0]] = info[1].replace('"', '')
            except IndexError:
                raise IndexError(info_list, info)

        if 'exon_number' in self.attributes:
            self.attributes['exon_number'] = int(self.attributes['exon_number'])
        assert 'gene_id', 'transcript_id' in self.attributes

        if 'nearest_ref' in self.attributes:
            self.nearest_ref = self.attributes['nearest_ref']
        if 'tss_id' in self.attributes:
            self.tss_id = self.attributes['tss_id']

        for tag in filter(lambda att: att not in
                          ('gene_id', 'transcript_id', 'nearest_ref', 'tss_id', 'class_code'),
                          self.attributes.keys()):
            self.__dict__[tag.lower()] = self.attributes[tag]

    def _format_attributes(self):

        """
        Private method to format the last field of the GTF line
        prior to printing.
        :return:
        """

        info_list = []
        assert 'gene_id', 'transcript_id' in self.attributes
        if isinstance(self.gene, list):
            gene = ",".join(self.gene)
        else:
            gene = self.gene
        self.attributes['gene_id'] = gene
        if self.attributes['transcript_id'] is None:
            self.attributes['transcript_id'] = self.transcript

        order = ['gene_id', 'transcript_id', 'exon_number', 'gene_name', 'transcript_name']

        for tag in order:
            if tag in self.attributes:
                if isinstance(self.attributes[tag], list):
                    val = ",".join(self.attributes[tag])
                else:
                    val = self.attributes[tag]
                info_list.append("{0} \"{1}\"".format(tag, val))

        for info in filter(lambda x: x not in order, self.attributes.keys()):
            if info == "Parent" and \
                            self.attributes[info] in (self.gene,
                                                      self.transcript,
                                                      self.parent):
                continue
            if info == "ID" and self.attributes[info] in (self.gene, self.transcript):
                continue

            if isinstance(self.attributes[info], list):
                val = ",".join(self.attributes[info])
            else:
                val = self.attributes[info]
            info_list.append("{0} \"{1}\"".format(info, val))
        attributes = ";".join(info_list) + ";"
        return attributes

    @property
    def name(self):
        """
        Alias for id.
        :rtype: str
        """
        return self.id

    @name.setter
    def name(self, *args):
        """
        Setter for name. The argument must be a string.
        :param args:
        :type args: list[(str)] | str
        """
        if isinstance(args[0], str):
            raise TypeError("Invalid value for name: {0}".format(args[0]))
        self.id = args[0]

    @property
    def is_transcript(self):
        """
        Flag. True if feature is "transcript" or contains "RNA", False in all other cases.
        :rtype : bool
        """
        if self.feature is None:
            return False
        if "transcript" == self.feature or "RNA" in self.feature:
            return True
        return False

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

        if self.is_transcript is True:
            return [self.gene]
        else:
            return [self.transcript]

    @parent.setter
    def parent(self, parent):
        """
        Setter for parent. Acceptable values
        :param parent: the new parent value
        :type parent: list | str
        """
        if isinstance(parent, str):
            parent = parent.split(",")
        self.attributes["Parent"] = parent

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

    @property
    def gene(self):
        """
        Return the "gene_id" field.
        :rtype : str | None
        """
        return self.attributes["gene_id"]

    @gene.setter
    def gene(self, gene):
        """
        Setter for gene.
        :param gene:
        :rtype : str

        """

        self.attributes["gene_id"] = self.__gene = gene
        if self.is_transcript:
            self.attributes["Parent"] = gene

    @property
    def transcript(self):
        """
        This property returns the "transcript_id" field of the GTF line.
        :rtype : str
        """
        return self.attributes["transcript_id"]

    @transcript.setter
    def transcript(self, transcript):
        """
        Setter for the transcript attribute. It also modifies the "transcript_id" field.
        :param transcript:
        :type transcript: str
        """

        self.attributes["transcript_id"] = self._transcript = transcript
        if self.is_transcript is True:
            self.attributes["ID"] = transcript
        else:
            self.attributes["Parent"] = [transcript]

    @property
    def is_parent(self):
        """
        True if we are looking at a transcript, False otherwise.
        :rtype : bool
        """
        if self.is_transcript is True:
            return True
        return False

    @property
    def is_derived(self):
        """
        Property. It checks whether there is a "Derives_from" attribute among the line attributes.
        :rtype bool
        """
        return "derives_from" in [x.lower() for x in self.attributes]

    @property
    def derived_from(self):
        """
        Boolean property. True if the GTF line has a "derives_from" tag,
        False otherwise.
        """
        if self.is_derived is False:
            return None
        else:
            key = list(filter(lambda x: x.lower() == "derives_from", self.attributes.keys()))[0]
            return self.attributes[key].split(",")

    @property
    def is_gene(self):
        """
        In a GTF this should always evaluate to False
        """
        return self.feature == "gene"


class GTF(Parser):
    """The parsing class."""

    def __init__(self, handle):
        super().__init__(handle)

    def __next__(self):
        line = self._handle.readline()
        if line == '':
            raise StopIteration
        return GtfLine(line)
