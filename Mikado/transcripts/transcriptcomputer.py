from collections import namedtuple
from .transcript import Transcript
from ..exceptions import InvalidCDS


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
    data_tuple = namedtuple("transcript_data", data_fields)

    def __init__(self, *args, **kwargs):
        kwargs.pop("accept_undefined_multi", None)
        kwargs.pop("trust_orf", None)
        super().__init__(accept_undefined_multi=True, trust_orf=True, *args, **kwargs)
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
