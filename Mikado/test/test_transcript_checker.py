from Mikado.loci.transcriptchecker import TranscriptChecker
from Mikado.loci.transcript import Transcript
from Mikado.parsers.GFF import GffLine
import unittest
import pyfaidx
import pkg_resources
import tempfile
import gzip
import os

# TODO: write this test!


class TChekerTester(unittest.TestCase):

    temp_genome = None

    @classmethod
    def setUpClass(cls):

        # Prepare the genome
        cls.temp_genome = tempfile.NamedTemporaryFile(mode="wb", suffix=".fa")
        with pkg_resources.resource_stream("Mikado.test", "chr5.fas.gz") as comp:
            cls.temp_genome.write(gzip.decompress(comp.read()))
        cls.temp_genome.flush()
        cls.fasta = pyfaidx.Fasta(cls.temp_genome.name)


    def setUp(self):

        # Prepare the model
        self.model_lines= """Chr5	tair10	transcript	26584797	26595528	100	+	.	ID=c58_g1_i3.mrna1.19;Parent=c58_g1_i3.path1.19;Name=c58_g1_i3.mrna1.19;gene_name=c58_g1_i3
    Chr5	tair10	exon	26584797	26584879	.	+	.	ID=c58_g1_i3.mrna1.19.exon1;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26585220	26585273	.	+	.	ID=c58_g1_i3.mrna1.19.exon2;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26585345	26585889	.	+	.	ID=c58_g1_i3.mrna1.19.exon3;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26585982	26586294	.	+	.	ID=c58_g1_i3.mrna1.19.exon4;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26586420	26586524	.	+	.	ID=c58_g1_i3.mrna1.19.exon5;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26586638	26586850	.	+	.	ID=c58_g1_i3.mrna1.19.exon6;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26586934	26586996	.	+	.	ID=c58_g1_i3.mrna1.19.exon7;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26587084	26587202	.	+	.	ID=c58_g1_i3.mrna1.19.exon8;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26587287	26587345	.	+	.	ID=c58_g1_i3.mrna1.19.exon9;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26587427	26587472	.	+	.	ID=c58_g1_i3.mrna1.19.exon10;Parent=c58_g1_i3.mrna1.19
    Chr5	tair10	exon	26595411	26595528	.	+	.	ID=c58_g1_i3.mrna1.19.exon11;Parent=c58_g1_i3.mrna1.19"""

        self.gff_lines = []
        for line in self.model_lines.split("\n"):
            line = line.rstrip().lstrip()
            line = GffLine(line)
            self.gff_lines.append(line)
    
        self.model = Transcript(self.gff_lines[0])
        self.model.add_exons(self.gff_lines[1:])
        self.model.finalize()
    
        self.exons = [self.fasta[line.chrom][line.start - 1:line.end] for line in self.gff_lines[1:]]

        self.assertEqual(sum([len(exon) for exon in self.exons]), 1718, self.exons)
        # We need the whole genomic fragment
        self.model_fasta = self.fasta["Chr5"][self.model.start -1:self.model.end]
        self.assertEqual(self.gff_lines[1].start, 26584797)
        self.assertEqual(self.gff_lines[1].end, 26584879)
        self.assertEqual(self.model.exons[0][0], self.gff_lines[1].start)
        self.assertEqual(self.model.exons[0][1], self.gff_lines[1].end)

    @classmethod
    def tearDownClass(cls):
        # Remove the genome
        if hasattr(cls.temp_genome, "close"):
            cls.temp_genome.close()
            os.remove("{}.fai".format(cls.temp_genome.name))

    def test_translation_table(self):

        self.assertEqual(TranscriptChecker.get_translation_table(),
                         {65: 84, 67: 71, 71: 67, 84: 65})

    def test_rev_complement(self):

        string = "AGTCGTGCAGNGTCGAAGTGCAACAGTGC"

        self.assertEqual(TranscriptChecker.rev_complement(string),
                         "GCACTGTTGCACTTCGACNCTGCACGACT")

    def test_init(self):

        tcheck = TranscriptChecker(self.model, self.model_fasta)
        self.assertEqual(tcheck.cdna_length, 1718)
        self.assertEqual(sorted(tcheck.exons), sorted([(exon.start, exon.end) for exon in self.exons]))
        self.assertEqual(tcheck.fasta_seq, self.model_fasta)

    def test_check_reverse_strand(self):
        
        self.model.strand = "-"
        tcheck = TranscriptChecker(self.model, self.model_fasta, strand_specific=False)
        tcheck.check_strand()
        self.assertEqual(tcheck.strand, "+")

    def test_check_strand_not_reversed(self):
        self.model.strand = "-"
        tcheck = TranscriptChecker(self.model, self.model_fasta, strand_specific=True)
        tcheck.check_strand()
        self.assertEqual(tcheck.strand, "-")
        self.assertTrue(tcheck.attributes["canonical_on_reverse_strand"])

    def test_monoexonic(self):

        exon = self.gff_lines[1]
        transcript_line = self.gff_lines[0]
        transcript_line.end = exon.end
        model = Transcript(transcript_line)
        model.add_exon(exon)
        model.finalize()
        fasta = self.fasta[model.chrom][model.start - 1: model.end]

        tcheck = TranscriptChecker(model.copy(), fasta, strand_specific=False)
        tcheck.check_strand()
        self.assertIsNone(tcheck.strand)

        tcheck = TranscriptChecker(model.copy(), fasta, strand_specific=True)
        tcheck.check_strand()
        self.assertEqual(tcheck.strand, "+")

        neg = model.copy()
        neg.strand = "-"

        tcheck = TranscriptChecker(neg.copy(), fasta, strand_specific=False)
        tcheck.check_strand()
        self.assertIsNone(tcheck.strand)

        tcheck = TranscriptChecker(neg.copy(), fasta, strand_specific=True)
        tcheck.check_strand()
        self.assertEqual(tcheck.strand, "-")


if __name__ == '__main__':
    unittest.main()
