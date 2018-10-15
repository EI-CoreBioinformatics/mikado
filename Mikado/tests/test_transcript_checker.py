import gzip
import os
import tempfile
import unittest
import pkg_resources
import pyfaidx
from Mikado.transcripts.transcriptchecker import TranscriptChecker
from Mikado.parsers.GFF import GffLine
from Mikado.parsers.GTF import GtfLine
from Mikado.transcripts.transcript import Transcript
from Mikado.exceptions import InvalidTranscript
from Mikado.utilities.log_utils import create_default_logger
import Bio.Seq


class TChekerTester(unittest.TestCase):

    temp_genome = None

    @classmethod
    def setUpClass(cls):

        # Prepare the genome
        cls.temp_genome = tempfile.NamedTemporaryFile(mode="wb", suffix=".fa")
        with pkg_resources.resource_stream("Mikado.tests", "chr5.fas.gz") as comp:
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
            cls.fasta.close()
            os.remove("{}.fai".format(cls.temp_genome.name))

    def test_translation_table(self):

        self.assertEqual(TranscriptChecker.get_translation_table(),
                         {65: 84, 67: 71, 71: 67, 84: 65})

    def test_rev_complement(self):

        string = "AGTCGTGCAGNGTCGAAGTGCAACAGTGC"

        self.assertEqual(TranscriptChecker.rev_complement(string),
                         "GCACTGTTGCACTTCGACNCTGCACGACT")

    def test_init(self):

        with self.assertRaises(ValueError):
            tcheck = TranscriptChecker(self.model, None)

        for wrong_splices in ["AGGT", None, 100]:
            with self.assertRaises(ValueError):
                tcheck = TranscriptChecker(self.model, self.model_fasta, canonical_splices=wrong_splices)

        tcheck = TranscriptChecker(self.model, self.model_fasta)
        self.assertEqual(tcheck.cdna_length, 1718)
        self.assertEqual(sorted(tcheck.exons), sorted([(exon.start, exon.end) for exon in self.exons]))
        self.assertEqual(tcheck.fasta_seq, self.model_fasta)
        transl_table = str.maketrans("ACGT", "TGCA")
        self.assertEqual(transl_table, tcheck.translation_table)

        with self.subTest(initializer=Bio.Seq.Seq):
            _ = TranscriptChecker(self.model, Bio.Seq.Seq(str(self.model_fasta)))

        with self.subTest(initializer=str):
            _ = TranscriptChecker(self.model, str(self.model_fasta))

        with self.subTest(initializer=pyfaidx.Sequence):
            _ = TranscriptChecker(self.model, pyfaidx.Sequence(seq=str(self.model_fasta), name=tcheck.id))

        # Now check initializing with a GFF/GTF line
        for out_format in ["gtf", "gff3"]:
            with self.subTest(out_format=out_format):
                line = self.model.format(out_format).split("\n")[0]
                try:
                    tcheck = TranscriptChecker(line, self.model_fasta)
                except ValueError as exc:
                    raise ValueError(line)
                # tcheck.finalize()
                # self.assertGreater(tcheck.cdna_length, 0)

    def test_check_reverse_strand(self):
        
        self.model.strand = "-"
        tcheck = TranscriptChecker(self.model, self.model_fasta, strand_specific=False)
        tcheck.check_strand()
        self.assertEqual(tcheck.strand, "+")

    def test_check_strand_not_reversed(self):
        self.model.strand = "-"
        tcheck = TranscriptChecker(self.model, self.model_fasta, strand_specific=True)
        tcheck.check_strand()
        self.assertEqual(tcheck.strand, "+")
        self.assertFalse(tcheck.attributes.get("canonical_on_reverse_strand", False))
        self.assertFalse(tcheck.suspicious_splicing)

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

    def test_monoexonic_cds(self):

        # Chr5	tair10	exon	26584797	26584879	.	+	.	ID=c58_g1_i3.mrna1.19.exon1;Parent=c58_g1_i3.mrna1.19
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                exon = self.gff_lines[1]
                transcript_line = self.gff_lines[0]
                transcript_line.end = exon.end
                transcript_line.strand = strand
                exon.strand = strand
                model = Transcript(transcript_line)
                model.add_exon(exon)
                model.add_exon((exon.start + 2, exon.end), feature="CDS")
                model.finalize()
                self.assertTrue(model.is_coding)
                fasta = self.fasta[model.chrom][model.start - 1: model.end]
                tcheck = TranscriptChecker(model.copy(), fasta, force_keep_cds=True, strand_specific=False)
                tcheck.check_strand()
                self.assertEqual(model.strand, strand)
                self.assertTrue(model.is_coding)

    def test_negative(self):

        gtf_lines = """Chr5	Cufflinks	transcript	26575364	26578163	1000	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";exon_number "1";FPKM "2.9700103727";conf_hi "3.260618";frac "0.732092";cov "81.895309";conf_lo "2.679403";
Chr5	Cufflinks	exon	26575364	26575410	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26575495	26575620	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26575711	26575797	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26575885	26575944	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26576035	26576134	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26576261	26577069	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26577163	26577288	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26577378	26577449	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";
Chr5	Cufflinks	exon	26577856	26578163	.	-	.	gene_id "cufflinks_star_at.23553";transcript_id "cufflinks_star_at.23553.1";"""

        gtf_lines = [GtfLine(line) for line in gtf_lines.split("\n")]

        self.assertEqual(len([_ for _ in gtf_lines if _.header]), 0)

        transcript = Transcript(gtf_lines[0])
        transcript.add_exons(gtf_lines[1:])
        transcript.finalize()
        fasta_seq = self.fasta[transcript.chrom][transcript.start - 1:transcript.end]

        tr_neg = transcript.copy()
        tchecker = TranscriptChecker(tr_neg, fasta_seq, strand_specific=False)
        self.assertEqual(tchecker.strand, "-")
        self.assertEqual(tchecker.fasta_seq, fasta_seq)
        tchecker.check_strand()
        self.assertEqual(tchecker.strand, "-")

        tr_neg = transcript.copy()
        tr_neg.strand = "+"
        for ss in (False, True):
            with self.subTest(ss=ss):
                tchecker = TranscriptChecker(tr_neg.copy(), fasta_seq, strand_specific=ss)
                tchecker.check_strand()
                self.assertEqual(tchecker.strand, "-")

    def test_suspicious(self):

        self.model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertTrue(self.model.suspicious_splicing)
        del self.model.attributes["mixed_splices"]
        self.assertFalse(self.model.suspicious_splicing)

        self.model.attributes["canonical_number"] = 0
        self.assertFalse(self.model.suspicious_splicing)

        del self.model.attributes["canonical_number"]

        self.model.attributes["canonical_on_reverse_strand"] = True
        self.assertTrue(self.model.suspicious_splicing)
        self.model.attributes["canonical_on_reverse_strand"] = False
        self.assertFalse(self.model.suspicious_splicing)
        self.model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertTrue(self.model.suspicious_splicing)

        del self.model.attributes["mixed_splices"]
        del self.model.attributes["canonical_on_reverse_strand"]
        self.model.attributes["canonical_number"] = 0
        self.assertFalse(self.model.suspicious_splicing)
        self.assertTrue(self.model.only_non_canonical_splicing)
        self.model.attributes["canonical_on_reverse_strand"] = True
        self.assertTrue(self.model.suspicious_splicing)
        self.assertTrue(self.model.only_non_canonical_splicing)
        del self.model.attributes["canonical_on_reverse_strand"]
        self.model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertTrue(self.model.suspicious_splicing)
        self.assertTrue(self.model.only_non_canonical_splicing)

    def test_reverse_with_cds_negative(self):

        model = Transcript()
        model.chrom, model.start, model.end = "Chr5", 1251, 5043
        model.id, model.parent, model.strand = "AT5G01010.1", "AT5G01010", "-"
        model.add_exons([(4765, 5043), (4552, 4679), (4335, 4467), (4102, 4258), (3927, 4005), (3762, 3802),
                         (3543, 3659), (3303, 3383), (2872, 2934), (2748, 2799), (2435, 2509), (1914, 1961),
                         (1745, 1780), (1572, 1646), (1251, 1459)])
        model.add_exons([(4765, 4924), (4552, 4679), (4335, 4467), (4102, 4258), (3927, 4005), (3762, 3802),
                         (3543, 3659), (3303, 3383), (2872, 2934), (2748, 2799), (2435, 2509), (1914, 1961),
                         (1745, 1780), (1572, 1646), (1388, 1459)],
                        features="CDS")
        model.finalize()
        model_fasta = self.fasta["Chr5"][model.start - 1:model.end]
        check_model = TranscriptChecker(model, model_fasta)
        check_model.check_strand()
        self.assertEqual(check_model.strand, "-")
        self.assertGreater(check_model.combined_cds_length, 0)
        model.unfinalize()
        model.strand = "+"
        check_model = TranscriptChecker(model, model_fasta)
        check_model.check_strand()
        self.assertEqual(check_model.strand, "-")
        self.assertFalse(check_model.is_coding)
        check_model = TranscriptChecker(model, model_fasta, force_keep_cds=True)
        # Check that if we want to keep the CDS, this will raise an error
        with self.assertRaises(InvalidTranscript):
            check_model.check_strand()

    def test_reverse_with_cds_positive(self):

        model = Transcript()
        model.chrom, model.start, model.end = "Chr5", 9930, 13235
        model.id, model.parent, model.strand = "AT5G01030.1", "AT5G01030", "+"
        model.add_exons([(9930, 10172), (10620, 12665), (12797, 13235)])
        model.add_exons([(10638, 12665), (12797, 13003)], features="CDS")
        model.finalize()
        model_fasta = self.fasta["Chr5"][model.start - 1:model.end]
        check_model = TranscriptChecker(model, model_fasta)
        check_model.check_strand()
        self.assertEqual(check_model.strand, "+")
        self.assertGreater(check_model.combined_cds_length, 0)
        model.unfinalize()
        model.strand = "-"
        check_model = TranscriptChecker(model, model_fasta)
        check_model.check_strand()
        self.assertEqual(check_model.strand, "+")
        self.assertFalse(check_model.is_coding)
        check_model = TranscriptChecker(model, model_fasta, force_keep_cds=True)
        with self.assertRaises(InvalidTranscript):
            check_model.check_strand()

    def test_monoexonic_suspicious(self):

        """A monoexonic transcript should never appear as a suspicious transcript, in terms of splicing."""

        exon = self.gff_lines[1]
        transcript_line = self.gff_lines[0]
        transcript_line.end = exon.end
        model = Transcript(transcript_line)
        model.add_exon(exon)
        model.finalize()

        model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertFalse(model.suspicious_splicing)
        del model.attributes["mixed_splices"]
        self.assertFalse(model.suspicious_splicing)
        
        model.attributes["canonical_number"] = 0
        self.assertFalse(model.suspicious_splicing)
        
        del model.attributes["canonical_number"]
        
        model.attributes["canonical_on_reverse_strand"] = True
        self.assertFalse(model.suspicious_splicing)
        model.attributes["canonical_on_reverse_strand"] = False
        self.assertFalse(model.suspicious_splicing)
        model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertFalse(model.suspicious_splicing)
        
        del model.attributes["mixed_splices"]
        del model.attributes["canonical_on_reverse_strand"]
        model.attributes["canonical_number"] = 0
        self.assertFalse(model.suspicious_splicing)
        self.assertFalse(model.only_non_canonical_splicing)
        model.attributes["canonical_on_reverse_strand"] = True
        self.assertFalse(model.suspicious_splicing)
        self.assertFalse(model.only_non_canonical_splicing)
        del model.attributes["canonical_on_reverse_strand"]
        model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertFalse(model.suspicious_splicing)
        self.assertFalse(model.only_non_canonical_splicing)

    def test_sequence_reversed(self):

        model = Transcript()
        model.chrom, model.start, model.end, model.strand = "Chr5", 1001, 1500, "+"
        model.add_exon((1001, 1500))
        model.id, model.parent = "foo.1", "foo"
        model.finalize()
        seq = str(self.fasta["Chr5"][1001-1:1500].seq)
        self.assertEqual(len(seq), len(model))
        model = TranscriptChecker(model, seq, strand_specific=True)
        model.reverse_strand()
        fasta = "".join(model.fasta.split("\n")[1:])
        self.assertEqual(model.strand, "-")
        self.assertEqual(fasta, TranscriptChecker.rev_complement(seq))


class StopCodonChecker(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.__genomefile__ = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".fa", prefix="prepare")

        with pkg_resources.resource_stream("Mikado.tests", "chr5.fas.gz") as _:
            cls.__genomefile__.write(gzip.decompress(_.read()))
        cls.__genomefile__.flush()
        cls.genome = pyfaidx.Fasta(cls.__genomefile__.name)
        cls.maxDiff = None

    @classmethod
    def tearDownClass(cls):
        cls.__genomefile__.close()
        cls.genome.close()
        os.remove(cls.__genomefile__.name)
        if os.path.exists(cls.__genomefile__.name + ".fai"):
            os.remove(cls.__genomefile__.name + ".fai")

    def test_positive_strand(self):

        gtf_lines = """chr1A	Self_CESAR/windows_chr1A.gp	transcript	265021906	265026255	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265021906	265021971	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	CDS	265021906	265021971	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265022056	265026255	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";
chr1A	Self_CESAR/windows_chr1A.gp	CDS	265022056	265026252	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";
chr1A	Self_CESAR/windows_chr1A.gp	start_codon	265021906	265021908	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	stop_codon	265026253	265026255	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";"""

        gtf_lines = [GtfLine(_) for _ in gtf_lines.split("\n")]
        t = Transcript(gtf_lines[0])
        t.add_exons(gtf_lines[1:])
        t.finalize()
        self.assertEqual(t.start, t.combined_cds_start)
        self.assertEqual(t.end, t.combined_cds_end)

    def test_positive_strand_split(self):

        gtf_lines="""chr1A	Self_CESAR/windows_chr1A.gp	transcript	265021906	265026355	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265021906	265021971	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	CDS	265021906	265021971	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265022056	265026253	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";
chr1A	Self_CESAR/windows_chr1A.gp	CDS	265022056	265026252	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265026354	265026355	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	start_codon	265021906	265021908	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	stop_codon	265026253	265026253	.	+	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";
chr1A	Self_CESAR/windows_chr1A.gp	stop_codon	265026354	265026355	.	+	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";"""

        gtf_lines = [GtfLine(_) for _ in gtf_lines.split("\n")]
        t = Transcript(gtf_lines[0])
        t.add_exons(gtf_lines[1:])
        t.finalize()
        self.assertEqual(t.start, t.combined_cds_start)
        self.assertEqual(t.end, t.combined_cds_end)

    def test_negative_strand(self):

        gtf_lines = """chr1A	Self_CESAR/windows_chr1A.gp	transcript	265021806	265026255	.	-	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265021806	265021807	.	-	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265021908	265021971	.	-	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	CDS	265021909	265021971	.	-	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	exon	265022056	265026255	.	-	.	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";
chr1A	Self_CESAR/windows_chr1A.gp	CDS	265022056	265026255	.	-	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";
chr1A	Self_CESAR/windows_chr1A.gp	stop_codon	265021806	265021807	.	-	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	stop_codon	265021908	265021908	.	-	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "1"; exon_id "TraesCS1A01G152900.1.1";
chr1A	Self_CESAR/windows_chr1A.gp	start_codon	265026253	265026255	.	-	0	gene_id "TraesCS1A01G152900.1"; transcript_id "TraesCS1A01G152900.1"; exon_number "2"; exon_id "TraesCS1A01G152900.1.2";"""

        gtf_lines = [GtfLine(_) for _ in gtf_lines.split("\n")]
        t = Transcript(gtf_lines[0])
        t.add_exons(gtf_lines[1:])
        t.finalize()
        self.assertEqual(t.start, t.combined_cds_end)
        self.assertEqual(t.end, t.combined_cds_start)

    def test_codon_finder_negative(self):
        gtf_lines = """Chr5	TAIR10	mRNA	5335	5769	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	CDS	5697	5766	.	-	0	transcript_id "AT5G01015.1"; gene_id "AT5G01015";;
Chr5	TAIR10	exon	5697	5769	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	CDS	5335	5576	.	-	1	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	exon	5335	5576	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";"""

        gtf_lines = [GtfLine(_) for _ in gtf_lines.split("\n")]
        t = Transcript(gtf_lines[0])
        t.add_exons(gtf_lines[1:])
        t.finalize()

        seq = self.genome[t.chrom][t.start - 1:t.end]
        logger = create_default_logger("test_codon_finder_negative", level="WARNING")
        tc = TranscriptChecker(t, seq, logger=logger)
        tc.finalize()
        tc.check_orf()
        self.assertTrue(tc.is_coding)
        self.assertIn("has_stop_codon", tc.attributes)
        self.assertIn("has_start_codon", tc.attributes)
        self.assertTrue(tc.has_stop_codon)
        self.assertFalse(tc.has_start_codon)

    def test_codon_finder_negative_2(self):
        gtf_lines = """Chr5	TAIR10	mRNA	5335	5769	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	CDS	5697	5766	.	-	0	transcript_id "AT5G01015.1"; gene_id "AT5G01015";;
Chr5	TAIR10	exon	5697	5769	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	CDS	5338	5576	.	-	1	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	exon	5335	5576	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";"""

        gtf_lines = [GtfLine(_) for _ in gtf_lines.split("\n")]
        t = Transcript(gtf_lines[0])
        t.add_exons(gtf_lines[1:])
        t.finalize()

        seq = self.genome[t.chrom][t.start - 1:t.end]
        logger = create_default_logger("test_codon_finder_negative_2", level="WARNING")
        self.assertTrue(t.has_start_codon)
        self.assertTrue(t.has_stop_codon)
        tc = TranscriptChecker(t, seq, logger=logger)
        tc.finalize()
        tc.check_orf()
        self.assertTrue(tc.is_coding)
        self.assertIn("has_stop_codon", tc.attributes)
        self.assertIn("has_start_codon", tc.attributes)
        self.assertFalse(tc.has_stop_codon)
        self.assertFalse(tc.has_start_codon)

    def test_codon_finder_negative_3(self):

        gtf_lines = """Chr5	TAIR10	mRNA	5335	5769	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	CDS	5697	5769	.	-	0	transcript_id "AT5G01015.1"; gene_id "AT5G01015";;
Chr5	TAIR10	exon	5697	5769	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	CDS	5335	5576	.	-	1	transcript_id "AT5G01015.1"; gene_id "AT5G01015";
Chr5	TAIR10	exon	5335	5576	.	-	.	transcript_id "AT5G01015.1"; gene_id "AT5G01015";"""

        gtf_lines = [GtfLine(_) for _ in gtf_lines.split("\n")]
        t = Transcript(gtf_lines[0])
        t.add_exons(gtf_lines[1:])
        t.finalize()

        seq = self.genome[t.chrom][t.start - 1:t.end]
        correct_seq = "".join("""ATGGAGTCTAGCTTGCATAGTGTGATTTTCTTAGGTTTGCTTGCGACGATTCTGGTTACG
ACCAATGGCCAAGGAGACGGGACGGGGCTAAATGCAGAAGAAATGTGGCCAGTGGAGGTG
GGGATGGAGTATAGAGTATGGAGGAGAAAGCTGATGACGCCATTGGAGCTGTGCTTGGAG
TGCAAATGCTGCTCCTCCACCACTTGTGCCACCATGCCTTGCTGTTTCGGCATCAATTGC
CAGCTTCCCAACAAGCCATTTGGCGTTTGTGCCTTTGTTCCCAAGTCATGCCATTGTAAT
TCTTGCTCCATTTGA""".split("\n"))
        logger = create_default_logger("test_codon_finder_negative_3", level="WARNING")
        tc = TranscriptChecker(t, seq, logger=logger)
        tc.finalize()
        correct_length = (5576 - 5335 + 1) + (5769 - 5697 + 1)
        self.assertEqual(correct_length, len(correct_seq), (correct_length, len(correct_seq)))
        self.assertEqual(tc.cdna_length, correct_length, (correct_length, tc.cdna_length))
        self.assertEqual(len(tc.cdna), tc.cdna_length)
        self.assertEqual(correct_seq, tc.cdna)

        tc.check_orf()
        tc_orfs = tc.find_overlapping_cds(tc.get_internal_orf_beds())
        self.assertEqual(1, len(tc_orfs))
        self.assertTrue(tc_orfs[0].has_stop_codon, (tc_orfs[0], tc_orfs[0].stop_codon))
        self.assertTrue(tc_orfs[0].has_start_codon, (tc_orfs[0], tc_orfs[0].start_codon))

        self.assertTrue(tc.is_coding)
        self.assertIn("has_stop_codon", tc.attributes)
        self.assertIn("has_start_codon", tc.attributes)
        self.assertTrue(tc.has_start_codon, tc.cdna)
        self.assertTrue(tc.has_stop_codon, tc.cdna)


if __name__ == '__main__':
    unittest.main()
