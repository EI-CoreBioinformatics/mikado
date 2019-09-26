import gzip
import os
import tempfile
import unittest
import pkg_resources
import pyfaidx
from ..transcripts.transcriptchecker import TranscriptChecker
import pysam
from ..parsers.GTF import GtfLine
from ..transcripts.transcript import Transcript
from ..exceptions import InvalidTranscript
from ..utilities.log_utils import create_default_logger
import Bio.Seq
from copy import copy


class TChekerTester(unittest.TestCase):

    temp_genome = None

    @classmethod
    def setUpClass(cls):

        # Prepare the genome
        cls.fasta = pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz")
        # cls.fasta_temp = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".fa.gz")
        # cls.fasta_temp.write(pkg_resources.resource_stream("Mikado.tests", "chr5.fas.gz").read())
        cls.fasta = pysam.FastaFile(cls.fasta)
        cls._model = Transcript()
        cls._model.chrom = "Chr5"
        cls._model.start, cls._model.end, cls._model.strand, cls._model.id, cls._model.parent = (
            26584797, 26595528, "+", "c58_g1_i3.mrna1.19", "c58_g1_i3.path1.19"
        )
        cls._model.add_exons([(26584797, 26584879),
                              (26585220, 26585273),
                              (26585345, 26585889),
                              (26585982, 26586294),
                              (26586420, 26586524),
                              (26586638, 26586850),
                              (26586934, 26586996),
                              (26587084, 26587202),
                              (26587287, 26587345),
                              (26587427, 26587472),
                              (26595411, 26595528)])
        cls._model.finalize()
        cls._exons = [cls.fasta.fetch(cls._model.chrom, line[0] - 1,  line[1]) for line in sorted(cls._model.exons)]
        # We need the whole genomic fragment
        cls._model_fasta = cls.fasta.fetch(cls._model.chrom, cls._model.start - 1, cls._model.end)

    def setUp(self):

        # Prepare the model

        self.model = self._model.copy()
        self.exons = self._exons.copy()

        self.assertEqual(sum([len(exon) for exon in self.exons]), 1718, self.exons)
        # We need the whole genomic fragment
        self.model_fasta = copy(self._model_fasta)

    @classmethod
    def tearDownClass(cls):
        # Remove the genome
        if hasattr(cls.temp_genome, "close"):
            cls.temp_genome.close()
            cls.fasta.close()
            os.remove("{}.fai".format(cls.temp_genome.name))

    def test_rev_complement(self):

        string = "AGTCGTGCAGNGTCGAAGTGCAACAGTGC"

        self.assertEqual(TranscriptChecker.rev_complement(string),
                         "GCACTGTTGCACTTCGACNCTGCACGACT")

        string = "agtcGTGCAGNGTCGAAGTGCAACAgtgc"

        self.assertEqual(TranscriptChecker.rev_complement(string),
                         "gcacTGTTGCACTTCGACNCTGCACgact")

    def test_init(self):

        with self.assertRaises(ValueError):
            tcheck = TranscriptChecker(self.model, None)

        for wrong_splices in ["AGGT", None, 100]:
            with self.assertRaises(ValueError):
                tcheck = TranscriptChecker(self.model, self.model_fasta, canonical_splices=wrong_splices)

        tcheck = TranscriptChecker(self.model, self.model_fasta)
        tcheck.finalize()
        self.assertEqual(tcheck.cdna_length, 1718)
        self.assertEqual(sorted(tcheck.exons), sorted([(exon[0], exon[1]) for exon in self.model.exons]))
        self.assertEqual(str(tcheck.fasta_seq.seq), self.model_fasta,
                         (type(tcheck.fasta_seq), type(self.model_fasta),
                          len(tcheck.fasta_seq), len(self.model_fasta)))

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
        logger = create_default_logger("test_check_reverse_strand", "DEBUG")
        tcheck.logger = logger
        self.assertFalse(tcheck.strand_specific)
        self.assertFalse(tcheck.is_reference)
        self.assertFalse(tcheck.lenient)
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

        self.model.unfinalize()
        for exon in sorted(self.model.exons)[1:]:
            self.model.remove_exon(exon)
        
        self.model.finalize()
        fasta = self.fasta[self.model.chrom][self.model.start - 1: self.model.end]
        tcheck = TranscriptChecker(self.model.copy(), fasta, strand_specific=False)
        tcheck.check_strand()
        self.assertIsNone(tcheck.strand)

        tcheck = TranscriptChecker(self.model.copy(), fasta, strand_specific=True)
        tcheck.check_strand()
        self.assertEqual(tcheck.strand, "+")

        neg = self.model.copy()
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
                model = self.model.copy()
                model.unfinalize()
                [model.remove_exon(exon) for exon in sorted(model.exons[1:])]
                exon = model.exons[0]
                model.add_exon((exon[0] + 2, exon[1]), feature="CDS")
                model.strand = strand
                model.finalize()
                self.assertTrue(model.is_coding)
                fasta = self.fasta[model.chrom][model.start - 1: model.end]
                tcheck = TranscriptChecker(model.copy(), fasta, force_keep_cds=True, strand_specific=False)
                tcheck.check_strand()
                self.assertTrue(tcheck.is_coding)
                self.assertEqual(tcheck.strand, strand)

    def test_negative(self):

        transcript = Transcript()
        transcript.chrom, transcript.start, transcript.end, transcript.strand = "Chr5", 26575364, 26578163, "-"
        transcript.id, transcript.parent = "cufflinks_star_at.23553.1", "cufflinks_star_at.23553"
        transcript.add_exons([(26575364, 26575410),
                              (26575495, 26575620),
                              (26575711, 26575797),
                              (26575885, 26575944),
                              (26576035, 26576134),
                              (26576261, 26577069),
                              (26577163, 26577288),
                              (26577378, 26577449),
                              (26577856, 26578163)])

        transcript.finalize()
        fasta_seq = self.fasta.fetch(transcript.chrom, transcript.start - 1, transcript.end)
        tr_neg = transcript.copy()
        tchecker = TranscriptChecker(tr_neg, fasta_seq, strand_specific=False)
        self.assertEqual(tchecker.strand, "-")
        self.assertEqual(str(tchecker.fasta_seq.seq), fasta_seq)
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
    
        self.model.unfinalize()
        [self.model.remove_exon(exon) for exon in sorted(self.model.exons)[1:]]
        self.model.finalize()

        self.model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertFalse(self.model.suspicious_splicing)
        del self.model.attributes["mixed_splices"]
        self.assertFalse(self.model.suspicious_splicing)
        
        self.model.attributes["canonical_number"] = 0
        self.assertFalse(self.model.suspicious_splicing)
        
        del self.model.attributes["canonical_number"]
        
        self.model.attributes["canonical_on_reverse_strand"] = True
        self.assertFalse(self.model.suspicious_splicing)
        self.model.attributes["canonical_on_reverse_strand"] = False
        self.assertFalse(self.model.suspicious_splicing)
        self.model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertFalse(self.model.suspicious_splicing)
        
        del self.model.attributes["mixed_splices"]
        del self.model.attributes["canonical_on_reverse_strand"]
        self.model.attributes["canonical_number"] = 0
        self.assertFalse(self.model.suspicious_splicing)
        self.assertFalse(self.model.only_non_canonical_splicing)
        self.model.attributes["canonical_on_reverse_strand"] = True
        self.assertFalse(self.model.suspicious_splicing)
        self.assertFalse(self.model.only_non_canonical_splicing)
        del self.model.attributes["canonical_on_reverse_strand"]
        self.model.attributes["mixed_splices"] = "6positive,1negative"
        self.assertFalse(self.model.suspicious_splicing)
        self.assertFalse(self.model.only_non_canonical_splicing)

    def test_sequence_reversed(self):

        model = Transcript()
        model.chrom, model.start, model.end, model.strand = "Chr5", 1001, 1500, "+"
        model.add_exon((1001, 1500))
        model.id, model.parent = "foo.1", "foo"
        model.finalize()
        seq = self.fasta.fetch("Chr5", 1001 - 1, 1500)
        self.assertEqual(len(seq), len(model))
        model = TranscriptChecker(model, seq, strand_specific=True)
        model.reverse_strand()
        fasta = "".join(model.fasta.split("\n")[1:])
        self.assertEqual(model.strand, "-")
        self.assertEqual(fasta, TranscriptChecker.rev_complement(seq))

    def test_reference_not_flipped(self):

        model = Transcript()
        model.chrom, model.start, model.end = "Chr5", 9930, 13235
        model.id, model.parent, model.strand = "AT5G01030.1", "AT5G01030", "-"
        model.add_exons([(9930, 10172), (10620, 12665), (12797, 13235)])
        model.add_exons([(10638, 12665), (12797, 13003)], features="CDS")
        model.finalize()
        model_fasta = self.fasta["Chr5"][model.start - 1:model.end]
        check_model = TranscriptChecker(model, model_fasta, is_reference=True)
        check_model.check_strand()
        self.assertEqual(check_model.strand, "-")
        self.assertGreater(check_model.combined_cds_length, 0)


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
