from .. import utilities, exceptions
from ..parsers.GFF import GffLine
from ..loci import Transcript
from ..loci.locus import Locus, expand_transcript
from ..parsers.bed12 import BED12
from ..subprograms.util.trim import trim_coding, trim_noncoding
import unittest
import pysam
import pkg_resources
from ..utilities.log_utils import create_null_logger, create_default_logger
from ..configuration.configurator import to_json
from pytest import mark


__author__ = 'Luca Venturini'


class TestTrimming(unittest.TestCase):

    def test_wrong_cds(self):

        transcript = Transcript()
        transcript.chrom = "15"
        transcript.source = "protein_coding"
        transcript.start = 47631264
        transcript.end = 48051999

        exons = [(47631264, 47631416),
                 (47704590, 47704669),
                 (47762671, 47762742),
                 (47893062, 47893093),
                 (47895572, 47895655),
                 (48051942, 48051999)]

        transcript.strand = "+"
        transcript.add_exons(exons)
        transcript.id = "ENST00000560636"
        transcript.parent = "ENSG00000137872"
        cds_line = "\t".join(["15",
                              "protein_coding",
                              "CDS",
                              "48051996",
                              "48051996",
                              ".",
                              "+",
                              "0",
                              "ID=ENST00000560636.cds1;Parent=ENST00000560636"])
        cds_line = GffLine(cds_line)
        transcript.add_exon(cds_line)
        logger = utilities.log_utils.create_null_logger()
        transcript.logger = logger
        with self.assertLogs("null", level="WARNING"):
            transcript.finalize()

        trimmed = trim_coding(transcript, logger, max_length=50)
        self.assertEqual(trimmed.start, 47631366)
        self.assertEqual(trimmed.end, 48051992)

    # @unittest.SkipTest
    def test_correct_cds(self):

        transcript = Transcript()
        transcript.chrom = "Chr1"
        transcript.source = "test"
        transcript.start = 10000
        transcript.end = 20000

        exons = [(10000, 11500),
                 (12000, 13000),
                 (15000, 18000),
                 (19000, 20000)]

        cds = [(11400, 11500),  # 101
               (12000, 13000),  # 1001 ==> 1102
               (15000, 17998)]  # 2998 == > 3090 (y)

        transcript.add_exons(exons)
        transcript.add_exons(cds, features="CDS")

        transcript.strand = "+"
        transcript.finalize()

        logger = utilities.log_utils.create_null_logger("correct_cds")

        copied = transcript.deepcopy()

        trimmed = trim_coding(copied, logger, max_length=50)
        self.assertEqual(trimmed.start, 11400)
        self.assertEqual(trimmed.end, 19050)

        copied = transcript.deepcopy()
        self.assertEqual(copied.start, 10000)
        trimmed = trim_coding(copied, logger, max_length=200)
        self.assertEqual(trimmed.start, 11300)
        self.assertEqual(trimmed.end, 19200)

    def test_noncoding(self):

        logger = utilities.log_utils.create_null_logger("correct_cds2")
        transcript = Transcript(logger=logger)
        transcript.chrom = "Chr1"
        transcript.source = "test"
        transcript.start = 10000
        transcript.end = 20000

        exons = [(10000, 11500),
                 (12000, 13000),
                 (15000, 18000),
                 (19000, 20000)]

        transcript.add_exons(exons)
        transcript.strand = "+"
        transcript.finalize()

        copied = transcript.deepcopy()

        trimmed = trim_noncoding(copied, max_length=50)
        self.assertEqual(trimmed.start, 11450)
        self.assertEqual(trimmed.end, 19050)

        copied = transcript.deepcopy()

        trimmed = trim_noncoding(copied, max_length=200)
        self.assertEqual(trimmed.start, 11300)
        self.assertEqual(trimmed.end, 19200)


class TestPadding(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

    def setUp(self):
        self.reference = "Chr5\t26574999\t26578625\tID=AT5G66600.3;coding=True;phase=0\t0\t-\t26575104\t26578315\t0\t11\t411,126,87,60,100,809,126,72,82,188,107\t0,495,711,885,1035,1261,2163,2378,2856,3239,3519"
        self.reference = Transcript(BED12(self.reference), source="TAIR10", is_reference=True)

    def test_basic_padding(self):

        logger = create_null_logger("test_basic_padding")
        logger.setLevel("INFO")
        template = self.reference.copy()
        template.id = "AT5G66600.3_exp"
        template.strip_cds()
        template.unfinalize()

        template.remove_exon((26575000, 26575410))  # First exon
        template.start = 26574650
        template.add_exon((26574970, 26575410))  # New exon, template at 5'
        template.add_exon((26574650, 26574820))  # New UTR exon

        template.remove_exon((26578519, 26578625))  # Last exon
        template.end = 26579700
        template.add_exon((26578519, 26578725))
        template.add_exon((26579325, 26579700))
        
        template.finalize()
        fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

        new5 = expand_transcript(self.reference, None, template, fai, logger)
        self.assertIn((26574970, 26575410), new5.exons)
        self.assertIn((26574650, 26574820), new5.exons)
        self.assertEqual(template.start, new5.start)
        self.assertEqual(self.reference.end, new5.end)

        new3 = expand_transcript(self.reference, template, None, fai, logger)
        self.assertIn((26578519, 26578725), new3.exons)
        self.assertIn((26579325, 26579700), new3.exons)
        self.assertEqual(self.reference.start, new3.start)
        self.assertEqual(template.end, new5.end)

        new53 = expand_transcript(self.reference, template, template, fai, logger)
        self.assertIn((26574970, 26575410), new53.exons)
        self.assertIn((26574650, 26574820), new53.exons)
        self.assertIn((26578519, 26578725), new53.exons)
        self.assertIn((26579325, 26579700), new53.exons)
        self.assertEqual(template.start, new53.start)
        self.assertEqual(template.end, new53.end)

    def test_locus_padding_equal_or_n(self):

        for num, exons_to_add in enumerate([
            ((26574970, 26575410), (26578519, 26578725)),
            ((26574970, 26575410), (26574650, 26574820),
             (26578519, 26578725), (26579325, 26579700))]):

            for num2, pad_transcripts in enumerate((False, True)):
                with self.subTest(exons_to_add=exons_to_add, pad_transcripts=pad_transcripts):
                    logger = create_null_logger("test_locus_padding_equal_or_n_" + str(num + num2 * 2))
                    logger.setLevel("DEBUG")
                    template = self.reference.copy()
                    del template.is_reference
                    template.id = "AT5G66600.3_exp"
                    template.unfinalize()
                    template.remove_exon((26575000, 26575410))  # First exon
                    template.remove_exon((26578519, 26578625))  # Last exon
                    template.start = min([_[0] for _ in exons_to_add])
                    template.end = max([_[1] for _ in exons_to_add])
                    template.add_exons(exons_to_add)  # New exon, template at 5'
                    template.finalize()
                    json_conf = to_json(None)
                    json_conf["reference"]["genome"] = self.fai
                    json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] = False
                    json_conf["pick"]["run_options"]["only_reference_update"] = True
                    locus = Locus(self.reference.copy(), logger=logger, json_conf=json_conf,
                                  pad_transcripts=pad_transcripts)
                    self.assertTrue(locus[self.reference.id].is_reference)
                    self.assertEqual(locus.perform_padding, pad_transcripts)
                    locus.add_transcript_to_locus(template)
                    if pad_transcripts is True:
                        self.assertIn(template.id, locus)
                    locus.finalize_alternative_splicing()
                    self.assertNotIn(template.id, locus)
                    if pad_transcripts is False:
                        self.assertEqual(locus[self.reference.id].start, self.reference.start)
                        self.assertEqual(locus[self.reference.id].end, self.reference.end)
                    else:
                        self.assertTrue(locus.perform_padding)
                        self.assertEqual(locus[self.reference.id].start, template.start,
                                         (locus[self.reference.id].exons[0], template.exons[0]))
                        self.assertEqual(locus[self.reference.id].end, template.end,
                                         (locus[self.reference.id].end, template.end))
                        self.assertNotIn(template.id, locus)

    @mark.triage
    def test_removal_after_padding(self):

        logger = create_default_logger("test_add_two_partials", "INFO")
        json_conf = to_json(None)
        json_conf["reference"]["genome"] = self.fai
        json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] = False
        json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] = True
        json_conf["pick"]["alternative_splicing"]["pad"] = True

        t1 = Transcript(BED12(
            "Chr5\t26584779\t26587869\tID=AT5G66610.1;coding=True;phase=0\t0\t+\t26585222\t26587755\t0\t11\t\
100,54,545,121,78,105,213,63,119,59,443\t0,440,565,1202,1437,1640,1858,2154,2304,2507,2647"
        ))

        t2_1 = Transcript(BED12(
            "Chr5\t26584773\t26586510\tID=AT5G66610.2_1;coding=True;phase=0\t0\t+\t26585222\t26586510\t0\t\
6\t177,54,545,121,78,85\t0,446,571,1208,1443,1652"
        ))
        t2_2 = Transcript(BED12(
            "Chr5\t26584873\t26587782\tID=AT5G66610.2_2;coding=True;phase=0\t0\t+\t26585222\t\
26587755\t0\t10\t77,54,545,121,78,99,213,63,119,496\t0,346,471,1108,1343,1552,1764,2060,2210,2413"
        ))

        t1.finalize()
        t2_1.finalize()
        t2_2.finalize()
        t1.is_reference = True

        locus = Locus(t1, logger=logger, json_conf=json_conf)
        locus.add_transcript_to_locus(t2_1, check_in_locus=False)
        locus.add_transcript_to_locus(t2_2, check_in_locus=False)
        self.assertTrue(locus.primary_transcript_id == t1.id)
        # locus.logger.setLevel("DEBUG")
        locus.finalize_alternative_splicing(_scores={t1.id: 20, t2_1.id: 15, t2_2.id: 10})

        self.assertIn(t1.id, locus.transcripts)
        if t2_1.id in locus.transcripts:
            from ..scales import Assigner
            import itertools
            for tid1, tid2 in itertools.combinations(locus.transcripts.keys(), 2):
                res, _ = Assigner.compare(locus[tid1], locus[tid2])
                print(tid1, tid2, res.ccode)
            self.assertNotIn(t2_1.id, locus.transcripts)

        self.assertIn(t2_2.id, locus.transcripts, "\n".join(tr.format("bed12") for tr in locus))
        self.assertTrue(locus[t2_2.id].attributes["padded"])
        self.assertEqual(locus[t2_2.id].start, t1.start,
                         ((locus[t2_2.id].start, t1.start, t2_1.start, t2_2.start),
                          (locus[t2_2.id].end, t1.end, t2_1.end, t2_2.end)),
                         )

    def test_add_two_partials(self):

        logger = create_null_logger("test_add_two_partials")
        logger.setLevel("INFO")
        json_conf = to_json(None)
        json_conf["reference"]["genome"] = self.fai
        json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] = False
        json_conf["pick"]["run_options"]["only_reference_update"] = True

        ref = Transcript(is_reference=True)
        ref.chrom, ref.strand, ref.id = "Chr5", "-", "AT5G66670.2"
        ref.add_exons([(26611258, 26612889)])
        ref.add_exons([(26611474, 26612700)], features=["CDS"])
        ref.finalize()
        self.assertTrue(ref.is_coding)

        # Chr5	TAIR10	mRNA	26611258	26612889	.	-	.	ID=AT5G66670.2;Parent=AT5G66670;Name=AT5G66670.2;index=1
        # Chr5	TAIR10	protein	26611474	26612700	.	-	.	ID=AT5G66670.2-Protein;Parent=AT5G66670.2;Name=AT5G66670.2;derives_from=AT5G66670.2
        # Chr5	TAIR10	three_prime_UTR	26611258	26611473	.	-	.	Parent=AT5G66670.2
        # Chr5	TAIR10	CDS	26611474	26612700	.	-	0	Parent=AT5G66670.2
        # Chr5	TAIR10	five_prime_UTR	26612701	26612889	.	-	.	Parent=AT5G66670.2
        # Chr5	TAIR10	exon	26611258	26612889	.	-	.	Parent=AT5G66670.2

        template1 = Transcript(is_reference=False)
        template1.chrom, template1.strand, template1.id = ref.chrom, ref.strand, ref.id + "_frag1"
        template1.add_exons(((26611116, 26611157), (26611258, 26612670)))
        template1.add_exons(((26611474, 26612670),), features=["CDS"])
        template1.finalize()
        self.assertTrue(template1.is_coding)

        template2 = Transcript(is_reference=False)
        template2.chrom, template2.strand, template2.id = ref.chrom, ref.strand, ref.id + "_frag2"
        template2.add_exons(((26611574, 26612889), (26613007, 26613403)))
        template2.add_exons(((26611574, 26612700),), features=["CDS"])
        template2.finalize()
        self.assertTrue(template2.is_coding)

        logger.setLevel("INFO")
        locus = Locus(ref, json_conf=json_conf, logger=logger, pad_transcripts=True)
        locus.add_transcript_to_locus(template1)
        locus.add_transcript_to_locus(template2)
        self.assertIn(template2.id, locus)
        # self.assertIn(template1.id, locus)
        # locus.logger.setLevel("DEBUG")
        # for tid in locus:
        #     locus[tid].logger.setLevel("DEBUG")
        locus.finalize_alternative_splicing()
        self.assertTrue(locus._finalized)
        self.assertNotIn(template1.id, locus, "\n" + str(locus))
        self.assertNotIn(template2.id, locus, "\n" + str(locus))
        self.assertEqual(locus[ref.id].end, template2.end,
                         ((locus[ref.id].end, ref.end, template2.end, template1.end),
                         (locus[ref.id].start, ref.start, template2.start, template1.start))
                         )


if __name__ == "__main__":
    unittest.main()
