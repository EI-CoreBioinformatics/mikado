import unittest
from Mikado.loci import Transcript
from Mikado.parsers.GFF import GffLine
from Mikado.utilities.log_utils import create_default_logger, create_null_logger


class PhaseChecker(unittest.TestCase):

    logger = create_default_logger("pcheck")
    logger.setLevel("DEBUG")

    def setUp(self):

        lines = """Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	mRNA	40282	46004	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960;Name=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2;aed=0.0;note=TRIAE_CS42_5DL_TGACv1_434051_AA1427960;confidence=High;has_start=True;has_stop=True;original_stop=True;protein_rank=P1;transcript_rank=T2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	40282	40933	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon1;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	three_prime_UTR	40282	40720	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.three_prime_UTR1;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	40721	40933	.	-	0	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS1;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	41018	41111	.	-	1	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS2;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	41018	41111	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon2;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	41227	41468	.	-	0	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS3;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	41227	41468	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon3;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	41673	41831	.	-	0	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS4;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	41673	41831	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon4;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	41946	42820	.	-	2	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS5;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	41946	42820	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon5;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	42905	42913	.	-	2	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS6;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	42905	42913	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon6;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	45373	45496	.	-	0	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS7;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	45373	45496	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon7;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	45600	45651	.	-	1	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS8;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	45600	45651	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon8;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	45726	45726	.	-	2	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS9;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	45726	45726	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon9;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	CDS	45875	45893	.	-	0	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.CDS10;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	exon	45875	46004	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.exon10;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2
Triticum_aestivum_CS42_TGACv1_scaffold_434051_5DL	TGACv1	five_prime_UTR	45894	46004	.	-	.	ID=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2.five_prime_UTR1;Parent=TRIAE_CS42_5DL_TGACv1_434051_AA1427960.2"""

        lines = [GffLine("\t".join(_.split())) for _ in lines.split("\n") if _]
        self.transcript = Transcript(lines[0], logger=self.logger)
        self.transcript.add_exons(lines[1:])
        self.correct_phases = {(40721, 40933): 2,
                               (41018, 41111): 0,
                               (41227, 41468): 2,
                               (41673, 41831): 2,
                               (41946, 42820): 1,
                               (42905, 42913): 1,
                               (45373, 45496): 2,
                               (45600, 45651): 0,
                               (45726, 45726): 2,
                               (45875, 45893): 0}

    @unittest.skip
    def test_check_phases(self):
        self.transcript.finalize()
        phases = dict((_[1], _[2]) for _ in self.transcript.internal_orfs[0]
                      if _[0] == "CDS")
        self.assertEqual(self.transcript.combined_cds_start, 45893)

        self.assertEqual(phases.keys(),
                         self.correct_phases.keys(),
                         list(zip(sorted(phases.keys()),
                                  sorted(self.correct_phases.keys()))))

        if self.correct_phases != phases:
            for key in sorted(phases.keys(), reverse=True):
                self.assertEqual(phases[key], self.correct_phases[key],
                                 (key, phases[key], self.correct_phases[key]))

        self.assertEqual(self.correct_phases,
                         phases,
                         (self.correct_phases, phases))

if __name__ == "__main__":
    unittest.main()
