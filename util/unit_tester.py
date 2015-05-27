import os,sys
import unittest
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects import *
from loci_objects import json_utils

class Tester(unittest.TestCase):
    def test_overlap(self):
        self.assertEqual(abstractlocus.abstractlocus.overlap(
                                                             (100,200),
                                                             (150,200)
                                                             ), 50)
    def test_gff3(self):
        gff=GFF.GFF3(open(os.path.join(os.path.dirname(__file__), "../sample_data/sample.gff3")))
        line=next(gff)
        self.assertEqual(line._line, "##gff-version 3\n")
        self.assertTrue(line.header)
        while line.header is True:
            line=next(gff)
        self.assertEqual(line.chrom, "Chr1")
        self.assertEqual(line.start, 101)
        self.assertEqual(line.id, "t0")
        gff.close()
        
    def test_bed(self):
        bed=bed12.BED12(open(os.path.join(os.path.dirname(__file__), "../sample_data/sample.bed")))
        line=next(bed)
        self.assertEqual(line.chrom, "t0")
        self.assertEqual(line.cdsStart, 21)
        bed.close()
        
    def test_transcript(self):
        gff=GFF.GFF3(open(os.path.join(os.path.dirname(__file__), "../sample_data/sample.gff3")))
        line=None
        tr=transcript.transcript()
        self.assertIs(tr.chrom, None)

        while line is None or line.header is True:
            line=next(gff)
        self.assertTrue(line.is_transcript)
        tr=transcript.transcript(line)
        self.assertEqual(line.chrom, tr.chrom)
        self.assertEqual(line.start, tr.start)
        self.assertEqual(line.id, tr.id)
        line=next(gff)
        tr.addExon(line)
        tr.finalize()
        self.assertTrue(tr.monoexonic)
        gff.close()

    def test_json(self):
        my_json = os.path.join(os.path.dirname(__file__), "../sample_data/scoring.json")
        json_utils.to_json(my_json)
        

unittest.main()