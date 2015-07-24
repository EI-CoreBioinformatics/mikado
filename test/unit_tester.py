import unittest
import os.path
from mikado_lib import json_utils
from mikado_lib import exceptions
from mikado_lib.parsers import GFF#,GTF, bed12
from mikado_lib.loci_objects import transcript, superlocus#,abstractlocus

class LocusTester(unittest.TestCase):
     
    def test_locus(self):
        '''Basic testing of the locus functionality.'''
 
        gff_transcript1="""Chr1\tfoo\ttranscript\t101\t200\t.\t+\t.\tID=t0
Chr1\tfoo\texon\t101\t200\t.\t+\t.\tID=t0:exon1;Parent=t0""".split("\n")
        gff_transcript1=[GFF.gffLine(x) for x in gff_transcript1]
        self.assertEqual(gff_transcript1[0].chrom, "Chr1", gff_transcript1[0])
        transcript1=transcript.transcript(gff_transcript1[0])
        for exon in gff_transcript1[1:]: transcript1.addExon(exon)
        transcript1.finalize()
        self.assertTrue(transcript1.monoexonic)
        self.assertEqual(transcript1.chrom, gff_transcript1[0].chrom)
         
        gff_transcript2="""Chr1\tfoo\ttranscript\t101\t600\t.\t+\t.\tID=t1
Chr1\tfoo\texon\t101\t200\t.\t+\t.\tID=t1:exon1;Parent=t1
Chr1\tfoo\texon\t301\t400\t.\t+\t.\tID=t1:exon2;Parent=t1
Chr1\tfoo\texon\t501\t600\t.\t+\t.\tID=t1:exon3;Parent=t1""".split("\n")
        gff_transcript2=[GFF.gffLine(x) for x in gff_transcript2]
        transcript2=transcript.transcript(gff_transcript2[0])
        for exon in gff_transcript2[1:-1]: transcript2.addExon(exon)
        #Test that a transcript cannot be finalized if the exons do not define the external boundaries
        with self.assertRaises(exceptions.InvalidTranscript):
            transcript2.finalize()
        transcript2.addExon(gff_transcript2[-1])
        transcript2.finalize()
        self.assertFalse(transcript2.monoexonic)
        self.assertEqual(transcript2.exon_num, len(gff_transcript2)-1)
        #Test that trying to modify a transcript after it has been finalized causes errors
        with self.assertRaises(exceptions.ModificationError):
            for exon in gff_transcript2[1:]: transcript2.addExon(exon)
        #Test that creating a superlocus without configuration fails
        with self.assertRaises(exceptions.NoJsonConfigError):
            _=superlocus.superlocus(transcript1)
        my_json =  os.path.join(os.path.dirname(__file__),  "../sample_data/configuration.yaml")
        my_json=json_utils.to_json(my_json)
        slocus=superlocus.superlocus(transcript1, json_dict=my_json)
        slocus.add_transcript_to_locus(transcript2)
        self.assertEqual(slocus.strand, transcript1.strand )
        self.assertEqual(slocus.start, min(transcript1.start,transcript2.start) )
        self.assertEqual(slocus.end, max(transcript1.end,transcript2.end) )
        slocus.define_subloci()
        self.assertEqual(len(slocus.subloci), 2)
        slocus.define_monosubloci()
        self.assertEqual(len(slocus.monosubloci), 2)
        slocus.define_loci()
        self.assertEqual(len(slocus.loci), 1)
        self.assertEqual(list(slocus.loci[0].transcripts.keys())[0], "t1")
        gff_transcript3="""Chr1\tfoo\ttranscript\t101\t200\t.\t-\t.\tID=tminus0
Chr1\tfoo\texon\t101\t200\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0""".split("\n")
        gff_transcript3=[GFF.gffLine(x) for x in gff_transcript3]
        transcript3=transcript.transcript(gff_transcript3[0])
        for exon in gff_transcript3[1:]: transcript3.addExon(exon)
        transcript3.finalize()
        minusuperlocus=superlocus.superlocus(transcript3, json_dict=my_json)
        minusuperlocus.define_loci()
        self.assertEqual(len(minusuperlocus.loci),1)
        self.assertTrue(transcript3.strand!=transcript1.strand)
        minuslocus=minusuperlocus.loci[0]
        pluslocus=slocus.loci[0]
        self.assertTrue(pluslocus.other_is_fragment(minuslocus))


unittest.main()