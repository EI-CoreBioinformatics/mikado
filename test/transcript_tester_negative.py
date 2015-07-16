import unittest
import re
import mikado_lib

class TranscriptTester(unittest.TestCase):

    tr_gff="""Chr1    TAIR10    mRNA    5928    8737    .    -    .    ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
Chr1    TAIR10    five_prime_UTR    8667    8737    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    8571    8666    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    8571    8737    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    8417    8464    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    8417    8464    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    8236    8325    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    8236    8325    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    7942    7987    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    7942    7987    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    7762    7835    .    -    2    Parent=AT1G01020.1;
Chr1    TAIR10    exon    7762    7835    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    7564    7649    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    7564    7649    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    7384    7450    .    -    1    Parent=AT1G01020.1;
Chr1    TAIR10    exon    7384    7450    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    7157    7232    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    7157    7232    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    6915    7069    .    -    2    Parent=AT1G01020.1;
Chr1    TAIR10    three_prime_UTR    6437    6914    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    exon    6437    7069    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    three_prime_UTR    5928    6263    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    exon    5928    6263    .    -    .    Parent=AT1G01020.1"""

    tr_lines = tr_gff.split("\n")
    for pos,line in enumerate(tr_lines):
        tr_lines[pos] = re.sub("\s+", "\t", line)
        assert len(tr_lines[pos].split("\t"))==9, line.split("\t")
        
    tr_gff_lines = [mikado_lib.parsers.GFF.gffLine(line) for line in tr_lines]
    
    for l in tr_gff_lines:
        assert l.header is False
#         print(l)
    
    def setUp(self):
        '''Basic creation test.'''
        
        self.tr = mikado_lib.loci_objects.transcript.transcript(self.tr_gff_lines[0])
        for line in self.tr_gff_lines[1:]:
            self.tr.addExon(line)
        self.tr.finalize()

        self.orf=mikado_lib.parsers.bed12.BED12()
        self.orf.chrom=self.tr.id
        self.orf.start = 0
        self.orf.end = self.tr.cdna_length
        self.orf.name = self.tr.id
        self.orf.strand="+"
        self.orf.score = 0
        self.orf.thickStart = self.tr.selected_start_distance_from_tss+1
        self.orf.thickEnd = self.tr.cdna_length - self.tr.selected_end_distance_from_tes
        self.orf.blockCount = 1
        self.orf.blockSize = self.tr.cdna_length
        self.orf.blockStarts = 0
        self.orf.has_start_codon = True
        self.orf.has_stop_codon = True
            

    def test_basics(self):
        
        self.assertEqual(self.tr.chrom, "Chr1")
        self.assertEqual(self.tr.strand, "-")
        self.assertEqual(self.tr.exon_num, 10)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 5928)
        self.assertEqual(self.tr.end, 8737)
        self.assertEqual(self.tr.exons,
                         [(5928,6263),(6437,7069),(7157,7232),(7384,7450),(7564,7649),(7762,7835),(7942,7987),(8236,8325),(8417,8464),(8571,8737)],
                         self.tr.exons)
        
    def test_cds(self):
        self.assertEqual(self.tr.combined_cds, self.tr.selected_cds)

        self.assertEqual(self.tr.combined_cds,
                         [(6915,7069),(7157,7232),(7384,7450),(7564,7649),(7762,7835),(7942,7987),(8236,8325),(8417,8464),(8571,8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8666)
        self.assertEqual(self.tr.selected_cds_end, 6915)
         
    def test_utr(self):
        self.assertEqual(self.tr.five_utr,[("UTR", 8667,8737)] )
        self.assertEqual(self.tr.three_utr,[("UTR",5928,6263 ), ("UTR", 6437,6914)] )
 
    def test_utr_metrics(self):
 
        '''Test for UTR exon num, start distance, etc.'''
 
        self.assertEqual(self.tr.five_utr_num, 1)
        self.assertEqual(self.tr.three_utr_num, 2)
        self.assertEqual(self.tr.five_utr_length, 8737+1-8667)
        self.assertEqual(self.tr.three_utr_length, 6263+1-5928 + 6914+1-6437 )
        self.assertEqual(self.tr.selected_start_distance_from_tss,8738-8667, self.tr.selected_end_distance_from_tes )
        self.assertEqual(self.tr.selected_end_distance_from_tes,6263+1-5928 + 6915-6437, self.tr.selected_end_distance_from_tes )
 
    def test_strip_cds(self):
 
        self.tr.strip_cds()
        self.assertEqual(self.tr.selected_cds_length, 0)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.selected_cds, [])
        self.assertEqual(self.tr.selected_cds_start, None)
        self.assertEqual(self.tr.selected_cds_end, None)
 
    def test_remove_utr(self):
        '''Test for CDS stripping. We remove the UTRs and verify that start/end have moved, no UTR is present, etc.'''
 
        self.tr.remove_utrs()
        self.assertEqual(self.tr.selected_cds_start, self.tr.end)
        self.assertEqual(self.tr.selected_cds_end, self.tr.start)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.combined_cds,
                         [(6915,7069),(7157,7232),(7384,7450),(7564,7649),(7762,7835),(7942,7987),(8236,8325),(8417,8464),(8571,8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_utr, [], self.tr.combined_utr)
 
    def test_load_orf(self):
         
        '''Test for loading a single ORF. We strip the CDS and reload it.'''
         
        self.tr.strip_cds()
        self.tr.load_orfs( [self.orf] )
        self.assertEqual(self.tr.combined_cds,
                         [(6915,7069),(7157,7232),(7384,7450),(7564,7649),(7762,7835),(7942,7987),(8236,8325),(8417,8464),(8571,8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8666)
        self.assertEqual(self.tr.selected_cds_end, 6915)
         
    def test_negative_orf(self):
        '''Test loading a negative strand ORF onto a multiexonic transcript. This should have no effect.'''
         
        self.orf.strand = "-"
        self.tr.strip_cds()
        self.tr.load_orfs( [self.orf] )
        self.assertEqual(self.tr.selected_cds_start, None)

unittest.main()        