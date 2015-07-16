import unittest
import re
import mikado_lib

class TranscriptTester(unittest.TestCase):

    tr_gff="""Chr2    TAIR10    mRNA    626842    628676    .    +    .    ID=AT2G02380.1;Parent=AT2G02380;Name=AT2G02380.1;Index=1
Chr2    TAIR10    exon    626842    626880    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    five_prime_UTR    626842    626877    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    626878    626880    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    626963    627059    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    626963    627059    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627137    627193    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627137    627193    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    627312    627397    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627312    627397    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    627488    627559    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627488    627559    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627696    627749    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627696    627749    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627840    627915    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627840    627915    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    628044    628105    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628044    628105    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    628182    628241    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628182    628241    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    628465    628676    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628465    628569    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    three_prime_UTR    628570    628676    .    +    .    Parent=AT2G02380.1"""

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
        
        self.assertEqual(self.tr.chrom, "Chr2")
        self.assertEqual(self.tr.strand, "+")
        self.assertEqual(self.tr.exon_num, 10)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 626842)
        self.assertEqual(self.tr.end, 628676)
        self.assertEqual(self.tr.exons,
                         [(626842,626880),(626963,627059),(627137,627193),(627312,627397),(627488,627559),(627696,627749),(627840,627915),(628044,628105),(628182,628241),(628465,628676)],
                         self.tr.exons)
        
    def test_cds(self):
        self.assertEqual(self.tr.combined_cds, self.tr.selected_cds)

        self.assertEqual(self.tr.combined_cds,
                         [(626878,626880),(626963,627059),(627137,627193),(627312,627397),(627488,627559),(627696,627749),(627840,627915),(628044,628105),(628182,628241),(628465,628569)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 626878)
        self.assertEqual(self.tr.selected_cds_end, 628569)
        
    def test_utr(self):
        self.assertEqual(self.tr.five_utr,[("UTR", 626842,626877)] )
        self.assertEqual(self.tr.three_utr,[("UTR", 628570,628676)] )

    def test_utr_metrics(self):

        '''Test for UTR exon num, start distance, etc.'''

        self.assertEqual(self.tr.five_utr_num, 1)
        self.assertEqual(self.tr.three_utr_num, 1)
        self.assertEqual(self.tr.five_utr_length, 626877+1-626842)
        self.assertEqual(self.tr.three_utr_length, 628676+1-628570)
        
        self.assertEqual(self.tr.selected_start_distance_from_tss,626878-626842, self.tr.selected_end_distance_from_tes )
        self.assertEqual(self.tr.selected_end_distance_from_tes,628676-628569, self.tr.selected_end_distance_from_tes )

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
        self.assertEqual(self.tr.selected_cds_start, self.tr.start)
        self.assertEqual(self.tr.selected_cds_end, self.tr.end)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.combined_cds,
                         [(626878,626880),(626963,627059),(627137,627193),(627312,627397),(627488,627559),(627696,627749),(627840,627915),(628044,628105),(628182,628241),(628465,628569)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_utr, [], self.tr.combined_utr)

    def test_load_orf(self):
        
        '''Test for loading a single ORF. We strip the CDS and reload it.'''
        
        self.tr.strip_cds()
        self.tr.load_orfs( [self.orf] )
        self.assertEqual(self.tr.combined_cds,
                         [(626878,626880),(626963,627059),(627137,627193),(627312,627397),(627488,627559),(627696,627749),(627840,627915),(628044,628105),(628182,628241),(628465,628569)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 626878)
        self.assertEqual(self.tr.selected_cds_end, 628569)
        
    def test_negative_orf(self):
        '''Test loading a negative strand ORF onto a multiexonic transcript. This should have no effect.'''
        
        self.orf.strand = "-"
        self.tr.strip_cds()
        self.tr.load_orfs( [self.orf] )
        self.assertEqual(self.tr.selected_cds_start, None)

unittest.main()        