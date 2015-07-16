import unittest
import mikado_lib

class test_checker(unittest.TestCase):
    
    def test_self(self):
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 2000
        query.strand = "+"
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.exons = [(100,300), (500, 1000), (1500, 2000)]
        query.finalize()
        
        result, _ = mikado_lib.scales.assigner.assigner.compare(query, query)
        self.assertEqual(result.ccode, ("=",))
        self.assertEqual(result.n_f1, (100,))
        self.assertEqual(result.n_prec, (100,))
        self.assertEqual(result.n_recall, (100,))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))

    def test_mono_self(self):
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 1000
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.strand = "+"
        query.exons = [(100,1000)]
        query.finalize()
        
        result, _ = mikado_lib.scales.assigner.assigner.compare( query, query)
        self.assertEqual(result.ccode, ("_",))
        self.assertEqual(result.n_f1, (100,))
        self.assertEqual(result.n_prec, (100,))
        self.assertEqual(result.n_recall, (100,))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        
    def test_equal(self):
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 2000
        query.strand = "+"
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.exons = [(100,300), (500, 1000), (1500, 2000)]
        query.finalize()
        
        target = mikado_lib.loci_objects.transcript.transcript()
        target.start = 200
        target.end = 1800
        target.strand = "+"
        target.chrom = "Chr1"
        target.id = "P1.1"
        target.parent = "P1"
        target.exons = [(200,300), (500, 1000), (1500, 1800)]
        target.finalize()
        
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("=",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual( result.n_prec[0],  100    , delta=0.1)
        self.assertAlmostEqual( result.n_recall[0],  100*(300-200+1 + 1000-500+1+1800-1500+1)/query.cdna_length, delta=0.1)
        
        target.strand = "-"
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("x",))
        
    def test_mono_equal(self):
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 1000
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.strand = "+"
        query.exons = [(100,1000)]
        query.finalize()
        
        target = mikado_lib.loci_objects.transcript.transcript()
        target.start = 105
        target.end = 995
        target.chrom = "Chr1"
        target.id = "G1.1"
        target.parent = "G1"
        target.strand = "+"
        target.exons = [(105,995)]
        target.finalize()
        
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("_",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual( result.n_prec[0],  100    , delta=0.1)
        self.assertAlmostEqual( result.n_recall[0],  100*(995-105+1)/query.cdna_length, delta=0.1)
        
        target.strand = "-"
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("x",))
        
    def test_mono_semiequal(self):
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 1000
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.strand = "+"
        query.exons = [(100,1000)]
        query.finalize()
        
        target = mikado_lib.loci_objects.transcript.transcript()
        target.start = 200
        target.end = 900
        target.chrom = "Chr1"
        target.id = "G1.1"
        target.parent = "G1"
        target.strand = "+"
        target.exons = [(200,900)]
        target.finalize()
        
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("c",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual( result.n_prec[0],  100    , delta=0.1)
        self.assertAlmostEqual( result.n_recall[0],  100*(900-200+1)/query.cdna_length, delta=0.1)
        
        target.strand = "-"
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("x",))
        
    def test_mono_overlap(self):
        
        '''Test that two monoexonic overlapping genes get a m'''
        
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 1000
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.strand = "+"
        query.exons = [(100,1000)]
        query.finalize()

        
        target = mikado_lib.loci_objects.transcript.transcript()
        target.start = 200
        target.end = 1200
        target.chrom = "Chr1"
        target.id = "G1.1"
        target.parent = "G1"
        target.strand = "+"
        target.exons = [(200,1200)]
        target.finalize()

        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("m",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual( result.n_prec[0],  100*801/target.cdna_length, delta=0.1)
        self.assertAlmostEqual( result.n_recall[0], 100*801/query.cdna_length, delta=0.1)
        
    def test_contained(self):
    
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 2000
        query.strand = "+"
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.exons = [(100,300), (500, 1000), (1500, 2000)]
        query.finalize()
        
        target = mikado_lib.loci_objects.transcript.transcript()
        target.start = 600
        target.end = 1800
        target.strand = "+"
        target.chrom = "Chr1"
        target.id = "P1.1"
        target.parent = "P1"
        target.exons = [(600, 1000), (1500, 1800)]
        target.finalize()
        
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("c",))
        self.assertAlmostEqual(result.j_f1[0], (100*(2/3),)[0], delta=0.1)
        self.assertAlmostEqual(result.j_prec, (100,))
        self.assertAlmostEqual(result.j_recall, (100*1/2,), delta=0.1)
        self.assertAlmostEqual( result.n_prec[0],  100    , delta=0.1)
#         self.assertAlmostEqual( result.n_recall[0],  100*(300-200+1 + 1000-500+1+1800-1500+1)/query.cdna_length, delta=0.1)
        
        target.strand = "-"
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("x",))
    
    def test_alternative(self):
    
        query = mikado_lib.loci_objects.transcript.transcript()
        query.start = 100
        query.end = 2000
        query.strand = "+"
        query.chrom = "Chr1"
        query.id = "G1.1"
        query.parent = "G1"
        query.exons = [(100,300), (500, 1000), (1500, 2000)]
        query.finalize()
        
        target = mikado_lib.loci_objects.transcript.transcript()
        target.start = 100
        target.end = 2000
        target.strand = "+"
        target.chrom = "Chr1"
        target.id = "P1.1"
        target.parent = "P1"
        target.exons = [ (100,300), (600, 1000), (1500, 2000)]
        target.finalize()
        
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], (100*((3/4)/1),)[0], delta=0.1)
        self.assertAlmostEqual(result.j_prec, (100*3/4,), delta=0.1)
        self.assertAlmostEqual(result.j_recall, (100*3/4,), delta=0.1)
#         self.assertAlmostEqual( result.n_prec[0],  100    , delta=0.1)
#         self.assertAlmostEqual( result.n_recall[0],  100*(300-200+1 + 1000-500+1+1800-1500+1)/query.cdna_length, delta=0.1)
        
        target.strand = "-"
        result, _ = mikado_lib.scales.assigner.assigner.compare( target, query)
        self.assertEqual(result.ccode, ("j",))

    
        
unittest.main() 