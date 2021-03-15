# coding: utf-8

"""
Very basic, all too basic test for some functionalities of locus-like classes.
"""
import dataclasses
import operator
import random
import re
import unittest
import os.path
import logging
from collections import namedtuple
from copy import deepcopy

import io
import marshmallow
import pkg_resources
import pytest
from numpy import arange
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from Mikado._transcripts.scoring_configuration import NumBoolEqualityFilter, ScoringFile, SizeFilter, RangeFilter, \
    InclusionFilter, MinMaxScore, TargetScore
from Mikado._transcripts.transcript_base import Metric
from Mikado.configuration import configurator, MikadoConfiguration, DaijinConfiguration
from Mikado import exceptions
from Mikado.parsers import GFF  # ,GTF, bed12
from Mikado.parsers.GTF import GtfLine
from Mikado.transcripts.transcript import Transcript
from Mikado.loci import Superlocus, Abstractlocus, Locus, Monosublocus, MonosublocusHolder, Sublocus, Excluded
from Mikado.loci.locus import pad_transcript
from Mikado.transcripts.reference_gene import Gene
from Mikado.utilities.log_utils import create_null_logger, create_default_logger
from Mikado.utilities.overlap import overlap
import itertools
from Mikado.utilities import Interval
from Mikado import loci
from Mikado.loci.abstractlocus import to_bool
import pickle
import inspect
from Mikado.parsers.bed12 import BED12
import pysam
from pytest import mark
from itertools import combinations_with_replacement
from Mikado.scales.assignment.assigner import Assigner
logging.getLogger("matplotlib").setLevel(logging.WARNING)


class OverlapTester(unittest.TestCase):

    def test_overlap(self):
        """
        Test for overlap function
        :return:
        """

        self.assertEqual(Abstractlocus.overlap((100, 200), (100, 200)),
                         100)
        self.assertEqual(Abstractlocus.overlap((100, 200), (100, 200)),
                         overlap((100, 200), (100, 200)))


class ToBoolTester(unittest.TestCase):

    def test_bool(self):
        for param in (1, True, "True", "tRue", "true", 1.0):
            with self.subTest(param=param):
                self.assertEqual(to_bool(param), True)
        for param in (0, False, "False", "faLse", "false", 0.0):
            with self.subTest(param=param):
                self.assertEqual(to_bool(param), False)
        for param in (2, 5.5, "hello", dict(), "Falsee", "Trrue"):
            with self.subTest(param=param), self.assertRaises((AttributeError, TypeError, ValueError)):
                to_bool(param)


class ExcludedTester(unittest.TestCase):

    def test_excluded(self):
        transcript = Transcript(BED12(
            "Chr5\t26603002\t26604376\tID=AT5G66650.1;coding=True;phase=0\t0\t-\t26603203\t26604257\t0\t2\t636,650\t0,724"
        ))
        transcript.finalize()
        excluded = Excluded(transcript)
        transcript2 = Transcript(BED12(
            "Chr5\t26611257\t26612891\tID=AT5G66670.1;coding=True;phase=0\t0\t-\t26611473\t26612700\t0\t2\t1470,46\t0,1588"
        ))
        comparison = Assigner.compare(transcript, transcript2)[0]
        self.assertEqual(0, comparison.n_f1[0])
        excluded.add_transcript_to_locus(transcript2)
        self.assertIn(transcript2.id, excluded)
        with self.assertRaises(NotImplementedError):
            str(excluded)
        with self.assertRaises(NotImplementedError):
            excluded.filter_and_calculate_scores()
        with self.assertRaises(NotImplementedError):
            excluded.define_monosubloci()
        with self.assertRaises(NotImplementedError):
            excluded.is_intersecting()

        excluded.remove_transcript_from_locus(transcript2.id)
        self.assertNotIn(transcript2.id, excluded)
        self.assertEqual(repr(excluded), "{name}\t{chrom}\t{start}\t{end}\t{strand}\t{transcripts}".format(
            name="excluded_transcripts",
            chrom=excluded.chrom, start=26603003, end=26604376, strand="-",
            transcripts="AT5G66650.1"))
        monosub = Monosublocus(transcript2)
        excluded.add_monosublocus(monosub)
        self.assertIn(transcript2.id, excluded)
        # test repr
        self.assertEqual(repr(excluded), "{name}\t{chrom}\t{start}\t{end}\t{strand}\t{transcripts}".format(
            name="excluded_transcripts",
            chrom=excluded.chrom, start=26603003, end=26612891, strand="-",
            transcripts="AT5G66650.1,AT5G66670.1"))


class SuperlocusTester(unittest.TestCase):

    def test_find_lost_transcripts(self):
        conf = MikadoConfiguration()
        # With this scoring, we are giving a score of 10 to transcripts with a CDS and a cDNA length >= 300.
        # We are giving a negative score to any transcript shorter than 300 bps, coding or not
        conf.scoring.scoring = {
            "cdna_length": MinMaxScore(rescaling="min", multiplier=-5,
                                       filter=SizeFilter(value=300, operator="le")),
            "combined_cds_length": MinMaxScore(rescaling="max", multiplier=10,
                                               filter=SizeFilter(value=300, operator="ge"))
        }
        conf.scoring.requirements.expression = ["cdna_length"]
        conf.scoring.requirements.parameters = {"cdna_length": SizeFilter(value=1, operator="ge")}
        conf.scoring.requirements._check_my_requirements()

        conf.scoring.cds_requirements.expression = ["cdna_length"]
        conf.scoring.cds_requirements.parameters = {"cdna_length": SizeFilter(value=1, operator="ge")}
        conf.scoring.cds_requirements._check_my_requirements()

        conf.pick.clustering.flank = 500

        t1 = Transcript(configuration=conf)
        t1.chrom, t1.start, t1.end, t1.strand, t1.id = "Chr1", 1101, 2000, "+", "valid"
        t1.add_exons([(1101, 1500), (1601, 2000)])
        t1.add_exons([(1201, 1500), (1601, 1900)], features="CDS")
        t1.finalize()
        self.assertGreater(t1.cdna_length, 300)
        self.assertGreater(t1.combined_cds_length, 300)

        t2 = Transcript(configuration=conf)
        t2.chrom, t2.start, t2.end, t2.strand, t2.id = "Chr1", 700, 900, "+", "invalid.1"
        t2.add_exons([(700, 900)])
        t2.finalize()

        t3 = Transcript(configuration=conf)
        t3.chrom, t3.start, t3.end, t3.strand, t3.id = "Chr1", 2100, 2310, "+", "invalid.2"
        t3.add_exons([(2100, 2310)])
        t3.add_exons([(2100, 2310)], features="CDS")
        t3.finalize()

        logger = create_default_logger("test_find_lost_transcripts")
        sl = Superlocus(t1, configuration=conf, flank=conf.pick.clustering.flank, logger=logger)
        sl.add_transcript_to_locus(t2)
        sl.add_transcript_to_locus(t3)
        self.assertTrue(sorted(sl.transcripts.keys()), sorted([t1.id, t2.id, t3.id]))
        sl.filter_and_calculate_scores()
        self.assertTrue(sorted(sl.transcripts.keys()), sorted([t1.id, t2.id, t3.id]))
        self.assertGreater(sl.scores[t1.id]["score"], 0, sl.scores[t1.id])
        self.assertLessEqual(sl.scores[t2.id]["score"], 0)
        self.assertLessEqual(sl.scores[t3.id]["score"], 0)

        locus = Locus(t1)
        sl.loci[locus.id] = locus
        sl._find_lost_transcripts()
        self.assertEqual(sl.lost_transcripts, dict())
        sl.configuration.pick.clustering.purge = False
        sl.logger.setLevel("DEBUG")
        sl._find_lost_transcripts()
        self.assertEqual(sorted(list(sl.lost_transcripts.keys())), sorted([t2.id, t3.id]))
        sl.loci = dict()
        sl.configuration.pick.clustering.purge = True
        sl.define_loci()
        self.assertEqual(len(sl.loci), 1)
        sl = Superlocus(t1, configuration=conf, flank=conf.pick.clustering.flank, logger=logger)
        sl.add_transcript_to_locus(t2)
        sl.add_transcript_to_locus(t3)
        sl.logger.setLevel("DEBUG")
        sl.configuration.pick.clustering.purge = False
        sl.define_loci()
        self.assertEqual(len(sl.loci), 3)

    def test_sl_is_intersecting(self):
        
        t1 = Transcript()
        t1.chrom, t1.start, t1.end, t1.strand, t1.id = "Chr1", 1101, 2000, "+", "multi.1"
        t1.add_exons([(1101, 1500), (1601, 2000)])
        t1.add_exons([(1201, 1500), (1601, 1900)], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id = "Chr1", 1001, 2200, "+", "multi.2"
        t2.add_exons([(1001, 1500), (1601, 2200)])
        t2.add_exons([(1201, 1500), (1601, 1900)], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.start, t3.end, t3.strand, t3.id = "Chr1", 901, 2200, "+", "multi.3"
        t3.add_exons([(901, 1500), (1601, 2200)])
        t3.finalize()
        self.assertEqual(t3.selected_cds_introns, set())

        t4 = Transcript()
        t4.chrom, t4.start, t4.end, t4.strand, t4.id = "Chr1", 901, 2700, "+", "multi.4"
        t4.add_exons([(901, 1500), (1601, 2100), (2201, 2700)])
        t4.add_exons([(1801, 2100), (2201, 2500)], features="CDS")
        t4.finalize()
        self.assertEqual(t4.selected_cds_introns, {(2101, 2200)}, (t4.selected_cds_introns, t4.is_coding))

        for prod in itertools.chain(itertools.combinations([t1, t2, t3], 2), [(t3, t4)]):
            self.assertTrue(Superlocus.is_intersecting(*prod, cds_only=False))
            self.assertTrue(Superlocus.is_intersecting(*prod, cds_only=True))

        for prod in [(t1, t4), (t2, t4)]:
            self.assertTrue(Superlocus.is_intersecting(*prod, cds_only=False))
            self.assertFalse(Superlocus.is_intersecting(*prod, cds_only=True))

        t5 = Transcript()
        t5.chrom, t5.start, t5.end, t5.strand, t5.id = "Chr1", 1101, 1630, "+", "mono.1"
        t5.add_exon((1101, 1630))
        t5.add_exon((1101, 1530), feature="CDS")
        t5.finalize()

        t6 = Transcript()
        t6.chrom, t6.start, t6.end, t6.strand, t6.id = "Chr1", 1101, 1530, "+", "mono.2"
        t6.add_exon((1101, 1530))
        t6.finalize()

        t7 = Transcript()
        t7.chrom, t7.start, t7.end, t7.strand, t7.id = "Chr1", 1400, 1800, "+", "mono.3"
        t7.add_exon((1400, 1800))
        t7.add_exon((1591, 1800), feature="CDS")
        t7.finalize()

        t8 = Transcript()
        t8.chrom, t8.start, t8.end, t8.strand, t8.id = "Chr1", 1801, 2110, "+", "mono.4"
        t8.add_exon((1801, 2110))
        t8.add_exon((1801, 2100), feature="CDS")
        t8.finalize()
        
        for mono in [t5, t6, t7]:
            self.assertFalse(Superlocus.is_intersecting(mono, t8, cds_only=False))
            self.assertFalse(Superlocus.is_intersecting(mono, t8, cds_only=True))

        for prod in itertools.combinations([t5, t6, t7], 2):
            self.assertTrue(Superlocus.is_intersecting(*prod, cds_only=False))
            # CDS-only is activated only when both are coding
            if prod[0].is_coding and prod[1].is_coding:
                self.assertFalse(Superlocus.is_intersecting(*prod, cds_only=True))
            else:
                self.assertTrue(Superlocus.is_intersecting(*prod, cds_only=True))

        # Now check that monoexonic and multi-exonic transcripts never intersect
        for prod in itertools.product([t1, t2, t3, t4], [t5, t6, t7, t8]):
            self.assertFalse(Superlocus.is_intersecting(*prod, cds_only=False))
            self.assertFalse(Superlocus.is_intersecting(*prod, cds_only=True))

    def test_define_graph(self):
        t1 = Transcript()
        t1.chrom, t1.start, t1.end, t1.strand, t1.id = "Chr1", 1101, 2000, "+", "multi.1"
        t1.add_exons([(1101, 1500), (1601, 2000)])
        t1.add_exons([(1201, 1500), (1601, 1900)], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id = "Chr1", 1001, 2200, "+", "multi.2"
        t2.add_exons([(1001, 1500), (1601, 2200)])
        t2.add_exons([(1201, 1500), (1601, 1900)], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.start, t3.end, t3.strand, t3.id = "Chr1", 901, 2200, "+", "multi.3"
        t3.add_exons([(901, 1500), (1601, 2200)])
        t3.finalize()
        self.assertEqual(t3.selected_cds_introns, set())

        t4 = Transcript()
        t4.chrom, t4.start, t4.end, t4.strand, t4.id = "Chr1", 901, 2700, "+", "multi.4"
        t4.add_exons([(901, 1500), (1601, 2100), (2201, 2700)])
        t4.add_exons([(1801, 2100), (2201, 2500)], features="CDS")
        t4.finalize()
        self.assertEqual(t4.selected_cds_introns, {(2101, 2200)}, (t4.selected_cds_introns, t4.is_coding))

        t5 = Transcript()
        t5.chrom, t5.start, t5.end, t5.strand, t5.id = "Chr1", 1101, 1630, "+", "mono.1"
        t5.add_exon((1101, 1630))
        t5.add_exon((1101, 1530), feature="CDS")
        t5.finalize()

        t6 = Transcript()
        t6.chrom, t6.start, t6.end, t6.strand, t6.id = "Chr1", 1101, 1530, "+", "mono.2"
        t6.add_exon((1101, 1530))
        t6.finalize()

        t7 = Transcript()
        t7.chrom, t7.start, t7.end, t7.strand, t7.id = "Chr1", 1400, 1800, "+", "mono.3"
        t7.add_exon((1400, 1800))
        t7.add_exon((1591, 1800), feature="CDS")
        t7.finalize()

        t8 = Transcript()
        t8.chrom, t8.start, t8.end, t8.strand, t8.id = "Chr1", 1801, 2110, "+", "mono.4"
        t8.add_exon((1801, 2110))
        t8.add_exon((1801, 2100), feature="CDS")
        t8.finalize()

        # This is monoexonic for the CDS *but* Multiexonic for the cDNA
        t9 = Transcript()
        t9.chrom, t9.start, t9.end, t9.strand, t9.id = "Chr1", 1121, 2000, "+", "mono_multi"
        t9.add_exons([(1121, 1500), (1601, 2000)])
        t9.add_exons([(1121, 1490)], features="CDS")
        t9.finalize()
        self.assertTrue(t9.is_coding)

        sl = Superlocus(t1, flank=2000)
        sl.configuration.pick.clustering.cds_only = False
        [sl.add_transcript_to_locus(t) for t in [t2, t3, t4, t5, t6, t7, t8, t9]]
        self.assertEqual(len(sl.transcripts), 9)
        graph = sl.define_graph()
        multi_ids = [t.id for t in [t1, t2, t3, t4, t9]]
        mono_ids = [t.id for t in [t5, t6, t7, t8]]
        self.assertTrue(all(edge in graph.edges() for edge in itertools.combinations(multi_ids, 2)))
        self.assertTrue(all(edge in graph.edges() for edge in itertools.combinations(mono_ids[:3], 2)))
        self.assertTrue(not any(edge in graph.edges() for edge in itertools.product(mono_ids[:3], [t8.id])))
        self.assertTrue(not any(edge in graph.edges() for edge in itertools.product(multi_ids,
                                                                                    mono_ids)))
        sl.configuration.pick.clustering.cds_only = True
        graph = sl.define_graph()
        self.assertFalse(all(edge in graph.edges() for edge in itertools.combinations(multi_ids, 2)))

        for comb in itertools.combinations(multi_ids, 2):
            if set.issubset(set(comb), {t1.id, t2.id, t3.id}):
                self.assertIn(comb, graph.edges())
            else:
                self.assertNotIn(comb, graph.edges())

        for comb in itertools.combinations(mono_ids, 2):
            if t5.id in comb and t6.id in comb:
                self.assertIn(comb, graph.edges())
            else:
                self.assertNotIn(comb, graph.edges())

        for comb in itertools.product(mono_ids, [t9.id]):
            if set.issubset(set(comb), {t5.id, t6.id, t9.id}):
                self.assertIn(comb, graph.edges())
            else:
                self.assertNotIn(comb, graph.edges())

        self.assertTrue(not any(edge in graph.edges() for edge in itertools.product([t.id for t in [t1, t2, t3, t4]],
                                                                                    [t9.id] + mono_ids)))


class RefGeneTester(unittest.TestCase):

    def setUp(self):

        # Chr5	26608315	26608866	ID=AT5G66658.1;coding=True;phase=0	0	+	26608328	26608553	0	1	551	0
        # Chr5	26609307	26610504	ID=AT5G66660.1;coding=True;phase=0	0	-	26609307	26610504	0	1	1197	0
        # Chr5	26611257	26612889	ID=AT5G66670.2;coding=True;phase=0	0	-	26611473	26612700	0	1	1632	0
        # Chr5	26611257	26612891	ID=AT5G66670.1;coding=True;phase=0	0	-	26611473	26612700	0	2	1470,460,1588

        self.tothers = Transcript(BED12(
            "Chr5	26608315	26608866	ID=AT5G66658.1;coding=True;phase=0	0	+	26608328	26608553	0	1	551	0"
        ))
        self.tothers.parent = "AT5G66658"
        self.tout = Transcript(BED12(
            "Chr5	26609307	26610504	ID=AT5G66660.1;coding=True;phase=0	0	-	26609307	26610504	0	1	1197	0"
        ))
        self.tout.parent = "AT5G66660"

        self.t1 = Transcript(BED12(
            "Chr5	26611257	26612889	ID=AT5G66670.2;coding=True;phase=0	0	-	26611473	26612700	0	1	1632	0"
        ))
        self.t1.parent = "AT5G66670"
        self.t2 = Transcript(BED12(
            "Chr5\t26611257\t26612891\tID=AT5G66670.1;coding=True;phase=0\t0\t-\t26611473\t26612700\t0\t2\t1470,46\t0,1588"
        ))
        self.t2.parent = self.t1.parent[:]

    def test_basic(self):
        gene = Gene(self.t1)
        self.assertEqual(gene.id, self.t1.parent[0])
        for attr in ["chrom", "source", "start", "end", "strand"]:
            self.assertEqual(getattr(gene, attr), getattr(self.t1, attr))
        gene.add(self.t2)
        self.assertIn(self.t2.id, gene)
        gene.finalize()
        self.assertEqual(gene.start, min(self.t1.start, self.t2.start))
        self.assertEqual(gene.end, max(self.t1.end, self.t2.end),
                         (self.t1.end, self.t2.end, max(self.t1.end, self.t2.end)))
        with self.assertRaises(TypeError):
            gene.logger = "foo"
        gene.logger = None
        self.assertEqual(gene.logger.name, "gene_{0}".format(gene.id))
        gene.id = None
        gene.logger = None
        self.assertEqual(gene.logger.name, "gene_generic")
        del gene.logger
        self.assertIs(gene.logger, None)

        new_gene = pickle.loads(pickle.dumps(gene))
        self.assertEqual(gene, new_gene)
        self.assertEqual(gene.transcripts, new_gene.transcripts)
        with self.assertRaises(ValueError):
            g = gene.format("foo")

    def test_less_than(self):

        g1 = Gene(self.tothers)
        g2 = Gene(self.tout)
        g3 = Gene(self.t1)
        self.assertLess(g1, g2)
        self.assertGreater(g3, g2)
        self.assertGreater(g3, g1)

    def test_deletion(self):
        gene = Gene(self.t1)
        self.assertEqual(gene.id, self.t1.parent[0])
        for attr in ["chrom", "source", "start", "end", "strand"]:
            self.assertEqual(getattr(gene, attr), getattr(self.t1, attr))
        gene.add(self.t2)
        self.assertIn(self.t2.id, gene)
        gene.finalize()
        gene.remove(self.t2.id)
        self.assertEqual((gene.start, gene.end), (self.t1.start, self.t1.end))
        gene.remove(self.t1.id)
        self.assertEqual((gene.start, gene.end), (None, None))

    def test_different_strand(self):
        gene = Gene(self.t1)
        with self.assertRaises(AssertionError):
            gene.add(self.tothers)
        self.tout.unfinalize()
        self.tout.chrom = "Chr2"
        self.tout.finalize()
        with self.assertRaises(AssertionError):
            gene.add(self.tout)

    def test_non_coding(self):
        gene = Gene(self.t1)
        self.t2.strip_cds()
        gene.add(self.t2)
        self.assertIn(self.t2.id, gene)
        gene.finalize()
        self.assertIn(self.t2.id, gene)

        gene = Gene(self.t1, only_coding=True)
        self.t2.strip_cds()
        gene.add(self.t2)
        self.assertIn(self.t2.id, gene)
        gene.finalize()
        self.assertNotIn(self.t2.id, gene)

    def test_properties(self):
        gene = Gene(self.t2)
        self.assertFalse(gene.has_monoexonic)
        gene.add(self.t1)
        self.assertTrue(gene.has_monoexonic)

        gene = Gene(self.t2)
        self.t1.unfinalize()
        exon = self.t1.exons[0]
        self.t1.remove_exon(exon)
        row = GtfLine("#")
        row.chrom, row.strand, row.start, row.end, row.feature, row.transcript, row.gene = (
            self.t1.chrom, self.t1.strand, exon[0], exon[1], "exon", self.t1.id, self.t1.parent[0]
        )
        gene.add(self.t1)
        gene.add_exon(row)
        gene.finalize()
        self.assertIn(exon, gene[self.t1.id].exons)
        self.assertEqual(gene.exons, {exon} | set(self.t2.exons))
        self.assertEqual(gene.introns, self.t2.introns)


class AbstractLocusTester(unittest.TestCase):

    @pytest.fixture(autouse=True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def setUp(self):
        gff_transcript1 = """Chr1\tfoo\ttranscript\t101\t400\t.\t+\t.\tID=t0
    Chr1\tfoo\texon\t101\t400\t.\t+\t.\tID=t0:exon1;Parent=t0
    Chr1\tfoo\tCDS\t101\t350\t.\t+\t.\tID=t0:exon1;Parent=t0""".split("\n")
        gff_transcript1 = [GFF.GffLine(x) for x in gff_transcript1]
        self.assertEqual(gff_transcript1[0].chrom, "Chr1", gff_transcript1[0])
        self.transcript1 = Transcript(gff_transcript1[0])
        for exon in gff_transcript1[1:]:
            self.transcript1.add_exon(exon)
        self.transcript1.finalize()

        gff_transcript2 = """Chr1\tfoo\ttranscript\t1001\t1400\t.\t+\t.\tID=t0
            Chr1\tfoo\texon\t1001\t1400\t.\t+\t.\tID=t0:exon1;Parent=t0
            Chr1\tfoo\tCDS\t1001\t1350\t.\t+\t.\tID=t0:exon1;Parent=t0""".split("\n")
        gff_transcript2 = [GFF.GffLine(x) for x in gff_transcript2]
        self.assertEqual(gff_transcript2[0].chrom, "Chr1", gff_transcript2[0])
        self.transcript2 = Transcript(gff_transcript2[0])
        for exon in gff_transcript2[1:]:
            self.transcript2.add_exon(exon)
        self.transcript2.finalize()

        self.assertTrue(self.transcript1.monoexonic)
        self.assertEqual(self.transcript1.chrom, gff_transcript1[0].chrom)

        self.assertTrue(self.transcript2.monoexonic)
        self.assertEqual(self.transcript2.chrom, gff_transcript2[0].chrom)
        self.configuration = configurator.load_and_validate_config(None)
        self.assertIsNotNone(self.configuration.scoring, self.configuration)
        self.transcript1.configuration = self.configuration
        self.transcript2.configuration = self.configuration
        self.assertEqual(self.transcript1.configuration.seed, self.transcript2.configuration.seed)

    def test_create_metrics_row(self):

        transcript = Transcript()
        transcript.chrom, transcript.start, transcript.end, transcript.strand = "Chr1", 101, 1000, "+"
        transcript.add_exons([(101, 500), (601, 1000)])
        transcript.add_exons([(201, 500), (601, 900)], features="CDS")
        transcript.id = "foo"
        transcript.finalize()
        transcript.is_reference = True
        for locus in [Locus, Sublocus, Superlocus]:
            for score, tpm in itertools.product([None, 10, 1.56789], [None, 10, 2.5689]):
                transcript.score = score
                transcript.external_scores["tpm"] = [tpm]
                instance = locus(transcript)
                # locus.scores[transcript.id] = score
                instance.get_metrics()
                instance.scores_calculated = True
                instance.scores[transcript.id] = {"score": score}
                self.assertTrue(instance.scores_calculated)
                self.assertEqual(transcript.score, score)
                row = instance._create_metrics_row(transcript.id, instance._metrics[transcript.id],
                                                   transcript)
                self.assertEqual(transcript.score, score)
                self.assertIsInstance(row, dict)
                self.assertIn("tid", row.keys())
                if locus == Locus:
                    self.assertEqual(row["alias"], transcript.id)
                else:
                    self.assertEqual(row["tid"], transcript.id)
                self.assertIn("external.tpm", row)
                if tpm is None:
                    self.assertEqual(row["external.tpm"], "NA")
                elif tpm == 10:
                    self.assertEqual(row["external.tpm"], 10)
                elif tpm == 2.5689:
                    self.assertEqual(row["external.tpm"], 2.57)
                if score is None:
                    self.assertEqual(row["score"], "NA", locus)
                elif score == 10:
                    self.assertEqual(row["score"], 10, locus)
                else:
                    self.assertEqual(row["score"], 1.57, locus)

    def test_removal(self):
        for cls in [Superlocus, Sublocus, Monosublocus, Locus]:
            obj = cls(self.transcript1, configuration=self.configuration)
            if cls != Monosublocus:
                obj.add_transcript_to_locus(self.transcript2, check_in_locus=False)
            self.assertEqual(obj.start, min(self.transcript2.start, self.transcript1.start))
            self.assertEqual(obj.chrom, "Chr1")
            obj._remove_all()
            self.assertEqual(obj.chrom, None)
            self.assertEqual(obj.start, float("Inf"))
            self.assertEqual(obj.end, float("-Inf"))
            self.assertEqual(len(obj.transcripts), 0)

    def test_invalid_conf(self):
        for cls in [Superlocus, Sublocus, Monosublocus, Locus]:
            for invalid in ([], ("hello",), 10, "I_do_not_exist.yaml"):
                with self.subTest(cls=cls, invalid=invalid), self.assertRaises(exceptions.InvalidConfiguration):
                    cls(self.transcript1, configuration=invalid)

    def test_not_implemented(self):
        with self.assertRaises(TypeError):
            _ = Abstractlocus(self.transcript1)

    def test_equality(self):
        for child1, child2 in combinations_with_replacement([Superlocus, Sublocus, Monosublocus, Locus], 2):
            obj1, obj2 = (child1(self.transcript1, configuration=self.configuration),
                          child2(self.transcript1, configuration=self.configuration))
            if child1 == child2:
                self.assertEqual(obj1, obj2)
            else:
                self.assertNotEqual(obj1, obj2)

    def test_less_than(self):

        for child in [Superlocus, Sublocus, Monosublocus, Locus]:
            child1, child2 = (child(self.transcript1), child(self.transcript2))
            self.assertLess(child1, child2)
            self.assertLessEqual(child1, child2)
            self.assertLessEqual(child1, child1)
            self.assertLessEqual(child2, child2)
            self.assertGreater(child2, child1)
            self.assertGreaterEqual(child2, child1)
            self.assertGreaterEqual(child2, child2)
            self.assertGreaterEqual(child1, child1)

    def test_invalid_logger(self):
        for locus in [Locus, Sublocus, Superlocus, Monosublocus, MonosublocusHolder]:
            l = locus(None)
            for invalid in [10, dict(), "hello", 5.0]:
                with self.subTest(invalid=invalid):
                    with self.assertRaises(TypeError):
                        l.logger = invalid
            l.logger = None
            self.assertIs(l.logger, create_null_logger())

    def test_invalid_source(self):
        for locus in [Locus, Sublocus, Superlocus, Monosublocus, MonosublocusHolder]:
            l = locus(None)
            for invalid in [10, dict(), b"hello", 5.0]:
                with self.subTest(invalid=invalid):
                    with self.assertRaises(AssertionError):
                        l.source = invalid
            l.source = None
            self.assertEqual(l.source, "Mikado")
            l.source = "Test"
            self.assertEqual(l.source, "Test")

    def test_invalid_stranded(self):
        for locus in [Locus, Sublocus, Superlocus, Monosublocus, MonosublocusHolder]:
            l = locus(None)
            for invalid in [10, dict(), b"hello", 5.0]:
                with self.subTest(invalid=invalid):
                    with self.assertRaises(ValueError):
                        l.stranded = invalid
            l.stranded = False
            self.assertFalse(l.stranded)
            l.stranded = True
            self.assertTrue(l.stranded)

    def test_invalid_flank(self):
        for locus in [Locus, Sublocus, Superlocus, Monosublocus, MonosublocusHolder]:
            l = locus(None)
            for invalid in [-10, dict(), b"hello", 5.0]:
                with self.subTest(invalid=invalid):
                    with self.assertRaises(TypeError):
                        l.flank = invalid
            for valid in [0, 100, 1000]:
                l.flank = valid
                self.assertEqual(l.flank, valid)

    def test_calculate_score_empty_locus(self):
        logger = create_default_logger("test_calculate_score_empty_locus", level="DEBUG")
        for locus in [Locus, Sublocus, Superlocus, Monosublocus, MonosublocusHolder]:
            l = locus(None, logger=logger)
            for metric in l.available_metrics:
                with self.assertLogs(logger, level="DEBUG") as cmo:
                    l._calculate_score(metric)
                self.assertTrue(any([re.search(r"No transcripts in.*\. Returning.", out) is not None for out in
                                     cmo.output]))

    def test_get_denominator(self):
        # _get_denominator(param: Union[MinMaxScore, TargetScore],
        #   use_raw: str, metrics: dict) -> (Union[float,int], Union[float, int, bool])
        min_score = MinMaxScore(rescaling="min", filter=None)
        metrics = {"foo.1": .5, "foo.2": .1, "foo.3": 1}
        self.assertEqual(Abstractlocus._get_denominator(min_score, use_raw=False, metrics=dict()), (1, None))
        self.assertEqual(Abstractlocus._get_denominator(min_score, use_raw=True, metrics=dict()), (1, None))
        self.assertEqual(Abstractlocus._get_denominator(min_score, use_raw=True, metrics=metrics),
                         (-1, None))
        self.assertEqual(Abstractlocus._get_denominator(min_score, use_raw=False, metrics=metrics),
                         (.9, None))
        max_score = MinMaxScore(rescaling="max", filter=None)

        self.assertEqual(Abstractlocus._get_denominator(max_score, use_raw=True, metrics=dict()),
                         (1, None))
        self.assertEqual(Abstractlocus._get_denominator(max_score, use_raw=False, metrics=dict()),
                         (1, None))
        self.assertEqual(Abstractlocus._get_denominator(max_score, use_raw=True, metrics=metrics),
                         (1, None))
        self.assertEqual(Abstractlocus._get_denominator(max_score, use_raw=False, metrics=metrics),
                         (.9, None))

        target_score = TargetScore(value=.5, rescaling="target", filter=None)
        self.assertEqual(Abstractlocus._get_denominator(target_score, use_raw=True, metrics=dict()),
                         (1, .5))
        self.assertEqual(Abstractlocus._get_denominator(target_score, use_raw=False, metrics=dict()),
                         (1, .5))
        self.assertEqual(Abstractlocus._get_denominator(target_score, use_raw=True, metrics=metrics),
                         (.5, .5))
        self.assertEqual(Abstractlocus._get_denominator(target_score, use_raw=False, metrics=metrics),
                         (.5, .5))

    def test_get_score_for_metric(self):

        # _get_score_for_metric(tid_metric, use_raw, target, denominator, param, min_val, max_val)

        for invalid in [10, b"10", dict(), 100.0]:
            with self.assertRaises(ValueError):
                Abstractlocus._get_score_for_metric(invalid, use_raw=True,
                                                    target=None, denominator=1, param="fraction",
                                                    min_val=0, max_val=1)

        max_score = MinMaxScore(rescaling="max", filter=None)
        min_score = MinMaxScore(rescaling="min", filter=None)
        target_score = TargetScore(rescaling="target", value=5, filter=None)

        for num, denominator, rescaling in itertools.product(
                arange(0, 1.05, .05), arange(.5, 10, .5), [max_score, min_score, target_score]):
            if rescaling.rescaling == "target":
                with self.assertRaises(ValueError) as exc:
                    Abstractlocus._get_score_for_metric(num, use_raw=True,
                                                        target=None, denominator=denominator, param=rescaling,
                                                        min_val=0, max_val=1)
            else:
                self.assertEqual(Abstractlocus._get_score_for_metric(num, use_raw=True,
                                            target=None, denominator=denominator, param=rescaling,
                                            min_val=0, max_val=1), round(num / denominator, 2))

        min_val = 0
        max_val = 1.05
        for num, denominator, multiplier in itertools.product(
                list(arange(0, 1.05, .05)) + [True, False], arange(.5, 10, .5),
                arange(.5, 10, .5)):
            max_score = MinMaxScore(rescaling="max", filter=None, multiplier=multiplier)
            min_score = MinMaxScore(rescaling="min", filter=None, multiplier=multiplier)
            max_result = Abstractlocus._get_score_for_metric(num, use_raw=False,
                                            target=None, denominator=denominator, param=max_score,
                                            min_val=min_val, max_val=max_val)
            self.assertEqual(max_result, round(multiplier * abs((num - min_val) / denominator), 2))
            min_result = Abstractlocus._get_score_for_metric(num, use_raw=False,
                                                             target=None, denominator=denominator, param=min_score,
                                                             min_val=min_val, max_val=max_val)
            self.assertEqual(min_result, round(multiplier * abs(1 - (num - min_val) / denominator), 2),
                             (num, denominator, multiplier))

        for num, denominator, target, multiplier in itertools.product(
                list(arange(0, 1.05, .05)) + [True, False],
                arange(.5, 10, .5), arange(.05, 1.05, .5), arange(.5, 10, .5)):
            num, denominator, target, multiplier = float(num), float(denominator), float(target), float(multiplier)
            target_score = TargetScore(rescaling="target", filter=None, multiplier=multiplier, value=target)
            target_result = Abstractlocus._get_score_for_metric(num, use_raw=False,
                                                                target=target, denominator=denominator,
                                                                param=target_score,
                                                                min_val=min_val, max_val=max_val)
            # score = 1 - abs(tid_metric - target) / denominator
            expected = round(multiplier * (1 - abs((num - target) / denominator)), 2)
            self.assertEqual(target_result, expected, (num, denominator, multiplier, target))

    def test_check_usable_raw(self):
        # _check_usable_raw(transcript, param, use_raw, rescaling, logger=create_null_logger())
        for invalid in [b"cdna_length", 10, None, dict()]:
            with self.assertRaises(TypeError) as exc:
                Abstractlocus._check_usable_raw(None, invalid, False, "max")
            self.assertIsNotNone(re.search(r"Invalid type of parameter", str(exc.exception)))
            with self.assertRaises(TypeError) as exc:
                Abstractlocus._check_usable_raw(invalid, "cdna_length", False, "max")
            self.assertIsNotNone(re.search(r"Invalid transcript type", str(exc.exception)))
        t = Transcript()
        t.chrom, t.start, t.end, t.strand, t.id = "Chr1", 101, 1000, "+", "foo"
        t.add_exons([(101, 500), (601, 1000)])
        t.add_exons([(201, 500), (601, 900)], features="CDS")
        t.finalize()
        t.external_scores.tpm = [10, False]
        t.external_scores.fraction = [.5, True]
        t.attributes["FPKM"] = 10
        logger = create_default_logger("test_check_usable_raw", level="WARNING")
        for param in t.get_available_metrics():
            for use_raw, rescaling in itertools.product([False, True], ["min", "max", "target"]):
                metr = getattr(Transcript, param)
                if not isinstance(metr, Metric):
                    continue
                usable_raw = metr.usable_raw
                if use_raw is True and (usable_raw is False or rescaling == "target"):
                    with self.assertLogs(logger, level="WARNING") as cmo:
                        found_use_raw = Abstractlocus._check_usable_raw(t, param, use_raw=use_raw, rescaling=rescaling,
                                                                  logger=logger)
                    self.assertTrue(any(re.search(r"Switching to False", out) for out in cmo.output))
                    self.assertFalse(found_use_raw)
                else:
                    found_use_raw = Abstractlocus._check_usable_raw(t, param, use_raw=use_raw, rescaling=rescaling,
                                                                    logger=logger)
                    self.assertEqual(use_raw, found_use_raw)

        # Now external things
        found_use_raw = Abstractlocus._check_usable_raw(t, "external.tpm",
                                                        use_raw=False, rescaling="target", logger=logger)
        self.assertFalse(found_use_raw)

        found_use_raw = Abstractlocus._check_usable_raw(t, "external.tpm",
                                                        use_raw=True, rescaling="max", logger=logger)
        self.assertFalse(found_use_raw)

        found_use_raw = Abstractlocus._check_usable_raw(t, "external.fraction",
                                                        use_raw=False, rescaling="target", logger=logger)
        self.assertFalse(found_use_raw)

        found_use_raw = Abstractlocus._check_usable_raw(t, "external.fraction",
                                                        use_raw=True, rescaling="max", logger=logger)
        self.assertTrue(found_use_raw)

        # Now attributes. Attributes are just returning use_raw
        found_use_raw = Abstractlocus._check_usable_raw(t, "attributes.FPKM",
                                                        use_raw=False, rescaling="max", logger=logger)
        self.assertFalse(found_use_raw)
        found_use_raw = Abstractlocus._check_usable_raw(t, "attributes.FPKM",
                                                        use_raw=True, rescaling="max", logger=logger)
        self.assertTrue(found_use_raw)

        found_use_raw = Abstractlocus._check_usable_raw(t, "attributes.FPKM",
                                                        use_raw=False, rescaling="target", logger=logger)
        self.assertFalse(found_use_raw)
        found_use_raw = Abstractlocus._check_usable_raw(t, "attributes.FPKM",
                                                        use_raw=True, rescaling="target", logger=logger)
        self.assertFalse(found_use_raw)

    def test_get_param_metrics(self):
        """"""
        
        param = "combined_cds_length"

        t1 = Transcript()
        t1.chrom, t1.start, t1.end, t1.strand, t1.id = "Chr1", 101, 1000, "+", "foo"
        # cDNA 800 bps, CDS 60 bps
        t1.add_exons([(101, 500), (601, 1000)])
        t1.add_exons([(471, 500), (601, 630)], features="CDS")
        t1.finalize()

        # cDNA 210 bps, CDS 210 bps
        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id = "Chr1", 451, 710, "+", "bar"
        t2.add_exons([(451, 500), (551, 710)])
        t2.add_exons([(451, 500), (551, 710)], features="CDS")
        t2.finalize()

        self.assertEqual(operator.attrgetter(param)(t1), t1.combined_cds_length)
        self.assertEqual(operator.attrgetter(param)(t2), t2.combined_cds_length)
        self.assertEqual(t1.combined_cds_length, 60)
        self.assertEqual(t2.combined_cds_length, 210)

        transcripts = {"foo": t1, "bar": t2}
        metrics = {t1.id: {param: t1.combined_cds_length},
                   t2.id: {param: t2.combined_cds_length}}

        param_conf = MinMaxScore(rescaling="max", filter=None)

        # First case: raise errors if the dictionary is not a dictionary ...
        for invalid in [None, 10, b"20", list()]:
            with self.assertRaises(AttributeError):
                Abstractlocus._get_param_metrics(invalid, metrics, param, param_conf)
            with self.assertRaises((AttributeError, TypeError)):
                Abstractlocus._get_param_metrics(transcripts, invalid, param, param_conf)
            with self.assertRaises((AttributeError, TypeError)):
                Abstractlocus._get_param_metrics(transcripts, metrics, invalid, param_conf)
            with self.assertRaises((AttributeError, TypeError)):
                Abstractlocus._get_param_metrics(transcripts, metrics, param, invalid)

        # Second case: no filter. Presume that the metrics come back.
        param_metrics, restored_metrics = Abstractlocus._get_param_metrics(transcripts, dict(), param, param_conf)
        self.assertEqual(restored_metrics, metrics)
        self.assertEqual({t1.id: t1.combined_cds_length, t2.id: t2.combined_cds_length}, param_metrics)
        param_metrics, _ = Abstractlocus._get_param_metrics(transcripts, metrics, param, param_conf)
        self.assertEqual({t1.id: t1.combined_cds_length, t2.id: t2.combined_cds_length},
                         param_metrics)
        # Third case. Now we have a filter but without a name, so we are presuming that it applies to
        # combined_cds_length.
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    value=150))
        param_metrics, _ = Abstractlocus._get_param_metrics(transcripts, metrics, param, param_conf)
        self.assertEqual(param_metrics, {t2.id: t2.combined_cds_length})
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    value=300))
        param_metrics, _ = Abstractlocus._get_param_metrics(transcripts, metrics, param, param_conf)
        self.assertEqual(param_metrics, {})

        # Fourth case, now we are filtering on a different parameter, cdna_length
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    metric="cdna_length",
                                                                    value=10))
        self.assertEqual(param_conf.filter.metric, "cdna_length")
        param_metrics, _ = Abstractlocus._get_param_metrics(transcripts, metrics, param, param_conf)
        self.assertEqual({t1.id: t1.combined_cds_length, t2.id: t2.combined_cds_length},
                         param_metrics)
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    metric="cdna_length",
                                                                    value=300))
        self.assertEqual(param_conf.filter.metric, "cdna_length")
        param_metrics, _ = Abstractlocus._get_param_metrics(transcripts, metrics, param, param_conf)
        self.assertEqual({t1.id: t1.combined_cds_length},
                         param_metrics)
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    metric="cdna_length",
                                                                    value=2000))
        param_metrics, _ = Abstractlocus._get_param_metrics(transcripts, metrics, param, param_conf)
        self.assertEqual({}, param_metrics)

        # Fifth case. Considering external and attributes

        t1.external_scores["tpm"] = [10, False]
        t1.external_scores["fraction"] = [.5, True]
        t1.attributes["FPKM"] = 10
        t2.external_scores["tpm"] = [5, False]
        t2.external_scores["fraction"] = [1, True]
        t2.attributes["FPKM"] = 50

        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    metric="attributes.FPKM",
                                                                    value=5))

        param_metrics, restored_metrics = Abstractlocus._get_param_metrics(transcripts, metrics,
                                                                           "external.tpm", param_conf)
        self.assertEqual(param_metrics, {t1.id: t1.external_scores.tpm[0],
                                         t2.id: t2.external_scores.tpm[0]})
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    metric="attributes.FPKM",
                                                                    value=50))

        param_metrics, restored_metrics = Abstractlocus._get_param_metrics(transcripts, metrics,
                                                                           "external.tpm", param_conf)
        self.assertEqual(param_metrics, {})
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="gt",
                                                                    metric=None,
                                                                    value=5))
        param_metrics, restored_metrics = Abstractlocus._get_param_metrics(transcripts, metrics,
                                                                           "external.tpm", param_conf)
        self.assertEqual(param_metrics, {t1.id: t1.external_scores.tpm[0]})
        param_conf = MinMaxScore(rescaling="max", filter=SizeFilter(operator="lt",
                                                                    metric="external.fraction",
                                                                    value=.6))
        param_metrics, restored_metrics = Abstractlocus._get_param_metrics(transcripts, metrics,
                                                                           "external.tpm", param_conf)
        self.assertEqual(param_metrics, {t1.id: t1.external_scores.tpm[0]})

    def test_score(self):
        t = Transcript()
        t.chrom, t.start, t.end, t.id, t.score = "Chr1", 101, 1000, "foo", None
        t.add_exons([(101, 1000)])
        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.id, t2.score = "Chr1", 101, 1000, "foo2", 100
        t2.add_exons([(101, 1000)])
        for locus in [Locus, Sublocus, Superlocus]:
            l = locus(None)
            self.assertIsNone(l.score)
            l = locus(t)
            self.assertEqual(l.score, 0)
            l.add_transcript_to_locus(t2, check_in_locus=False)
            self.assertEqual(l.score, 100)

            for invalid in [10, -5, dict(), None]:
                with self.assertRaises(ValueError):
                    l._use_transcript_scores = invalid

            for valid in [False, True]:
                l._use_transcript_scores = valid
                self.assertEqual(valid, l._use_transcript_scores)
                self.assertEqual(valid, l.scores_calculated)

    def test_serialisation(self):
        engine = create_engine("sqlite://")
        maker = sessionmaker(bind=engine)
        session = maker()
        for child in [Superlocus, Sublocus, Monosublocus, Locus, Excluded]:
            with self.subTest(child=child):
                child1 = child(self.transcript1)
                child1.session = session
                child1.sessionmaker = sessionmaker
                # Check compiled in dictionary
                self.assertIsInstance(child1.configuration, (MikadoConfiguration, DaijinConfiguration))
                obj = pickle.dumps(child1)
                nobj = pickle.loads(obj)
                self.assertEqual(child1, nobj)
                obj2 = child1.as_dict()
                self.assertIsInstance(obj2, dict)
                obj3 = child1.as_dict()
                self.assertEqual(obj2, obj3)
                nobj2 = child(None)
                nobj2.load_dict(obj2)
                self.assertEqual(child1, nobj2)
                # Now try to inject a daijin configuration
                d = DaijinConfiguration()
                obj2["json_conf"] = dataclasses.asdict(d)
                assert isinstance(obj2["json_conf"], dict)
                nobj2b = child(None)
                nobj2b.load_dict(obj2, load_configuration=True)
                self.assertIsInstance(nobj2b.configuration, DaijinConfiguration)
                # Now avoid loading the configuration altogether
                nobj2b = child(None, configuration=MikadoConfiguration())
                seed = nobj2b.configuration.seed = child1.configuration.seed + 1
                obj2["json_conf"] = dataclasses.asdict(d)
                nobj2b.load_dict(obj2, load_configuration=False)
                self.assertEqual(nobj2b.configuration.seed, seed)
                self.assertIsInstance(nobj2b.configuration, MikadoConfiguration)
                # Finally let's validate that an invalid configuration would be found out
                obj2["json_conf"] = {"I am an invalid": "dictionary"}
                with self.assertRaises(marshmallow.exceptions.MarshmallowError):
                    nobj2b = child(None)
                    nobj2b.load_dict(obj2, load_configuration=True)

    def test_slocus_dicts(self):

        self.assertEqual(self.transcript1.configuration.seed, self.transcript2.configuration.seed)
        locus = Superlocus(self.transcript1)
        locus.add_transcript_to_locus(self.transcript2, check_in_locus=False)
        locus.subloci = [Sublocus(self.transcript1)]
        l = Locus(self.transcript1)
        locus.loci = {l.id: l}
        ml = MonosublocusHolder(Monosublocus(self.transcript1))
        locus.monoholders = [ml]
        locus.excluded = Excluded(self.transcript2, configuration=locus.configuration)
        conf = locus.configuration.copy()
        _without = locus.as_dict(with_subloci=False, with_monoholders=False)
        self.assertEqual(_without["subloci"], [])
        self.assertEqual(_without["monoholders"], [])
        self.assertEqual(_without["excluded"], locus.excluded.as_dict())
        self.assertEqual(_without["loci"], {l.id: l.as_dict()})
        _with = locus.as_dict(with_subloci=True, with_monoholders=True)
        self.assertIsNotNone(_with["json_conf"]["seed"])
        self.assertEqual(_with["json_conf"]["seed"], conf.seed)
        self.assertEqual(_with["subloci"], [locus.subloci[0].as_dict()])
        self.assertEqual(_with["monoholders"], [ml.as_dict()])
        self.assertEqual(conf.seed, locus.configuration.seed)
        self.assertEqual(conf.seed, self.transcript2.configuration.seed)
        excl = Excluded(self.transcript2, configuration=conf)
        self.assertEqual(excl.configuration.seed, locus.configuration.seed)
        self.assertEqual(_with["excluded"], Excluded(self.transcript2, configuration=conf).as_dict())
        self.assertEqual(_with["loci"], {l.id: l.as_dict()})
        self.assertIsInstance(_with["json_conf"], dict)
        # Now test the reloading
        # def load_dict(self, state, print_subloci=True, print_monoloci=True, load_transcripts=True,
        #                   load_configuration=True):
        conf.threads = 10
        self.assertNotEqual(conf, locus.configuration)
        num = 0
        for print_subloci in (True, False):
            for print_monoloci in (True, False):
                for load_transcripts in (True, False):
                    for load_configuration in (True, False):
                        with self.subTest(
                            print_subloci=print_subloci,
                            print_monoloci=print_monoloci,
                            load_transcripts=load_transcripts,
                            load_configuration=load_configuration
                        ):
                            num += 1
                            test = Superlocus(None, configuration=conf)
                            self.assertIsInstance(_with["json_conf"], dict, num)
                            test.load_dict(
                                deepcopy(_with), print_subloci=print_subloci, print_monoloci=print_monoloci,
                                load_transcripts=load_transcripts, load_configuration=load_configuration)
                            if load_configuration is False:
                                self.assertEqual(conf, test.configuration)
                            else:
                                self.assertEqual(locus.configuration, test.configuration)
                            if load_transcripts is True:
                                self.assertEqual(locus.transcripts, test.transcripts)
                            else:
                                self.assertEqual(test.transcripts, dict((tid, transcript.as_dict())
                                                                        for tid, transcript in
                                                                        locus.transcripts.items()))
                            if print_subloci is True:
                                self.assertEqual(locus.excluded, test.excluded)
                                self.assertEqual(locus.subloci, test.subloci)
                            else:
                                self.assertEqual(test.subloci, [])
                            if print_monoloci is True:
                                self.assertEqual(test.monoholders, locus.monoholders)
                            else:
                                self.assertEqual(test.monoholders, [])

    def test_evaluate(self):
        """Test to verify the abstractlocus.evaluate static method"""
        # param: Union[str, int, bool, float],
        #                  conf: Union[SizeFilter, InclusionFilter, NumBoolEqualityFilter, RangeFilter]

        # Params: eq, ne, gt, ge, lt, le, in, not in, within, not within

        # eq, ne
        for val in (2, 3.0, True):
            for oval in (2, 3.0, True, 4.0, 5, False):
                for op in ["ne", "eq"]:
                    filt = NumBoolEqualityFilter.Schema().load({"value": val, "operator": op})
                    for obj in [Superlocus, Locus, Sublocus, Monosublocus, MonosublocusHolder]:
                        result = obj.evaluate(oval, filt)
                        if (val == oval and op == "eq") or (val != oval and op == "ne"):
                            self.assertTrue(result)
                        else:
                            self.assertFalse(result)
        # ge, lt, le, gt
        for val in (2, 3.0, 4):
            for oval in (2, 3, 4, 5):
                for op in ["gt", "ge", "lt", "le"]:
                    filt = SizeFilter.Schema().load({"value": val, "operator": op})
                    for obj in [Superlocus, Locus, Sublocus, Monosublocus, MonosublocusHolder]:
                        result = obj.evaluate(oval, filt)
                        if op == "lt":
                            self.assertEqual(oval < val, result, (oval, val, op))
                        elif op == "le":
                            self.assertEqual(oval <= val, result, (oval, val, op))
                        elif op == "ge":
                            self.assertEqual(oval >= val, result, (oval, val, op))
                        else:  # "gt"
                            self.assertEqual(oval > val, result, (oval, val, op))

        # within, not within
        within = RangeFilter.Schema().load({"operator": "within", "value": [20, 10]})
        without = RangeFilter.Schema().load({"operator": "not within", "value": [20, 10]})
        for obj in [Superlocus, Locus, Sublocus, Monosublocus, MonosublocusHolder]:
            self.assertTrue(obj.evaluate(15, within))
            self.assertTrue(obj.evaluate(5, without))
            self.assertFalse(obj.evaluate(15, without))
            self.assertFalse(obj.evaluate(5, within))

        # in, not in
        within = InclusionFilter.Schema().load({"operator": "in", "value": ["valid_one", "valid_two", 10]})
        without = InclusionFilter.Schema().load({"operator": "not in", "value": ["valid_one", "valid_two", 10]})
        for obj in [Superlocus, Locus, Sublocus, Monosublocus, MonosublocusHolder]:
            self.assertTrue(obj.evaluate("valid_one", within))
            self.assertTrue(obj.evaluate("valid_three", without))
            self.assertTrue(obj.evaluate(10, within))
            self.assertTrue(obj.evaluate(15, without))
            self.assertFalse(obj.evaluate("valid_one", without))
            self.assertFalse(obj.evaluate("valid_three", within))
            self.assertFalse(obj.evaluate(10, without))
            self.assertFalse(obj.evaluate(15, within))

        # Mock to demonstrate that ValueError is raised
        _mock = namedtuple("mock", ["operator", "value"])
        mock = _mock("foo", "hello")
        for obj in [Superlocus, Locus, Sublocus, Monosublocus, MonosublocusHolder]:
            with self.assertRaises(ValueError):
                obj.evaluate(10, mock)

    def test_calculate_metrics(self):
        logger = create_default_logger("test_calculate_metrics", level="DEBUG")
        stream = io.StringIO()
        logger.addHandler(logging.StreamHandler(stream))
        for locus_class in [Locus, Sublocus, Superlocus, Monosublocus]:
            locus = locus_class(self.transcript1, logger=logger)
            self.assertFalse(locus.metrics_calculated)
            with self.assertLogs(logger, level="DEBUG"):
                locus.calculate_metrics(self.transcript1.id)
            self.assertFalse(locus.metrics_calculated)
            _ = stream.getvalue()
            with self._caplog.at_level(logging.DEBUG, logger=logger.name):
                locus.get_metrics()
            self.assertGreater(len(stream.getvalue()), 0, self._caplog.text)
            self._caplog.clear()
            self.assertTrue(locus.metrics_calculated)
            stream.close()
            locus.logger.removeHandler(stream)
            stream = io.StringIO()
            locus.logger.addHandler(logging.StreamHandler(stream))
            with self._caplog.at_level(logging.DEBUG, logger=logger.name) as cmo:
                locus.get_metrics()
            logged = stream.getvalue()
            self.assertEqual(len(logged), 0, logged)
            self._caplog.clear()
        stream.close()

    def test_in_locus(self):

        for child in [Superlocus, Sublocus, Monosublocus, Locus]:
            child1 = child(self.transcript1)
            self.assertTrue(child.in_locus(child1, self.transcript1))
            self.assertTrue(Abstractlocus.in_locus(child1, self.transcript1))
            self.assertFalse(child.in_locus(child1, self.transcript2))
            self.assertTrue(child.in_locus(child1, self.transcript2, flank=abs(self.transcript1.end - self.transcript2.end)))
            self.assertFalse(Abstractlocus.in_locus(child1, self.transcript2))
            self.assertTrue(Abstractlocus.in_locus(child1, self.transcript2,
                                                   flank=abs(self.transcript1.end - self.transcript2.end)))

            with self.assertRaises(TypeError):
                child1.in_locus(child1, child1)

            with self.assertRaises(TypeError):
                Abstractlocus.in_locus(child1, child1)

        # Check that we have a suitable error
        with self.assertRaises(TypeError):
            Abstractlocus.in_locus(self.transcript1, self.transcript2)

    def test_load_scores(self):

        scores = {self.transcript1.id: 10}
        empty_scores = dict()
        false_scores = set()

        for child in [Superlocus, Sublocus, Monosublocus, Locus]:
            child1 = child(self.transcript1)
            with self.assertRaises(ValueError):
                child1._load_scores(false_scores)
            with self.assertRaises(KeyError):
                child1._load_scores(empty_scores)
            child1._load_scores(scores)
            self.assertEqual(child1.scores[self.transcript1.id], 10)
            self.assertTrue(child1.scores_calculated)
            self.assertTrue(child1.metrics_calculated)

    def test_evaluate_overlap(self):

        for child in [Superlocus, Sublocus, Monosublocus, Locus]:
            child1 = child(self.transcript1)
            self.assertFalse(child1._evaluate_transcript_overlap(self.transcript1, self.transcript2)[0])
            self.assertTrue(child1._evaluate_transcript_overlap(self.transcript1, self.transcript1)[0])

    def test_invalid_sublocus(self):

        with self.assertRaises((OSError, FileNotFoundError, exceptions.InvalidConfiguration)):
            _ = Sublocus(self.transcript1, configuration="test")

    def test_sublocus_from_sublocus(self):

        s = Sublocus(self.transcript1)
        s2 = Sublocus(s)
        self.assertFalse(s.fixed_size)
        self.assertTrue(s2.fixed_size)
        for attr in ["parent", "chrom", "start", "end", "strand", "attributes"]:
            self.assertEqual(getattr(s, attr), getattr(s2, attr))

    def test_pickling_unpickling(self):

        bed_line = "Chr5\t26585506\t26586850\tID=c58_g1_i2.mrna1.35;coding=False\t99.0\t+\t26585506\t26585507\t0\t5\t383,121,78,105,213\t0,475,710,913,1131"
        conf = MikadoConfiguration()
        conf.prepare.files.source_score = {"at": 5, "tr": -1, "pb": 1, "st": 0}
        t = Transcript(bed_line, source="tr", configuration=conf)
        t.finalize()
        for metrics in t.get_available_metrics():
            try:
                rtype = operator.attrgetter("{metric}.rtype")(Transcript)
            except AttributeError:
                continue
            if rtype == "float":
                value = random.random()
            elif rtype == "bool":
                value = random.choice([True, False])
            elif rtype == "int":
                value = random.randint(0, 1000)
            else:
                continue
            try:
                setattr(t, metrics, value)
            except AttributeError:
                continue
        for loctype in Sublocus, Monosublocus, Superlocus, Locus:
            locus = loctype(t)
            nocus = pickle.loads(pickle.dumps(locus))
            for metric in t.get_available_metrics():
                original, new = getattr(t, metric), getattr(locus[t.id], metric)
                self.assertEqual(original, new, (metric, original, new))


class LocusTester(unittest.TestCase):

    logger = create_null_logger(inspect.getframeinfo(inspect.currentframe())[2])
    logger_name = logger.name

    @classmethod
    def setUpClass(cls):

        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

    def setUp(self):

        gff_transcript1 = """Chr1\tfoo\ttranscript\t101\t400\t.\t+\t.\tID=t0
Chr1\tfoo\texon\t101\t400\t.\t+\t.\tID=t0:exon1;Parent=t0
Chr1\tfoo\tCDS\t101\t350\t.\t+\t.\tID=t0:exon1;Parent=t0""".split("\n")
        gff_transcript1 = [GFF.GffLine(x) for x in gff_transcript1]
        self.assertEqual(gff_transcript1[0].chrom, "Chr1", gff_transcript1[0])
        self.transcript1 = Transcript(gff_transcript1[0])
        for exon in gff_transcript1[1:]:
            self.transcript1.add_exon(exon)
        self.transcript1.finalize()
        self.assertTrue(self.transcript1.monoexonic)
        self.assertEqual(self.transcript1.chrom, gff_transcript1[0].chrom)

        gff_transcript2 = """Chr1\tfoo\ttranscript\t101\t600\t.\t+\t.\tID=t1
Chr1\tfoo\texon\t101\t200\t.\t+\t.\tID=t1:exon1;Parent=t1
Chr1\tfoo\texon\t301\t400\t.\t+\t.\tID=t1:exon2;Parent=t1
Chr1\tfoo\texon\t501\t600\t.\t+\t.\tID=t1:exon3;Parent=t1""".split("\n")
        gff_transcript2 = [GFF.GffLine(x) for x in gff_transcript2]
        self.transcript2 = Transcript(gff_transcript2[0], logger=self.logger)

        for exon in gff_transcript2[1:-1]:
            self.transcript2.add_exon(exon)
        # Test that a transcript cannot be finalized if
        # the exons do not define the external boundaries
        with self.assertLogs(logger=self.logger_name, level="WARNING") as _:
            self.transcript2.finalize()
        with self.assertRaises(exceptions.ModificationError):
            self.transcript2.add_exon(gff_transcript2[-1])

        self.transcript2.finalized = False
        self.transcript2.start = 101
        self.transcript2.end = 600
        self.transcript2.add_exon(gff_transcript2[-1])
        self.transcript2.finalize()
        self.assertFalse(self.transcript2.monoexonic)
        self.assertEqual(self.transcript2.exon_num, len(gff_transcript2) - 1)
        # Test that trying to modify a transcript after it has been finalized causes errors
        with self.assertRaises(exceptions.ModificationError):
            for exon in gff_transcript2[1:]:
                self.transcript2.add_exon(exon)
        # # Test that creating a superlocus without configuration fails
        # with self.assertRaises(exceptions.NoJsonConfigError):
        #     _ = Superlocus(self.transcript1)
        self.my_json = os.path.join(os.path.dirname(__file__), "configuration.yaml")

        self.my_json = configurator.load_and_validate_config(self.my_json)
        self.my_json.reference.genome = self.fai.filename.decode()
        self.assertIsInstance(self.my_json.scoring, ScoringFile, self.my_json.scoring)

    def test_locus(self):
        """Basic testing of the Locus functionality."""

        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2])
        logger.setLevel("CRITICAL")
        logger.info("Started")
        self.transcript1.logger = logger
        self.transcript2.logger = logger
        self.assertTrue(self.transcript1.monoexonic)
        slocus = Superlocus(self.transcript1,
                            configuration=self.my_json, logger=logger)
        slocus.add_transcript_to_locus(self.transcript2)
        self.assertEqual(len(slocus.transcripts), 2)
        self.assertEqual(slocus.strand, self.transcript1.strand)
        self.assertEqual(slocus.start, min(self.transcript1.start, self.transcript2.start))
        self.assertEqual(slocus.end, max(self.transcript1.end, self.transcript2.end))
        logger.info(slocus.transcripts)
        slocus.define_subloci(check_requirements=False)
        self.assertEqual(len(slocus.transcripts), 2)
        logger.info(slocus.subloci)
        logger.info(slocus.transcripts)
        self.assertEqual(len(slocus.transcripts), 2)
        self.assertEqual(len(slocus.subloci), 2)
        slocus.define_monosubloci(check_requirements=False)
        self.assertEqual(len(slocus.monosubloci), 2)
        slocus.calculate_mono_metrics()
        self.assertEqual(len(slocus.monoholders), 1)
        slocus.logger.setLevel("DEBUG")
        slocus.define_loci(check_requirements=False)
        self.assertEqual(len(slocus.loci), 1)
        # self.assertFalse(slocus["t0"].is_coding, slocus["t0"].format("gtf"))
        self.assertFalse(slocus["t1"].is_coding, slocus["t1"].format("gtf"))
        self.assertEqual(sorted(list(slocus.loci[
                                  list(slocus.loci.keys())[0]].transcripts.keys())), ["t0"])
        gff_transcript3 = """Chr1\tfoo\ttranscript\t101\t1000\t.\t-\t.\tID=tminus0
Chr1\tfoo\texon\t101\t600\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0
Chr1\tfoo\tCDS\t201\t500\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0
Chr1\tfoo\texon\t801\t1000\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0""".split("\n")
        gff_transcript3 = [GFF.GffLine(x) for x in gff_transcript3]
        transcript3 = Transcript(gff_transcript3[0])
        for exon in gff_transcript3[1:]:
                transcript3.add_exon(exon)
        transcript3.finalize()
        self.assertGreater(transcript3.combined_cds_length, 0)
        self.my_json.pick.clustering.purge = True
        logger.setLevel("WARNING")
        minusuperlocus = Superlocus(transcript3, configuration=self.my_json)
        minusuperlocus.logger = logger
        minusuperlocus.define_subloci()
        self.assertGreater(len(minusuperlocus.subloci), 0)
        minusuperlocus.calculate_mono_metrics()
        self.assertGreater(len(minusuperlocus.monoholders), 0)
        minusuperlocus.define_loci()
        self.assertEqual(len(minusuperlocus.loci), 1)
        self.assertTrue(transcript3.strand != self.transcript1.strand)

    def test_cannot_add(self):

        for strand, stranded, ostrand in itertools.product(("+", "-", None), (True, False), ("+", "-", None)):
            with self.subTest(strand=strand, stranded=stranded, ostrand=ostrand):
                t1 = "1\t100\t2000\tID=T1;coding=False\t0\t{strand}\t100\t2000\t0\t1\t1900\t0".format(
                    strand=strand if strand else ".", )
                t1 = Transcript(BED12(t1))
                t2 = "1\t105\t2300\tID=T2;coding=False\t0\t{strand}\t105\t2300\t0\t1\t2195\t0".format(
                    strand=strand if strand else ".")
                t2 = Transcript(BED12(t2))
                sl = loci.Superlocus(t1, stranded=stranded, configuration=None)
                self.assertIn(t1.id, sl)
                if not stranded or t2.strand == t1.strand:
                    sl.add_transcript_to_locus(t2)
                    self.assertIn(t2.id, sl)
                else:
                    with self.assertRaises(exceptions.NotInLocusError):
                        sl.add_transcript_to_locus(t2)

        with self.subTest():
            t1 = "1\t100\t2000\tID=T1;coding=False\t0\t+\t100\t2000\t0\t1\t1900\t0"
            t1 = Transcript(BED12(t1))
            t1.finalize()
            t2 = t1.copy()
            t2.unfinalize()
            t2.chrom = "2"
            t2.id = "T2"
            t2.finalize()
            sl = loci.Superlocus(t1, stranded=False)
            with self.assertRaises(exceptions.NotInLocusError):
                sl.add_transcript_to_locus(t2)

        st1 = "1\t100\t2000\tID=T1;coding=False\t0\t+\t100\t2000\t0\t1\t1900\t0"
        t1 = Transcript(BED12(st1))
        t1.finalize()
        t2 = BED12(st1)
        t2.start += 10000
        t2.end += 10000
        t2.thick_start += 10000
        t2.thick_end += 1000
        t2 = Transcript(t2)
        for flank in [0, 1000, 10000, 20000]:
            with self.subTest(flank=flank):
                sl = Superlocus(t1, flank=flank)
                if flank < 10000:
                    with self.assertRaises(exceptions.NotInLocusError):
                        sl.add_transcript_to_locus(t2)
                else:
                    sl.add_transcript_to_locus(t2)
                    self.assertIn(t2.id, sl)

    def test_empty_locus(self):

        st1 = "1\t100\t2000\tID=T1;coding=False\t0\t+\t100\t2000\t0\t1\t1900\t0"
        t1 = Transcript(BED12(st1))
        t1.finalize()
        sl = Superlocus(t1)
        sl.check_configuration()
        sl.remove_transcript_from_locus(t1.id)
        _ = sl.segmenttree

    def test_verified_introns(self):

        """This method will check that during run-time, the verified introns are considered at
        the Superlocus level, not at the Sublocus one."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = "1", "+", "t1"
        t1.start, t1.end = 100, 1000
        t1.add_exons([(100, 200), (300, 500), (600, 1000)])
        t1.finalize()
        t1.verified_introns.add((201, 299))
        t1.verified_introns.add((501, 599))

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = "1", "+", "t2"
        t2.start, t2.end = 100, 1000
        t2.add_exons([(100, 200), (300, 1000)])
        t2.finalize()
        t2.verified_introns.add((201, 299))

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = "1", "+", "t3"
        t3.start, t3.end = 100, 1000
        t3.add_exons([(100, 250), (300, 505), (600, 1000)])
        t3.finalize()

        jconf = configurator.load_and_validate_config(None)

        loc = Superlocus(t1, configuration=jconf)
        loc.add_transcript_to_locus(t2)
        loc.add_transcript_to_locus(t3)

        loc.define_subloci()

        self.assertEqual(t1.proportion_verified_introns, 1)
        self.assertEqual(t1.proportion_verified_introns_inlocus, 1)

        self.assertEqual(t2.proportion_verified_introns, 1)
        self.assertEqual(t2.proportion_verified_introns_inlocus, 0.5)

        self.assertEqual(t3.proportion_verified_introns, 0)
        self.assertEqual(t3.proportion_verified_introns_inlocus, 0)

    def test_boolean_requirement(self):

        logger = create_null_logger()
        logger.setLevel("DEBUG")
        logger.info("Started")

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = "1", "+", "t1"
        t1.start, t1.end = 100, 1000
        t1.add_exons([(100, 200), (300, 500), (600, 1000)])
        t1.finalize()
        t1.verified_introns.add((201, 299))
        t1.verified_introns.add((501, 599))

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = "1", "+", "t2"
        t2.start, t2.end = 100, 1000
        t2.add_exons([(100, 200), (300, 1000)])
        t2.finalize()
        t2.verified_introns.add((201, 299))

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = "1", "+", "t3"
        t3.start, t3.end = 100, 1000
        t3.add_exons([(100, 250), (300, 505), (600, 1000)])
        t3.finalize()

        jconf = configurator.load_and_validate_config(None)
        log = create_default_logger("tester", level="DEBUG")

        jconf.scoring.requirements.parameters = dict()
        jconf.scoring.requirements.expression = ["suspicious_splicing"]
        jconf.scoring.requirements.parameters["suspicious_splicing"] = NumBoolEqualityFilter.Schema().load(
            {"operator": "ne", "name": "suspicious_splicing", "value": True})

        jconf.pick.alternative_splicing.report = False
        jconf.scoring.check(minimal_orf_length=jconf.pick.orf_loading.minimal_orf_length)
        # Necessary to make sure that the externally-specified requirements are taken in
        self.assertEqual(jconf.scoring.requirements._expression, "evaluated[\"suspicious_splicing\"]")

        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2])
        for suspicious in (False, True):
            with self.subTest(suspicious=suspicious):
                logger.setLevel("WARNING")
                loc = Superlocus(t1, configuration=jconf, logger=logger)
                t2.attributes["canonical_on_reverse_strand"] = suspicious
                loc.add_transcript_to_locus(t2)
                loc.add_transcript_to_locus(t3)
                self.assertEqual(len(loc.transcripts), 3)
                # loc.logger.setLevel("DEBUG")
                loc.define_subloci()
                self.assertEqual(len(loc.transcripts), 3 if not suspicious else 2)

    def test_reducing_methods_one(self):

        t1 = Transcript()
        t1.chrom, t1.start, t1.end, t1.strand, t1.id, = "Chr5", 1000, 2000, "+", "t1"
        t1.add_exons([(1000, 1200), (1500, 2000)])
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id, = "Chr5", 999, 2001, "+", "t2"
        t2.add_exons([(999, 1200), (1500, 2001)])
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.start, t3.end, t3.strand, t3.id, = "Chr5", 999, 2002, "+", "t3"
        t3.add_exons([(999, 1200), (1500, 2002)])
        t3.finalize()

        logger = create_default_logger("test_reducing_methods_one", level="INFO")
        locus = Superlocus(t1, logger=logger)
        locus.add_transcript_to_locus(t2)
        locus.add_transcript_to_locus(t3)
        self.assertIn("t1", locus)
        self.assertIn("t2", locus)
        self.assertIn("t3", locus)
        locus.logger.setLevel(level="DEBUG")
        transcript_graph = locus.define_graph()
        self.assertIn("t1", transcript_graph.neighbors("t2"))
        locus.reduce_method_one(transcript_graph)
        self.assertNotIn("t1", locus)
        self.assertNotIn("t2", locus)
        self.assertIn("t3", locus)

    def test_reducing_methods_one_2(self):

        t1 = Transcript()
        t1.chrom, t1.start, t1.end, t1.strand, t1.id, = "Chr5", 1000, 2000, "+", "t1"
        t1.add_exons([(1000, 1200), (1500, 2000)])
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id, = "Chr5", 999, 2001, "+", "t2"
        t2.add_exons([(999, 1200), (1500, 2001)])
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.start, t3.end, t3.strand, t3.id, = "Chr5", 999, 2002, "+", "t3"
        t3.add_exons([(999, 1200), (1500, 2002)])
        t3.finalize()
        logger = create_default_logger("test_reducing_methods_one_2", level="WARNING")  # level="DEBUG")
        locus = Superlocus(t1, logger=logger)
        locus.add_transcript_to_locus(t2)
        locus.add_transcript_to_locus(t3)
        locus._complex_limit = (3, 10**4)
        locus.define_subloci(check_requirements=False)
        self.assertNotIn("t1", locus)
        self.assertNotIn("t2", locus)
        self.assertIn("t3", locus)

    def test_reducing_methods_one_3(self):

        t1 = Transcript()
        t1.chrom, t1.start, t1.end, t1.strand, t1.id, = "Chr5", 1000, 2000, "+", "t1"
        t1.add_exons([(1000, 2000)])
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id, = "Chr5", 999, 2001, "+", "t2"
        t2.add_exons([(999, 2001)])
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.start, t3.end, t3.strand, t3.id, = "Chr5", 999, 2002, "+", "t3"
        t3.add_exons([(999, 2002)])
        t3.finalize()
        logger = create_default_logger("test_reducing_methods_one_2", level="DEBUG")
        locus = Superlocus(t1, logger=logger)
        locus.add_transcript_to_locus(t2)
        locus.add_transcript_to_locus(t3)
        locus._complex_limit = (3, 10**4)
        locus.define_subloci(check_requirements=False)
        self.assertNotIn("t1", locus)
        self.assertNotIn("t2", locus)
        self.assertIn("t3", locus)

    def test_reducing_methods_one_4(self):

        t1 = Transcript(is_reference=True)
        t1.chrom, t1.start, t1.end, t1.strand, t1.id, = "Chr5", 1000, 2000, "+", "t1"
        t1.add_exons([(1000, 2000)])
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id, = "Chr5", 999, 2001, "+", "t2"
        t2.add_exons([(999, 2001)])
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.start, t3.end, t3.strand, t3.id, = "Chr5", 999, 2002, "+", "t3"
        t3.add_exons([(999, 2002)])
        t3.finalize()
        logger = create_default_logger("test_reducing_methods_one_2", level="DEBUG")
        locus = Superlocus(t1, logger=logger)
        locus.add_transcript_to_locus(t2)
        locus.add_transcript_to_locus(t3)
        locus._complex_limit = (3, 10**4)
        locus.define_subloci(check_requirements=False)
        self.assertIn("t1", locus)
        self.assertNotIn("t2", locus)
        self.assertIn("t3", locus)

    def test_reducing_methods_one_and_two(self):

        for is_ref in (True, False):
            with self.subTest(is_reference=is_ref):
                t1 = Transcript(is_reference=is_ref)
                t1.chrom, t1.start, t1.end, t1.strand, t1.id, = "Chr5", 1000, 2500, "+", "t1"
                t1.add_exons([(1000, 1500), (1800, 2000), (2200, 2500)])
                t1.finalize()

                t2 = Transcript()
                t2.chrom, t2.start, t2.end, t2.strand, t2.id, = "Chr5", 999, 2501, "+", "t2"
                t2.add_exons([(999, 1500), (1800, 2000), (2200, 2501)])
                t2.finalize()

                t3 = Transcript()
                t3.chrom, t3.start, t3.end, t3.strand, t3.id, = "Chr5", 999, 1999, "+", "t3"
                t3.add_exons([(999, 1500), (1800, 1999)])
                t3.finalize()

                logger = create_default_logger("test_reducing_methods_one_2", level="WARNING")
                locus = Superlocus(t1, logger=logger)
                locus.add_transcript_to_locus(t2)
                locus.add_transcript_to_locus(t3)
                locus._complex_limit = (1, 10**4)
                locus.logger.setLevel("DEBUG")
                locus.define_subloci(check_requirements=False)
                if is_ref:
                    self.assertIn("t1", locus)
                else:
                    self.assertNotIn("t1", locus)
                self.assertIn("t2", locus)
                self.assertNotIn("t3", locus)

    def test_reducing_methods_two(self):

        for is_ref in (True, False):
            with self.subTest(is_ref=is_ref):
                t1 = Transcript(is_reference=is_ref)
                t1.chrom, t1.start, t1.end, t1.strand, t1.id, = "Chr5", 1000, 2500, "+", "t1"
                t1.add_exons([(1000, 1500), (1800, 2000), (2200, 2500)])
                t1.finalize()

                t2 = Transcript()
                t2.chrom, t2.start, t2.end, t2.strand, t2.id, = "Chr5", 999, 2501, "+", "t2"
                t2.add_exons([(999, 1500), (1800, 2000), (2200, 2501)])
                t2.finalize()

                t3 = Transcript()
                t3.chrom, t3.start, t3.end, t3.strand, t3.id, = "Chr5", 999, 1999, "+", "t3"
                t3.add_exons([(999, 1500), (1800, 1999)])
                t3.finalize()

                logger = create_default_logger("test_reducing_methods_one_2", level="DEBUG")
                locus = Superlocus(t1, logger=logger)
                locus.add_transcript_to_locus(t2)
                locus.add_transcript_to_locus(t3)
                graph, edges = locus.reduce_method_one(locus.define_graph())
                if is_ref:
                    self.assertIn("t1", graph.nodes)
                else:
                    self.assertNotIn("t1", graph.nodes)
                self.assertIn("t2", locus)
                self.assertIn("t3", locus)
                locus = Superlocus(t1, logger=logger)
                locus.add_transcript_to_locus(t2)
                locus.add_transcript_to_locus(t3)
                _ = locus.reduce_method_two(locus.define_graph())
                if is_ref:
                    self.assertIn("t1", locus)
                else:
                    self.assertNotIn("t1", locus)
                self.assertIn("t2", locus)
                self.assertNotIn("t3", locus)

    def test_choose_best(self):
        conf = MikadoConfiguration()
        conf.seed = 10
        conf.prepare.files.labels = ["Foo", "Bar", "Baz"]
        conf.prepare.files.reference = [False, True, True]
        t1 = Transcript(source="Foo")
        t1.chrom, t1.strand, t1.start, t1.end, t1.id = "Chr1", "+", 101, 1000, "Foo.1"
        t1.add_exons([(101, 1000)])
        t1.finalize()

        t2 = Transcript(source="Bar")
        t2.chrom, t2.strand, t2.start, t2.end, t2.id = "Chr1", "+", 101, 1000, "Foo.2"
        t2.add_exons([(101, 1000)])
        t2.finalize()

        t3 = Transcript(source="Baz")
        t3.chrom, t3.strand, t3.start, t3.end, t3.id = "Chr1", "+", 101, 1000, "Foo.3"
        t3.add_exons([(101, 1000)])
        t3.source = "Baz"
        t3.finalize()

        reference_sources = {source for source, is_reference in
                             zip(conf.prepare.files.labels,
                                 conf.prepare.files.reference) if is_reference is True}
        self.assertEqual(reference_sources, {"Bar", "Baz"})
        self.assertEqual(t1.original_source, "Foo")
        self.assertEqual(t2.original_source, "Bar")
        self.assertEqual(t3.original_source, "Baz")
        self.assertTrue(t3.original_source in reference_sources and
                        t2.original_source in reference_sources)

        for locus_type in [Sublocus, Superlocus, Locus]:
            locus = locus_type(t1, configuration=conf)
            locus.configuration.pick.run_options.reference_update = False
            locus.configuration.pick.run_options.only_reference_update = False
            locus.configuration.prepare.files.reference = [False, True, True]
            locus.add_transcript_to_locus(t2, check_in_locus=False)
            locus.add_transcript_to_locus(t3, check_in_locus=False)
            locus["Foo.1"].score = 10
            locus["Foo.2"].score = 5
            locus["Foo.3"].score = 1
            self.assertEqual(locus.choose_best(locus.transcripts), "Foo.1")
            locus.configuration.pick.run_options.reference_update = True
            self.assertTrue(locus.reference_update)
            self.assertEqual(locus.choose_best(locus.transcripts), "Foo.2")
            locus.configuration.pick.run_options.reference_update = False
            locus.configuration.pick.run_options.only_reference_update = True
            self.assertEqual(locus.choose_best(locus.transcripts), "Foo.2")
            locus.configuration.prepare.files.reference = [False, False, True]
            self.assertEqual(locus.choose_best(locus.transcripts), "Foo.3")
            locus.configuration.prepare.files.reference = [False, False, False]
            self.assertEqual(locus.choose_best(locus.transcripts), "Foo.1")

    def test_find_retained_introns_monoexonic(self):
        t1 = Transcript(source="Foo")
        t1.chrom, t1.strand, t1.start, t1.end, t1.id = "Chr1", "+", 101, 1000, "Foo.1"
        t1.add_exons([(101, 1000)])
        t1.finalize()
        for locus_type in [Sublocus, Superlocus, Locus]:
            t1.retained_introns = []
            locus = locus_type(t1)
            locus.find_retained_introns(locus[t1.id])
            self.assertEqual(locus[t1.id].retained_introns, tuple([]))

    def test_skip_evaluate_overlap(self):
        t1 = Transcript(source="Foo")
        t1.chrom, t1.strand, t1.start, t1.end, t1.id = "Chr1", "+", 101, 1000, "Foo.1"
        t1.add_exons([(101, 1000)])
        t1.finalize()

        t2 = Transcript(source="Foo")
        t2.chrom, t2.strand, t2.start, t2.end, t2.id = "Chr2", "+", 101, 1000, "Foo.2"
        t2.add_exons([(101, 1000)])
        t2.finalize()

        t3 = Transcript(source="Foo")
        t3.chrom, t3.strand, t3.start, t3.end, t3.id = "Chr1", "+", 1101, 2000, "Foo.3"
        t3.add_exons([(1101, 2000)])
        t3.finalize()

        t4 = Transcript(source="Foo")
        t4.chrom, t4.strand, t4.start, t4.end, t4.id = "Chr1", "+", 999, 2000, "Foo.4"
        t4.add_exons([(999, 2000)])
        t4.finalize()

        t1.is_reference = True
        t2.is_reference = True
        t3.is_reference = True

        for locus_type in [Sublocus, Superlocus, Locus]:
            t4.is_reference = False
            for t in (t2, t3, t4):
                intersecting, reason = locus_type._evaluate_transcript_overlap(
                    t1, t, min_cdna_overlap=1, min_cds_overlap=1, check_references=False
                )
                self.assertFalse(intersecting, reason)
            t4.is_reference = True
            intersecting, reason = locus_type._evaluate_transcript_overlap(
                t1, t4, min_cdna_overlap=1, min_cds_overlap=1, check_references=False
            )
            self.assertTrue(intersecting, reason)


class ASeventsTester(unittest.TestCase):

    logger = create_null_logger("ASevents")

    def setUp(self):
        
        self.conf = configurator.load_and_validate_config(None)
        self.conf.pick.alternative_splicing.report = True
        self.conf.pick.alternative_splicing.valid_ccodes = ["j", "J", "o"]
        self.conf.pick.alternative_splicing.redundant_ccodes = ["c", "=", "_", "m"]
        self.conf.pick.alternative_splicing.only_confirmed_introns = False
        self.conf.pick.alternative_splicing.min_score_perc = 0.5
        self.conf.pick.alternative_splicing.keep_retained_introns = True
        self.conf.pick.alternative_splicing.min_cdna_overlap = 0.2
        self.conf.pick.alternative_splicing.min_cds_overlap = 0.2
        self.conf.pick.alternative_splicing.max_isoforms = 3
        self.conf.pick.alternative_splicing.cds_only = False
    
        self.t1 = Transcript()
        self.t1.chrom = "Chr1"
        self.t1.strand = "+"
        self.t1.score = 20
        self.t1.id = "G1.1"
        self.t1.parent = "G1"
        self.t1.start = 101
        self.t1.end = 1500
        
        self.t1.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1500)],
                          "exon")
        self.t1.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                          "CDS")
        self.t1.finalize()
        
        self.locus = Locus(self.t1)
        self.locus.logger = self.logger
        self.locus.configuration = self.conf

    def test_not_intersecting(self):

        # This one is contained and should be rejected
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G1.1"
        t2.parent = "G1"
        t2.start = 601
        t2.end = 1420
        t2.add_exons([(601, 700), (1001, 1300), (1401, 1420)],
                     "exon")
        t2.add_exons([(601, 700), (1001, 1300), (1401, 1420)],
                     "CDS")
        t2.finalize()

        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (False, "c"))

    def test_non_coding_primary(self):
        
        t1 = self.t1.copy()
        t1.strip_cds()
        locus = Locus(t1, configuration=self.conf, logger=self.logger)
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G1.1"
        t2.parent = "G1"
        t2.start = 101
        t2.end = 1460
        t2.add_exons([(101, 500), (601, 700), (1001, 1270), (1401, 1460)],
                          "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1270), (1401, 1440)],
                          "CDS")
        t2.finalize()
        self.assertEqual(locus.is_alternative_splicing(t2), (False, "NA", None))

    def test_valid_as(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600
        
        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()        
        self.locus.configuration.pick.alternative_splicing.cds_only = False

        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (True, "J"))

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)
        
    def test_redundant_as(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600
        
        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")

        t2.finalize()

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)        
        
        t3 = Transcript()
        t3.chrom = "Chr1"
        t3.strand = "+"
        t3.score = 20
        t3.id = "G3.1"
        t3.parent = "G3"
        t3.start = 201
        t3.end = 1630
        
        t3.add_exons([(201, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1630)],
                     "exon")
        t3.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t3.finalize()

        self.assertEqual(self.locus.is_alternative_splicing(t3)[:2], (False, "J"))
        self.locus.add_transcript_to_locus(t3)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)

    def test_non_redundant_as(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        # self.t1.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1500)],
        #                   "exon")
        # self.t1.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
        #                   "CDS")

        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        # self.locus.add_transcript_to_locus(t2)
        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (True, "J"))
        self.locus.configuration.pick.clustering.cds_only = True
        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (True, "J"))
        self.locus.configuration.pick.alternative_splicing.cds_only = True
        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (False, "="))

    def test_redundant_cds_non_redundant_cdna(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 670), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 670), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)

        t3 = Transcript()
        t3.chrom = "Chr1"
        t3.strand = "+"
        t3.score = 20
        t3.id = "G3.1"
        t3.parent = "G3"
        t3.start = 201
        t3.end = 1630

        t3.add_exons([(201, 500), (601, 670), (1001, 1300), (1401, 1460), (1601, 1630)],
                     "exon")
        t3.add_exons([(401, 500), (601, 670), (1001, 1300), (1401, 1440)],
                     "CDS")
        t3.logger = self.logger
        t3.finalize()

        self.assertEqual(self.locus.is_alternative_splicing(t3)[:2], (True, "j"))
        self.locus.configuration.pick.clustering.cds_only = True
        self.assertEqual(self.locus.is_alternative_splicing(t3)[:2], (True, "j"))
        self.locus.configuration.pick.alternative_splicing.cds_only = True
        self.assertEqual(self.locus.is_alternative_splicing(t3)[:2], (False, "j"))
        self.locus.configuration.pick.clustering.cds_only = False
        self.locus.configuration.pick.alternative_splicing.cds_only = False
        self.locus.add_transcript_to_locus(t3)
        self.assertEqual(len(self.locus.transcripts), 3, self.locus.transcripts)
        self.locus.remove_transcript_from_locus(t3.id)
        self.locus.configuration.pick.clustering.cds_only = True
        self.locus.add_transcript_to_locus(t3)
        self.assertEqual(len(self.locus.transcripts), 3, self.locus.transcripts)
        self.locus.remove_transcript_from_locus(t3.id)
        self.locus.configuration.pick.alternative_splicing.cds_only = True
        self.locus.add_transcript_to_locus(t3)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)

    def test_lowscore(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)


class MonoHolderTester(unittest.TestCase):

    logger = create_default_logger("MonoHolderTester")

    def setUp(self):

        self.conf = dict()

        self.t1 = Transcript()
        self.t1.chrom = "Chr1"
        self.t1.strand = "+"
        self.t1.score = 20
        self.t1.id = "G1.1"
        self.t1.parent = "G1"
        self.t1.start = 101
        self.t1.end = 1500

        self.t1.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1500)],
                          "exon")
        self.t1.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                          "CDS")
        self.t1.finalize()
        self.assertIs(self.t1.is_coding, True)
        self.assertTrue(hasattr(self.t1, "configuration"))

    def testCdsOverlap(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_intronMatch(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 700), (1001, 1320), (1451, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1320), (1451, 1460), (1501, 1510)],
                     "CDS")
        t2.finalize()

        self.assertTrue(self.t1.is_coding)
        self.assertTrue(t2.is_coding)

        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2, logger=self.logger))
        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2, cds_only=True, logger=self.logger))

    def test_intronOverlap(self):

        self.t1.strip_cds()
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1470
        t2.add_exons([(101, 510), (601, 700), (960, 1350), (1420, 1470)])

        t2.finalize()
        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_intron_contained_in_exon(self):

        """Here the intron is completely contained within an exon. Returns true."""

        self.t1.strip_cds()
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1250
        t2.end = 2000
        t2.add_exons([(1250, 1560), (1800, 2000)])
        t2.finalize()
        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_intron_not_contained_in_exon(self):

        self.t1.strip_cds()
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1400
        t2.end = 3000
        t2.add_exons([(1400, 1560), (2800, 3000)])
        t2.finalize()

        logger = create_default_logger("test_intron_not_contained_in_exon")

        for min_cdna_overlap in (0.01, 1):
            with self.subTest(min_cdna_overlap=min_cdna_overlap):
                self.assertIs(MonosublocusHolder.is_intersecting(
                    self.t1, t2,
                    logger=logger,
                    cds_only=False,
                    min_cdna_overlap=min_cdna_overlap,
                    min_cds_overlap=min_cdna_overlap), (min_cdna_overlap < 0.28))

    def test_noCDSOverlap(self):

        self.t1.strip_cds()
        self.assertEqual(self.t1.combined_cds_introns, set())
        self.t1.finalized = False
        self.t1.add_exons([(401, 500), (601, 700), (1001, 1100)],
                          "CDS")
        self.t1.finalize()

        t2 = Transcript()
        t2.logger = self.logger
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1470
        t2.add_exons([(101, 510), (601, 700), (960, 1350), (1421, 1470)])
        t2.add_exons([(1201, 1350), (1421, 1450)], "CDS")
        t2.finalize()

        self.assertTrue(self.t1.is_coding)
        self.assertTrue(t2.is_coding)

        self.assertGreaterEqual(0,
                                overlap(
                                    (self.t1.combined_cds_start, self.t1.combined_cds_end),
                                    (t2.combined_cds_start, t2.combined_cds_end)),
                                [(self.t1.combined_cds_start, self.t1.combined_cds_end),
                                 (t2.combined_cds_start, t2.combined_cds_end)])

        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2, logger=self.logger))
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2, cds_only=True, logger=self.logger))

    def test_only_CDS_overlap(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1250
        t2.end = 2000
        t2.add_exons([(1250, 1560), (1801, 2000)])
        t2.add_exons([(1402, 1560), (1801, 1851)], "CDS")
        t2.finalize()
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2], level="WARNING")

        for min_cds_overlap in [0.05, 0.1, 0.15, 0.2, 0.5]:
            with self.subTest(min_cds_overlap=min_cds_overlap):
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1, t2,
                                                                 cds_only=True,
                                                                 logger=logger,
                                                                 min_cds_overlap=min_cds_overlap,
                                                                 min_cdna_overlap=0.01),
                              (min_cds_overlap <= 0.19),
                              (self.t1.internal_orfs, t2.internal_orfs))

        t2.strip_cds()
        t2.finalized = False
        t2.add_exons([(1461, 1560), (1801, 1850)], "CDS")
        t2.finalize()
        self.assertGreater(len(t2.introns), 0)
        self.assertGreater(len(t2.combined_cds_introns), 0)
        # No CDS overlap this time, but cDNA overlap.
        for cds_only in (True, False):
            with self.subTest(cds_only=cds_only):
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1,
                                                            t2,
                                                            cds_only=cds_only,
                                                            logger=logger), not cds_only)

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1350
        t2.end = 3850
        t2.add_exons([(1350, 1560), (2801, 3850)])
        t2.add_exons([(1402, 1560), (2801, 3850)], "CDS")
        # logger.setLevel("DEBUG")
        t2.logger = logger
        t2.finalize()
        self.assertTrue(t2.is_coding)
        for min_overlap in [0.1, 0.2, 0.3, 0.5]:
            with self.subTest(min_overlap=min_overlap):
                cds_overlap = 0
                for frame in range(3):
                    cds_overlap += len(set.intersection(
                        self.t1.frames[frame], t2.frames[frame]
                    ))

                self.assertIs(MonosublocusHolder.is_intersecting(self.t1, t2,
                                                                 cds_only=False,
                                                                 min_cds_overlap=0.07,
                                                                 min_cdna_overlap=min_overlap,
                                                                 logger=logger), (min_overlap <= 0.12),
                              ((t2.selected_internal_orf_cds, self.t1.selected_internal_orf_cds),
                               cds_overlap))
        self.assertTrue(t2.is_coding)

        for min_overlap in [0.01, 0.05, 0.1, 0.2]:
            with self.subTest(min_overlap=min_overlap):
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1,
                                                                 t2,
                                                                 cds_only=True,
                                                                 min_cds_overlap=min_overlap,
                                                                 min_cdna_overlap=min_overlap,
                                                                 logger=logger), (min_overlap <= 0.07))

    def test_frame_compatibility(self):

        """Check that the phase method functions"""
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2])
        for phase in [0, 1, 2]:
            with self.subTest(phase=phase):
                t2 = Transcript()
                t2.chrom = "Chr1"
                t2.strand = "+"
                t2.score = 1
                t2.id = "G2.1"
                t2.parent = "G2"
                t2.start = 1350 + phase
                t2.end = 3850 + phase
                t2.add_exons([(t2.start, 1560), (2801, t2.end)])
                t2.add_exons([(1402 + phase, 1560), (2801, 3850 + phase)], "CDS")
                self.assertIs(t2.is_coding, True)
                self.assertIsInstance(self.t1, Transcript)
                self.assertIsInstance(t2, Transcript)
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1,
                                                                 t2,
                                                                 cds_only=True,
                                                                 min_cds_overlap=0.05,
                                                                 min_cdna_overlap=0.05,
                                                                 logger=logger), (phase == 0))

        self.t1.unfinalize()
        self.t1.strand = "-"
        self.t1.phases = {}  # Necessary for the correct recalculation of phases!
        self.t1.logger = logger
        self.t1.finalize()
        self.assertIs(self.t1.is_coding, True, "Something went wrong in finalising T1")
        for phase in [0, 1, 2]:
            with self.subTest(phase=phase):
                t2 = Transcript()
                t2.chrom = "Chr1"
                t2.strand = "-"
                t2.score = 1
                t2.id = "G2.1"
                t2.parent = "G2"
                t2.start = 1350 + phase
                t2.end = 3850 + phase
                t2.add_exons([(t2.start, 1560), (2801, t2.end)])
                t2.add_exons([(1402 + phase, 1560), (2801, 3850 + phase)], "CDS")
                self.assertIs(t2.is_coding, True)
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1,
                                                                 t2,
                                                                 cds_only=True,
                                                                 min_cds_overlap=0.05,
                                                                 min_cdna_overlap=0.05,
                                                                 logger=logger), (phase == 0))

    def test_no_overlap(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1600
        t2.end = 2000
        t2.add_exons([(1600, 1700), (1801, 2000)])
        t2.add_exons([(1661, 1700), (1801, 1850)], "CDS")
        t2.finalize()
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_sameness(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G1.1"
        t2.parent = "G1"
        t2.start = 1250
        t2.end = 2000
        t2.add_exons([(1250, 1560), (1801, 2000)])
        t2.add_exons([(1401, 1560), (1801, 1850)], "CDS")
        t2.finalize()
        # This fails because they have the same ID
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_holder_clustering(self):

        """This test has been added starting from the annotation of IWGSC.
        It verifies that in a complex locus we create the holders correctly."""

        chrom, strand = "chr7A", "+"
        transcripts = dict()
        transcripts["TA_PGSB_v1_dez2016_mRNA_662095"] = [
            [(711041145, 711041431), (711041641, 711042803), (711059154, 711059935)],
            36.17, "sublocus:chr3A+:711041145-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662100"] = [
            [(711056723, 711056806), (711056870, 711057549), (711057994, 711059935)],
            49.8, "sublocus:chr3A+:711040605-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662101"] = [
            [(711056723, 711057549), (711057991, 711059935)],
            48.02, "sublocus:chr3A+:711056723-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662106"] = [
            [(711056723, 711057549), (711057995, 711058007)],
            39.51, "sublocus:chr3A+:711056723-711058007.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662109"] = [
            [(711056723, 711057141), (711057213, 711057237)],
            35.85, "sublocus:chr3A+:711056723-711057237.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662111"] = [
            [(711056723, 711057549), (711057994, 711059935)],
            49.97, "sublocus:chr3A+:711040605-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662116"] = [
            [(711056723, 711057610)],
            30.7, "sublocus:chr3A+:711056723-711057610.mono"
        ]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662121"] = [
            [(711058325, 711059913), (711060068, 711060089)],
            36.48, "sublocus:chr3A+:711058325-711060089.multi"
        ]

        superlocus = None
        configuration = configurator.load_and_validate_config(None)

        for transcript in transcripts:
            tr = Transcript()
            exons, score, sublocus = transcripts[transcript]
            tr.chrom, tr.strand, tr.id = chrom, strand, transcript
            tr.add_exons(exons, features="exon")
            tr.add_exons(exons, features="CDS")
            tr.finalize()
            subl = Sublocus(tr, configuration=configuration)
            if superlocus is None:
                superlocus = Superlocus(tr, configuration=configuration)
            else:
                superlocus.add_transcript_to_locus(tr, check_in_locus=False)
            superlocus.subloci.append(subl)
            superlocus.scores[transcript] = score

        superlocus.subloci_defined = True
        self.assertEqual(len(superlocus.subloci), len(transcripts))
        superlocus.logger = create_default_logger("test_holder_clustering", level="WARNING")
        self.assertFalse(superlocus.configuration.pick.clustering.simple_overlap_for_monoexonic)
        superlocus.configuration.scoring.requirements = None
        superlocus.calculate_mono_metrics()
        self.assertEqual(len(superlocus.monoholders), 1,
                         "\n".join([", ".join(list(_.transcripts.keys())) for _ in superlocus.monoholders]))

    def test_alternative_splicing_monoexonic_not_enough_overlap(self):

        """This test verifies that while we can cluster together the transcripts at the holder stage,
        if the overlap is not enough they will fail to be recognised as valid AS events."""

        jconf = configurator.load_and_validate_config(None)

        t1, t2 = Transcript(), Transcript()
        t1.chrom, t2.chrom = "1", "1"
        t1.strand, t2.strand = "+", "+"
        t1.add_exons([(1260208, 1260482), (1262216, 1262412), (1262621, 1263851)])
        t1.add_exons([(1262291, 1262412), (1262621, 1263143)], features="CDS")
        t1.id = "cls-0-sta-combined-0_1.27.12"

        t2.add_exons([(1262486, 1264276)])
        t2.add_exons([(1263571, 1264236)], features="CDS")
        t2.id = "trn-0-sta-combined-0_1_TRINITY_GG_1373_c0_g2_i1.path1"

        self.assertTrue(MonosublocusHolder.is_intersecting(
            t1, t2, cds_only=False,
            min_cdna_overlap=jconf.pick.alternative_splicing.min_cdna_overlap,
            min_cds_overlap=jconf.pick.alternative_splicing.min_cds_overlap,
            simple_overlap_for_monoexonic=True))
        self.assertFalse(MonosublocusHolder.is_intersecting(
            t1, t2, cds_only=False,
            min_cdna_overlap=jconf.pick.clustering.min_cdna_overlap,
            min_cds_overlap=jconf.pick.clustering.min_cds_overlap,
            simple_overlap_for_monoexonic=False))

        for simple in (True, False):
            with self.subTest(simple=simple):
                jconf.pick.clustering.simple_overlap_for_monoexonic = simple
                self.logger.setLevel("DEBUG")
                slocus = Superlocus(t1, configuration=jconf, logger=self.logger)
                slocus.add_transcript_to_locus(t2)
                locus = Locus(t1, configuration=jconf)
                slocus.loci[locus.id] = locus
                slocus.define_alternative_splicing()
                self.assertEqual(len(slocus.loci[locus.id].transcripts), 1)


class TestLocus(unittest.TestCase):

    """
    This unit test is focused on the locus definition and alternative splicings.
    """

    logger = create_default_logger("tester")

    @classmethod
    def setUpClass(cls):
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

    def setUp(self):

        """Set up for the unit test."""

        # Mock dictionary to be used for the alternative splicing checks
        self.configuration = configurator.load_and_validate_config(None)
        self.configuration.pick.alternative_splicing.report = True
        self.configuration.pick.alternative_splicing.pad = True
        self.configuration.pick.alternative_splicing.max_isoforms = 3
        self.configuration.pick.alternative_splicing.keep_retained_introns = False
        self.configuration.pick.alternative_splicing.min_cds_overlap = 0
        self.configuration.pick.alternative_splicing.min_cdna_overlap = 0
        self.configuration.pick.alternative_splicing.min_score_perc = 0.1
        self.configuration.pick.alternative_splicing.valid_ccodes = ["j", "G", "g"]
        self.configuration.pick.alternative_splicing.redundant_ccodes = ["c", "=", "_", "m", "n"]
        self.configuration.pick.alternative_splicing.only_confirmed_introns = False

        # self.configuration = configurator.check_json(self.configuration)

        t1 = """Chr1\tfoo\ttranscript\t1001\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t2101\t2500\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t2101\t2500\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t2801\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t2801\t2902\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        """

        t1lines = [GtfLine(line) for line in t1.split("\n") if line]

        self.assertEqual(t1lines[0].transcript, "Chr1.1.1", t1lines[0].transcript)
        self.assertEqual(t1lines[0].parent, ["Chr1.1"], t1lines[0].parent)
        self.assertEqual(t1lines[0].gene, "Chr1.1", t1lines[0].parent)
        self.assertEqual(t1lines[0].id, "Chr1.1.1", t1lines[0].id)

        self.t1 = loci.Transcript(t1lines[0])
        for exon in t1lines[1:]:
            if exon.header:
                continue
            self.t1.add_exon(exon)
        self.t1.score = 20
        self.t1.finalize()

        # Just a fragment of the best
        t1_contained = """Chr1\tfoo\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\tCDS\t1701\t1902\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        """

        t1_contained_lines = [GtfLine(line) for line in t1_contained.split("\n") if line]
        self.t1_contained = loci.Transcript(t1_contained_lines[0])
        for exon in t1_contained_lines[1:]:
            if exon.header:
                continue
            self.t1_contained.add_exon(exon)
        self.t1_contained.score = 15
        self.t1_contained.finalize()

        # Valid AS event
        t1_as = """Chr1\tfoo\ttranscript\t1001\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t2101\t2400\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t2101\t2400\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t2801\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t2801\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        """

        t1_as_lines = [GtfLine(line) for line in t1_as.split("\n") if line]
        self.t1_as = loci.Transcript(t1_as_lines[0])
        for exon in t1_as_lines[1:]:
            if exon.header:
                continue
            self.t1_as.add_exon(exon)
        self.t1_as.score = 19
        self.t1_as.finalize()

        # Retained intron AS event
        t1_retained = """Chr1\tfoo\ttranscript\t1001\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t2101\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t2101\t2472\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        """

        t1_retained_lines = [GtfLine(line) for line in t1_retained.split("\n") if line]
        self.t1_retained = loci.Transcript(t1_retained_lines[0])
        for exon in t1_retained_lines[1:]:
            if exon.header:
                continue
            self.t1_retained.add_exon(exon)
        self.t1_retained.score = 10
        self.t1_retained.finalize()

        self.logger.setLevel(logging.WARNING)
        self.configuration.reference.genome = self.fai.filename.decode()

    def test_validity(self):
        """
        First test of validity to ensure the CCodes are as expected.
        :return:
        """

        # The fragment should have a c assigned
        result, _ = Assigner.compare(self.t1_contained, self.t1)
        self.assertEqual(result.ccode[0], "c")

        # The valid AS should have a j assigned
        result, _ = Assigner.compare(self.t1_as, self.t1)
        self.assertEqual(result.ccode[0], "j")

        # The retained intron AS should have a j assigned
        result, _ = Assigner.compare(self.t1_retained, self.t1)
        self.assertEqual(result.ccode[0], "j", result.ccode)

    def testCreate(self):

        """
        Test the creation of the locus
        :return:
        """

        locus = loci.Locus(self.t1, logger=self.logger)
        locus.configuration = self.configuration
        self.assertEqual(len(locus.transcripts), 1)

    def test_exclude_contained(self):

        """Test that we exclude a transcript with a contained class code (c)"""

        locus = loci.Locus(self.t1, logger=self.logger)
        locus.configuration = self.configuration
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_contained)
        self.assertEqual(len(locus.transcripts), 1)

    def test_add_contained(self):

        """Test that we add a transcript with a contained class code (c) if
        we explicitly ask for it"""

        locus = loci.Locus(self.t1, logger=self.logger)
        locus.configuration = self.configuration.copy()
        locus._add_to_alternative_splicing_codes("c")
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_contained)
        self.assertEqual(len(locus.transcripts), 2)

    def test_addValid(self):

        """Test that we can successfully add a transcript to the locus if
        it passes the muster."""

        locus = loci.Locus(self.t1, logger=self.logger)
        locus.configuration = self.configuration
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 2)

    def test_excludeValid(self):

        """Test that a usually valid AS is excluded when:
        - we ask for no more than one AS event
        - we exclude its ccode (j)
        - we ask for perfect (100%) CDS overlap
        """

        configuration = self.configuration.copy()
        configuration.pick.alternative_splicing.max_isoforms = 3
        configuration.pick.alternative_splicing.valid_ccodes = ["n", "O", "h"]

        locus = loci.Locus(self.t1, logger=self.logger, configuration=configuration)
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 1)

        locus._add_to_alternative_splicing_codes("j")
        locus.configuration.pick.alternative_splicing.min_cds_overlap = 1

        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 1)

    def test_exclude_opposite_strand(self):

        candidate = self.t1_as
        candidate.reverse_strand()
        logger = self.logger
        # logger.setLevel(logging.DEBUG)
        locus = loci.Locus(self.t1, logger=logger)
        locus.configuration = self.configuration.copy()
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(candidate)
        self.assertEqual(len(locus.transcripts), 1)

    def test_serialisation(self):

        """Check that the main types can be serialised correctly."""

        candidate = self.t1
        pickle.dumps(candidate)

        configuration = configurator.load_and_validate_config(None)

        for obj in Superlocus, Sublocus, Locus:
            with self.subTest(obj=obj):
                locus = obj(candidate, configuration=configuration)
                dumped = pickle.dumps(locus)
                recovered = pickle.loads(dumped)
                self.assertIsInstance(recovered, obj)
                self.assertEqual(recovered, locus)

    def test_double_orf(self):

        t = Transcript()
        t.add_exons([(101, 1000), (1101, 1200), (2001, 2900)])
        t.id = "t1"
        t.strand = "+"

        orf1 = BED12()
        orf1.transcriptomic = True
        orf1.chrom = t.id
        orf1.start = 1
        orf1.end = sum([_[1] - _[0] + 1 for _ in t.exons])
        orf1.strand = "+"
        orf1.name = "t1.orf1"
        orf1.block_sizes = (900,)
        orf1.thick_start = 1
        orf1.thick_end = 900
        orf1.block_starts = (1,)
        orf1.block_count = 1

        orf2 = BED12()
        orf2.transcriptomic = True
        orf2.strand = "+"
        orf2.chrom = t.id
        orf2.start = 1
        orf2.end = sum([_[1] - _[0] + 1 for _ in t.exons])
        orf2.name = "t1.orf2"
        orf2.block_sizes = (900,)
        orf2.thick_start = 1001
        orf2.thick_end = 1900
        orf2.block_starts = (1,)
        orf2.block_count = 1

        self.assertFalse(orf1.invalid)
        self.assertFalse(orf2.invalid)

        t.load_orfs([orf1, orf2])
        self.assertEqual(t.number_internal_orfs, 2)

        locus = Locus(t)
        locus.filter_and_calculate_scores()
        self.assertTrue(list(locus.scores.keys()), [t.id])
        rows = list(locus.print_scores())
        self.assertEqual(len(rows), 1, rows)
        self.assertEqual(rows[0]["tid"], t.id, rows[0])

    def test_remove_AS_overlapping(self):

        logger = create_null_logger(inspect.getframeinfo(inspect.currentframe())[2],
                                    level="WARNING")

        t1, t2, t1_1, t2_1 = Transcript(), Transcript(), Transcript(), Transcript()
        t1.chrom = t2.chrom = t1_1.chrom = t2_1.chrom = "1"
        t1.id, t2.id, t1_1.id, t2_1.id = "t1", "t2", "t1_1", "t2_1"
        t1.strand = t2.strand = t1_1.strand = t2_1.strand = "+"
        t1.add_exons([(101, 500), (801, 1000)])
        t1.add_exons([(101, 500), (801, 1000)], features="CDS")
        t1_1.add_exons([(101, 500), (903, 1100), (1301, 1550)])
        t1_1.add_exons([(101, 500), (903, 1100), (1301, 1550)], features="CDS")
        t2.add_exons([(1601, 1800), (1901, 2000)])
        t2.add_exons([(1601, 1800), (1901, 2000)], features="CDS")
        t2_1.add_exons([(1351, 1550), (1651, 1851), (1901, 2000)])
        t2_1.add_exons([(1351, 1550), (1651, 1851), (1901, 2000)], features="CDS")

        for tr in [t1, t2, t1_1, t2_1]:
            with self.subTest(tr=tr):
                tr.finalize()
                self.assertGreater(tr.combined_cds_length, 0, tr.id)

        conf = configurator.load_and_validate_config(None)
        conf.pick.alternative_splicing.valid_ccodes = ["j", "J", "g", "G"]
        conf.pick.alternative_splicing.only_confirmed_introns = False

        conf.as_requirements = {"_expression": "cdna_length",
                                   "expression": "evaluated['cdna_length']",
                                   "parameters": {
                                       "cdna_length": {"operator": "gt", "value": 0, "name": "cdna_length"}
                                   }}
        conf.pick.alternative_splicing.pad = False
        conf.pick.alternative_splicing.min_cds_overlap = 0.75

        with self.subTest():
            superlocus_one = Superlocus(t1, configuration=conf)
            superlocus_one.add_transcript_to_locus(t1_1)

            locus_one = Locus(t1, configuration=conf)
            locus_one.logger = logger
            superlocus_one.loci[locus_one.id] = locus_one
            superlocus_one.loci_defined = True
            with self.assertLogs(logger=logger, level="DEBUG") as cm:
                superlocus_one.logger = logger
                superlocus_one.define_alternative_splicing()
                self.assertEqual(len(superlocus_one.loci), 1)
                locus_id = [_ for _ in superlocus_one.loci.keys() if
                            t1.id in superlocus_one.loci[_].transcripts][0]
                self.assertEqual(len(superlocus_one.loci[locus_id].transcripts), 2,
                                 cm.output)

        with self.subTest():
            superlocus_two = Superlocus(t2, configuration=conf)
            superlocus_two.add_transcript_to_locus(t2_1)
            locus_two = Locus(t2, configuration=conf)
            superlocus_two.loci[locus_two.id] = locus_two
            superlocus_two.loci_defined = True
            superlocus_two.logger = logger
            superlocus_two.define_alternative_splicing()
            self.assertEqual(len(superlocus_two.loci), 1)
            locus_id = [_ for _ in superlocus_two.loci.keys() if
                        t2.id in superlocus_two.loci[_].transcripts][0]
            self.assertEqual(len(superlocus_two.loci[locus_id].transcripts), 2)
        with self.subTest():
            superlocus = Superlocus(t1, configuration=conf)
            superlocus.add_transcript_to_locus(t2, check_in_locus=False)
            superlocus.add_transcript_to_locus(t1_1)
            superlocus.add_transcript_to_locus(t2_1)
            locus_one = Locus(t1_1, configuration=conf, logger=logger)
            locus_two = Locus(t2, configuration=conf, logger=logger)
            superlocus.loci[locus_one.id] = locus_one
            superlocus.loci[locus_two.id] = locus_two
            self.assertEqual(len(superlocus.loci[locus_one.id].transcripts), 1)
            self.assertEqual(len(superlocus.loci[locus_two.id].transcripts), 1)
            superlocus.loci_defined = True
            with self.assertLogs(logger=logger, level="DEBUG") as cm:
                self.assertEqual(len(superlocus.loci), 2)
                superlocus.define_alternative_splicing()
                locus_one_id = [_ for _ in superlocus.loci.keys() if
                                t1_1.id in superlocus.loci[_].transcripts][0]
                locus_two_id = [_ for _ in superlocus.loci.keys() if
                                t2.id in superlocus.loci[_].transcripts][0]
                self.assertNotEqual(locus_one_id, locus_two_id)
                self.assertEqual(len(superlocus.loci), 2)
                self.assertEqual(len(superlocus.loci[locus_two_id].transcripts), 1,
                                 (cm.output, superlocus.loci[locus_one_id].transcripts.keys()))
                self.assertEqual(len(superlocus.loci[locus_one_id].transcripts), 1,
                                 (cm.output, superlocus.loci[locus_one_id].transcripts.keys()))


class EmptySuperlocus(unittest.TestCase):

    def test_empty(self):

        logger = create_null_logger()
        logger.setLevel("WARNING")
        with self.assertLogs(logger=logger, level="WARNING"):
            _ = Superlocus(transcript_instance=None)


class WrongSplitting(unittest.TestCase):

    def test_split(self):

        t1 = Transcript(BED12("Chr1\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"))
        t2 = Transcript(BED12("Chr1\t100\t1000\tID=t2;coding=False\t0\t-\t100\t1000\t0\t1\t900\t0"))
        sl = Superlocus(t1, stranded=False)
        sl.add_transcript_to_locus(t2)
        splitted = list(sl.split_strands())
        self.assertEqual(len(splitted), 2)
        self.assertIsInstance(splitted[0], Superlocus)
        self.assertIsInstance(splitted[1], Superlocus)
        self.assertTrue(splitted[0].stranded)
        self.assertTrue(splitted[1].stranded)

    def test_invalid_split(self):
        t1 = Transcript(BED12("Chr1\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"))
        t2 = Transcript(BED12("Chr1\t100\t1000\tID=t2;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"))

        logger = create_default_logger("test_invalid_split", level="WARNING")
        with self.assertLogs(logger=logger, level="WARNING") as cm:
            sl = Superlocus(t1, stranded=True, logger=logger)
            sl.add_transcript_to_locus(t2)
            splitted = list(sl.split_strands())

        self.assertEqual(splitted[0], sl)
        self.assertEqual(len(splitted), 1)
        self.assertIn("WARNING:test_invalid_split:Trying to split by strand a stranded Locus, {}!".format(sl.id),
                      cm.output, cm.output)


class WrongLoadingAndIntersecting(unittest.TestCase):

    def test_wrong_intersecting(self):
        t1 = Transcript(BED12("Chr1\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"))
        sl = Superlocus(t1, stranded=True)

        with self.subTest():
            self.assertFalse(sl.is_intersecting(t1, t1))
        t2 = Transcript(BED12("Chr1\t100\t1000\tID=t1;coding=False\t0\t-\t100\t1000\t0\t1\t900\t0"))

        with self.subTest():
            self.assertTrue(sl.is_intersecting(t1, t2))

    def test_coding_intersecting(self):
        t1 = Transcript(BED12("Chr1\t100\t1000\tID=t1;coding=True\t0\t+\t200\t500\t0\t1\t900\t0"))
        sl = Superlocus(t1, stranded=True)
        t2 = Transcript(BED12("Chr1\t100\t1000\tID=t2;coding=True\t0\t+\t600\t900\t0\t1\t900\t0"))
        t3 = Transcript(BED12("Chr1\t100\t1000\tID=t3;coding=True\t0\t+\t300\t600\t0\t1\t900\t0"))
        t1.finalize()
        t2.finalize()
        t3.finalize()
        self.assertTrue(t1.is_coding)
        self.assertTrue(t2.is_coding)
        self.assertTrue(t3.is_coding)
        self.assertNotEqual(t1, t2)
        self.assertNotEqual(t1, t3)

        with self.subTest():
            self.assertTrue(sl.is_intersecting(t1, t2, cds_only=False))
            self.assertFalse(sl.is_intersecting(t1, t2, cds_only=True))
        with self.subTest():
            self.assertTrue(sl.is_intersecting(t1, t3, cds_only=False))
            self.assertTrue(sl.is_intersecting(t1, t3, cds_only=True))


class RetainedIntronTester(unittest.TestCase):

    def setUp(self):
        self.my_json = os.path.join(os.path.dirname(__file__), "configuration.yaml")
        self.my_json = configurator.load_and_validate_config(self.my_json)

    def test_5utr_vs_3utr(self):

        t1 = Transcript()
        t1.chrom, t1.id = "Chr1", "t1"
        t1.add_exons([(101, 300), (501, 700), (901, 1000), (1301, 1500), (1801, 2000)])
        #                 100         200         100          200         60
        t1.add_exons([(201, 300), (501, 700), (901, 1000), (1301, 1500), (1801, 1860)], features="CDS")

        t2 = Transcript()
        t2.chrom, t2.id = "Chr1", "t2"
        t2.add_exons([(101, 700), (901, 1000), (1301, 1500), (1801, 2000)])
        #                 40         200         60
        t2.add_exons([(961, 1000), (1301, 1500), (1801, 1860)], features="CDS")

        for strand, cds_disr in zip(["+", "-"], [True, False]):
            with self.subTest(strand=strand, cds_disr=cds_disr):
                t1.unfinalize()
                t2.unfinalize()
                t1.strand, t2.strand = strand, strand
                t1.finalize()
                t2.finalize()
                s = Superlocus(t1)
                s.add_transcript_to_locus(t2, check_in_locus=False)
                s.find_retained_introns(s["t1"])
                s.find_retained_introns(s["t2"])
                self.assertTrue(s["t2"].retained_introns, [(101, 700)])
                self.assertEqual(s["t2"].cds_disrupted_by_ri, cds_disr)

    def test_ri_in_monoexonic(self):

        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                logger = create_default_logger("test_ri_in_monoexonic_{strand}".format(**locals()),
                                               level="INFO")
                t = Transcript()
                t.chrom, t.strand, t.id, t.parent = "Chr1", strand, "t1", "gene"
                t.add_exons([(101, 200), (301, 600), (901, 1200), (1501, 1800)])
                t.add_exons([(171, 200), (301, 600), (901, 1200), (1501, 1530)], features="CDS")

                t2 = Transcript()
                t2.chrom, t2.strand, t2.id, t2.parent = "Chr1", strand, "t2", "gene"
                t2.add_exons([(101, 650)])
                t2.add_exons([(251, 580)], features="CDS")

                t3 = Transcript()
                t3.chrom, t3.strand, t3.id, t3.parent = "Chr1", strand, "t3", "gene"
                t3.add_exons([(151, 470)])
                t3.add_exons([(251, 470)], features="CDS")

                t.finalize()
                t2.finalize()
                t3.finalize()
                sl = Superlocus(t, logger=logger)
                sl.add_transcript_to_locus(t2)
                sl.add_transcript_to_locus(t3)
                sl.logger.setLevel("DEBUG")
                sl.find_retained_introns(t2)
                sl.find_retained_introns(t3)
                self.assertEqual(t2.retained_introns, (t2.exons[0],))
                self.assertEqual(t3.retained_introns, (t3.exons[0],))
                self.assertTrue(t2.cds_disrupted_by_ri)
                self.assertTrue(t3.cds_disrupted_by_ri)

    def test_cds_disr_five_utr(self):

        lines = """Chr1	100	2682	ID=test_0;coding=True;phase=0	0	+	497	2474	0	7	208,201,41,164,106,170,715	0,351,780,1075,1439,1616,1867
Chr1	100	2682	ID=test_1;coding=True;phase=0	0	+	497	2474	0	6	552,41,164,106,170,715	0,780,1075,1439,1616,1867
Chr1	100	2682	ID=test_2;coding=True;phase=0	0	+	1228	2474	0	5	821,164,106,170,715	0,1075,1439,1616,1867
Chr1	100	2682	ID=test_3;coding=True;phase=0	0	+	497	2474	0	7	234,201,41,164,106,170,715	0,351,780,1075,1439,1616,1867
"""
        beds = [BED12(line) for line in lines.split("\n")]
        transcripts = dict((bed.id, Transcript(bed)) for bed in beds if not bed.header)
        sl = Superlocus(transcripts["test_0"])
        for tid in transcripts:
            sl.add_transcript_to_locus(transcripts[tid], check_in_locus=False)

        sl.filter_and_calculate_scores()
        for tid in sl:
            if tid == "test_2":
                self.assertTrue(sl[tid].cds_disrupted_by_ri)
                self.assertEqual(sl[tid].retained_introns, ((101, 921),))
            else:
                self.assertFalse(sl[tid].cds_disrupted_by_ri)

    def test_real_retained_pos(self):

        """Here we verify that a real retained intron is called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1180)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        for pred, retained in [(t2, True), (t3, False)]:
            with self.subTest(pred=pred, retained=retained):
                sup = Superlocus(t1, configuration=self.my_json)
                sup.add_transcript_to_locus(pred)
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained, (pred.id, retained))

    def test_false_ri(self):

        t1 = [[(11, 100), (301, 400), (451, 500)],
              [(71, 100), (301, 400), (451, 470)]]
        t2 = [[(11, 150), (301, 400), (451, 500)],
              [(121, 150), (301, 400), (451, 470)]]

        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                tr1 = Transcript()
                tr1.add_exons(t1[0])
                tr1.add_exons(t1[0], features="CDS")
                tr1.strand = strand
                tr1.id = "t1"
                tr1.finalize()
                self.assertTrue(tr1.is_coding)
                tr2 = Transcript()
                tr2.add_exons(t2[0])
                tr2.add_exons(t2[0], features="CDS")
                tr2.strand = strand
                tr2.id = "t2"
                self.assertTrue(tr2.is_coding)
                logger = create_default_logger("test_false_ri", level="DEBUG")
                sup = Superlocus(tr1, logger=logger)
                sup.add_transcript_to_locus(tr2)
                sup.filter_and_calculate_scores(check_requirements=False)
                self.assertFalse(sup["t1"].cds_disrupted_by_ri)
                self.assertEqual(sup["t1"].retained_intron_num, 0)
                self.assertFalse(sup["t2"].cds_disrupted_by_ri)
                self.assertEqual(sup["t2"].retained_intron_num, 0)

    def test_false_positive_retained_introns(self):

        bed1 = "\t".join(str(_) for _ in
                         ["Chr1", 33677, 37762,
                          "ID=mikado.Chr1G2.1;coding=True;phase=0;alias=trinity_c0_g1_i1.mrna1.2",
                          19.0, "-", 33991, 37061, 0, 9,
                          "650,1074,81,234,62,112,181,26,194",	"0,723,1889,2052,2946,3132,3345,3695,3891"])
        bed2 = "\t".join(str(_) for _ in
                         ["Chr1", 33677, 37762,
                          "ID=mikado.Chr1G2.2;coding=True;phase=0;alias=stringtie_Stringtie_TH.27.2",
                          19.0, "-", 33991, 37061, 0, 9,
                          "650,1071,81,234,62,112,181,26,194", "0,723,1889,2052,2946,3132,3345,3695,3891"])
        bed3 = "\t".join(str(_) for _ in
                         ["Chr1", 33677, 37762,
                          "ID=mikado.Chr1G2.3;coding=True;phase=0;alias=stringtie_Stringtie_TH.27.5",
                          11.0, "-", 34833, 37061, 0, 10,
                          "650,455,545,81,234,62,112,181,26,194", "0,723,1252,1889,2052,2946,3132,3345,3695,3891"])
        bed4 = "\t".join(str(_) for _ in
                         ["Chr1", 33677, 37762,
                          "ID=mikado.Chr1G1.2;coding=True;phase=0;alias=stringtie_Stringtie_TH.27.6",
                          15.15,
                          "-", 34833, 37061, 0, 11,
                          "181,347,455,545,81,234,62,112,181,26,194",
                          "0,303,723,1252,1889,2052,2946,3132,3345,3695,3891"])

        beds = [BED12(line) for line in [bed1, bed2, bed3, bed4]]
        transcripts = [Transcript(bed) for bed in beds]
        [transcript.finalize() for transcript in transcripts]

        sstart = 0
        for transcript in transcripts:
            print(transcript.id, end="\t")
            print(transcript.start - sstart, transcript.end - sstart, transcript.combined_cds_start - sstart,
                  transcript.combined_cds_end - sstart, end="\t", sep="\t")
            [print((_[0] - sstart, _[1] - sstart), end="\t") for _ in sorted(transcript.introns)]
            print("\n")

        self.assertTrue([transcript.is_coding for transcript in transcripts])

        sup = Superlocus(transcripts[0], stranded=True)
        for transcript in transcripts[1:]:
            sup.add_transcript_to_locus(transcript, check_in_locus=False)

        [sup.find_retained_introns(transcript) for transcript in sup.transcripts.values()]
        for transcript in sup.transcripts.values():
            self.assertFalse(transcript.cds_disrupted_by_ri, (transcript.id,
                                                              transcript.retained_introns,
                                                              transcript.cds_disrupted_by_ri))

    def test_retained_pos_truncated(self):
        """Here we verify that a real retained intron is called as such,
        even when the transcript is truncated."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1420)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 1420)

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        logger = create_default_logger("test_retained_pos_truncated")
        for pred, retained in [(t2, True), (t3, False)]:
            with self.subTest(pred=pred, retained=retained):
                logger.setLevel("WARNING")
                sup = Superlocus(t1, configuration=self.my_json, logger=logger)
                sup.add_transcript_to_locus(pred)
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained, (pred.id, retained, pred.retained_introns))
                # Now check that things function also after unpickling
                unpickled_t1 = pickle.loads(pickle.dumps(t1))
                unpickled_other = pickle.loads(pickle.dumps(pred))
                logger.setLevel("WARNING")
                sup = Superlocus(unpickled_t1, configuration=self.my_json, logger=logger)
                sup.add_transcript_to_locus(unpickled_other)
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained)

    def test_real_retained_pos_truncated(self):
        """Here we verify that a real retained intron is *NOT* called as such when
        the transcript is truncated and we elect not to investigate the 3' end."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1420)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 1420)

        logger = create_default_logger("test_real_retained_pos_truncated")
        sup = Superlocus(t1, configuration=self.my_json, logger=logger)
        sup.add_transcript_to_locus(t2)

        sup.logger.setLevel("DEBUG")
        sup.find_retained_introns(t2)
        self.assertEqual(sup.transcripts["t2"].retained_introns, ((1201, 1420),))
        self.assertTrue(sup.transcripts["t2"].cds_disrupted_by_ri)

    def test_real_retained_neg_truncated(self):
        """Here we verify that a real retained intron is called as such,
        even when the transcript is truncated."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(601, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 601)

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "-", "t3"
        t3.add_exons([(551, 580), (801, 1000), (1201, 1300), (1501, 1800)])
        t3.add_exons([(551, 580),
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t3.finalize()
        self.assertEqual(t3.combined_cds_end, 551)

        for pred, retained in [(t2, True), (t3, False)]:
            with self.subTest(pred=pred, retained=retained):
                sup = Superlocus(t1, configuration=self.my_json)
                sup.add_transcript_to_locus(pred)
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained, (pred.id, pred.retained_introns))

    def test_real_retained_neg_truncated_2(self):
        """Here we verify that a real retained intron is *NOT* called as such when
        the transcript is truncated and we elect not to investigate the 3' end."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(601, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 601)

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)
        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_introns, ((601, 1000),))
        self.assertTrue(sup.transcripts["t2"].cds_disrupted_by_ri)

    def test_real_retained_pos_noCDS(self):
        """Here we verify that a real retained intron is called as such, even when the transcript lacks a CDS"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.finalize()

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.logger = create_default_logger("test_real_retained_pos_noCDS", level="DEBUG")
        # sup.logger.setLevel("DEBUG")
        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_introns, ((1201, 1600),))

    def test_not_retained_pos(self):

        """Here we verify that a false retained intron is not called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1540),  # 340
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        for pred in [t2, t3]:
            with self.subTest(pred=pred):
                sup = Superlocus(t1, configuration=self.my_json)
                sup.add_transcript_to_locus(pred)
                sup.find_retained_introns(pred)
                if pred == t2:
                    self.assertEqual(sup.transcripts[pred.id].retained_intron_num, 1)
                else:
                    self.assertEqual(sup.transcripts[pred.id].retained_intron_num, 0)
                self.assertFalse(sup.transcripts[pred.id].cds_disrupted_by_ri)
                unpickled_t1 = pickle.loads(pickle.dumps(t1))
                unpickled_other = pickle.loads(pickle.dumps(pred))
                sup = Superlocus(unpickled_t1, configuration=self.my_json)
                sup.add_transcript_to_locus(unpickled_other)
                sup.find_retained_introns(unpickled_other)
                if pred == t2:
                    self.assertEqual(sup.transcripts[unpickled_other.id].retained_intron_num, 1)
                else:
                    self.assertEqual(sup.transcripts[unpickled_other.id].retained_intron_num, 0)
                self.assertFalse(sup.transcripts[unpickled_other.id].cds_disrupted_by_ri)

    def test_neg_retained_example(self):

        t1 = Transcript()
        t1.chrom = "Chr1"
        t1.id = "t1"
        t1.add_exons([(3168512, 3168869),(3168954, 3169234),(3169327, 3169471),
                      (3169589, 3170045),(3170575, 3170687),(3170753, 3170803)])
        t1.strand = "-"
        t1.add_exons(
            [(3168568, 3168869), (3168954, 3169234), (3169327, 3169471), (3169589, 3170045), (3170575, 3170682)],
            features="CDS"
        )
        t1.finalize()

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.id = "t2"
        t2.strand = "-"
        t2.add_exons(
            [(3168483, 3168869),(3168954, 3169234),(3169327, 3169471),(3169589, 3170816)]
        )
        t2.add_exons(
            [(3168568, 3168869),(3168954, 3169234),(3169327, 3169471),(3169589, 3170192)],
        features="CDS")

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)
        sup.find_retained_introns(t2)
        self.assertGreater(sup.transcripts[t2.id].retained_intron_num, 0)

    def test_real_retained_neg(self):
        """Here we verify that a real retained intron is called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(401, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (771, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        with self.subTest():
            sup = Superlocus(t1, configuration=self.my_json)
            sup.add_transcript_to_locus(t2)
            sup.logger = create_default_logger("test_real_retained_neg",
                                               level="DEBUG")
            sup.find_retained_introns(t2)
            self.assertEqual(sup.transcripts["t2"].retained_introns, ((401, 1000),))

        with self.subTest():
            unpickled_t1 = pickle.loads(pickle.dumps(t1))
            unpickled_other = pickle.loads(pickle.dumps(t2))
            sup = Superlocus(unpickled_t1, configuration=self.my_json)
            sup.add_transcript_to_locus(unpickled_other)
            sup.find_retained_introns(unpickled_other)
            self.assertEqual(sup.transcripts["t2"].retained_introns, ((401, 1000),))

    def test_real_retained_neg(self):

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (771, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(401, 1000), (1201, 1300), (1501, 1800)])
        t3.add_exons([(831, 1000),  # 200
                      (1201, 1300),
                      (1501, 1530)
                      ], features="CDS")
        t3.finalize()

        graph = Abstractlocus._calculate_graph([t1, t2, t3])
        exons = set.union(*[set(_.exons) for _ in [t1, t2, t3]])
        introns = set.union(*[_.introns for _ in [t1, t2, t3]])

        segmenttree = Abstractlocus._calculate_segment_tree(exons, introns)
        logger=create_default_logger("test_real_retained_neg", level="WARNING")
        # logger.setLevel("DEBUG")
        self.assertTrue(
            Abstractlocus._is_exon_retained((401, 1000),
                                            "-",
                                            segmenttree,
                                            graph,
                                            [Interval(401, 830)],
                                            introns=set.union(
                                                t1.introns, t2.introns, t3.introns
                                            ),
                                            cds_introns=set.union(
                                                t1.combined_cds_introns,
                                                t2.combined_cds_introns,
                                                t3.combined_cds_introns
                                            ),
                                            internal_splices={1000},
                                            logger=logger))

        for alt, num_retained, cds_disrupted in zip([t2, t3], [1, 1], [True, True]):
            unpickled_t1 = pickle.loads(pickle.dumps(t1))
            unpickled_alt = pickle.loads(pickle.dumps(alt))
            with self.subTest(alt=alt):
                sup = Superlocus(t1, configuration=self.my_json, logger=logger)
                logger.setLevel("DEBUG")
                sup.find_retained_introns(alt)
                self.assertEqual(alt.retained_intron_num, num_retained,
                                 (alt.id, alt.retained_introns))
                self.assertEqual(alt.cds_disrupted_by_ri, cds_disrupted,
                                 (alt.id, alt.retained_introns))
                logger.setLevel("WARNING")
            with self.subTest(alt=alt):
                sup = Superlocus(unpickled_t1, configuration=self.my_json)
                sup.find_retained_introns(unpickled_alt)
                self.assertEqual(unpickled_alt.retained_intron_num, num_retained,
                                 unpickled_alt.retained_introns)

    def test_consider_cds_only(self):

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(831, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "-", "t3"
        t3.add_exons([(401, 1000), (1201, 1300), (1501, 1800)])
        t3.add_exons([(831, 1000),  # 200
                      (1201, 1300),
                      (1501, 1530)
                      ], features="CDS")
        t3.finalize()

        graph = Abstractlocus._calculate_graph([t1, t3])
        exons = set.union(*[set(_.combined_cds) for _ in [t1, t3]])
        introns = set.union(*[_.combined_cds_introns for _ in [t1, t3]])

        logger = create_default_logger("test_consider_cds_only", level="DEBUG")
        logger.debug("Exons: %s", exons)
        logger.debug("Introns: %s", introns)

        segmenttree = Abstractlocus._calculate_segment_tree(exons, introns)
        logger.setLevel("DEBUG")
        self.assertEqual((False, False),
                         Abstractlocus._is_exon_retained(
                             (401, 1000), "-", segmenttree, graph, [Interval(401, 830)],
                             introns=set.union(t1.introns, t3.introns),
                             cds_introns=set.union(t1.combined_cds_introns, t3.combined_cds_introns),
                             internal_splices={1000},
                             logger=logger))

    def test_issue_255(self):

        beds = ['Calendula_officinalis_EIV1.1_Contig_0000010\t81356\t84036\tID=CoDisc.hq_lq_UnnamedSample_HQ_transcript/78575.mrna1;coding=True;phase=0\t20.82\t-\t81605\t83927\t0\t12\t315,15,211,44,178,185,62,88,93,96,44,179\t0,391,518,818,960,1216,1486,1640,1819,2169,2365,2501',
            'Calendula_officinalis_EIV1.1_Contig_0000010\t81356\t84036\tID=CoDisc_scallop_CoDisc_sca.200.0.3;coding=True;phase=1\t7.99\t-\t82161\t84036\t0\t9\t862,178,185,62,88,93,96,44,179\t0,960,1216,1486,1640,1819,2169,2365,2501',
            'Calendula_officinalis_EIV1.1_Contig_0000010\t81356\t84036\tID=CoDisc_scallop_CoDisc_sca.200.0.0;coding=True;phase=0\t15.07\t-\t81605\t83927\t0\t11\t315,211,44,178,185,62,88,93,96,44,179\t0,518,818,960,1216,1486,1640,1819,2169,2365,2501',
            'Calendula_officinalis_EIV1.1_Contig_0000010\t81356\t84036\tID=CoRay_scallop_CoRay_sca.180.0.2;coding=True;phase=0\t16.010000228881836\t-\t81714\t83927\t0\t11\t406,211,44,178,185,62,88,93,96,44,179\t0,518,818,960,1216,1486,1640,1819,2169,2365,2501',
            'Calendula_officinalis_EIV1.1_Contig_0000010\t81356\t84036\tID=CoDisc_scallop_CoDisc_sca.200.0.4;coding=True;phase=0\t9.46\t-\t81838\t83927\t0\t10\t729,44,178,185,62,88,93,96,44,179\t0,818,960,1216,1486,1640,1819,2169,2365,2501']

        beds = [BED12(bed) for bed in beds]
        assert all([bed.header is False for bed in beds])
        transcripts = dict((bed.id, Transcript(bed)) for bed in beds)
        [transcripts[tid].finalize() for tid in transcripts]
        self.assertEqual(len(transcripts), 5)
        conf = configurator.load_and_validate_config(None)
        conf.pick.alternative_splicing.keep_cds_disrupted_by_ri = False
        conf.pick.alternative_splicing.pad = False
        logger = create_default_logger("test_issue_255", level="WARNING")

        scores = dict((tid, transcripts[tid].score) for tid in transcripts)
        scores = list(sorted(scores.items(), key=operator.itemgetter(1), reverse=True))
        self.assertEqual(len(scores), 5)

        locus = Locus(transcripts[scores[0][0]], configuration=conf, logger=logger)
        for tid, score in scores[1:]:
            locus.add_transcript_to_locus(transcripts[tid], check_in_locus=False)

        self.assertEqual(len(locus.transcripts), len(beds), scores)
        locus.logger.setLevel("DEBUG")
        locus.find_retained_introns(locus["CoDisc_scallop_CoDisc_sca.200.0.4"])
        self.assertTrue(locus["CoDisc_scallop_CoDisc_sca.200.0.4"].retained_intron_num > 0)
        self.assertTrue(locus["CoDisc_scallop_CoDisc_sca.200.0.4"].cds_disrupted_by_ri)
        locus.find_retained_introns(locus["CoDisc_scallop_CoDisc_sca.200.0.3"])
        self.assertTrue(locus["CoDisc_scallop_CoDisc_sca.200.0.3"].retained_intron_num > 0)
        self.assertTrue(locus["CoDisc_scallop_CoDisc_sca.200.0.3"].cds_disrupted_by_ri)

    def test_not_retained_neg(self):
        """Here we verify that a false retained intron is not called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(301, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (471, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)

        self.assertEqual(t2.cds_tree.find(301, 1000),
                         [Interval(471, 1000)])

        self.assertEqual(Abstractlocus._exon_to_be_considered((301, 1000), t2),
                         (True, [(301, 470)], {1000}),
                         Abstractlocus._exon_to_be_considered((301, 1000), t2))

        graph = Abstractlocus._calculate_graph([t1, t2])

        segmenttree = Abstractlocus._calculate_segment_tree(set.union(set(t1.exons), set(t2.exons)),
                                                            set.union(t1.introns, t2.introns))

        self.assertEqual(
            (True, False),
            Abstractlocus._is_exon_retained((301, 1000),
                                            "-",
                                            segmenttree, graph, [(301, 470)],
                                            introns=set.union(t1.introns, t2.introns),
                                            internal_splices={1000},
                                            cds_introns=set.union(t1.combined_cds_introns, t2.combined_cds_introns),
                                            logger=create_default_logger("test_not_retained_neg",
                                                                         level="DEBUG")))

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 1,
                         sup.transcripts["t2"].retained_introns)
        self.assertEqual(sup.transcripts["t2"].retained_introns, ((301, 1000),),
                         sup.transcripts["t2"].retained_introns)
        self.assertFalse(sup.transcripts["t2"].cds_disrupted_by_ri)

    def test_exon_switching_pos(self):

        """Checking that an exon switching is treated correctly as a NON-retained intron. Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)

    def test_exon_switching_pos_noCDS(self):
        """Checking that an exon switching is treated correctly as a NON-retained intron even when the CDS is absent.
        Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        # t2.add_exons([(201, 500),  # 300
        #               (801, 1000),  # 200
        #               (1201, 1300),  # 100
        #               (1501, 1530)  # 30
        #               ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0,
                         sup.transcripts["t2"].retained_introns)

    def test_exon_switching_neg(self):
        """Checking that an exon switching is treated correctly as a NON-retained intron. Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (2201, 2300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (2201, 2300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(101, 500), (1701, 2000), (2201, 2300), (2501, 2800)])
        t2.add_exons([
                      (1801, 2000),  # 200
                      (2201, 2300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(len(t2.cds_tree), len(t2.combined_cds) + len(t2.combined_cds_introns))
        self.assertEqual(len(t2.cds_tree), 5)

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)

    def test_exon_switching_neg_noCDS(self):
        """Checking that an exon switching is treated correctly as a NON-retained intron even when the CDS is absent.
        Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        # t2.add_exons([(201, 500),  # 300
        #               (801, 1000),  # 200
        #               (1201, 1300),  # 100
        #               (1501, 1530)  # 30
        #               ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, configuration=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)

    def test_neg_delayed_cds(self):

        t1 = Transcript()
        t1.chrom = "Chr1"
        t1.start, t1.end, t1.strand, t1.id = 47498, 49247, "-", "cls-0-hst-combined-0_Chr1.7.0"
        t1.add_exons([(47498, 47982), (48075, 48852), (48936, 49247)])
        t1.add_exons([(47705, 47982), (48075, 48852), (48936, 49166)], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.start, t2.end, t2.strand, t2.id = 47485, 49285, "-", "scl-1-hst-combined-0_gene.13.0.1"
        t2.add_exons([(47485, 47982), (48075, 49285)])
        t2.add_exons([(47705, 47982), (48075, 48813)], features="CDS")
        t2.finalize()

        logger = create_default_logger("test_neg_delayed_cds", level="WARNING")
        sup = Superlocus(t1, configuration=self.my_json, logger=logger)
        sup.add_transcript_to_locus(t2)
        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts[t2.id].retained_intron_num, 1, sup.combined_cds_exons)

    def test_mixed_strands(self):

        """Verify that no retained intron is called if the strands are mixed."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (771, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, configuration=self.my_json, stranded=False)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)

    def test_extensive_cases(self):
        """This method implements the case described by D. Swarbreck in issue #255"""

        basis = {
            "t1": [[(11, 100), (301, 400), (451, 600), (751, 900), (1041, 1200), (1501, 1600)],
                   [(551, 600), (751, 900), (1041, 1080)]],
            # Difference: partially within 5'UTR intron. NOT RI
            "t2": [[(11, 180), (301, 400), (451, 600), (751, 900), (1041, 1200), (1501, 1600)],
                   [(551, 600), (751, 900), (1041, 1080)]],
            # Difference: RI in the UTR. NOT CDS disrupted. Same CDS start and end
            "t3": [[(11, 400), (451, 600), (751, 900), (1041, 1200), (1501, 1600)],
                   [(551, 600), (751, 900), (1041, 1080)]],
            # Difference: RI in the CDS. NOT CDS disrupted. Different CDS start and end
            "t4": [[(11, 100), (301, 400), (451, 600), (751, 1200), (1501, 1600)],
                   [(551, 600), (751, 1200), (1540, 1600)]],
            # Difference: RI in the first CDS exon, leading to premature stop
            "t5": [[(11, 100), (301, 400), (451, 680)],
                   [(551, 640)]],
            # Difference: Start of the transcript within first intron. Same CDS
            "t6": [[(201, 400), (451, 600), (751, 900), (1041, 1200), (1501, 1600)],
                   [(551, 600), (751, 900), (1041, 1080)]],
            # Difference: novel exon within first intron
            "t7": [[(201, 250), (301, 400), (451, 600), (751, 900), (1041, 1200), (1501, 1600)],
                   [(551, 600), (751, 900), (1041, 1080)]],
            # Difference: same UTR, last exon ending within a CDS intron. CDS disrupted
            "t8": [[(11, 100), (301, 400), (451, 600), (751, 950)],
                   [(551, 600), (751, 940)]],
            # Difference: AS in the second-to-last exon
            "t9": [[(11, 100), (301, 400), (451, 600), (751, 980), (1501, 1600)],
                   [(551, 600), (751, 970)]],
            # Difference: AS in the last exon
            "t10": [[(11, 100), (301, 400), (451, 600), (751, 900), (1041, 1200), (1401, 1600)],
                    [(551, 600), (751, 900), (1041, 1080)]]
        }

        res = {"t1": ((), False), "t2": ((), False), "t3": (((11, 400),), False), "t4": ((), False),
               "t5": (((451, 680),), True), "t6": (((201, 400),), False), "t7": ((), False), "t8": (((751, 950),), True),
               "t9": ((), False), "t10": ((), False)}

        logger = create_default_logger("test_extensive_cases", "WARNING")
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                transcripts = dict()
                for tid in basis:
                    transcripts[tid] = Transcript()
                    transcripts[tid].chrom, transcripts[tid].strand = "Chr1", strand
                    transcripts[tid].add_exons(basis[tid][0])
                    transcripts[tid].add_exons(basis[tid][1], features="CDS")
                    transcripts[tid].id = tid
                    transcripts[tid].finalize()
                    self.assertTrue(transcripts[tid].is_coding, tid)
                sup = Superlocus(transcripts["t1"])
                sup.logger = logger
                sup.configuration.pick.clustering.purge = False
                self.assertIn("t1", sup.transcripts.keys())
                for tid in transcripts:
                    if tid == "t1":
                        continue
                    elif tid == "t5" and strand == "-":
                        sup.logger.setLevel("DEBUG")
                    else:
                        logger.setLevel("WARNING")
                    sup.add_transcript_to_locus(transcripts[tid], check_in_locus=True)
                    self.assertIn(tid, sup.transcripts.keys())

                sup.filter_and_calculate_scores()
                sup.logger.setLevel("WARNING")
                for tid in sorted(res.keys()):
                    self.assertIn(tid, sup.transcripts.keys())
                    err_message = [strand,
                                   tid,
                                   res[tid],
                                   (sup.transcripts[tid].retained_introns, sup.transcripts[tid].cds_disrupted_by_ri)]
                    self.assertEqual(
                        (sup.transcripts[tid].retained_introns, sup.transcripts[tid].cds_disrupted_by_ri),
                        res[tid],
                        err_message)

    def test_mixed_intron(self):
        """This test verifies that a retained intron with one boundary from one transcript and another from a second
        is correctly marked."""

        basis = {
            "t1": [[(101, 300), (701, 900), (1101, 1300), (1801, 2000)],
                   [(271, 300), (701, 900), (1101, 1300), (1801, 1850)]],  # 30 + 200 + 200 + 60 = 480
            "t2": [[(401, 500), (661, 900), (1101, 1500), (1800, 2000)],
                   [(451, 500), (661, 900), (1101, 1410)]],  # 50 + 240 + 320 = 290 + 310 = 600
            "t3": [[(401, 500), (661, 1300), (1800, 2000)],
                   [(451, 500), (661, 1030)]]  # 50 + 380 = 420
        }

        res = {"t1": ((), False), "t2": ((), False), "t3": (((661, 1300),), True)}
        logger = create_default_logger("test_mixed_intron", "WARNING")
        for strand in ("+", "-"):
            transcripts = dict()
            for tid in basis:
                transcripts[tid] = Transcript()
                transcripts[tid].chrom, transcripts[tid].strand = "Chr1", strand
                transcripts[tid].add_exons(basis[tid][0])
                transcripts[tid].add_exons(basis[tid][1], features="CDS")
                transcripts[tid].id = tid
                transcripts[tid].finalize()
                self.assertTrue(transcripts[tid].is_coding, tid)
            sup = Superlocus(transcripts["t1"])
            self.assertIn("t1", sup.transcripts.keys())
            [sup.add_transcript_to_locus(transcripts[tid]) for tid in transcripts if tid != "t1"]
            sup.logger = logger
            sup.filter_and_calculate_scores()
            for tid in sorted(res.keys()):
                self.assertEqual((sup.transcripts[tid].retained_introns, sup.transcripts[tid].cds_disrupted_by_ri),
                                 res[tid],
                                 [tid, res[tid],
                                  (sup.transcripts[tid].retained_introns, sup.transcripts[tid].cds_disrupted_by_ri)])


class PicklingTest(unittest.TestCase):

    def setUp(self):
        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1180)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        self.t1, self.t2, self.t3 = t1, t2, t3
        self.configuration = configurator.load_and_validate_config(None)

    def test_transcript_pickling(self):

        for transcript in [self.t1, self.t2, self.t3]:
            with self.subTest(transcript=transcript):
                pickled = pickle.dumps(transcript)
                unpickled = pickle.loads(pickled)
                self.assertEqual(transcript, unpickled)
                self.assertEqual(len(transcript.combined_cds) + len(transcript.combined_cds_introns),
                                 len(unpickled.cds_tree))
                self.assertEqual(len(transcript.segmenttree), len(unpickled.segmenttree))

    def test_locus_unpickling(self):

        for transcript in [self.t1, self.t2, self.t3]:
            for (loc_type, loc_name) in [(_, _.__name__) for _ in (Superlocus, Sublocus, Monosublocus, Locus)]:
                with self.subTest(transcript=transcript, loc_type=loc_type, loc_name=loc_name):
                    loc = loc_type(transcript, configuration=self.configuration)
                    pickled = pickle.dumps(transcript)
                    unpickled = pickle.loads(pickled)
                    self.assertEqual(transcript, unpickled)


class PaddingTester(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

    @staticmethod
    def load_from_bed(manager, resource):
        transcripts = dict()
        with pkg_resources.resource_stream(manager, resource) as bed:
            for line in bed:
                line = line.decode()
                line = BED12(line, coding=True)
                line.coding = True
                transcript = Transcript(line)
                assert transcript.start > 0
                assert transcript.end > 0
                assert transcript.is_coding, transcript.format("bed12")
                transcript.finalize()
                transcript.verified_introns = transcript.introns
                transcript.parent = "{}.gene".format(transcript.id)
                transcripts[transcript.id] = transcript
        return transcripts

    def test_pad_utr(self):
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2], level="WARNING")
        transcripts = self.load_from_bed("Mikado.tests", "pad_utr.bed12")
        locus = Locus(transcripts["mikado.Chr5G2.1"], logger=logger)
        # [locus._add_to_alternative_splicing_codes(code) for code in ("J", "j", "G", "h")]
        self.assertIn("J", locus.configuration.pick.alternative_splicing.valid_ccodes)
        locus.configuration.reference.genome = self.fai.filename.decode()
        # We need to pad
        locus.configuration.pick.alternative_splicing.pad = True
        locus.configuration.pick.alternative_splicing.ts_distance = 10000
        locus.configuration.pick.alternative_splicing.ts_max_splices = 10
        locus.configuration.pick.alternative_splicing.ts_distance = 1000
        locus.configuration.pick.alternative_splicing.ts_max_splices = 10
        locus.configuration.pick.alternative_splicing.only_confirmed_introns = False
        locus.configuration.pick.alternative_splicing.min_cdna_overlap = 0.1

        locus.logger = create_default_logger("test_pad_utr", level="DEBUG")
        with self.assertLogs("test_pad_utr", "DEBUG") as cmo:
            locus.add_transcript_to_locus(transcripts["mikado.Chr5G2.2"])

        self.assertIn("mikado.Chr5G2.1", locus.transcripts, cmo.output)
        self.assertIn("mikado.Chr5G2.2", locus.transcripts, cmo.output)
        locus.pad_transcripts()
        for tid in locus:
            self.assertEqual(locus[tid].end, locus.end, tid)

    def test_remove_redundant(self):
        """Test to verify that the routine to remove redundant transcripts after padding functions as desired"""
        tmult = Transcript()
        tmult.chrom, tmult.start, tmult.end, tmult.strand, tmult.id = "Chr1", 101, 2000, "+", "foo.1"
        tmult.add_exons([(101, 500), (801, 1200), (1501, 2000)])
        tmult.add_exons([(201, 500), (801, 1200), (1501, 1700)], features="CDS")
        tmult.finalize()
        tmult_a = tmult.copy()
        tmult_a.id = "foo.1a"

        tmono = Transcript()
        tmono.chrom, tmono.start, tmono.end, tmono.strand, tmono.id = "Chr1", 601, 1800, "+", "bar.1"
        tmono.add_exons([(601, 1800)])
        # Note that the CDS *must be in-frame* with tmult
        tmono.add_exons([(711, 1580)], features="CDS")
        tmono.finalize()
        tmono_a = tmono.copy()
        tmono_a.id = "bar.1a"

        conf = MikadoConfiguration()
        comp = Assigner.compare(tmono, tmult)
        rev_comp = Assigner.compare(tmult, tmono)

        # Make sure that the transcripts can coexist
        conf.pick.alternative_splicing.valid_ccodes.extend([comp[0].ccode[0], rev_comp[0].ccode[0]])
        conf.pick.alternative_splicing.min_cds_overlap = 0.1
        conf.pick.alternative_splicing.min_cdna_overlap = 0.1
        conf.pick.alternative_splicing.min_score_perc = 0.01
        conf.seed = 10

        # Now for the real tests
        logger = create_default_logger("test_remove_redundant", level="DEBUG")
        for primary, is_reference, score in itertools.product([tmult.id, tmono.id], [False, True], [5, 10, 1]):
            with self.subTest():
                if primary == tmult.id:
                    locus = Locus(tmult, configuration=conf)
                    secondary = tmono.id
                    alt_secondary = tmono_a.id
                    locus.add_transcript_to_locus(tmono, check_in_locus=False)
                else:
                    locus = Locus(tmono, configuration=conf)
                    secondary = tmult.id
                    alt_secondary = tmult_a.id
                    locus.add_transcript_to_locus(tmult, check_in_locus=False)

                locus.primary_transcript_id = primary
                random.seed(conf.seed)
                locus.add_transcript_to_locus(tmono_a, check_in_locus=False)
                locus.add_transcript_to_locus(tmult_a, check_in_locus=False)
                self.assertIn(tmult.id, locus, (primary, is_reference, score))
                self.assertIn(tmono.id, locus, (primary, is_reference, score))
                self.assertIn(tmult_a.id, locus, (primary, is_reference, score))
                self.assertIn(tmono_a.id, locus, (primary, is_reference, score))

                locus[tmono.id].score = score
                locus[tmono_a.id].score = 5
                locus[tmult.id].score = score
                locus[tmult_a.id].score = 5
                locus[tmult.id].is_reference = locus[tmono.id].is_reference = is_reference

                locus.logger = logger
                locus._remove_redundant_after_padding()
                self.assertTrue(len(locus.transcripts), 2)
                if is_reference is False:
                    if score == 5:
                        self.assertEqual(locus.primary_transcript_id, primary)
                        self.assertTrue(secondary in locus or alt_secondary in locus)
                    elif score == 10:
                        self.assertEqual(sorted(locus.transcripts.keys()), sorted([tmult.id, tmono.id]),
                                         (primary, is_reference, score))
                    elif score == 1:
                        self.assertEqual(sorted(locus.transcripts.keys()), sorted([alt_secondary, primary]),
                                         (primary, is_reference, score))
                else:
                    self.assertEqual(sorted(locus.transcripts.keys()), sorted([primary, secondary]),
                                     (primary, is_reference, score))
                # Check that calling again does nothing
                current = list(sorted(locus.transcripts.keys()))[:]
                locus._remove_redundant_after_padding()
                self.assertTrue(current == list(sorted(locus.transcripts.keys())))

        # Now check that a single transcript in a locus causes the function to return immediately
        locus = Locus(tmult, configuration=conf, logger=logger)
        with self.assertLogs(logger.name, level="DEBUG") as cmo:
            locus._remove_redundant_after_padding()

        self.assertTrue(any([re.search(r"only has one transcript, no redundancy removal needed.", output) is not None
                             for output in cmo.output]))

    def test_ad_three_prime(self):
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2], level="WARNING")
        transcripts = self.load_from_bed("Mikado.tests", "pad_three_neg.bed12")
        locus = Locus(transcripts["mikado.Chr5G486.1"], logger=logger)
        [locus._add_to_alternative_splicing_codes(code) for code in ("J", "j", "G", "h")]
        self.assertIn("J", locus.configuration.pick.alternative_splicing.valid_ccodes)
        locus.configuration.reference.genome = self.fai.filename.decode()
        locus.add_transcript_to_locus(transcripts["mikado.Chr5G486.2"])
        # We need to pad
        locus.configuration.pick.alternative_splicing.pad = True
        locus.configuration.pick.alternative_splicing.ts_distance = 10000
        locus.configuration.pick.alternative_splicing.ts_max_splices = 10
        locus.configuration.pick.alternative_splicing.ts_distance = 1000
        locus.configuration.pick.alternative_splicing.ts_max_splices = 10
        locus.configuration.pick.alternative_splicing.only_confirmed_introns = False
        locus.configuration.pick.alternative_splicing.min_cdna_overlap = 0.1
        self.assertIn("mikado.Chr5G486.1", locus.transcripts.keys())
        self.assertIn("mikado.Chr5G486.2", locus.transcripts.keys())
        locus.logger.setLevel("DEBUG")
        locus.pad_transcripts()
        for tid in locus:
            self.assertEqual(locus[tid].start, locus.start, tid)

    def test_one_off(self):
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2], level="WARNING")
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                logger.setLevel("WARNING")
                t1 = Transcript()
                t1.chrom, t1.strand, t1.start, t1.end, t1.id, t1.parent = ["Chr5", strand, 100, 1000, "t1", "loc"]
                t1.add_exons([(100, 200), (300, 500), (700, 800), (900, 1000)])
                t1.finalize()
                loc = Locus(t1, logger=logger)
                loc.configuration.reference.genome = self.fai.filename.decode()
                # We need these to be padded
                loc.configuration.pick.alternative_splicing.ts_distance = 1000
                loc.configuration.pick.alternative_splicing.ts_max_splices = 10
                loc.configuration.pick.alternative_splicing.only_confirmed_introns = False
                loc.configuration.pick.alternative_splicing.min_cdna_overlap = 0.1
                t2 = Transcript()
                t2.chrom, t2.strand, t2.start, t2.end, t2.id, t2.parent = ["Chr5", strand, 299, 1000, "t2", "loc"]
                t2.add_exons([(299, 400), (700, 800), (900, 1000)])
                t2.finalize()
                loc.add_transcript_to_locus(t2)
                self.assertIn(t2.id, loc)
                t3 = Transcript()
                t3.chrom, t3.strand, t3.start, t3.end, t3.id, t3.parent = ["Chr5", strand, 100, 801, "t3", "loc"]
                t3.add_exons([(100, 150), (350, 500), (700, 801)])
                t3.finalize()
                loc.add_transcript_to_locus(t3)
                self.assertIn(t3.id, loc)
                t4 = Transcript()
                t4.chrom, t4.strand, t4.start, t4.end, t4.id, t4.parent = ["Chr5", strand, 300, 1000, "t4", "loc"]
                t4.add_exons([(300, 320), (600, 800), (900, 1000)])
                t4.finalize()
                self.assertGreaterEqual(t4.cdna_length, 300)
                loc.add_transcript_to_locus(t4)
                self.assertIn(t4.id, loc)
                t5 = Transcript()
                t5.chrom, t5.strand, t5.start, t5.end, t5.id, t5.parent = ["Chr5", strand, 100, 800, "t5", "loc"]
                t5.add_exons([(100, 140), (360, 650), (700, 800)])
                t5.finalize()
                loc.add_transcript_to_locus(t5)
                loc._load_scores({"t1": 20, "t2": 10, "t3": 10, "t4": 15, "t5": 15})
                self.assertIn(t5.id, loc)
                self.assertIn(t1.id, loc)
                self.assertEqual(loc["t1"].score, 20)
                loc.pad_transcripts()
                self.assertEqual(loc["t4"].start, 100)
                self.assertEqual(loc["t5"].end, 1000)
                self.assertEqual(loc["t2"].start, 299)
                self.assertEqual(loc["t3"].end, 801)

    @mark.slow
    def test_complete_padding(self):

        transcripts = self.load_from_bed("Mikado.tests", "complete_padding.bed12")
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2], level="WARNING")
        for idx in range(1, 5):
            self.assertIn('AT5G01030.{}'.format(idx), transcripts.keys(), transcripts.keys())
        cds_coordinates = dict()
        genome = pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz")
        # Distance of .4 to .1 or .2: 600
        # Distance of .3 to .1: 270
        # Distance of .3 to .2: 384

        for pad_distance, max_splice, coding, best in itertools.product(
                (300, 400, 601),
                (0, 1, 2, 3),
                (False, True,),
                tuple(["AT5G01030.1", "AT5G01030.2"])):
            with self.subTest(pad_distance=pad_distance, max_splice=max_splice, coding=coding, best=best):
                primary = transcripts[best].copy()
                if coding is False:
                    primary.strip_cds()
                locus = loci.Locus(primary)
                locus.configuration.reference.genome = genome
                for t in transcripts:
                    if t == locus.primary_transcript_id:
                        continue
                    trans = transcripts[t].copy()
                    if coding is False:
                        trans.strip_cds()
                    locus.add_transcript_to_locus(trans)

                # Now add the scores
                scores = {best: 15}
                for tid in transcripts.keys():
                    if tid not in ("AT5G01030.1", "AT5G01030.2"):
                        scores[tid] = 9
                    else:
                        scores[tid] = 10

                locus._load_scores(scores=scores)
                locus._load_scores(scores=scores)

                if coding:
                    cds_coordinates = dict()
                    for transcript in locus:
                        cds_coordinates[transcript] = (
                            locus[transcript].combined_cds_start, locus[transcript].combined_cds_end)

                locus.logger = logger
                locus.configuration.pick.alternative_splicing.ts_distance = pad_distance
                locus.configuration.pick.alternative_splicing.ts_max_splices = max_splice
                locus.logger.setLevel("DEBUG")
                locus.pad_transcripts()
                locus.logger.setLevel("WARNING")
                self.assertEqual(transcripts["AT5G01030.2"].start, 9869)
                self.assertEqual(locus[best].start, 9869)
                self.assertIn(best, locus)
                if max_splice < 2 or pad_distance <= 250:
                    with self.assertLogs(logger, "DEBUG") as cm:
                        locus.logger.setLevel("DEBUG")
                        share = locus._share_five_prime(transcripts["AT5G01030.3"],
                                                        transcripts["AT5G01030.1"])
                        self.assertEqual(False, share, cm.output)

                    self.assertEqual(locus["AT5G01030.3"].start, transcripts["AT5G01030.3"].start,
                                     (pad_distance, max_splice, coding, best))
                    self.assertNotIn("padded", locus["AT5G01030.3"].attributes,
                                     (pad_distance, max_splice, coding, best))
                else:
                    self.assertEqual(locus["AT5G01030.3"].start, transcripts["AT5G01030.2"].start,
                                     (locus["AT5G01030.3"].start, pad_distance, max_splice, coding, best))
                    self.assertTrue(locus["AT5G01030.3"].attributes.get("padded", False),
                                    (pad_distance, max_splice, coding, best))
                    self.assertEqual(locus["AT5G01030.3"].exons,
                                     [(9869, 10172), (10574, 12665), (12803, 13235)])

                    if coding:
                        self.assertEqual(locus["AT5G01030.3"].combined_cds_start,
                                         transcripts["AT5G01030.2"].combined_cds_start)
                # self.assertFalse(locus[best].attributes.get("padded", False))
                if max_splice < 2 or pad_distance < 600:
                    self.assertEqual(locus["AT5G01030.4"].end, transcripts["AT5G01030.4"].end)
                    self.assertNotIn("padded", locus["AT5G01030.4"].attributes,
                                     (pad_distance, max_splice, coding, best))
                else:
                    self.assertEqual(locus["AT5G01030.4"].end, transcripts["AT5G01030.2"].end,
                                     (pad_distance, max_splice, coding, best))
                    self.assertTrue(locus["AT5G01030.4"].attributes.get("padded", False))
                    if coding:
                        self.assertTrue(locus["AT5G01030.4"].is_coding,
                                        locus["AT5G01030.4"])
                        self.assertTrue(transcripts["AT5G01030.2"].is_coding)
                        self.assertTrue(locus["AT5G01030.2"].is_coding)
                        self.assertEqual(locus["AT5G01030.4"].combined_cds_end,
                                         # transcripts["AT5G01030.2"].combined_cds_end,
                                         locus["AT5G01030.2"].combined_cds_end,
                                         "\n".join([locus["AT5G01030.2"].format("bed12"),
                                                    locus["AT5G01030.4"].format("bed12")]))

    @mark.slow
    def test_negative_padding(self):
        genome = pkg_resources.resource_filename("Mikado.tests", "neg_pad.fa")
        transcripts = self.load_from_bed("Mikado.tests", "neg_pad.bed12")
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2], level="WARNING")
        self.assertIn('Human_coding_ENSP00000371111.2.m1', transcripts.keys(), transcripts.keys())
        locus = loci.Locus(transcripts['Human_coding_ENSP00000371111.2.m1'], logger=logger)
        locus.configuration.reference.genome = genome
        for t in transcripts:
            if t == locus.primary_transcript_id:
                continue
            locus.add_transcript_to_locus(transcripts[t])

        self.assertEqual(transcripts['Human_coding_ENSP00000371111.2.m1'].combined_cds_end, 1646)
        self.assertEqual(transcripts['Human_coding_ENSP00000371111.2.m1'].combined_cds_start, 33976)
        self.assertEqual(transcripts['Human_coding_ENSP00000371111.2.m1'].combined_cds_end,
                         transcripts['Human_coding_ENSP00000371111.2.m1'].start)
        self.assertEqual(transcripts['Human_coding_ENSP00000371111.2.m1'].combined_cds_start,
                         transcripts['Human_coding_ENSP00000371111.2.m1'].end)

        cds_coordinates = dict()
        for transcript in locus:
            cds_coordinates[transcript] = (locus[transcript].combined_cds_start, locus[transcript].combined_cds_end)

        corr = {1: "Human_coding_ENSP00000371111.2.m1",  # 1645	33976
                2: "Mikado_gold_mikado.0G230.1",  # 1	34063
                3: "ACOCA10068_run2_woRNA_ACOCA10068_r3_0032600.1"  # 1032	34095
                }

        for pad_distance, max_splice in zip((130, 700, 1500, 2000), (1, )):
            with self.subTest(pad_distance=pad_distance, max_splice=max_splice):
                logger = create_default_logger("logger", level="WARNING")
                locus.logger = logger
                locus.configuration.pick.alternative_splicing.ts_distance = pad_distance
                locus.configuration.pick.alternative_splicing.ts_max_splices = max_splice
                locus.logger.setLevel("DEBUG")
                with self.assertLogs(logger, level="DEBUG") as pado:
                    locus.pad_transcripts()

                for tid in corr:
                    self.assertIn(corr[tid], locus.transcripts, corr[tid])

                for transcript in locus:
                    self.assertGreater(locus[transcript].combined_cds_length, 0, transcript)
                    self.assertEqual(locus[transcript].combined_cds_start, cds_coordinates[transcript][0])
                    self.assertEqual(locus[transcript].combined_cds_end, cds_coordinates[transcript][1])
                if pad_distance > 720:  # Ends must be uniform
                    self.assertEqual(locus[corr[1]].end, locus[corr[3]].end,
                                     ([locus[corr[_]].end for _ in range(1, 4)],
                                     locus._share_extreme(transcripts[corr[1]],
                                                          transcripts[corr[2]],
                                                          three_prime=False))
                                     )
                    self.assertEqual(locus[corr[1]].end, locus[corr[2]].end,
                                     ([locus[corr[_]].end for _ in range(1, 4)],
                                      locus._share_extreme(transcripts[corr[1]],
                                                          transcripts[corr[2]],
                                                          three_prime=False))
                                     )

                elif pad_distance < 20:
                    self.assertNotEqual(locus[corr[1]].end, locus[corr[3]].end)
                    self.assertNotEqual(locus[corr[1]].end, locus[corr[2]].end)
                    self.assertNotEqual(locus[corr[2]].end, locus[corr[3]].end)

                if pad_distance >= (abs(transcripts[corr[1]].start - transcripts[corr[2]].start)):
                    self.assertEqual(locus[corr[1]].start,
                                     locus[corr[2]].start)
                    self.assertEqual(locus[corr[1]].start,
                                     locus[corr[3]].start)
                else:

                    self.assertNotEqual(locus[corr[1]].start, locus[corr[2]].start,
                                        (abs(transcripts[corr[1]].start - transcripts[corr[2]].start),
                                         pad_distance,
                                         locus._share_extreme(transcripts[corr[1]], transcripts[corr[2]],
                                                              three_prime=True)
                                        ))

                if pad_distance >= (abs(transcripts[corr[1]].start - transcripts[corr[3]].start)):
                    self.assertEqual(locus[corr[3]].start,
                                     locus[corr[1]].start)
                elif pad_distance <= 130:
                    with self.assertLogs(logger, "DEBUG") as cm:
                        five_graph = locus.define_graph(objects=transcripts,
                                                        inters=locus._share_extreme, three_prime=False)
                        three_graph = locus.define_graph(objects=transcripts,
                                                        inters=locus._share_extreme, three_prime=True)
                        print(five_graph.edges)
                        print(three_graph.edges)

                        boundaries = locus._find_communities_boundaries(five_graph, three_graph)
                        print(boundaries)

                        # locus.logger.setLevel("DEBUG")
                        shared = locus._share_five_prime(transcripts[corr[2]], transcripts[corr[3]])
                        self.assertTrue(shared is False, cm.output)
                        shared = locus._share_five_prime(transcripts[corr[1]], transcripts[corr[3]])
                        self.assertTrue(shared is False, cm.output)

                    # self.assertEqual()
                    self.assertNotEqual(locus[corr[3]].start, locus[corr[1]].start,
                                        pado.output)

    @mark.slow
    def test_padding(self):
        genome = pkg_resources.resource_filename("Mikado.tests", "padding_test.fa")
        transcripts = self.load_from_bed("Mikado.tests", "padding_test.bed12")

        ids = ["mikado.44G2.{}".format(_) for _ in range(1, 6)]

        params = {
            "mikado.44G2.1": (sum([exon[1] + 1 - max(exon[0], transcripts["mikado.44G2.2"].end)
                                   for exon in transcripts["mikado.44G2.1"].exons
                                   if exon[1] > transcripts["mikado.44G2.2"].end]), 1),
            "mikado.44G2.5": (sum([exon[1] + 1 - max(exon[0], transcripts["mikado.44G2.2"].end)
                                   for exon in transcripts["mikado.44G2.5"].exons
                                   if exon[1] > transcripts["mikado.44G2.2"].end]), 4)
        }

        print(params)

        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2], level="INFO")
        for pad_distance, max_splice, coding, best in itertools.product((200, 1000, 1200, 5000), (1, 1, 5),
                                                                        (True, False),
                                                                        ("mikado.44G2.1", "mikado.44G2.5")):
            with self.subTest(pad_distance=pad_distance, max_splice=max_splice, coding=coding, best=best):
                primary = transcripts[best].copy()
                if coding is False:
                    primary.strip_cds()
                locus = loci.Locus(primary)
                locus.configuration.reference.genome = genome
                for t in transcripts:
                    if t == locus.primary_transcript_id:
                        continue
                    trans = transcripts[t].copy()
                    if coding is False:
                        trans.strip_cds()
                    locus.add_transcript_to_locus(trans)

                # Now add the scores
                scores = {best: 15}
                for tid in ids:
                    if tid != best:
                        scores[tid] = 10

                locus._load_scores(scores=scores)

                if coding:
                    cds_coordinates = dict()
                    for transcript in locus:
                        cds_coordinates[transcript] = (
                        locus[transcript].combined_cds_start, locus[transcript].combined_cds_end)

                locus.logger = logger
                locus.configuration.pick.alternative_splicing.ts_distance = pad_distance
                locus.configuration.pick.alternative_splicing.ts_max_splices = max_splice
                locus.pad_transcripts()

                # The .1 transcript can NEVER be expanded, it ends within an intron.
                self.assertFalse(locus[best].attributes.get("padded", False))
                # self.assertFalse(locus["mikado.44G2.1"].attributes.get("padded", False))

                if params[best][0] <= pad_distance and params[best][1] <= max_splice:
                    for trans in ids:
                        if trans in params.keys():
                            continue
                        self.assertTrue(locus[trans].attributes.get("padded", False),
                                        (locus[trans].id, best, locus[trans].end, pad_distance, max_splice,
                                         params[best],
                                         {item for item in locus[trans].attributes.items() if "ts" in item[0]}))
                        self.assertEqual(locus[trans].end, locus[best].end,
                                         (locus[trans].id, best, locus[trans].end, pad_distance, max_splice,
                                          params[best],
                                          {item for item in locus[trans].attributes.items() if "ts" in item[0]}
                                          ))
                else:
                    for trans in ids:
                        if trans in params.keys():
                            continue
                        self.assertFalse(locus[trans].attributes.get("padded", False),
                                         ((locus[trans].id, locus[trans].start, locus[trans].end,
                                          transcripts[trans].start, transcripts[trans].end),
                                          best,
                                          pad_distance, max_splice,
                                          params[best],
                                          {item for item in locus[trans].attributes.items() if "ts" in item[0]}))

                self.assertEqual(locus["mikado.44G2.3"].end, locus["mikado.44G2.2"].end)
                self.assertEqual(locus["mikado.44G2.4"].end, locus["mikado.44G2.2"].end)

                if coding:
                    for transcript in locus:
                        self.assertGreater(locus[transcript].combined_cds_length, 0)
                        self.assertEqual(locus[transcript].combined_cds_start, cds_coordinates[transcript][0])
                        self.assertEqual(locus[transcript].combined_cds_end, cds_coordinates[transcript][1])

    @mark.slow
    def test_phasing(self):

        transcripts = self.load_from_bed("Mikado.tests", "phasing_padding.bed12")
        # We have to test that the CDS is reconstructed correctly even when considering the phasing
        genome = self.fai
        logger = create_null_logger("test_phasing", level="INFO")
        for phase in (0, 1, 2):
            with self.subTest(phase=phase):
                locus = Locus(transcripts["AT5G01030.1"], logger=logger)
                locus.configuration.reference.genome = genome
                other = transcripts["AT5G01030.2"].deepcopy()
                self.assertNotEqual(other.start, locus["AT5G01030.1"].start)
                other.unfinalize()
                other.start += (3 - phase) % 3
                other.remove_exon((10644, 12665))
                other.add_exon((other.start, 12665))
                other.add_exon((other.start, 12665), feature="CDS", phase=phase)
                with self.assertLogs(other.logger, level="DEBUG") as cmo:
                    other.finalize()
                self.assertTrue(other.is_coding, cmo.output)
                self.assertEqual(other.phases[(other.start, 12665)], phase)
                self.assertEqual(other.combined_cds_start, other.start)
                locus.add_transcript_to_locus(other)
                self.assertIn(other.id, locus.transcripts)
                locus.logger.setLevel("DEBUG")
                locus.pad_transcripts()
                # self.assertEqual("", locus[other.id].format("bed12"))
                self.assertEqual(locus[other.id].start, transcripts["AT5G01030.1"].start, phase)
                self.assertEqual(locus[other.id].combined_cds_start, transcripts["AT5G01030.1"].combined_cds_start)

    def test_pad_monoexonic(self):

        transcript = Transcript()
        transcript.chrom, transcript.strand, transcript.id = "Chr5", "+", "mono.1"
        transcript.add_exons([(2001, 3000)])
        transcript.finalize()
        backup = transcript.deepcopy()

        template_one = Transcript()
        template_one.chrom, template_one.strand, template_one.id = "Chr5", "+", "multi.1"
        template_one.add_exons([(1931, 2500), (2701, 3500)])
        template_one.finalize()
        logger = create_null_logger("test_pad_monoexonic")

        for case in range(3):
            with self.subTest(case=case):
                transcript = backup.deepcopy()
                start = template_one if case % 2 == 0 else False
                end = template_one if case > 0 else False

                expanded_one = pad_transcript(transcript,
                                              transcript.deepcopy(),
                                              start, end, self.fai, logger=logger)
                if start:
                    self.assertEqual(expanded_one.start, template_one.start)
                else:
                    self.assertEqual(expanded_one.start, transcript.start)
                if end:
                    self.assertEqual(expanded_one.end, template_one.end)
                else:
                    self.assertEqual(expanded_one.end, transcript.end)

        # Now monoexonic template
        template_two = Transcript()
        template_two.chrom, template_two.strand, template_two.id = "Chr5", "+", "multi.1"
        template_two.add_exons([(1931, 3500)])
        template_two.finalize()

        for case in range(3):
            with self.subTest(case=case):
                transcript = backup.deepcopy()
                start = template_two if case % 2 == 0 else False
                end = template_two if case > 0 else False

                expanded_one = pad_transcript(transcript,
                                              transcript.deepcopy(),
                                              start, end, self.fai, logger=logger)
                if start:
                    self.assertEqual(expanded_one.start, template_two.start)
                else:
                    self.assertEqual(expanded_one.start, transcript.start)
                if end:
                    self.assertEqual(expanded_one.end, template_two.end)
                else:
                    self.assertEqual(expanded_one.end, transcript.end)

        # Now monoexonic template
        template_three = Transcript()
        template_three.chrom, template_three.strand, template_three.id = "Chr5", "+", "multi.1"
        template_three.add_exons([(1501, 1700), (1931, 3500), (4001, 5000)])
        template_three.finalize()

        for case in range(3):
            with self.subTest(case=case):
                transcript = backup.deepcopy()
                start = template_three if case % 2 == 0 else False
                end = template_three if case > 0 else False

                expanded_one = pad_transcript(transcript,
                                              transcript.deepcopy(),
                                              start, end, self.fai, logger=logger)
                if start:
                    self.assertEqual(expanded_one.start, start.start)
                    self.assertIn((1501, 1700), expanded_one.exons)
                else:
                    self.assertEqual(expanded_one.start, transcript.start)
                    self.assertNotIn((1501, 1700), expanded_one.exons)
                if end:
                    self.assertEqual(expanded_one.end, end.end)
                    self.assertIn((4001, 5000), expanded_one.exons)
                else:
                    self.assertEqual(expanded_one.end, transcript.end)
                    self.assertNotIn((4001, 5000), expanded_one.exons)

    def test_pad_multiexonic(self):

        transcript = Transcript()
        transcript.chrom, transcript.strand, transcript.id = "Chr5", "+", "mono.1"
        transcript.add_exons([(2001, 2400), (2800, 3000)])
        transcript.finalize()
        backup = transcript.deepcopy()

        template_one = Transcript()
        template_one.chrom, template_one.strand, template_one.id = "Chr5", "+", "multi.1"
        template_one.add_exons([(1931, 2500), (2701, 3500)])
        template_one.finalize()
        logger = create_null_logger("test_pad_monoexonic")

        for case in range(3):
            with self.subTest(case=case):
                transcript = backup.deepcopy()
                start = template_one if case % 2 == 0 else False
                end = template_one if case > 0 else False

                expanded_one = pad_transcript(transcript,
                                              transcript.deepcopy(),
                                              start, end, self.fai, logger=logger)
                if start:
                    self.assertEqual(expanded_one.start, template_one.start)
                else:
                    self.assertEqual(expanded_one.start, backup.start)
                if end:
                    self.assertEqual(expanded_one.end, template_one.end)
                else:
                    self.assertEqual(expanded_one.end, backup.end)

        # Now monoexonic template
        template_two = Transcript()
        template_two.chrom, template_two.strand, template_two.id = "Chr5", "+", "multi.1"
        template_two.add_exons([(1931, 3500)])
        template_two.finalize()

        for case in range(3):
            with self.subTest(case=case):
                transcript = backup.deepcopy()
                start = template_two if case % 2 == 0 else False
                end = template_two if case > 0 else False

                expanded_one = pad_transcript(transcript,
                                              transcript.deepcopy(),
                                              start, end, self.fai, logger=logger)
                if start:
                    self.assertEqual(expanded_one.start, template_two.start)
                else:
                    self.assertEqual(expanded_one.start, backup.start)
                if end:
                    self.assertEqual(expanded_one.end, template_two.end)
                else:
                    self.assertEqual(expanded_one.end, transcript.end)

        # Now monoexonic template
        template_three = Transcript()
        template_three.chrom, template_three.strand, template_three.id = "Chr5", "+", "multi.1"
        template_three.add_exons([(1501, 1700), (1931, 3500), (4001, 5000)])
        template_three.finalize()

        for case in range(3):
            with self.subTest(case=case):
                transcript = backup.deepcopy()
                start = template_three if case % 2 == 0 else False
                end = template_three if case > 0 else False

                expanded_one = pad_transcript(transcript,
                                              transcript.deepcopy(),
                                              start, end, self.fai, logger=logger)
                if start:
                    self.assertEqual(expanded_one.start, start.start)
                    self.assertIn((1501, 1700), expanded_one.exons)
                else:
                    self.assertEqual(expanded_one.start, backup.start)
                    self.assertNotIn((1501, 1700), expanded_one.exons)
                if end:
                    self.assertEqual(expanded_one.end, end.end)
                    self.assertIn((4001, 5000), expanded_one.exons)
                else:
                    self.assertEqual(expanded_one.end, backup.end)
                    self.assertNotIn((4001, 5000), expanded_one.exons)

    def test_expand_multi_end(self):

        transcript = Transcript()
        transcript.chrom, transcript.strand, transcript.id = "Chr5", "-", "multi.1"
        transcript.add_exons([
            (12751486, 12751579),
            (12751669, 12751808),
            (12751895, 12752032),
            (12752078, 12752839)])
        transcript.finalize()

        template = Transcript()
        template.chrom, template.strand, template.id = "Chr5", "-", "template"
        template.add_exons([
            (12751151, 12751579),
            (12751669, 12751808),
            (12751895, 12752839),  # This exon terminates exactly as the last exon of the transcript ends
            (12752974, 12753102)
        ])
        template.finalize()

        logger = create_null_logger("test_expand_multi_end")

        # Now let us expand on both ends

        with self.subTest():
            expanded = pad_transcript(transcript, transcript.deepcopy(), template, template,
                                      fai=self.fai, logger=logger)
            self.assertEqual(expanded.exons,
                             [(12751151, 12751579),
                              (12751669, 12751808), (12751895, 12752032),
                              (12752078, 12752839), (12752974, 12753102)])

    def test_expand_both_sides(self):

        transcript = Transcript()
        transcript.chrom, transcript.strand, transcript.id = "Chr5", "+", "test"
        transcript.add_exons([(100053, 100220), (100640, 101832)])
        transcript.finalize()

        template = Transcript()
        template.chrom, template.strand, template.id = "Chr5", "+", "template"
        template.add_exons([(99726, 100031), (100657, 102000)])
        template.finalize()

        with self.subTest():
            backup = transcript.deepcopy()
            logger = create_null_logger()
            pad_transcript(transcript, transcript.deepcopy(), template, template, self.fai, logger=logger)
            self.assertEqual(
                transcript.exons,
                [(99726, 100220), (100640, 102000)]

            )

    def test_no_expansion(self):
        transcript = Transcript()
        transcript.chrom, transcript.strand, transcript.id = "Chr5", "+", "test"
        transcript.add_exons([(100053, 100220), (100657, 101832)])
        transcript.finalize()
        backup = transcript.deepcopy()
        logger = create_null_logger(level="DEBUG")
        with self.assertLogs(logger=logger, level="DEBUG") as cm:

            pad_transcript(transcript, backup, None, None, self.fai, logger)
            self.assertEqual(transcript, backup)

        self.assertIn("DEBUG:null:test does not need to be expanded, exiting", cm.output)

    def test_edge_expansion(self):

        transcript = Transcript()
        transcript.id, transcript.chrom, transcript.strand = "test", "Chr5", "+"
        transcript.add_exons([(194892, 195337), (195406, 195511),
                              (195609, 195694), (195788, 195841),
                              (195982, 196098), (196207, 196255),
                              (196356, 196505), (196664, 196725),
                              (197652, 197987)])
        transcript.finalize()
        backup = transcript.deepcopy()

        start_transcript = Transcript()
        start_transcript.id, start_transcript.chrom, start_transcript.strand = "template", "Chr5", "+"
        start_transcript.add_exons([(194741, 194891), (195179, 195337), (195406, 195511),
                                    (195609, 195694), (195788, 195841), (195982, 196098),
                                    (196207, 196255), (196356, 196505), (196664, 196725),
                                    (196848, 196943)])

        logger = create_null_logger()
        with self.assertLogs(logger=logger, level="DEBUG"):
            pad_transcript(transcript, backup, start_transcript, False, self.fai, logger)
            self.assertNotEqual(transcript, backup)
            self.assertEqual(transcript.exons,
                             [(194741, 195337), (195406, 195511),
                              (195609, 195694), (195788, 195841),
                              (195982, 196098), (196207, 196255),
                              (196356, 196505), (196664, 196725),
                              (197652, 197987)]
                             )

    def test_swap_single(self):

        transcript = Transcript()
        transcript.id, transcript.chrom, transcript.strand = "test", "Chr5", "+"
        transcript.add_exons([(101, 1000)])
        transcript.finalize()
        new = transcript.deepcopy()

        locus = Locus(transcript)
        self.assertEqual(locus.primary_transcript, transcript)
        self.assertEqual(len(locus.exons), 1)

        # False swap
        locus._swap_transcript(transcript, transcript)

        new.unfinalize()
        new.remove_exon((101, 1000))
        new.start, new.end = 51, 1200
        new.add_exons([(51, 200), (501, 1200)])
        new.finalize()
        self.assertEqual(transcript.id, new.id)

        locus._swap_transcript(transcript, new)
        self.assertEqual(len(locus.exons), 2)
        self.assertEqual(locus.exons, set(new.exons))
        self.assertEqual(locus.primary_transcript, new)

        new2 = transcript.deepcopy()
        new2.id = "test2"

        with self.assertRaises(KeyError):
            locus._swap_transcript(transcript, new2)

    def test_swap_see_metrics(self):

        transcript = Transcript()
        transcript.id, transcript.chrom, transcript.strand = "test", "Chr5", "+"
        transcript.add_exons([(101, 1000), (1201, 1500)])
        transcript.finalize()
        new = transcript.deepcopy()
        locus = Locus(transcript)
        locus.configuration.pick.alternative_splicing.only_confirmed_introns = False
        second = Transcript()
        second.id, second.chrom, second.strand = "test2", "Chr5", "+"
        second.add_exons([(101, 1000), (1301, 1600)])
        second.finalize()
        locus.add_transcript_to_locus(second)
        self.assertEqual(len(locus.transcripts), 2)
        locus.filter_and_calculate_scores()
        self.assertAlmostEqual(locus[second.id].exon_fraction, 2/3, places=3)
        self.assertAlmostEqual(locus[transcript.id].exon_fraction, 2 / 3, places=3)
        self.assertEqual(locus.primary_transcript, transcript)

        new.unfinalize()
        new.remove_exon((101, 1000))
        new.start, new.end = 51, 1500
        new.add_exons([(51, 200), (501, 1000)])
        new.finalize()
        self.assertEqual(transcript.id, new.id)
        self.assertEqual(locus.primary_transcript_id, transcript.id)

        locus._swap_transcript(transcript, new)
        self.assertEqual(locus.primary_transcript, new)
        self.assertEqual(locus.exons, {(51, 200), (501, 1000), (101, 1000), (1301, 1600), (1201, 1500)})
        locus.filter_and_calculate_scores()
        self.assertAlmostEqual(locus[second.id].exon_fraction, 2 / 5, places=3)
        self.assertAlmostEqual(locus[transcript.id].exon_fraction, 3 / 5, places=3)


if __name__ == '__main__':
    unittest.main(verbosity=2)
