import re
import unittest

from .. import create_default_logger
from .._transcripts.scoring_configuration import ScoringFile, MinMaxScore, SizeFilter, RangeFilter
from ..loci import Transcript, Superlocus
import pkg_resources
import os
from ..configuration.configurator import load_and_validate_config
from ..serializers.external import ExternalSerializer
import pandas as pd


class ExternalTester(unittest.TestCase):

    def setUp(self):
        self.conf = load_and_validate_config(None)
        self.transcript = Transcript()
        self.transcript.chrom = "15"
        self.transcript.source = "protein_coding"
        self.transcript.start = 47631264
        self.transcript.end = 48051999

        exons = [(47631264, 47631416),
                 (47704590, 47704669),
                 (47762671, 47762742),
                 (47893062, 47893093),
                 (47895572, 47895655),
                 (48051942, 48051999)]

        self.transcript.strand = "+"
        self.transcript.add_exons(exons)
        self.transcript.id = "ENST00000560636"
        self.transcript.parent = "ENSG00000137872"
        self.transcript2 = self.transcript.copy()
        self.transcript2.id = "ENST00000560637"
        # self.assertIn("scoring", self.conf)
        
        # This is for the is_fragment test
        self.other_transcript = Transcript()
        self.other_transcript.start = 48053001
        self.other_transcript.end = 48053301

        exons = [(48053001, 48053301)]

        self.other_transcript.strand = "+"
        self.other_transcript.add_exons(exons)
        self.other_transcript.id = "ENST00000560646"
        self.other_transcript.parent = "ENSG00000137873"

    def test_copying(self):
        self.transcript.external_scores.update({"test": 0, "test1": 1})
        self.assertEqual(self.transcript.external_scores.test, 0)
        self.assertEqual(self.transcript.external_scores.test1, 1)
        transcript = self.transcript.deepcopy()
        self.assertEqual(transcript.external_scores.test, 0)
        self.assertEqual(transcript.external_scores.test1, 1)

    def test_no_handle(self):
        logger = create_default_logger("test_no_handle", level="DEBUG")
        with self.assertLogs(logger.name, level="WARNING") as cmo:
            ExternalSerializer(logger=logger, configuration=load_and_validate_config(None),
                               handle=None)
        self.assertTrue(any([record.msg == "No input file specified. Exiting." for record in cmo.records]))

    def test_no_fasta(self):
        conf = load_and_validate_config(None)
        conf.serialise.files.transcripts = "foo"
        logger = create_default_logger("test_no_fasta", level="DEBUG")
        with self.assertRaises(AssertionError), self.assertLogs(logger.name, level="CRITICAL") as cmo:
            ExternalSerializer(logger=logger, configuration=load_and_validate_config(None),
                               handle="trinity.bed12")
        self.assertTrue(any([record.msg == """I cannot find the mikado prepared FASTA file with the transcripts to analyse.
Please run mikado serialise in the folder with the correct files, and/or modify the configuration or the command line options."""
                             for record in cmo.records]))

    def test_invalid_tsv(self):
        conf = load_and_validate_config(None)
        conf.serialise.files.transcripts = pkg_resources.resource_filename("Mikado.tests",
                                                                           os.path.join("test_external",
                                                                                        "mikado_prepared.fasta"))
        self.assertTrue(os.path.exists(conf.serialise.files.transcripts),
                        conf.serialise.files.transcripts)
        logger = create_default_logger("test_no_fasta", level="DEBUG")
        h5 = pkg_resources.resource_filename("Mikado.tests", os.path.join("test_external", "abundance.h5"))
        with self.assertRaises((pd.errors.ParserError, UnicodeDecodeError)):
            ExternalSerializer(logger=logger, configuration=conf, handle=h5)
        wrong_tsv_one = pkg_resources.resource_filename("Mikado.tests", os.path.join("test_external", "wrong_1.tsv"))
        with self.assertRaises((ValueError, pd.errors.ParserError)), \
                self.assertLogs(logger.name, level="CRITICAL") as cmo:
            ExternalSerializer(logger=logger, configuration=conf, handle=wrong_tsv_one)
        self.assertTrue(any([re.search(r"Invalid input file", str(record.msg)) is not None for record in cmo.records]))
        wrong_tsv_two = pkg_resources.resource_filename("Mikado.tests", os.path.join("test_external", "wrong_2.tsv"))
        with self.assertRaises((ValueError, pd.errors.ParserError)), \
                self.assertLogs(logger.name, level="CRITICAL") as cmo:
            ExternalSerializer(logger=logger, configuration=conf, handle=wrong_tsv_two)
        self.assertTrue(any([re.search(r"Invalid input file", str(record.msg)) is not None for record in cmo.records]))

    def test_real(self):
        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        self.assertTrue(hasattr(self.conf.scoring.requirements, "parameters"))
        self.assertTrue(hasattr(checked_conf.scoring.requirements, "parameters"))
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float"})
        # self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        logger = create_default_logger(name="test_real", level="WARNING")
        sup = Superlocus(self.transcript, configuration=checked_conf, logger=logger)
        self.assertIn("attributes.tpm", sup.configuration.scoring.scoring)
        sup.get_metrics()
        self.assertIn("attributes.tpm", sup._metrics[self.transcript.id])
        self.assertEqual(sup._metrics[self.transcript.id]["attributes.tpm"], 10)

        sup2 = Superlocus(self.transcript2, configuration=checked_conf)
        sup2.get_metrics()
        self.assertIn("attributes.tpm", sup2._metrics[self.transcript2.id])
        self.assertEqual(sup2._metrics[self.transcript2.id]["attributes.tpm"], 0)

    def test_multiplier(self):

        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4})

        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        tid = self.transcript.id
        self.assertIn(tid, sup.transcripts)
        sup.filter_and_calculate_scores(check_requirements=False)
        self.assertEqual(sup.scores[tid]['attributes.tpm'], 4)

    def test_default_attribute_score(self):
        self.transcript.attributes["foo"] = True

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.foo"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "rtype": "bool", "default": False})

        self.assertIn('attributes.foo', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        sup.get_metrics()
        self.assertIn("attributes.foo", sup._metrics[self.transcript.id])
        self.assertEqual(sup._metrics[self.transcript.id]["attributes.foo"], True)

        sup2 = Superlocus(self.transcript2, configuration=checked_conf)
        sup2.get_metrics()
        self.assertIn("attributes.foo", sup2._metrics[self.transcript2.id])
        self.assertEqual(sup2._metrics[self.transcript2.id]["attributes.foo"], False)

    def test_error_attribute(self):
        self.transcript.attributes["tpm"] = "10a"
        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4})
        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)
        sup = Superlocus(self.transcript, configuration=checked_conf)
        with self.assertRaises(ValueError):
            sup.get_metrics()

    def test_attribute_use_raw(self):

        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4,
             'use_raw': True})

        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        tid = self.transcript.id
        self.assertIn(tid, sup.transcripts)
        with self.assertRaises(ValueError):
            sup.get_metrics()

    def test_attribute_use_raw_percentage(self):

        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4,
             'use_raw': True, 'percentage': True})

        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        tid = self.transcript.id
        self.assertIn(tid, sup.transcripts)
        sup.get_metrics()
        self.assertIn("attributes.tpm", sup._metrics[self.transcript.id])
        self.assertEqual(sup._metrics[self.transcript.id]["attributes.tpm"], 0.1)
        sup.define_loci()

    def test_attribute_use_as_alternative_splicing(self):
        from Mikado.loci import Locus
        checked_conf = self.conf.copy()
        self.transcript.attributes['cov'] = 10
        loci = Locus(self.transcript, configuration=checked_conf)
        loci.configuration.scoring.as_requirements.parameters['attributes.cov'] = SizeFilter(operator='ge', value=3,
                                                                                             metric=None,
                                                                                             name='attributes.cov')
        loci.configuration.scoring.as_requirements.expression = [
            loci.configuration.scoring.as_requirements.expression[0] + ' and attributes.cov']
        loci.configuration.scoring.as_requirements._check_my_requirements()
        self.assertTrue(loci._check_as_requirements(self.transcript))

    def test_attribute_use_as_alternative_splicing_fail_filter(self):
        from Mikado.loci import Locus
        checked_conf = self.conf.copy()
        self.transcript.attributes['cov'] = 10
        loci = Locus(self.transcript, configuration=checked_conf)
        loci.configuration.scoring.as_requirements.parameters['attributes.cov'] = SizeFilter(operator='ge', value=15,
                                                                                             metric=None,
                                                                                             name='attributes.cov')
        loci.configuration.scoring.as_requirements.expression = [
            loci.configuration.scoring.as_requirements.expression[0] + ' and attributes.cov']
        loci.configuration.scoring.as_requirements._check_my_requirements()
        self.assertFalse(loci._check_as_requirements(self.transcript))

    def test_attributes_use_as_cds_requirements(self):
        from Mikado.loci import Locus
        checked_conf = self.conf.copy()
        checked_conf.scoring.cds_requirements.expression = ["attributes.something"]
        checked_conf.scoring.cds_requirements.parameters = {"attributes.something": SizeFilter(value=1, operator="ge")}
        checked_conf.scoring.cds_requirements._check_my_requirements()
        self.transcript.attributes['something'] = 2
        loci = Locus(self.transcript, configuration=checked_conf)
        not_passing = loci._check_not_passing(section_name="cds_requirements")
        self.assertSetEqual(not_passing, set())

    def test_attributes_use_as_fragment_filter(self):
        from Mikado.loci import Locus
        checked_conf = self.conf.copy()
        checked_conf.pick.fragments.max_distance = checked_conf.pick.clustering.flank = 5000
        checked_conf.scoring.not_fragmentary.expression = ["attributes.something"]
        checked_conf.scoring.not_fragmentary.parameters = {"attributes.something": SizeFilter(value=1, operator="ge")}
        checked_conf.scoring.not_fragmentary._check_my_requirements()
        self.transcript.attributes['something'] = 2
        logger = create_default_logger("test_attributes_use_as_fragment_filter", level="DEBUG")
        loci = Locus(self.transcript, configuration=checked_conf, logger=logger)
        self.other_transcript.attributes["something"] = 0.5
        fragment = Locus(self.other_transcript, configuration=checked_conf, logger=logger)
        self.assertEqual(fragment.other_is_fragment(loci), (False, None))
        self.assertTrue(loci.other_is_fragment(fragment)[0])
        del self.other_transcript.attributes["something"]
        fragment = Locus(self.other_transcript, configuration=checked_conf, logger=logger)
        self.assertEqual(fragment.other_is_fragment(loci), (False, None))
        self.assertTrue(loci.other_is_fragment(fragment)[0])

    def test_attributes_range_use_as_fragment_filter(self):
        from Mikado.loci import Locus
        checked_conf = self.conf.copy()
        checked_conf.pick.fragments.max_distance = checked_conf.pick.clustering.flank = 5000
        checked_conf.scoring.not_fragmentary.expression = ["attributes.something"]
        checked_conf.scoring.not_fragmentary.parameters = {"attributes.something": RangeFilter(value=[10, 50],
                                                                                               operator="within")}
        checked_conf.scoring.not_fragmentary._check_my_requirements()
        self.transcript.attributes['something'] = 11
        logger = create_default_logger("test_attributes_range_as_fragment_filter", level="DEBUG")
        loci = Locus(self.transcript, configuration=checked_conf, logger=logger)
        self.other_transcript.attributes["something"] = 0.5
        fragment = Locus(self.other_transcript, configuration=checked_conf, logger=logger)
        self.assertEqual(fragment.other_is_fragment(loci), (False, None))
        self.assertTrue(loci.other_is_fragment(fragment)[0])
        del self.other_transcript.attributes["something"]
        fragment = Locus(self.other_transcript, configuration=checked_conf, logger=logger)
        self.assertEqual(fragment.other_is_fragment(loci), (False, None))
        self.assertTrue(loci.other_is_fragment(fragment)[0])
        # Now change the default.
        checked_conf.scoring.not_fragmentary.parameters = {"attributes.something": RangeFilter(value=[10, 50],
                                                                                               operator="within",
                                                                                               default=20)}
        checked_conf.scoring.not_fragmentary._check_my_requirements()
        not_fragment = Locus(self.other_transcript, configuration=checked_conf, logger=logger)
        self.assertNotIn("something", not_fragment.transcripts[self.other_transcript.id].attributes)
        # Not a fragment any more, the default of 15 is within the range.
        self.assertFalse(loci.other_is_fragment(fragment)[0])

    def test_print_metrics_with_attributes_from_scoring(self):
        from Mikado.loci import Locus, Sublocus
        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring = {"attributes.something": MinMaxScore(default=0, rescaling="max", filter=None)}
        checked_conf.scoring.cds_requirements.expression = ["cdna_length"]
        checked_conf.scoring.cds_requirements.parameters = {"cdna_length": SizeFilter(operator="ge", value=1)}
        checked_conf.scoring.requirements.expression = ["cdna_length"]
        checked_conf.scoring.requirements.parameters = {"cdna_length": SizeFilter(operator="ge", value=1)}
        checked_conf.scoring.check(checked_conf.pick.orf_loading.minimal_orf_length)
        for lclass in [Locus, Sublocus]:
            logger = create_default_logger(f"test_print_metrics_with_attributes_{lclass.name}", level="DEBUG")
            locus = lclass(self.transcript, configuration=checked_conf, logger=logger)
            self.assertIn(self.transcript.id, locus.transcripts)
            rows = list(locus.print_metrics())
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertNotIn("something", locus.transcripts[self.transcript.id].attributes)
            self.assertIn("attributes.something", row.keys())
            self.assertEqual(row["attributes.something"], checked_conf.scoring.scoring["attributes.something"].default)
            self.transcript.attributes["something"] = 5
            locus = Locus(self.transcript, configuration=checked_conf)
            rows = list(locus.print_metrics())
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertIn("something", locus.transcripts[self.transcript.id].attributes)
            self.assertIn("attributes.something", row.keys())
            self.assertEqual(row["attributes.something"], 5)
            del self.transcript.attributes["something"]

    def test_print_metrics_with_attributes_from_requirements(self):
        from Mikado.loci import Locus, Sublocus
        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring = {"attributes.something": MinMaxScore(default=0, rescaling="max", filter=None)}
        checked_conf.scoring.cds_requirements.expression = ["attributes.cds"]
        checked_conf.scoring.cds_requirements.parameters = {"attributes.cds": SizeFilter(operator="ge", value=1,
                                                                                         default=1)}
        checked_conf.scoring.requirements.expression = ["attributes.req"]
        checked_conf.scoring.requirements.parameters = {"attributes.req": SizeFilter(operator="ge", value=1,
                                                                                     default=1)}
        checked_conf.scoring.as_requirements.expression = ["attributes.as"]
        checked_conf.scoring.as_requirements.parameters = {"attributes.as": SizeFilter(operator="ge", value=1,
                                                                                       default=1)}
        checked_conf.scoring.not_fragmentary.expression = ["attributes.frag"]
        checked_conf.scoring.not_fragmentary.parameters = {"attributes.frag": SizeFilter(operator="ge", value=1,
                                                                                         default=1)}
        sections = {"something": checked_conf.scoring.scoring,
                    "req": checked_conf.scoring.requirements.parameters,
                    "as": checked_conf.scoring.as_requirements.parameters,
                    "cds": checked_conf.scoring.cds_requirements.parameters,
                    "frag": checked_conf.scoring.not_fragmentary.parameters}

        checked_conf.scoring.check(checked_conf.pick.orf_loading.minimal_orf_length)
        for lclass in [Locus, Sublocus]:
            logger = create_default_logger(f"test_print_metrics_with_attributes_{lclass.name}", level="DEBUG")
            locus = lclass(self.transcript, configuration=checked_conf, logger=logger)
            self.assertIn(self.transcript.id, locus.transcripts)
            rows = list(locus.print_metrics())
            self.assertEqual(len(rows), 1)
            row = rows[0]
            for key, section in sections.items():
                self.assertNotIn(key, locus.transcripts[self.transcript.id].attributes)
                self.assertIn(f"attributes.{key}", row.keys())
                self.assertEqual(row[f"attributes.{key}"],
                                 section[f"attributes.{key}"].default)
            for key in sections:
                self.transcript.attributes[key] = 5
            locus = Locus(self.transcript, configuration=checked_conf)
            rows = list(locus.print_metrics())
            self.assertEqual(len(rows), 1)
            row = rows[0]
            for key, section in sections.items():
                self.assertIn(key, locus.transcripts[self.transcript.id].attributes)
                self.assertIn(f"attributes.{key}", row.keys())
                self.assertEqual(row["attributes.something"], 5)
            # Reset the object
            for key in sections:
                del self.transcript.attributes[key]
