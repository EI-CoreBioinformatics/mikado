from .. import configuration
import unittest
import os
import pkg_resources
import glob


__author__ = 'Luca Venturini'


class TestScoring(unittest.TestCase):

    def setUp(self):
        pass

    def test_missing_both(self):

        scor_conf = configuration.MikadoConfiguration()
        scor_conf.requirements = dict()
        scor_conf.requirements["expression"] = ["combined_cds_length"]
        scor_conf.requirements["parameters"] = dict()
        scor_conf.requirements["parameters"]["combined_cds_length"] = dict()
        scor_conf.requirements["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf.requirements["parameters"]["combined_cds_length"]["value"] = 0

        configuration.configurator.check_all_requirements(scor_conf)
        self.assertGreater(len(scor_conf.as_requirements), 0, scor_conf)
        self.assertGreater(len(scor_conf.not_fragmentary), 0, scor_conf)

    def test_missing_as(self):

        scor_conf = configuration.MikadoConfiguration()
        scor_conf.requirements = dict()
        scor_conf.requirements["expression"] = ["cdna_length"]
        scor_conf.requirements["parameters"] = dict()
        scor_conf.requirements["parameters"]["cdna_length"] = dict()
        scor_conf.requirements["parameters"]["cdna_length"]["operator"] = "ge"
        scor_conf.requirements["parameters"]["cdna_length"]["value"] = 0

        scor_conf.not_fragmentary = dict()
        scor_conf.not_fragmentary["expression"] = ["combined_cds_length"]
        scor_conf.not_fragmentary["parameters"] = dict()
        scor_conf.not_fragmentary["parameters"]["combined_cds_length"] = dict()
        scor_conf.not_fragmentary["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf.not_fragmentary["parameters"]["combined_cds_length"]["value"] = 0

        configuration.configurator.check_all_requirements(scor_conf)
        self.assertGreater(len(scor_conf.as_requirements), 0, scor_conf)
        self.assertGreater(len(scor_conf.not_fragmentary), 0, scor_conf)

        self.assertEqual(scor_conf.as_requirements, scor_conf.requirements)
        self.assertNotEqual(scor_conf.not_fragmentary, scor_conf.requirements)

    def test_missing_nf(self):
        scor_conf = configuration.MikadoConfiguration()
        scor_conf.requirements = dict()
        scor_conf.requirements["expression"] = ["cdna_length"]
        scor_conf.requirements["parameters"] = dict()
        scor_conf.requirements["parameters"]["cdna_length"] = dict()
        scor_conf.requirements["parameters"]["cdna_length"]["operator"] = "ge"
        scor_conf.requirements["parameters"]["cdna_length"]["value"] = 0

        scor_conf.as_requirements = dict()
        scor_conf.as_requirements["expression"] = ["combined_cds_length"]
        scor_conf.as_requirements["parameters"] = dict()
        scor_conf.as_requirements["parameters"]["combined_cds_length"] = dict()
        scor_conf.as_requirements["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf.as_requirements["parameters"]["combined_cds_length"]["value"] = 0

        configuration.configurator.check_all_requirements(scor_conf)
        self.assertGreater(len(scor_conf.as_requirements), 0, scor_conf)
        self.assertGreater(len(scor_conf.not_fragmentary), 0, scor_conf)

        self.assertEqual(scor_conf.not_fragmentary, scor_conf.requirements)
        self.assertNotEqual(scor_conf.as_requirements, scor_conf.requirements)

    def test_all_in(self):
        scor_conf = configuration.MikadoConfiguration()
        scor_conf.requirements = dict()
        scor_conf.requirements["expression"] = ["cdna_length"]
        scor_conf.requirements["parameters"] = dict()
        scor_conf.requirements["parameters"]["cdna_length"] = dict()
        scor_conf.requirements["parameters"]["cdna_length"]["operator"] = "ge"
        scor_conf.requirements["parameters"]["cdna_length"]["value"] = 0

        scor_conf.as_requirements = dict()
        scor_conf.as_requirements["expression"] = ["combined_cds_length"]
        scor_conf.as_requirements["parameters"] = dict()
        scor_conf.as_requirements["parameters"]["combined_cds_length"] = dict()
        scor_conf.as_requirements["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf.as_requirements["parameters"]["combined_cds_length"]["value"] = 0

        scor_conf.not_fragmentary = dict()
        scor_conf.not_fragmentary["expression"] = ["selected_cds_length"]
        scor_conf.not_fragmentary["parameters"] = dict()
        scor_conf.not_fragmentary["parameters"]["selected_cds_length"] = dict()
        scor_conf.not_fragmentary["parameters"]["selected_cds_length"]["operator"] = "ge"
        scor_conf.not_fragmentary["parameters"]["selected_cds_length"]["value"] = 0

        configuration.configurator.check_all_requirements(scor_conf)
        self.assertGreater(len(scor_conf.as_requirements), 0, scor_conf)
        self.assertGreater(len(scor_conf.not_fragmentary), 0, scor_conf)

        self.assertNotEqual(scor_conf.not_fragmentary, scor_conf.requirements)
        self.assertNotEqual(scor_conf.as_requirements, scor_conf.requirements)
        self.assertNotEqual(scor_conf.as_requirements, scor_conf.not_fragmentary)

    def test_available_scoring_files(self):

        for scorer in glob.iglob(os.path.join(pkg_resources.resource_filename("Mikado.configuration", "scoring_files"),
                                 "**", "*yaml"), recursive=True):
            conf = configuration.MikadoConfiguration()
            conf.pick.scoring_file = scorer
            conf.filename = None
            self.assertGreater(len(configuration.configurator.check_json(conf).scoring), 0, scorer)

    def test_analyse_twice(self):

        """Verify that parsing and reparsing the configuration dictionary will not cause any errors."""

        scor_conf = configuration.configurator.to_json(None)
        self.assertGreater(len(scor_conf.as_requirements), 0)
        scor_conf = configuration.configurator.check_json(scor_conf)
        self.assertGreater(len(scor_conf.as_requirements), 0)


if __name__ == "__main__":
    unittest.main()
