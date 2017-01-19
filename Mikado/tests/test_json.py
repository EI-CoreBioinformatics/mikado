import Mikado.utilities
import Mikado.exceptions
import Mikado.configuration
# from Mikado.parsers.GFF import GffLine
# from Mikado.loci import Transcript
# from Mikado.subprograms.util.trim import trim_coding, trim_noncoding
import unittest
# import tempfile
import os
import pkg_resources

__author__ = 'Luca Venturini'


class TestScoring(unittest.TestCase):

    def setUp(self):
        pass

    def test_missing_both(self):

        scor_conf = dict()
        scor_conf["requirements"] = dict()
        scor_conf["requirements"]["expression"] = ["combined_cds_length"]
        scor_conf["requirements"]["parameters"] = dict()
        scor_conf["requirements"]["parameters"]["combined_cds_length"] = dict()
        scor_conf["requirements"]["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf["requirements"]["parameters"]["combined_cds_length"]["value"] = 0

        Mikado.configuration.configurator.check_all_requirements(scor_conf)
        self.assertIn("as_requirements", scor_conf, scor_conf)
        self.assertIn("not_fragmentary", scor_conf, scor_conf)

    def test_missing_as(self):

        scor_conf = dict()
        scor_conf["requirements"] = dict()
        scor_conf["requirements"]["expression"] = ["cdna_length"]
        scor_conf["requirements"]["parameters"] = dict()
        scor_conf["requirements"]["parameters"]["cdna_length"] = dict()
        scor_conf["requirements"]["parameters"]["cdna_length"]["operator"] = "ge"
        scor_conf["requirements"]["parameters"]["cdna_length"]["value"] = 0

        scor_conf["not_fragmentary"] = dict()
        scor_conf["not_fragmentary"]["expression"] = ["combined_cds_length"]
        scor_conf["not_fragmentary"]["parameters"] = dict()
        scor_conf["not_fragmentary"]["parameters"]["combined_cds_length"] = dict()
        scor_conf["not_fragmentary"]["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf["not_fragmentary"]["parameters"]["combined_cds_length"]["value"] = 0

        Mikado.configuration.configurator.check_all_requirements(scor_conf)
        self.assertIn("as_requirements", scor_conf, scor_conf)
        self.assertIn("not_fragmentary", scor_conf, scor_conf)

        self.assertEqual(scor_conf["as_requirements"], scor_conf["requirements"])
        self.assertNotEqual(scor_conf["not_fragmentary"], scor_conf["requirements"])

    def test_missing_nf(self):
        scor_conf = dict()
        scor_conf["requirements"] = dict()
        scor_conf["requirements"]["expression"] = ["cdna_length"]
        scor_conf["requirements"]["parameters"] = dict()
        scor_conf["requirements"]["parameters"]["cdna_length"] = dict()
        scor_conf["requirements"]["parameters"]["cdna_length"]["operator"] = "ge"
        scor_conf["requirements"]["parameters"]["cdna_length"]["value"] = 0

        scor_conf["as_requirements"] = dict()
        scor_conf["as_requirements"]["expression"] = ["combined_cds_length"]
        scor_conf["as_requirements"]["parameters"] = dict()
        scor_conf["as_requirements"]["parameters"]["combined_cds_length"] = dict()
        scor_conf["as_requirements"]["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf["as_requirements"]["parameters"]["combined_cds_length"]["value"] = 0

        Mikado.configuration.configurator.check_all_requirements(scor_conf)
        self.assertIn("as_requirements", scor_conf, scor_conf)
        self.assertIn("not_fragmentary", scor_conf, scor_conf)

        self.assertEqual(scor_conf["not_fragmentary"], scor_conf["requirements"])
        self.assertNotEqual(scor_conf["as_requirements"], scor_conf["requirements"])

    def test_all_in(self):
        scor_conf = dict()
        scor_conf["requirements"] = dict()
        scor_conf["requirements"]["expression"] = ["cdna_length"]
        scor_conf["requirements"]["parameters"] = dict()
        scor_conf["requirements"]["parameters"]["cdna_length"] = dict()
        scor_conf["requirements"]["parameters"]["cdna_length"]["operator"] = "ge"
        scor_conf["requirements"]["parameters"]["cdna_length"]["value"] = 0

        scor_conf["as_requirements"] = dict()
        scor_conf["as_requirements"]["expression"] = ["combined_cds_length"]
        scor_conf["as_requirements"]["parameters"] = dict()
        scor_conf["as_requirements"]["parameters"]["combined_cds_length"] = dict()
        scor_conf["as_requirements"]["parameters"]["combined_cds_length"]["operator"] = "ge"
        scor_conf["as_requirements"]["parameters"]["combined_cds_length"]["value"] = 0

        scor_conf["not_fragmentary"] = dict()
        scor_conf["not_fragmentary"]["expression"] = ["selected_cds_length"]
        scor_conf["not_fragmentary"]["parameters"] = dict()
        scor_conf["not_fragmentary"]["parameters"]["selected_cds_length"] = dict()
        scor_conf["not_fragmentary"]["parameters"]["selected_cds_length"]["operator"] = "ge"
        scor_conf["not_fragmentary"]["parameters"]["selected_cds_length"]["value"] = 0

        Mikado.configuration.configurator.check_all_requirements(scor_conf)
        self.assertIn("as_requirements", scor_conf, scor_conf)
        self.assertIn("not_fragmentary", scor_conf, scor_conf)

        self.assertNotEqual(scor_conf["not_fragmentary"], scor_conf["requirements"])
        self.assertNotEqual(scor_conf["as_requirements"], scor_conf["requirements"])
        self.assertNotEqual(scor_conf["as_requirements"], scor_conf["not_fragmentary"])

    def test_available_scoring_files(self):

        for scorer in pkg_resources.resource_listdir("Mikado.configuration", "scoring_files"):
            conf = dict()
            conf["pick"] = dict()
            conf["pick"]["scoring_file"] = scorer
            conf["filename"] = os.path.join(os.getcwd(), "foo.yaml")
            self.assertIsInstance(Mikado.configuration.configurator.check_json(conf),
                                  dict, scorer)

    def test_analyse_twice(self):

        """Verify that parsing and reparsing the configuration dictionary will not cause any errors."""

        scor_conf = Mikado.configuration.configurator.to_json(None)
        self.assertIn("as_requirements", scor_conf)
        scor_conf = Mikado.configuration.configurator.check_json(scor_conf)
        self.assertIn("as_requirements", scor_conf)


if __name__ == "__main__":
    unittest.main()
