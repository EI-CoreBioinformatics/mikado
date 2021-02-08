import unittest
from ..configuration import configurator
from .._transcripts.scoring_configuration import SizeFilter
from ..loci.abstractlocus import Abstractlocus
import pickle


class TestConfigurator(unittest.TestCase):
    def setUp(self) -> None:
        self.configuration = configurator.load_and_validate_config(None)

    def test_pickling(self):
        reloaded = pickle.loads(pickle.dumps(self.configuration))
        self.assertEqual(reloaded, self.configuration)

    def test_eval_compiled(self):
        self.configuration.scoring.requirements.expression = ["cdna_length"]
        self.configuration.scoring.requirements.parameters = dict()
        self.configuration.scoring.requirements.parameters["cdna_length"] = SizeFilter(
            name="cdna_length", value=50, operator="ge")
        self.configuration.scoring.requirements.parameters["selected_cds_length"] = SizeFilter(
            name="cdna_length", value=200, operator="ge")
        self.configuration.scoring.requirements._check_my_requirements()
        values = {"cdna_length": 100, "selected_cds_length": 150}
        evaluated = dict()
        evaluated["cdna_length"] = Abstractlocus.evaluate(
            values["cdna_length"],
            self.configuration.scoring.requirements.parameters["cdna_length"])
        self.assertTrue(evaluated["cdna_length"])
        evaluated["selected_cds_length"] = Abstractlocus.evaluate(
            values["selected_cds_length"],
            self.configuration.scoring.requirements.parameters["selected_cds_length"])
        self.configuration.scoring.requirements._check_my_requirements()
        self.assertFalse(evaluated["selected_cds_length"])
        self.assertTrue(eval(self.configuration.scoring.requirements.compiled))
        # Now we are changing the contents of the expression
        self.configuration.scoring.requirements.expression = ["selected_cds_length"]
        self.configuration.scoring.requirements._check_my_requirements()
        self.assertFalse(eval(self.configuration.scoring.requirements.compiled))
