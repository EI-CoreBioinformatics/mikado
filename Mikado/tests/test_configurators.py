import shutil
import unittest
from dataclasses import asdict

import marshmallow
import tempfile

import pkg_resources

from .. import create_default_logger
from ..configuration import configurator
from .._transcripts.scoring_configuration import SizeFilter, TargetScore
from ..configuration.configuration import *
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

    def test_target_score(self):
        """The order of the Union matters for loading!"""

        for value in (0.3, True, False, 10, 100.4):
            with self.subTest(value=value, msg=f"Testing value {value}"):
                schema = {"rescaling": "target", "value": value}
                loaded = TargetScore.Schema().load(schema)
                self.assertEqual(loaded.value, schema["value"])
                if isinstance(value, bool):
                    self.assertIsInstance(loaded.value, bool)
                else:
                    self.assertIsInstance(loaded.value, float)

        for value in ("hello", dict(), []):
            with self.subTest(value=value, msg=f"Testing value {value}"):
                schema = {"rescaling": "target", "value": value}
                with self.assertRaises(marshmallow.ValidationError):
                    TargetScore.Schema().load(schema)

    def test_codon_loading(self):
        for value in [0, "0", "1", 1, "Standard", "Flatworm Mitochondrial"]:
            with self.subTest(value=value, msg=f"Testing value {value}"):
                schema = {"codon_table": value, "max_regression": 0.2, "substitution_matrix": "blosum62"}
                loaded = SerialiseConfiguration.Schema().load(schema)
                if value in [0, "0", "1", 1]:
                    self.assertIsInstance(loaded.codon_table, int)
                else:
                    self.assertIsInstance(loaded.codon_table, str)

        for value in ["8", 8, "Universal"]:
            with self.subTest(value=value, msg=f"Testing value {value}"):
                schema = {"codon_table": value, "max_regression": 0.2, "substitution_matrix": "blosum62"}
                with self.assertRaises(marshmallow.ValidationError):
                    SerialiseConfiguration.Schema().load(schema)

    def test_load_plant_scoring(self):
        mammalian = MikadoConfiguration(pick=PickConfiguration(scoring_file="mammalian.yaml"))
        plant = MikadoConfiguration(pick=PickConfiguration(scoring_file="plant.yaml"))
        insect = MikadoConfiguration(pick=PickConfiguration(scoring_file="HISTORIC/insects.yaml"))

        self.assertTrue(mammalian.scoring.requirements.expression == plant.scoring.requirements.expression
                        != insect.scoring.requirements.expression)
        self.assertTrue(mammalian.scoring.requirements.compiled == plant.scoring.requirements.compiled
                        != insect.scoring.requirements.compiled)
        self.assertTrue(mammalian.scoring.requirements.parameters != plant.scoring.requirements.parameters)

    def test_load_invalid_scoring(self):
        erroneous = tempfile.NamedTemporaryFile(suffix=".yaml", mode="wt")
        plant = MikadoConfiguration(pick=PickConfiguration(scoring_file="plant.yaml"))
        err_scoring = copy.deepcopy(plant.scoring)
        key = next(iter(err_scoring.scoring.keys()))
        err_scoring.scoring[key].rescaling = "invalid"
        key = next(iter(err_scoring.requirements.parameters.keys()))
        del err_scoring.requirements.parameters[key]
        yaml.dump(asdict(err_scoring), erroneous)
        erroneous.flush()
        plant.pick.scoring_file = erroneous.name
        logger = create_default_logger("test_load_invalid_scoring", level="DEBUG")
        with self.assertRaises(InvalidConfiguration):
            plant.load_scoring(logger=logger)
        current = os.getcwd()
        os.chdir(os.path.dirname(erroneous.name))
        plant.scoring_file = os.path.basename(erroneous.name)
        self.assertFalse(os.path.exists("plant.yaml"))
        os.link(erroneous.name, "plant.yaml")
        plant._loaded_scoring = None
        plant.pick.scoring_file = "plant.yaml"
        with self.assertRaises(InvalidConfiguration):
            plant.load_scoring(logger=logger)
        os.remove(os.path.join(tempfile.gettempdir(), "plant.yaml"))
        plant._loaded_scoring = None
        plant.load_scoring(logger=logger)
        self.assertEqual(plant._loaded_scoring,
                         pkg_resources.resource_filename("Mikado.configuration", os.path.join("scoring_files",
                                                                                              "plant.yaml")))
        os.chdir(current)
