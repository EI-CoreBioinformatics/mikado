import copy
import dataclasses
from dataclasses import field
from marshmallow import validate, ValidationError
from marshmallow_dataclass import dataclass, Optional
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration, create_null_logger
from Mikado._transcripts.scoring_configuration import ScoringFile
import os
import yaml
import toml
import json
from ..exceptions import InvalidJson
from pkg_resources import resource_filename
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader


@dataclass
class ReferenceConfiguration:
    genome: str = field(default="")
    genome_fai: str = field(default="")
    transcriptome: str = field(default="")


@dataclass
class MikadoConfiguration:
    """
    Configuration properties for Mikado.
    """

    scoring: Optional[ScoringFile] = field(default=None)
    threads: int = field(default=1, metadata={
        "metadata": {"description": "Threads to be used per process"},
        "validate": validate.Range(min=1),
        "required": True
    })
    seed: int = field(default=0, metadata={
        "metadata": {"description": "Random number generator seed, to ensure reproducibility across runs"},
        "validate": validate.Range(min=0, max=2 ** 32 - 1)
    })
    multiprocessing_method: Optional[str] = field(default="spawn", metadata={
        "metadata": {"description": "Which method (fork, spawn, forkserver) Mikado should use for multiprocessing"},
        "validate": validate.OneOf(["spawn", "fork", "fork-server"])
    })
    log_settings: LoggingConfiguration = field(default_factory=LoggingConfiguration, metadata={
                "metadata": {"description": "Settings related to the verbosity of logs"}
    })
    db_settings: Optional[DBConfiguration] = field(default_factory=DBConfiguration, metadata={
                "metadata": {"description": "Settings related to DB connection"}
    })
    serialise: SerialiseConfiguration = field(default_factory=SerialiseConfiguration, metadata={
                "metadata": {"description": "Settings related to data serialisation"}
    })
    prepare: PrepareConfiguration = field(default_factory=PrepareConfiguration, metadata={
                "metadata": {"description": "Settings related to the input data preparation"},
    })
    pick: PickConfiguration = field(default_factory=PickConfiguration, metadata={
                "metadata": {"description": "Settings related to the Mikado pick stage"},
    })
    reference: Optional[ReferenceConfiguration] = field(default_factory=ReferenceConfiguration, metadata={
            })

    _loaded_scoring = None
    filename: Optional[str] = field(default=None)

    def __post_init__(self):
        self.check()

    def copy(self):
        return copy.deepcopy(self)

    def check(self):
        if self.scoring is None or not hasattr(self.scoring.requirements, "parameters"):
            self.load_scoring()
        self.scoring.check(minimal_orf_length=self.pick.orf_loading.minimal_orf_length)
        self.Schema().validate(dataclasses.asdict(self))

    def load_scoring(self, logger=None):
        """
        The purpose of this section is the following:
        - check that the scoring file exists somewhere different from the system folder. If it does, check whether it is
          a valid file.
        - If the local file is broken / wrong / inexistent, use the system one.
        - load the scoring file into the configuration object.

        :param configuration: configuration object to check
        :param logger: logger to use
        """

        if logger is None:
            logger = create_null_logger("check_scoring")
        if self.pick.scoring_file is None:
            if self._loaded_scoring != self.pick.scoring_file:            
                logger.debug("Resetting the scoring to its previous value")
                self.pick.scoring_file = self._loaded_scoring
        elif self._loaded_scoring != self.pick.scoring_file:
            logger.debug("Overwriting the scoring self using '%s' as scoring file", self.pick.scoring_file)
            self.scoring_file = None
        else:
            logger.debug("Restarting")

        options = [os.path.abspath(self.pick.scoring_file),
                   os.path.abspath(os.path.join(os.path.dirname(self.filename or ""),
                                                self.pick.scoring_file)),
                   os.path.abspath(os.path.join(resource_filename("Mikado.configuration", "scoring_files"),
                                                self.pick.scoring_file))]

        if self.filename is not None:
            options.append(os.path.join(os.path.dirname(os.path.abspath(self.filename)),
                                        os.path.basename(self.pick.scoring_file)))
            if not os.path.isabs(self.pick.scoring_file):
                options.append(os.path.join(os.path.dirname(os.path.abspath(self.filename)),
                                            self.pick.scoring_file))

        found = False

        for option in options:
            if not os.path.exists(option):
                continue
            if option.endswith(("yaml", "json", "toml")):
                with open(option) as scoring_file:
                    if option.endswith("yaml"):
                        scoring = yaml.load(scoring_file, Loader=yLoader)
                    elif option.endswith("toml"):
                        scoring = toml.load(scoring_file)
                    else:
                        scoring = json.loads(scoring_file.read())
                    if not isinstance(scoring, dict):
                        continue
                    try:
                        checked = ScoringFile.Schema().load(scoring)
                        checked.check(minimal_orf_length=self.pick.orf_loading.minimal_orf_length)
                        self._loaded_scoring = self.pick.scoring_file
                        # self = check_scoring(self)
                        # self = check_all_requirements(self)
                    except (InvalidJson, ValidationError) as exc:
                        logger.debug("Invalid option: %s", option)
                        logger.warning(exc)
                        continue
                    # self.scoring = dataclasses.asdict(checked.scoring)
                    self.scoring = checked
                    found = True
                    self.pick.scoring_file = option
            if found is True:
                logger.info("Found the correct option: %s", option)
                self.pick.scoring_file = option
                break
        if not found:
            raise InvalidJson("No scoring configuration file found. Options: {}".format(",".join(options)))
        return
