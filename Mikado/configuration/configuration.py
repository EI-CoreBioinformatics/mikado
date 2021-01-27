from marshmallow_dataclass import dataclass, Optional
import copy
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from marshmallow import validate
from dataclasses import field


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

    threads: int = field(default=1, metadata={
        "description": "Threads to be used per process",
        "validate": validate.Range(min=1),
    })
    seed: int = field(default=0, metadata={
        "description": "Random number generator seed, to ensure reproducibility across runs",
        "validate": validate.Range(min=0, max=2**32 - 1),
    })
    multiprocessing_method: Optional[str] = field(default="spawn", metadata={
        "description": "Which method (fork, spawn, forkserver) Mikado should use for multiprocessing",
        "validate": validate.OneOf(["spawn", "fork", "fork-server"])
    })
    log_settings: LoggingConfiguration = field(default_factory=LoggingConfiguration, metadata={
        "name": "log_settings",
        "description": "Settings related to the verbosity of logs"
    })
    db_settings: Optional[DBConfiguration] = field(default_factory=DBConfiguration, metadata={
        "name": "db_settings",
        "description": "Settings related to DB connection"
    })
    serialise: SerialiseConfiguration = field(default_factory=SerialiseConfiguration, metadata={
        "name": "serialise",
        "description": "Settings related to data serialisation"
    })
    prepare: PrepareConfiguration = field(default_factory=PrepareConfiguration, metadata={
        "name": "prepare",
        "description": "Settings related to the input data preparation",
    })
    pick: PickConfiguration = field(default_factory=PickConfiguration, metadata={
        "name": "pick",
        "description": "Settings related to the Mikado pick stage",
    })
    reference: Optional[ReferenceConfiguration] = field(default_factory=ReferenceConfiguration, metadata={
        "name": "reference",
    })

    # These fields are loaded *from the scoring configuration*
    scoring: Optional[dict] = field(default=None)
    cds_requirements: Optional[dict] = field(default=None)
    as_requirements: Optional[dict] = field(default=None)
    requirements: Optional[dict] = field(default=None)
    not_fragmentary: Optional[dict] = field(default=None)
    filename: Optional[str] = field(default=None)

    def copy(self):
        return copy.deepcopy(self)
