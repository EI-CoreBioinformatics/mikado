from marshmallow_dataclass import dataclass, List, Optional
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
        "description": "Settings related to the verbosity of logs"
    })
    db_settings: DBConfiguration = field(default_factory=DBConfiguration, metadata={
        "description": "Settings related to DB connection"
    })
    serialise: SerialiseConfiguration = field(default_factory=SerialiseConfiguration, metadata={
        "description": "Settings related to data serialisation"
    })
    prepare: PrepareConfiguration = field(default_factory=PrepareConfiguration, metadata={
        "description": "Settings related to the input data preparation",
    })
    pick: PickConfiguration = field(default_factory=PickConfiguration, metadata={
        "description": "Settings related to the Mikado pick stage",
    })
    reference: ReferenceConfiguration = field(default_factory=ReferenceConfiguration, metadata={
        "description": "Settings related to the reference genome"
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
