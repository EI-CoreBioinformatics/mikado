import copy
from dataclasses import field

from marshmallow import validate
from marshmallow_dataclass import dataclass, Optional

from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration


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
        "metadata": {"description": "Threads to be used per process"},
        "validate": validate.Range(min=1),
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

    # These fields are loaded *from the scoring configuration*
    scoring: Optional[dict] = field(default=None)
    cds_requirements: Optional[dict] = field(default=None)
    as_requirements: Optional[dict] = field(default=None)
    requirements: Optional[dict] = field(default=None)
    not_fragmentary: Optional[dict] = field(default=None)
    filename: Optional[str] = field(default=None)

    def copy(self):
        return copy.deepcopy(self)
