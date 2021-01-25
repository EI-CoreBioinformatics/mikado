from dataclasses import dataclass, field
import copy
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from dacite import Optional


@dataclass
class ReferenceConfiguration:
    genome: str = ""
    genome_fai: str = ""
    transcriptome: str = ""


@dataclass
class MikadoConfiguration:
    filename: Optional[str] = None
    threads: int = 1
    seed: int = 0
    multiprocessing_method: str = "spawn"
    log_settings: LoggingConfiguration = field(default_factory=LoggingConfiguration)
    db_settings: DBConfiguration = field(default_factory=DBConfiguration)
    serialise: SerialiseConfiguration = field(default_factory=SerialiseConfiguration)
    prepare: PrepareConfiguration = field(default_factory=PrepareConfiguration)
    pick: PickConfiguration = field(default_factory=PickConfiguration)
    reference: ReferenceConfiguration = field(default_factory=ReferenceConfiguration)

    # These fields are loaded *from the scoring configuration*
    scoring: Optional[dict] = None
    cds_requirements: Optional[dict] = None
    as_requirements: Optional[dict] = None
    requirements: Optional[dict] = None
    not_fragmentary: Optional[dict] = None

    def copy(self):
        return copy.deepcopy(self)
