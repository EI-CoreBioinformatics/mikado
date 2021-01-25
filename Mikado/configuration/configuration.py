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
    threads: int = field(default=1, metadata={"validate": validate.Range(min=1)})
    seed: int = field(default=0, metadata={"validate": validate.Range(min=0, max=2**32 - 1)})
    multiprocessing_method: Optional[str] = field(
        default="spawn",
        metadata={"validate": validate.OneOf(["spawn", "fork", "fork-server"])}) 
    log_settings: LoggingConfiguration = field(default_factory=LoggingConfiguration)
    db_settings: DBConfiguration = field(default_factory=DBConfiguration)
    serialise: SerialiseConfiguration = field(default_factory=SerialiseConfiguration)
    prepare: PrepareConfiguration = field(default_factory=PrepareConfiguration)
    pick: PickConfiguration = field(default_factory=PickConfiguration)
    reference: ReferenceConfiguration = field(default_factory=ReferenceConfiguration)

    # These fields are loaded *from the scoring configuration*
    scoring: Optional[dict] = field(default=None)
    cds_requirements: Optional[dict] = field(default=None)
    as_requirements: Optional[dict] = field(default=None)
    requirements: Optional[dict] = field(default=None)
    not_fragmentary: Optional[dict] = field(default=None)
    filename: Optional[str] = field(default=None)

    def copy(self):
        return copy.deepcopy(self)
