from dataclasses import dataclass
import copy
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration


@dataclass
class ReferenceConfiguration:
    genome: str = ""
    genome_fai: str = ""
    transcriptome: str = ""


@dataclass
class MikadoConfiguration:
    filename: str = None
    threads: int = 1
    seed: int = 0
    multiprocessing_method: str = "spawn"
    log_settings: LoggingConfiguration = LoggingConfiguration()
    db_settings: DBConfiguration = DBConfiguration()
    serialise: SerialiseConfiguration = SerialiseConfiguration()
    prepare: PrepareConfiguration = PrepareConfiguration()
    pick: PickConfiguration = PickConfiguration()
    reference: ReferenceConfiguration = ReferenceConfiguration()

    def copy(self):
        return copy.deepcopy(self)
