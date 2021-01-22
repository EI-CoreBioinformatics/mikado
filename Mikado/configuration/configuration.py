from dataclasses import dataclass

from Mikado.picking.configuration import PickConfiguration
from Mikado.preparation.configuration import PrepareConfiguration
from Mikado.serializers.configuration import SerialiseConfiguration
from Mikado.utilities.dbutils import DBConfiguration
from Mikado.utilities.log_utils import LoggingConfiguration


@dataclass
class ReferenceConfiguration:
    genome: str = ""
    genome_fai: str = ""
    transcriptome: str = ""


@dataclass
class MikadoConfiguration:
    filename: str
    seed: int = 0
    multiprocessing_method: str = "spawn"
    log_settings: LoggingConfiguration = LoggingConfiguration()
    db_settings: DBConfiguration = DBConfiguration()
    serialise: SerialiseConfiguration = SerialiseConfiguration()
    prepare: PrepareConfiguration = PrepareConfiguration()
    pick: PickConfiguration = PickConfiguration()
    reference: ReferenceConfiguration = ReferenceConfiguration()
