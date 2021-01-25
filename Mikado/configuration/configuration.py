from dataclasses import dataclass, field
import copy
from typing import Optional

from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from marshmallow import Schema, fields


@dataclass
class ReferenceConfiguration(Schema):
    genome: str = fields.Str(missing="")
    genome_fai: str = fields.Str(missing="")
    transcriptome: str = fields.Str(missing="")


@dataclass
class MikadoConfiguration(Schema):
    filename: str = fields.Str()
    threads: int = fields.Int(missing=1)
    seed: int = fields.Int(missing=0)
    multiprocessing_method: str = fields.Str(missing="spawn")
    log_settings: LoggingConfiguration = fields.Nested(LoggingConfiguration)
    db_settings: DBConfiguration = fields.Nested(DBConfiguration)
    serialise: SerialiseConfiguration = fields.Nested(SerialiseConfiguration)
    prepare: PrepareConfiguration = fields.Nested(PrepareConfiguration)
    pick: PickConfiguration = fields.Nested(PickConfiguration)
    reference: ReferenceConfiguration = fields.Nested(ReferenceConfiguration)

    # These fields are loaded *from the scoring configuration*
    scoring: Optional[dict] = fields.Dict(missing=None)
    cds_requirements: Optional[dict] = fields.Dict(missing=None)
    as_requirements: Optional[dict] = fields.Dict(missing=None)
    requirements: Optional[dict] = fields.Dict(missing=None)
    not_fragmentary: Optional[dict] = fields.Dict(missing=None)

    def copy(self):
        return copy.deepcopy(self)
