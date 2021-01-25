from dataclasses import field
from marshmallow_dataclass import dataclass
from marshmallow import validate
from typing import List
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from .configuration import ReferenceConfiguration
from typing import List, Optional
from marshmallow import Schema, fields
import copy


@dataclass
class ShortReads:
    r1: List[str] = field(default_factory=lambda: [])
    r2: List[str] = field(default_factory=lambda: [])
    samples: List[str] = field(default_factory=lambda: [])
    strandedness: List[str] = field(default_factory=lambda: [])


class LongReads(Schema):
    files: list = fields.List(fields.Str(missing=""), missing=lambda: [""])
    samples: list = fields.List(fields.Str(missing=""), missing=lambda: [""])
    strandedness: list = fields.List(fields.Str(missing=""), missing=lambda: [""])
    skip_split: bool = fields.Bool(missing=True)


class OrfCalling(Schema):
    min_protein_len: int = fields.Int(missing=30)
    execute: bool = fields.Bool(missing=True)


class BlastX(Schema):
    prot_db: list = fields.List(fields.Str(missing=""), missing=lambda: [""])
    evalue: float = fields.Float(missing=1e-7)
    max_target_seqs: int = fields.Int(missing=10)
    chunks: int = fields.Int(missing=10)


class AlignMethods(Schema):
    tophat: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    hisat: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    star: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    gsnap: list = fields.List(fields.Str(missing=''), missing=lambda: [""])


class LongReadAlign(Schema):
    star: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    gmap: list = fields.List(fields.Str(missing=''), missing=lambda: [""])


class AsmMethods(Schema):
    stringtie: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    cufflinks: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    class2: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    scallop: list = fields.List(fields.Str(missing=''), missing=lambda: [""])
    trinity: list = fields.List(fields.Str(missing=''), missing=lambda: [""])


class ProgramLoader(Schema):
    tophat: str = fields.Str(missing='')
    gmap: str = fields.Str(missing='')
    star: str = fields.Str(missing='')
    hisat: str = fields.Str(missing='')
    samtools: str = fields.Str(missing='')
    cufflinks: str = fields.Str(missing='')
    scallop: str = fields.Str(missing='')
    trinity: str = fields.Str(missing='')
    stringtie: str = fields.Str(missing='')
    class2: str = fields.Str(missing='')
    transdecoder: str = fields.Str(missing='')
    prodigal: str = fields.Str(missing='')
    portcullis: str = fields.Str(missing='')
    mikado: str = fields.Str(missing='')
    diamond: str = fields.Str(missing='')
    blast: str = fields.Str(missing='')


class Portcullis(Schema):
    do: bool = fields.Bool(missing=True)
    canonical_juncs: str = fields.Str(missing="C,S")


class AlnIndex(Schema):
    star: str = fields.Str(missing="")


class DaijinMikadoConfiguration(Schema):
    use_diamond: bool = fields.Bool(missing=True)
    use_prodigal: bool = fields.Bool(missing=True)
    modes: List[str] = fields.List(fields.Str(missing=''), missing=lambda: ["stringent"])


class DaijinConfiguration(Schema):
    # Daijin specific sub-modules
    short_reads: ShortReads = fields.Nested(ShortReads)
    long_reads: LongReads = fields.Nested(LongReads)
    orf_calling: OrfCalling = fields.Nested(OrfCalling)
    blastx: BlastX = fields.Nested(BlastX)
    aln_index: AlnIndex = fields.Nested(AlnIndex)
    align_methods: AlignMethods = fields.Nested(AlignMethods)
    long_read_align_methods: LongReadAlign = fields.Nested(LongReadAlign)
    asm_methods: AsmMethods = fields.Nested(AsmMethods)
    portcullis: Portcullis = fields.Nested(Portcullis)
    mikado: DaijinMikadoConfiguration = fields.Nested(DaijinMikadoConfiguration)
    #
    filename: Optional[str] = fields.Str(missing=None)
    seed: int = fields.Int(missing=0)
    multiprocessing_method: str = fields.Str(missing="spawn")
    log_settings: LoggingConfiguration = fields.Nested(LoggingConfiguration)
    db_settings: DBConfiguration = fields.Nested(DBConfiguration)
    serialise: SerialiseConfiguration = fields.Nested(SerialiseConfiguration)
    prepare: PrepareConfiguration = fields.Nested(PrepareConfiguration)
    pick: PickConfiguration = fields.Nested(PickConfiguration)
    reference: ReferenceConfiguration = fields.Nested(ReferenceConfiguration)
    name: str = fields.Str(missing="Daijin")
    out_dir: str = fields.Str(missing="daijin")
    threads: int = fields.Int(missing=4)
    scheduler: str = fields.Str(missing=None)
    #
    # # These fields are loaded *from the scoring configuration*
    scoring: Optional[dict] = fields.Dict(missing=None)
    cds_requirements: Optional[dict] = fields.Dict(missing=None)
    as_requirements: Optional[dict] = fields.Dict(missing=None)
    requirements: Optional[dict] = fields.Dict(missing=None)
    not_fragmentary: Optional[dict] = fields.Dict(missing=None)
    load: ProgramLoader = fields.Nested(ProgramLoader)

    def copy(self):
        return copy.deepcopy(self)
