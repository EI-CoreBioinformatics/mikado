from dataclasses import dataclass, field
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from .configuration import ReferenceConfiguration
from typing import List, Optional
from marshmallow import Schema, fields
import copy


class ShortReads(Schema):
    r1: list = fields.List(fields.Str(), missing=list)
    r2: list = fields.List(fields.Str(), missing=list)
    samples: list = fields.List(fields.Str(), missing=list)
    strandedness: list = fields.List(fields.Str(), missing=list)


class LongReads(Schema):
    files: list = fields.List(fields.Str(), missing=list)
    samples: list = fields.List(fields.Str(), missing=list)
    strandedness: list = fields.List(fields.Str(), missing=list)
    skip_split: bool = True


class OrfCalling(Schema):
    min_protein_len: int = 30
    execute: bool = True


class BlastX(Schema):
    prot_db: list = fields.List(fields.Str(), missing=list)
    evalue: float = 1e-7
    max_target_seqs: int = 10
    chunks: int = 10


class AlignMethods(Schema):
    tophat: list = fields.List(fields.Str(), missing=list)
    hisat: list = fields.List(fields.Str(), missing=list)
    star: list = fields.List(fields.Str(), missing=list)
    gsnap: list = fields.List(fields.Str(), missing=list)


class LongReadAlign(Schema):
    star: list = fields.List(fields.Str(), missing=list)
    gmap: list = fields.List(fields.Str(), missing=list)


class AsmMethods(Schema):
    stringtie: list = fields.List(fields.Str(), missing=list)
    cufflinks: list = fields.List(fields.Str(), missing=list)
    class2: list = fields.List(fields.Str(), missing=list)
    scallop: list = fields.List(fields.Str(), missing=list)
    trinity: list = fields.List(fields.Str(), missing=list)


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
    do: bool = True
    canonical_juncs: str = fields.Str(missing="C,S")


class AlnIndex(Schema):
    star: str = fields.Str(missing="")


class DaijinMikadoConfiguration(Schema):
    use_diamond: bool = True
    use_prodigal: bool = True
    modes: List[str] = fields.List(fields.Str(), missing=lambda: ["stringent"])


class DaijinConfiguration(Schema):
    filename: str = fields.Str(missing=None)
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

    # Daijin specific sub-modules
    short_reads = fields.Nested(ShortReads)
    long_reads: LongReads = fields.Nested(LongReads)
    orf_calling: OrfCalling = fields.Nested(OrfCalling)
    blastx: BlastX = fields.Nested(BlastX)
    aln_index: AlnIndex = fields.Nested(AlnIndex)
    align_methods: AlignMethods = fields.Nested(AlignMethods)
    long_read_align_methods: LongReadAlign = fields.Nested(LongReadAlign)
    asm_methods: AsmMethods = fields.Nested(AsmMethods)
    load: ProgramLoader = fields.Nested(ProgramLoader)
    portcullis: Portcullis = fields.Nested(Portcullis)

    mikado: DaijinMikadoConfiguration = fields.Nested(DaijinMikadoConfiguration)

    # These fields are loaded *from the scoring configuration*
    scoring: Optional[dict] = fields.Dict(missing=None)
    cds_requirements: Optional[dict] = fields.Dict(missing=None)
    as_requirements: Optional[dict] = fields.Dict(missing=None)
    requirements: Optional[dict] = fields.Dict(missing=None)
    not_fragmentary: Optional[dict] = fields.Dict(missing=None)

    def copy(self):
        return copy.deepcopy(self)
