from dataclasses import dataclass, field
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from .configuration import ReferenceConfiguration
from typing import Union, List
import copy


@dataclass
class ShortReads:

    r1: list = field(default_factory=lambda: [])
    r2: list = field(default_factory=lambda: [])
    samples: list = field(default_factory=lambda: [])
    strandedness: list = field(default_factory=lambda: [])


@dataclass
class LongReads:
    files: list = field(default_factory=lambda: [])
    samples: list = field(default_factory=lambda: [])
    strandedness: list = field(default_factory=lambda: [])
    skip_split: bool = True


@dataclass
class OrfCalling:
    min_protein_len: int = 30
    execute: bool = True


@dataclass
class BlastX:
    prot_db: list = field(default_factory=lambda: [])
    evalue: float = 1e-7
    max_target_seqs: int = 10
    chunks: int = 10


@dataclass
class AlignMethods:
    tophat: list = field(default_factory=lambda: [])
    hisat: list = field(default_factory=lambda: [])
    star: list = field(default_factory=lambda: [])
    gsnap: list = field(default_factory=lambda: [])


@dataclass
class LongReadAlign:
    star: list = field(default_factory=lambda: [])
    gmap: list = field(default_factory=lambda: [])


@dataclass
class AsmMethods:
    stringtie: list = field(default_factory=lambda: [])
    cufflinks: list = field(default_factory=lambda: [])
    class2: list = field(default_factory=lambda: [])
    scallop: list = field(default_factory=lambda: [])
    trinity: list = field(default_factory=lambda: [])


@dataclass
class ProgramLoader:
    tophat: str = ''
    gmap: str = ''
    star: str = ''
    hisat: str = ''
    samtools: str = ''
    cufflinks: str = ''
    scallop: str = ''
    trinity: str = ''
    stringtie: str = ''
    class2: str = ''
    transdecoder: str = ''
    prodigal: str = ''
    portcullis: str = ''
    mikado: str = ''
    diamond: str = ''
    blast: str = ''


@dataclass
class Portcullis:
    do: bool = True
    canonical_juncs: str = "C,S"


@dataclass
class AlnIndex:
    star: str = ""


@dataclass
class DaijinMikadoConfiguration:
    use_diamond: bool = True
    use_prodigal: bool = True
    modes: List[str] = field(default_factory=lambda: ["stringent"])


@dataclass
class DaijinConfiguration:
    filename: Union[str, None] = None
    seed: int = 0
    multiprocessing_method: str = "spawn"
    log_settings: LoggingConfiguration = field(default_factory=LoggingConfiguration)
    db_settings: DBConfiguration = field(default_factory=DBConfiguration)
    serialise: SerialiseConfiguration = field(default_factory=SerialiseConfiguration)
    prepare: PrepareConfiguration = field(default_factory=PrepareConfiguration)
    pick: PickConfiguration = field(default_factory=PickConfiguration)
    reference: ReferenceConfiguration = field(default_factory=ReferenceConfiguration)
    name: str = "Daijin"
    out_dir: str = "daijin"
    threads: int = 4
    scheduler: str = None

    # Daijin specific sub-modules
    short_reads: ShortReads = field(default_factory=ShortReads)
    long_reads: LongReads = field(default_factory=LongReads)
    orf_calling: OrfCalling = field(default_factory=OrfCalling)
    blastx: BlastX = field(default_factory=BlastX)
    aln_index: AlnIndex = field(default_factory=AlnIndex)
    align_methods: AlignMethods = field(default_factory=AlignMethods)
    long_read_align_methods: LongReadAlign = field(default_factory=LongReadAlign)
    asm_methods: AsmMethods = field(default_factory=AsmMethods)
    load: ProgramLoader = field(default_factory=ProgramLoader)
    portcullis: Portcullis = field(default_factory=Portcullis)

    mikado: DaijinMikadoConfiguration = field(default_factory=DaijinMikadoConfiguration)

    # These fields are loaded *from the scoring configuration*
    scoring: dict = None
    cds_requirements: dict = None
    as_requirements: dict = None
    requirements: dict = None
    not_fragmentary: dict = None

    def copy(self):
        return copy.deepcopy(self)
