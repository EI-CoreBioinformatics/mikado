from dataclasses import field
from marshmallow_dataclass import dataclass, Optional, List
from marshmallow import validate
from .picking_config import PickConfiguration
from .prepare_config import PrepareConfiguration
from .serialise_config import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from .configuration import ReferenceConfiguration
import copy


@dataclass
class ShortReads:
    r1: List[str] = field(default_factory=lambda: [])
    r2: List[str] = field(default_factory=lambda: [])
    samples: List[str] = field(default_factory=lambda: [])
    strandedness: List[str] = field(default_factory=lambda: [])


@dataclass
class LongReads:
    files: List[str] = field(default_factory=lambda: [])
    samples: List[str] = field(default_factory=lambda: [])
    strandedness: List[str] = field(default_factory=lambda: [])
    skip_split: bool = field(default=True)


@dataclass
class OrfCalling:
    min_protein_len: int = field(default=30, metadata={"validate": validate.Range(min=0)})
    execute: bool = field(default=True)


@dataclass
class BlastX:
    prot_db: Optional[List[str]] = field(default_factory=lambda: [""])
    evalue: float = field(default=1e-7, metadata={"validate": validate})
    max_target_seqs: int = field(default=10, metadata={"validate": validate.Range(min=1)})
    chunks: int = field(default=10, metadata={"validate": validate.Range(min=1)})


@dataclass
class AlignMethods:
    tophat: Optional[List[str]] = field(default=None)
    hisat: Optional[List[str]] = field(default=None)
    star: Optional[List[str]] = field(default=None)
    gsnap: Optional[List[str]] = field(default=None)


@dataclass
class LongReadAlign:
    star: Optional[List[str]] = field(default=None)
    gmap: Optional[List[str]] = field(default=None)


@dataclass
class AsmMethods:
    stringtie: Optional[List[str]] = field(default=None)
    cufflinks: Optional[List[str]] = field(default=None)
    class2: Optional[List[str]] = field(default=None)
    scallop: Optional[List[str]] = field(default=None)
    trinity: Optional[List[str]] = field(default=None)


@dataclass
class ProgramLoader:
    tophat: str = field(default='')
    gmap: str = field(default='')
    star: str = field(default='')
    hisat: str = field(default='')
    samtools: str = field(default='')
    cufflinks: str = field(default='')
    scallop: str = field(default='')
    trinity: str = field(default='')
    stringtie: str = field(default='')
    class2: str = field(default='')
    transdecoder: str = field(default='')
    prodigal: str = field(default='')
    portcullis: str = field(default='')
    mikado: str = field(default='')
    diamond: str = field(default='')
    blast: str = field(default='')


@dataclass
class Portcullis:
    do: bool = field(default=True)
    canonical_juncs: str = field(default="C,S")


@dataclass
class AlnIndex:
    star: str = field(default="")


@dataclass
class DaijinMikadoConfiguration:
    use_diamond: bool = field(default=True)
    use_prodigal: bool = field(default=True)
    modes: List[str] = field(default_factory=lambda: ["stringent"])


@dataclass
class DaijinConfiguration:
    # Daijin specific sub-modules
    short_reads: ShortReads = field(default_factory=ShortReads)
    long_reads: LongReads = field(default_factory=LongReads)
    orf_calling: OrfCalling = field(default_factory=OrfCalling)
    blastx: BlastX = field(default_factory=BlastX)
    aln_index: AlnIndex = field(default_factory=AlnIndex)
    align_methods: AlignMethods = field(default_factory=AlignMethods)
    long_read_align_methods: LongReadAlign = field(default_factory=LongReadAlign)
    asm_methods: AsmMethods = field(default_factory=AsmMethods)
    portcullis: Portcullis = field(default_factory=Portcullis)
    mikado: DaijinMikadoConfiguration = field(default_factory=DaijinMikadoConfiguration)
    #
    filename: Optional[str] = field(default=None)
    seed: int = field(default=0, metadata={"validate": validate.Range(min=0, max=2**32 - 1)})
    multiprocessing_method: str = field(default="spawn",
                                        metadata={"validate": validate.OneOf(["spawn", "fork", "fork-server"])})
    log_settings: LoggingConfiguration = field(default_factory=LoggingConfiguration)
    db_settings: DBConfiguration = field(default_factory=DBConfiguration)
    serialise: SerialiseConfiguration = field(default_factory=SerialiseConfiguration)
    prepare: PrepareConfiguration = field(default_factory=PrepareConfiguration)
    pick: PickConfiguration = field(default_factory=PickConfiguration)
    reference: ReferenceConfiguration = field(default_factory=ReferenceConfiguration)
    name: str = field(default="Daijin")
    out_dir: str = field(default="Daijin")
    threads: int = field(default=4, metadata={"validate": validate.Range(min=1)})
    scheduler: Optional[str] = field(default=None,
                                     metadata={"validate": validate.OneOf(["SLURM", "LSF", "local", "PBS", ""])})

    # # These fields are loaded *from the scoring configuration*
    # TODO test this
    scoring: Optional[dict] = field(default_factory=dict, default=None)
    cds_requirements: Optional[dict] = field(default_factory=dict, default=None)
    as_requirements: Optional[dict] = field(default_factory=dict, default=None)
    requirements: Optional[dict] = field(default_factory=dict, default=None)
    not_fragmentary: Optional[dict] = field(default_factory=dict, default=None)
    load: ProgramLoader = field(default_factory=ProgramLoader)

    def copy(self):
        return copy.deepcopy(self)
