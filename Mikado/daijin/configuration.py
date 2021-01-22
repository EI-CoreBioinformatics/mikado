from dataclasses import dataclass, field
from ..picking.configuration import PickConfiguration
from ..preparation.configuration import PrepareConfiguration
from ..serializers.configuration import SerialiseConfiguration
from ..utilities.dbutils import DBConfiguration
from ..utilities.log_utils import LoggingConfiguration
from ..configuration.configuration import ReferenceConfiguration


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
    canonical_juncs: list = field(default_factory=lambda: ["C", "S"])


class AlnIndex:
    star: str = ""


@dataclass
class DaijinConfiguration:
    filename: str = None
    seed: int = 0
    multiprocessing_method: str = "spawn"
    log_settings: LoggingConfiguration = LoggingConfiguration()
    db_settings: DBConfiguration = DBConfiguration()
    serialise: SerialiseConfiguration = SerialiseConfiguration()
    prepare: PrepareConfiguration = PrepareConfiguration()
    pick: PickConfiguration = PickConfiguration()
    reference: ReferenceConfiguration = ReferenceConfiguration()
    name: str = "Daijin"
    out_dir: str = "daijin"
    threads: int = 4
    scheduler: str = None

    # Daijin specific sub-modules
    short_reads: ShortReads = ShortReads()
    long_reads: LongReads = LongReads()
    orf_calling: OrfCalling = OrfCalling()
    blastx: BlastX = BlastX()
    aln_index: AlnIndex = AlnIndex()
    align_methods: AlignMethods = AlignMethods()
    long_read_align: LongReadAlign = LongReadAlign()
    asm_methods: AsmMethods = AsmMethods()
    load: ProgramLoader = ProgramLoader()
    portcullis: Portcullis = Portcullis()