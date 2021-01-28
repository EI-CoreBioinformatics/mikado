from dataclasses import field
from marshmallow_dataclass import dataclass, Optional, List
from marshmallow import validate
from .configuration import MikadoConfiguration
import copy


@dataclass
class ShortReads:
    r1: List[str] = field(default_factory=lambda: [], metadata={
        "metadata": {"description": "Array of left read files."},
    })
    r2: List[str] = field(default_factory=lambda: [], metadata={
        "metadata": {"description": "Array of right read files. It must be of the same length of r1; if one or more of the samples are single-end reads, add an empty string."},
    })
    samples: List[str] = field(default_factory=lambda: [], metadata={
        "metadata": {"description": "Array of the sample names. It must be of the same length of r1."},
    })
    strandedness: List[str] = field(default_factory=lambda: [], metadata={
        "metadata": {"description": "Array of strand-specificity of the samples. It must be of the same length of r1. Valid values: fr-firststrand, fr-secondstrand, fr-unstranded."},
    })


@dataclass
class LongReads:
    files: List[str] = field(default_factory=lambda: [])
    samples: List[str] = field(default_factory=lambda: [])
    strandedness: List[str] = field(default_factory=lambda: [])
    skip_split: bool = field(default=True)


@dataclass
class OrfCalling:
    min_protein_len: int = field(default=30, metadata={
        "metadata": {"description": "minimum length of called proteins (in AAs). Default: 30 (90 nts)"},
        "validate": validate.Range(min=0)
    })
    execute: bool = field(default=True, metadata={
        "description": "boolean flag. Default: true, ie execute the ORF calling."
    })


@dataclass
class BlastX:
    prot_db: Optional[List[str]] = field(default_factory=lambda: [], metadata={
        "metadata": {"description": "FASTA file(s) to be used as database for the homology search. Multiple files can be provided and they will be merged into a single database before running."},
    })
    evalue: float = field(default=1e-7, metadata={
        "metadata": {"description": "Maximum e-value. Default 1e-7"},
        "validate": validate.Range(min=0)
    })
    max_target_seqs: int = field(default=10, metadata={
        "metadata": {"description": "Maximum number of targets that the homology will report. Default: 10"},
        "validate": validate.Range(min=1)
    })
    chunks: int = field(default=10, metadata={
        "metadata": {"description": "Number of chunks to divide the search into. Must be equal or greater than the number of processes."},
        "validate": validate.Range(min=1)
    })


@dataclass
class AlignMethods:
    tophat: List[str] = field(default_factory=lambda: [])
    hisat: List[str] = field(default_factory=lambda: [])
    star: List[str] = field(default_factory=lambda: [])
    gsnap: List[str] = field(default_factory=lambda: [])


@dataclass
class LongReadAlign:
    star: List[str] = field(default_factory=lambda: [])
    gmap: List[str] = field(default_factory=lambda: [])


@dataclass
class AsmMethods:
    stringtie: List[str] = field(default_factory=lambda: [])
    cufflinks: List[str] = field(default_factory=lambda: [])
    class2: List[str] = field(default_factory=lambda: [])
    scallop: List[str] = field(default_factory=lambda: [])
    trinity: List[str] = field(default_factory=lambda: [])


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
    use_diamond: bool = field(default=True, metadata={
        "metadata": {"description": "use DIAMOND instead of NCBI BLASTX. Default: true"},
    })
    use_prodigal: bool = field(default=True, metadata={
        "description": "Use Prodigal instead of TransDecoder for ORF calling. Default: true"
    })
    modes: List[str] = field(default_factory=lambda: ["stringent"], metadata={
        "metadata": {"description": "which mode(s) to run Mikado into. Default: permissive (split multiple ORF models unless there is strong BLAST evidence against the decision)."},
    })


@dataclass
class TGGConfiguration:
    max_mem: int = field(default=6000, metadata={
        "metadata": {"description": "Maximum memory to be used for the assembly. Default: 6000Mb"},
        "validate": validate.Range(min=1000)
    })
    npaths: int = field(default=0, metadata={
        "description": "Number of alignments per sequence, using GMAP. Default: 0 (one alignment per sequence, exclude chimeric)."
    })
    identity: float = field(default=0.95, metadata={
        "metadata": {"description": "minimum identity for any alignment. Default: 95%"},
        "validate": validate.Range(min=0.0, max=1.0)
    })
    coverage: float = field(default=0.70, metadata={
        "metadata": {"description": "minimum coverage for any alignment. Default: 70%"},
        "validate": validate.Range(min=0.0, max=1.0)
    })


@dataclass
class DaijinConfiguration(MikadoConfiguration):

    """
    Configuration properties for Daijin. This is an extended Mikado configuration file and can be given directly to
Mikado itself.
    """

    load: ProgramLoader = field(default_factory=ProgramLoader, metadata={
        "name": "load",
        "metadata": {"description": "Commands to use to load/select the versions of the programs to use. Leave an empty string if no loading is necessary."},
    })
    portcullis: Portcullis = field(default_factory=Portcullis, metadata={
        "name": "portcullis",
        "metadata": {"description": "Options related to portcullis"},
    })
    aln_index: AlnIndex = field(default_factory=AlnIndex, metadata={
        "name": "aln_index",
        "metadata": {"description": "Options related to indexing."},
    })
    long_reads: LongReads = field(default_factory=LongReads, metadata={
        "name": "long_reads",
        "metadata": {"description": "Parameters related to long reads to use for the assemblies."},
    })
    short_reads: ShortReads = field(default_factory=ShortReads, metadata={
        "name": "short_reads",
        "metadata": {"description": "Parameters related to the reads to use for the assemblies."},
    })
    name: str = field(default="Daijin", metadata={
        "description": "Name to be used for the project"
    })
    out_dir: str = field(default="Daijin", metadata={
        "description": "Output directory for the project"
    })
    scheduler: Optional[str] = field(default=None, metadata={
        "metadata": {"description": "Scheduler to be used for the project. Set to null if you plan to use DRMAA or are using a local machine."},
        "validate": validate.OneOf(["SLURM", "LSF", "local", "PBS", ""])
    })
    align_methods: AlignMethods = field(default_factory=AlignMethods, metadata={
        "name": "align_methods",
        "metadata": {"description": "Aligners to use. Each aligner can be invoked multiple times: the per-aligner list includes the extra command line arguments to be passed to the program"},
    })
    long_read_align_methods: LongReadAlign = field(default_factory=LongReadAlign, metadata={
        "name": "long_read_align_methods",
        "metadata": {"description": "Aligners for long reads to use. Each aligner can be invoked multiple times: the per-aligner list includes the extra command line arguments to be passed to the program"},
    })
    asm_methods: AsmMethods = field(default_factory=AsmMethods, metadata={
        "name": "asm_methods",
        "metadata": {"description": "Short-read assemblers to use. Each assembler can be invoked multiple times: the per-aligner list includes the extra command line arguments to be passed to the program"},
    })
    orf_calling: OrfCalling = field(default_factory=OrfCalling, metadata={
        "name": "orf_calling",
        "metadata": {"description": "Parameters related to the ORF calling:"},
    })
    blastx: BlastX = field(default_factory=BlastX, metadata={
        "name": "blastx",
        "description": "Parameters related to how DIAMOND/BLASTX will be run"
    })
    mikado: DaijinMikadoConfiguration = field(default_factory=DaijinMikadoConfiguration, metadata={
        "metadata": {"description": "Parameters related to the Mikado execution."},
        "name": "mikado",
    })
    tgg: TGGConfiguration = field(default_factory=TGGConfiguration, metadata={
        "metadata": {"description": "Options related to genome-guided Trinity."},
        "name": "tgg",
    })
