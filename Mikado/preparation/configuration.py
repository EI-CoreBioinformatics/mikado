from dataclasses import dataclass, field


@dataclass
class PrepareFilesConfiguration:
    output_dir: str = "."
    out: str = "mikado_prepared.gtf"
    out_fasta: str = "mikado_prepared.fasta"
    log: str = "prepare.log"
    gff: list = field(default_factory=list)
    labels: list = field(default_factory=list)
    strand_specific_assemblies: list = field(default_factory=list)
    reference: list = field(default_factory=list)
    exclude_redundant: list = field(default_factory=list)
    strip_cds: list = field(default_factory=list)
    source_score: dict = field(default_factory=dict)


@dataclass
class PrepareConfiguration:
    exclude_redundant: bool = False
    minimum_cdna_length: int = 200
    max_intron_length: int = 1000000
    strip_cds: bool = False
    strip_faulty_cds: bool = False
    single: bool = False
    lenient: bool = False
    strand_specific: bool = False
    canonical: list = field(default_factory=list)
    files: PrepareFilesConfiguration = PrepareFilesConfiguration()
