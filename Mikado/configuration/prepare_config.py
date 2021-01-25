from marshmallow import validate
from dataclasses import field
from typing import List
from marshmallow_dataclass import dataclass


@dataclass
class PrepareFilesConfiguration:
    output_dir: str = field(default="")
    out: str = field(default="mikado_prepared.gtf", metadata={"required": True})
    out_fasta: str = field(default="mikado_prepared.fasta", metadata={"required": True})
    log: str = field(default="prepare.log")
    gff: List[str] = field(default_factory=lambda: [""], metadata={"required": True})
    labels: List[str] = field(default_factory=lambda: [""])
    strand_specific_assemblies: list = field(default_factory=lambda: [""])
    reference: List[str] = field(default_factory=lambda: [""])
    exclude_redundant: List[str] = field(default_factory=lambda: [""])
    strip_cds: List[str] = field(default_factory=lambda: [""])
    source_score: dict = field(default_factory=dict)


@dataclass
class PrepareConfiguration:
    exclude_redundant: bool = field(default=False)
    minimum_cdna_length: int = field(default=200, metadata={"validate": validate.Range(min=1), "required": True})
    max_intron_length: int = field(default=1000000, metadata={"validate": validate.Range(min=20), "required": True})
    strip_cds: bool = field(default=False)
    strip_faulty_cds: bool = field(default=False)
    single: bool = field(default=False)
    lenient: bool = field(default=False)
    strand_specific: bool = field(default=False, metadata={"required": True})
    canonical: List[str] = field(default_factory=lambda: [["GT", "AG"], ["GC", "AG"], ["AT", "AC"]])
    files: PrepareFilesConfiguration = field(default_factory=PrepareFilesConfiguration, metadata={"required": True})
