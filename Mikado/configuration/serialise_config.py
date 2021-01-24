from dataclasses import dataclass, field
from typing import Union


@dataclass
class FilesConfiguration:
    junctions: list = field(default_factory=list)
    xml: list = field(default_factory=list)
    blast_loading_debug: bool = False
    external_scores: str = ""
    orfs: list = field(default_factory=list)
    transcripts: str = "mikado_prepared.fasta"
    log: str = "serialise.log"
    blast_targets: list = field(default_factory=list)
    output_dir: str = "../serializers"


@dataclass
class SerialiseConfiguration:
    files: FilesConfiguration = field(default_factory=FilesConfiguration)
    substitution_matrix: str = "blosum62"
    blast_flavour: str = "blastx"
    codon_table: Union[str,int] = 0
    max_objects: int = 10000000
    max_regression: float = 0.2
    start_adjustment: bool = True
    max_target_seqs: int = 100000
    force: bool = False
    single_thread: bool = False
