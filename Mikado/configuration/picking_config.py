from marshmallow_dataclass import dataclass, Optional
from dataclasses import field
from marshmallow import validate
from typing import List


@dataclass
class AlternativeSplicingConfiguration:
    report: bool = field(default=True)
    cds_only: bool = field(default=False)
    min_cds_overlap: float = field(default=0.5, metadata={"validate": validate.Range(min=0, max=1)})
    min_cdna_overlap: float = field(default=0.6, metadata={"validate": validate.Range(min=0, max=1)})
    keep_retained_introns: bool = field(default=True)
    keep_cds_disrupted_by_ri: bool = field(default=False)
    max_isoforms: int = field(default=10, metadata={"validate": validate.Range(min=1)})
    valid_ccodes: list = field(default_factory=lambda: ["j", "J", "G", "h"],
         metadata={"validate": validate.ContainsOnly(["j", "e", "o", "h", "J", "C", "g", "G", "=", "n", "_"])})
    redundant_ccodes: list = field(default_factory=lambda: ["c", "m", "_", "=", "n"],
        metadata={"validate": validate.ContainsOnly(
            ["j", "n", "O", "e", "o", "h", "J", "C", "c", "m", "mo", "=", "_", "x", "p", "P", "X", "I", "i"])})
    min_score_perc: float = field(default=0.5, metadata={"validate": validate.Range(min=0, max=1)})
    only_confirmed_introns: bool = field(default=True)
    ts_distance: int = field(default=2000, metadata={"validate": validate.Range(min=0)})
    pad: bool = field(default=True)
    ts_max_splices: int = field(default=2, metadata={"validate": validate.Range(min=0)})


@dataclass
class OutputFormatConfiguration:
    source: str = field(default="Mikado")
    id_prefix: str = field(default="mikado")
    report_all_orfs: bool = field(default=False)


@dataclass
class OrfLoadingConfiguration:
    minimal_secondary_orf_length: int = field(default=200, metadata={"validate": validate.Range(min=0)})
    minimal_orf_length: int = field(default=50, metadata={"validate": validate.Range(min=0)})
    strand_specific: bool = field(default=True)


@dataclass
class BlastParamsConfiguration:
    evalue: float = field(default=1e-06, metadata={"validate": validate.Range(min=0)})
    hsp_evalue: float = field(default=1e-06, metadata={"validate": validate.Range(min=0)})
    leniency: str = field(default="STRINGENT",
                          metadata={"validate": validate.OneOf(["STRINGENT", "LENIENT", "PERMISSIVE"])})
    max_target_seqs: int = field(default=3, metadata={"validate": validate.Range(min=1)})
    minimal_hsp_overlap: float = field(default=0.5, metadata={"validate": validate.Range(min=0, max=1)})
    min_overlap_duplication: float = field(default=0.8, metadata={"validate": validate.Range(min=0, max=1)})


@dataclass
class ChimeraSplitConfiguration:
    blast_check: bool = field(default=True)
    execute: bool = field(default=True)
    skip: List[bool] = field(default_factory=lambda: [], metadata={"validate": validate.Length(min=0)})
    blast_params: BlastParamsConfiguration = field(default_factory=BlastParamsConfiguration)


@dataclass
class RunOptionsConfiguration:
    shm: bool = field(default=False)
    exclude_cds: bool = field(default=False)
    # TODO this must be at most 2
    intron_range: List[int] = field(default_factory=lambda: [60, 10000],
                               metadata={"validate": validate.Length(min=2, max=2)})
    reference_update: bool = field(default=False)
    only_reference_update: bool = field(default=False)
    check_references: bool = field(default=False)
    single_thread: bool = field(default=False)


@dataclass
class ClusteringConfiguration:
    cds_only: bool = field(default=False)
    min_cds_overlap: float = field(default=0.2, metadata={"validate": validate.Range(min=0, max=1)})
    min_cdna_overlap: float = field(default=0.2, metadata={"validate": validate.Range(min=0, max=1)})
    purge: bool = field(default=True)
    flank: int = field(default=200, metadata={"validate": validate.Range(min=0)})
    simple_overlap_for_monoexonic: bool = field(default=False)


@dataclass
class FragmentsConfiguration:
    remove: bool = field(default=True)
    max_distance: int = field(default=2000, metadata={"validate": validate.Range(min=0)})
    valid_class_codes: List[str] = field(
        default_factory=lambda: ["p", "P", "x", "X", "i", "m", "_", "e", "o"],
        metadata={"validate": validate.ContainsOnly(["p", "P", "i", "I", "ri", "rI", "x", "X", "m", "_", "e", "o"])})


@dataclass
class FilesConfiguration:
    output_dir: str = field(default="")
    input: str = field(default="mikado_prepared.gtf")
    loci_out: str = field(default="mikado.loci.gff3")
    subloci_out: Optional[str] = field(default=None)
    monoloci_out: Optional[str] = field(default=None)
    log: str = field(default="pick.log")


@dataclass
class PickConfiguration:
    scoring_file: str = field(default="plant.yaml")
    alternative_splicing: AlternativeSplicingConfiguration = field(default_factory=AlternativeSplicingConfiguration)
    output_format: OutputFormatConfiguration = field(default_factory=OutputFormatConfiguration)
    orf_loading: OrfLoadingConfiguration = field(default_factory=OrfLoadingConfiguration)
    chimera_split: ChimeraSplitConfiguration = field(default_factory=ChimeraSplitConfiguration)
    run_options: RunOptionsConfiguration = field(default_factory=RunOptionsConfiguration)
    clustering: ClusteringConfiguration = field(default_factory=ClusteringConfiguration)
    fragments: FragmentsConfiguration = field(default_factory=FragmentsConfiguration)
    files: FilesConfiguration = field(default_factory=FilesConfiguration)
