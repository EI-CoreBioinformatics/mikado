from marshmallow_dataclass import dataclass
from marshmallow import Schema, fields


@dataclass
class AlternativeSplicingConfiguration(Schema):
    report: bool = fields.Bool(missing=True)
    cds_only: bool = fields.Bool(missing=False)
    min_cds_overlap: float = fields.Float(missing=0.5)
    min_cdna_overlap: float = fields.Float(missing=0.6)
    keep_retained_introns: bool = fields.Bool(missing=True)
    keep_cds_disrupted_by_ri: bool = fields.Bool(missing=False)
    max_isoforms: int = fields.Int(missing=10)
    valid_ccodes: list = fields.List(fields.Str(), missing=lambda: ["j", "J", "G", "h"])
    redundant_ccodes: list = fields.List(fields.Str(), missing=lambda: ["c", "m", "_", "=", "n"])
    min_score_perc: float = fields.Float(missing=0.5)
    only_confirmed_introns: bool = fields.Bool(missing=True)
    ts_distance: int = fields.Int(missing=2000)
    pad: bool = fields.Bool(missing=True)
    ts_max_splices: int = fields.Int(missing=2)


@dataclass
class OutputFormatConfiguration(Schema):
    source: str = fields.Str(missing="Mikado")
    id_prefix: str = fields.Str(missing="mikado")
    report_all_orfs: bool = fields.Bool(missing=False)


@dataclass
class OrfLoadingConfiguration(Schema):
    minimal_secondary_orf_length: int = fields.Int(missing=200)
    minimal_orf_length: int = fields.Int(missing=50)
    strand_specific: bool = fields.Bool(missing=True)


@dataclass
class BlastParamsConfiguration(Schema):
    evalue: float = fields.Float(missing=1e-06)
    hsp_evalue: float = fields.Float(missing=1e-06)
    leniency: str = fields.Str(missing="STRINGENT")
    max_target_seqs: int = fields.Int(missing=3)
    minimal_hsp_overlap: float = fields.Float(missing=0.8)
    min_overlap_duplication: float = fields.Float(missing=0.8)


@dataclass
class ChimeraSplitConfiguration(Schema):
    blast_check: bool = fields.Bool(missing=True)
    execute: bool = fields.Bool(missing=True)
    skip: list = fields.List(fields.Str(), missing=list)
    blast_params: BlastParamsConfiguration = fields.Nested(BlastParamsConfiguration)


@dataclass
class RunOptionsConfiguration(Schema):
    shm: bool = fields.Bool(missing=False)
    exclude_cds: bool = fields.Bool(missing=False)
    # TODO this must be at most 2
    intron_range: list = fields.List(fields.Int(), missing=lambda: [60, 10000])
    reference_update: bool = fields.Bool(missing=False)
    only_reference_update: bool = fields.Bool(missing=False)
    check_references: bool = fields.Bool(missing=False)
    single_thread: bool = fields.Bool(missing=False)


@dataclass
class ClusteringConfiguration(Schema):
    cds_only: bool = fields.Bool(missing=False)
    min_cds_overlap: float = fields.Float(missing=0.2)
    min_cdna_overlap: float = fields.Float(missing=0.2)
    purge: bool = fields.Bool(missing=True)
    flank: int = fields.Int(missing=200)
    simple_overlap_for_monoexonic: bool = fields.Bool(missing=False)


@dataclass
class FragmentsConfiguration(Schema):
    remove: bool = fields.Bool(missing=True)
    max_distance: int = fields.Int(missing=2000)
    valid_class_codes: list = fields.List(fields.Str(), missing=lambda: ["p", "P", "x", "X", "i", "m", "_", "e", "o"])


@dataclass
class FilesConfiguration(Schema):
    output_dir: str = fields.Str(missing="")
    input: str = fields.Str(missing="mikado_prepared.gtf")
    loci_out: str = fields.Str(missing="mikado.loci.gff3")
    subloci_out: str = fields.Str(missing="")
    monoloci_out: str = fields.Str(missing="")
    log: str = fields.Str(missing="pick.log")


@dataclass
class PickConfiguration(Schema):
    scoring_file: str = fields.Str(missing="plant.yaml")
    alternative_splicing: AlternativeSplicingConfiguration = fields.Nested(AlternativeSplicingConfiguration)
    output_format: OutputFormatConfiguration = fields.Nested(OutputFormatConfiguration)
    orf_loading: OrfLoadingConfiguration = fields.Nested(OrfLoadingConfiguration)
    chimera_split: ChimeraSplitConfiguration = fields.Nested(ChimeraSplitConfiguration)
    run_options: RunOptionsConfiguration = fields.Nested(RunOptionsConfiguration)
    clustering: ClusteringConfiguration = fields.Nested(ClusteringConfiguration)
    fragments: FragmentsConfiguration = fields.Nested(FragmentsConfiguration)
    files: FilesConfiguration = fields.Nested(FilesConfiguration)
