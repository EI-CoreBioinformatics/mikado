from dataclasses import dataclass, field

from marshmallow import Schema, fields


class AlternativeSplicingConfiguration(Schema):
    report: fields.Bool(missing=True)
    cds_only: fields.Bool(missing=False)
    min_cds_overlap: float = fields.Float(missing=0.5)
    min_cdna_overlap: float = fields.Float(missing=0.6)
    keep_retained_introns: fields.Bool(missing=True)
    keep_cds_disrupted_by_ri: fields.Bool(missing=False)
    max_isoforms: int = 10
    valid_ccodes: list = field(default_factory=lambda: ["j", "J", "G", "h"])
    redundant_ccodes: list = field(default_factory=lambda: ["c", "m", "_", "=", "n"])
    min_score_perc: float = fields.Float(missing=0.5)
    only_confirmed_introns: fields.Bool(missing=True)
    ts_distance: int = 2000
    pad: fields.Bool(missing=True)
    ts_max_splices: int = 2


class OutputFormatConfiguration(Schema):
    source: fields.Str(missing="Mikado")
    id_prefix: str = fields.Str(missing="mikado")
    report_all_orfs: fields.Bool(missing=False)


class OrfLoadingConfiguration(Schema):
    minimal_secondary_orf_length: int = 200
    minimal_orf_length: int = 50
    strand_specific: bool = fields.Bool(missing=True)


class BlastParamsConfiguration(Schema):
    evalue: float = fields.Float(missing=1e-06)
    hsp_evalue: float = fields.Float(missing=1e-06)
    leniency: str = fields.Str(missing="STRINGENT")
    max_target_seqs: int = fields.Int(missing=3)
    minimal_hsp_overlap: float = fields.Float(missing=0.8)
    min_overlap_duplication: float = fields.Float(missing=0.8)


class ChimeraSplitConfiguration(Schema):
    blast_check: fields.Bool(missing=True)
    execute: fields.Bool(missing=True)
    skip: list = fields.List(fields.Str(), missing=list)
    blast_params: BlastParamsConfiguration = fields.Nested(BlastParamsConfiguration)


class RunOptionsConfiguration(Schema):
    shm: fields.Bool(missing=False)
    exclude_cds: fields.Bool(missing=False)
    intron_range: list = field(default_factory=lambda: [60, 10000])
    reference_update: fields.Bool(missing=False)
    only_reference_update: fields.Bool(missing=False)
    check_references: fields.Bool(missing=False)
    single_thread: fields.Bool(missing=False)


class ClusteringConfiguration(Schema):
    cds_only: fields.Bool(missing=False)
    min_cds_overlap: float = fields.Float(missing=0.2)
    min_cdna_overlap: float = fields.Float(missing=0.2)
    purge: fields.Bool(missing=True)
    flank: int = fields.Int(missing=200)
    simple_overlap_for_monoexonic: fields.Bool(missing=False)


class FragmentsConfiguration(Schema):
    remove: fields.Bool(missing=True)
    max_distance: int = 2000
    valid_class_codes: list = field(default_factory=lambda: ["p", "P", "x", "X", "i", "m", "_", "e", "o"])


class FilesConfiguration(Schema):
    output_dir: str = fields.Str(missing="../picking")
    input: str = fields.Str(missing="mikado_prepared.gtf")
    loci_out: str = fields.Str(missing="mikado.loci.gff3")
    subloci_out: str = fields.Str(missing="")
    monoloci_out: str = fields.Str(missing="")
    log: str = fields.Str(missing="pick.log")


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
