from dataclasses import dataclass, field


@dataclass
class AlternativeSplicingConfiguration:
    report: bool = True
    cds_only: bool = False
    min_cds_overlap: float = 0.5
    min_cdna_overlap: float = 0.6
    keep_retained_introns: bool = True
    keep_cds_disrupted_by_ri: bool = False
    max_isoforms: int = 10
    valid_ccodes: list = field(default=["j", "J", "G", "h"])
    redundant_ccodes: list = field(default=["c", "m", "_", "=", "n"])
    min_score_perc: float = 0.5
    only_confirmed_introns: bool = True
    ts_distance: int = 2000
    pad: bool = True
    ts_max_splices: int = 2


@dataclass
class OutputFormatConfiguration:
    source: str = "Mikado"
    id_prefix: str = "mikado"
    report_all_orfs: bool = False


@dataclass
class OrfLoadingConfiguration:
    minimal_secondary_orf_length: int = 200
    minimal_orf_length: int = 50
    strand_specific: bool = True


@dataclass
class BlastParamsConfiguration:
    evalue: float = 1e-06
    hsp_evalue: float = 1e-06
    leniency: str = "STRINGENT"
    max_target_seqs: int = 3
    minimal_hsp_overlap: float = 0.8
    min_overlap_duplication: float = 0.8


@dataclass
class ChimeraSplitConfiguration:
    blast_check: bool = True
    execute: bool = True
    skip: list = field(default_factory=list)
    blast_params: BlastParamsConfiguration = BlastParamsConfiguration()


@dataclass
class RunOptionsConfiguration:
    shm: bool = False
    exclude_cds: bool = False
    intron_range: list = field(default=[60, 10000])
    reference_update: bool = False
    only_reference_update: bool = False
    check_references: bool = False
    single_thread: bool = False


@dataclass
class ClusteringConfiguration:
    cds_only: bool = False
    min_cds_overlap: float = 0.2
    min_cdna_overlap: float = 0.2
    purge: bool = True
    flank: int = 200
    simple_overlap_for_monoexonic: bool = False


@dataclass
class FragmentsConfiguration:
    remove: bool = True
    max_distance: int = 2000
    valid_ccodes: list = field(default=["p", "P", "x", "X", "i", "m", "_", "e", "o"])


@dataclass
class FilesConfiguration:
    output_dir: str = "."
    input: str = "mikado_prepared.gtf"
    loci_out: str = "mikado.loci.gff3"
    subloci_out: str = ""
    monoloci_out: str = ""
    log: str = "pick.log"


@dataclass
class PickConfiguration:
    scoring_file: str = "plant.yaml"
    alternative_splicing: AlternativeSplicingConfiguration = AlternativeSplicingConfiguration()
    output_format: OutputFormatConfiguration = OutputFormatConfiguration()
    orf_loading: OrfLoadingConfiguration = OrfLoadingConfiguration()
    chimera_split: ChimeraSplitConfiguration = ChimeraSplitConfiguration()
    run_options: RunOptionsConfiguration = RunOptionsConfiguration()
    clustering_options: ClusteringConfiguration = ClusteringConfiguration()
    fragments: FragmentsConfiguration = FragmentsConfiguration()
    files: FilesConfiguration = FilesConfiguration()
