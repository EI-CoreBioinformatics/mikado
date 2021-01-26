from marshmallow_dataclass import dataclass, Optional
from dataclasses import field
from marshmallow import validate
from typing import List


@dataclass
class AlternativeSplicingConfiguration:
    report: bool = field(default=True, metadata={
        "description": "Boolean flag, about whether Mikado should find and report AS events or simply define one transcript per locus."
    })
    cds_only: bool = field(default=False, metadata={
        "description": "Only consider the coding part of transcripts to determine the validity of the AS event.",
    })
    min_cds_overlap: float = field(default=0.5, metadata={
        "description": "Minimum CDS overlap threshold (for coding transcripts) to be cleared for two transcripts to be considered AS events of each other.",
        "validate": validate.Range(min=0, max=1),
    })
    min_cdna_overlap: float = field(default=0.6, metadata={
        "description": "Minimum cDNA overlap threshold to be cleared for two transcripts to be considered AS events of each other.",
        "validate": validate.Range(min=0, max=1)
    })
    keep_retained_introns: bool = field(default=True, metadata={
        "description": "Keep or discard AS events with a retained intron. Default: true",
    })
    keep_cds_disrupted_by_ri: bool = field(default=False, metadata={
        "description": "Keep or discard AS events with their CDS disrupted by a retained intron event, ie either having their stop codon or ending with a truncated CDS within the intron of another transcript. Default: false.",
    })
    max_isoforms: int = field(default=10, metadata={
        "description": "Maximum number of isoforms to report per locus. Default: 10.",
        "validate": validate.Range(min=1)
    })
    valid_ccodes: list = field(default_factory=lambda: ["j", "J", "G", "h"], metadata={
        "description": "AS event class codes considered as valid AS events. Valid codes are in categories 'Alternative splicing', 'Extension' (with junction F1 lower than 100%), and 'Overlap' (exluding m). Please run 'mikado util class_codes' or refer to the online documentation for an explanation of each code.",
        "validate": validate.ContainsOnly(["j", "e", "o", "h", "J", "C", "g", "G", "=", "n", "_"])})
    redundant_ccodes: list = field(default_factory=lambda: ["c", "m", "_", "=", "n"], metadata={
        "description": "AS event class codes considered as a duplicate of other transcripts in the locus. Please run 'mikado util class_codes' or refer to the online documentation for an explanation of each code.",
        "validate": validate.ContainsOnly(
            ["j", "n", "O", "e", "o", "h", "J", "C", "c", "m", "mo", "=", "_", "x", "p", "P", "X", "I", "i"])
    })
    min_score_perc: float = field(default=0.5, metadata={
        "description": "Minimum percentage of the score associated to an AS event *compared to the primary transcript* for the AS event to be considered valid. Default: 0.5, or 50%.",
        "validate": validate.Range(min=0, max=1)
    })
    only_confirmed_introns: bool = field(default=True, metadata={
        "description": "Boolean flag. If set to true (default), Mikado will only report AS events whose introns *not in common with the primary transcript* are verified by the junctions provided to serialise (usually Portcullis reliable junctions).",
    })
    ts_distance: int = field(default=2000, metadata={
        "description": "When padding, this value indicates how many bps can be added to *either* of the 5' or 3' section of the transcript, excluding introns.",
        "validate": validate.Range(min=0)
    })
    pad: bool = field(default=True, metadata={
        "description": "Boolean flag. If set to true, Mikado will pad transcripts. Please refer to the online documentation.",
    })
    ts_max_splices: int = field(default=2, metadata={
        "description": "When padding, this value indicates the maximum number of splicing junctions that can be added to *either* of the 5' or 3' section of the transcript.",
        "validate": validate.Range(min=0)
    })


@dataclass
class OutputFormatConfiguration:
    source: str = field(default="Mikado", metadata={
        "description": "Prefix for the source field in the mikado output.",
    })
    id_prefix: str = field(default="mikado", metadata={
        "description": "Prefix for the ID of the genes/transcripts in the output",
    })
    report_all_orfs: bool = field(default=False, metadata={
        "description": "Boolean switch. If set to true, Mikado will report all ORFs associated with a transcript in the final loci file.",
    })


@dataclass
class OrfLoadingConfiguration:
    minimal_secondary_orf_length: int = field(default=200, metadata={
        "description": "Minimum length of a *secondary* ORF to be loaded after the first, in bp. Default: 200 bps",
        "validate": validate.Range(min=0)
    })
    minimal_orf_length: int = field(default=50, metadata={
        "description": "Minimum length in bps of an ORF to be loaded, as the primary ORF, onto a transcript. Default: 50 bps",
        "validate": validate.Range(min=0)
    })
    strand_specific: bool = field(default=True, metadata={
        "description": "Boolean flag. If set to true, monoexonic transcripts with an available ORF on the opposite strand will still not be reversed.",
    })


@dataclass
class BlastParamsConfiguration:
    evalue: float = field(default=1e-06, metadata={
        "description": "Minimum evalue for the whole hit. Default: 1e-6",
        "validate": validate.Range(min=0)
    })
    hsp_evalue: float = field(default=1e-06, metadata={
        "description": "Minimum evalue for any HSP hit (some might be discarded even if the whole hit is valid). Default: 1e-6",
        "validate": validate.Range(min=0)
    })
    leniency: str = field(default="STRINGENT", metadata={
        "description": "One of 'STRINGENT', 'LENIENT', 'PERMISSIVE'. Please refer to the online documentation for details. Default: STRINGENT",
        "validate": validate.OneOf(["STRINGENT", "LENIENT", "PERMISSIVE"])})
    max_target_seqs: int = field(default=3, metadata={
        "description": "Maximum number of hits to consider. Default: 3",
        "validate": validate.Range(min=1)
    })
    minimal_hsp_overlap: float = field(default=0.5, metadata={
        "description": "Minimum overlap of the ORF with the HSP (*not* reciprocal). Default: 0.8, i.e. 80%",
        "validate": validate.Range(min=0, max=1)
    })
    min_overlap_duplication: float = field(default=0.8, metadata={
        "description": "min_overlap_duplication: minimum overlap (in %) for two ORFs to consider them as target duplications. This means that if two ORFs have no HSPs in common, but the coverage of their disjoint HSPs covers more than this percentage of the length of the *target*, they represent most probably a duplicated gene.",
        "validate": validate.Range(min=0, max=1)
    })


@dataclass
class ChimeraSplitConfiguration:
    blast_check: bool = field(default=True, metadata={
        "description": "Whether to use BLAST information to take a decision. See blast_params for details.",
    })
    execute: bool = field(default=True, metadata={
        "description": "Whether to split multi-ORF transcripts at all. Boolean.",
    })
    skip: List[bool] = field(default_factory=lambda: [], metadata={
        "description": "Input sources for which Mikado will skip the splitting, e.g. ultra-reliable full cDNA sequences.",
        "validate": validate.Length(min=0)})
    blast_params: BlastParamsConfiguration = field(default_factory=BlastParamsConfiguration, metadata={
        "description": "Parameters for the BLAST check prior to splitting.",
    })



@dataclass
class RunOptionsConfiguration:
    # From swagger-marshmallow-codegen (Unique validator)
    class Unique(validate.Validator):
        message = "{input} is Not unique"

        def __init__(self, error=None):
            self.error = error

        def _repr_args(self):
            return ""

        def _format_error(self, value):
            return self.message.format(input=value)

        def __call__(self, value):
            if len(value) != len(set(value)):
                raise v.ValidationError(self._format_error(value))
            return value

    shm: bool = field(default=False, metadata={
        "description": "boolean flag. If set and the DB is sqlite, it will be copied onto the /dev/shm faux partition, for a potentially faster execution.",
    })
    exclude_cds: bool = field(default=False, metadata={
        "description": "boolean flag. If set, the CDS information will not be printed in Mikado output. Default: false",
    })
    intron_range: List[int] = field(default_factory=lambda: [60, 10000], metadata={
        "description": "A range where most of the introns (99%) should fall into. Transcripts with too many introns larger or smaller than what is defined in this range will be penalised in the scoring. Default: [60, 900]",
        "validate": [validate.Length(min=2, max=2), Unique]
    })
    reference_update: bool = field(default=False, metadata={
        "description": "Boolean flag. If set, Mikado will run in reference-update mode, see documentation.",
    })
    only_reference_update: bool = field(default=False, metadata={
        "description": "Boolean flag. If set, Mikado will run in reference-update mode, see documentation. Additionally, Mikado will ignore any locus where there is not at least one reference transcript.",
    })
    check_references: bool = field(default=False, metadata={
        "description": "boolean flag. If set to true, transcripts marked as reference will still be checked for compliance with the requirements established in the scoring file.",
    })
    single_thread: bool = field(default=False, metadata={
        "description": "Boolean flag. If set, multithreading will be disabled - useful for profiling and debugging.",
    })


@dataclass
class ClusteringConfiguration:
    cds_only: bool = field(default=False, metadata={
        "description": "Boolean, it specifies whether to cluster transcripts only according to their CDS (if present).",
    })
    min_cds_overlap: float = field(default=0.2, metadata={
        "description": "Minimal CDS overlap for the second clustering, in percentage between 0 and 1. Default: 0.2, or 20%",
        "validate": validate.Range(min=0, max=1)
    })
    min_cdna_overlap: float = field(default=0.2, metadata={
        "description": "Minimal cDNA overlap for the second clustering, in percentage between 0 and 1. Default: 0.2, or 20%.",
        "validate": validate.Range(min=0, max=1)
    })
    purge: bool = field(default=True, metadata={
        "description": "Boolean, it specifies whether to remove transcripts which fail the minimum requirements check, or if instead to just assign them a score of 0 (potentially retaining them in the final output).",
    })
    flank: int = field(default=200, metadata={
        "description": "Maximum distance for transcripts to be clustered within the same superlocus.",
        "validate": validate.Range(min=0)})
    simple_overlap_for_monoexonic: bool = field(default=False, metadata={
        "description": "boolean. Disabled by default. If set to true, then any overlap, even minimal, will suffice to incude monoexonic transcripts in a locus.",
    })


@dataclass
class FragmentsConfiguration:
    remove: bool = field(default=True, metadata={
        "description": "boolean. Whether to remove fragments or leave them, properly tagged, in the output file. Default: remove them."
    })
    max_distance: int = field(default=2000, metadata={
        "description": "Maximum distance of a putative fragment from a valid gene, for it to be considered by this filter.",
        "validate": validate.Range(min=0)
    })
    valid_class_codes: List[str] = field(
        default_factory=lambda: ["p", "P", "x", "X", "i", "m", "_", "e", "o"], metadata={
            "description": "Which class codes will be considered as fragments. Default: (p, P, x, X, i, m, _). Choices: '_' plus any class code with category 'Intronic', 'Fragment', or 'Overlap'. Please refer to the online documentation or run 'mikado util class_codes for details.",
            "validate": validate.ContainsOnly(["p", "P", "i", "I", "ri", "rI", "x", "X", "m", "_", "e", "o"])
    })


@dataclass
class FilesConfiguration:
    output_dir: str = field(default="", metadata={
        "description": "Output directory for mikado pick",
    })
    input: str = field(default="mikado_prepared.gtf", metadata={
        "description": "Input GTF/GFF3/BED12 file. Default: mikado_prepared.gtf",
    })
    loci_out: str = field(default="mikado.loci.gff3", metadata={
        "description": "Main output GFF3 file from Mikado pick. Default: mikado.loci.gff3",
    })
    subloci_out: Optional[str] = field(default=None, metadata={
        "description": "Optional GFF file with the intermediate subloci. Default: no output",
    })
    monoloci_out: Optional[str] = field(default=None, metadata={
        "description": "optional GFF file with the intermediate monoloci. Default: no output",
    })
    log: str = field(default="pick.log", metadata={
        "description": "Log file for mikado pick."
    })


@dataclass
class PickConfiguration:
    scoring_file: str = field(default="plant.yaml", metadata={
        "description": "Scoring file to be used by Mikado.",
    })
    alternative_splicing: AlternativeSplicingConfiguration = field(default_factory=AlternativeSplicingConfiguration, metadata={
        "description": "Parameters related to how Mikado will select and report alternative splicing events.",
    })
    output_format: OutputFormatConfiguration = field(default_factory=OutputFormatConfiguration, metadata={
        "description": "Parameters related to the output format.",
    })
    orf_loading: OrfLoadingConfiguration = field(default_factory=OrfLoadingConfiguration, metadata={
        "description": "Parameters related to ORF loading.",
    })
    chimera_split: ChimeraSplitConfiguration = field(default_factory=ChimeraSplitConfiguration, metadata={
        "description": "Parameters related to the splitting of transcripts in the presence of two or more ORFs.",
    })
    run_options: RunOptionsConfiguration = field(default_factory=RunOptionsConfiguration, metadata={
        "description": "Generic run options for Mikado pick.",
    })
    clustering: ClusteringConfiguration = field(default_factory=ClusteringConfiguration, metadata={
        "description": "Parameters related to the clustering of transcripts into loci.",
    })
    fragments: FragmentsConfiguration = field(default_factory=FragmentsConfiguration, metadata={
        "description": "Parameters related to the handling of fragments."
    })
    files: FilesConfiguration = field(default_factory=FilesConfiguration, metadata={
        "description": "Input and output files for Mikado pick.",
    })
