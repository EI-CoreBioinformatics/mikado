from dataclasses import field
from typing import List

from marshmallow import validate
from marshmallow_dataclass import dataclass


@dataclass
class PrepareFilesConfiguration:
    output_dir: str = field(default="", metadata={
        "name": "output_dir",
        "description": "Output folder."
    })
    out: str = field(default="mikado_prepared.gtf", metadata={
        "name": "out",
        "description": "Output GTF file.",
        "required": True,
    })
    out_fasta: str = field(default="mikado_prepared.fasta", metadata={
        "name": "out_fasta",
        "default": "mikado_prepared.fasta",
        "required": True
    })
    log: str = field(default="prepare.log", metadata={
        "name": "log",
        "description": "Log file."
    })
    gff: List[str] = field(default_factory=lambda: [], metadata={
        "name": "gff",
        "description": "List of the input files.",
        "required": True
    })
    labels: List[str] = field(default_factory=lambda: [], metadata={
        "name": "labels",
        "description": "List of labels associated with each input GFF. This list must *either* be empty *or* be the same length as the 'gff' content above.",
    })
    strand_specific_assemblies: List[str] = field(default_factory=lambda: [], metadata={
        "name": "strand_specific_assemblies",
        "description": "List of input assemblies to be considered as strand specific. Any 'reference' input is automatically marked as strand-specific.",
    })
    reference: List[bool] = field(default_factory=lambda: [], metadata={
        "name": "reference",
        "description": "List of input assemblies to be considered as of 'reference' quality. Transcripts from this list will be excluded only if they have glaring mistakes eg. an incorrect CDS.",
    })
    exclude_redundant: List[bool] = field(default_factory=lambda: [], metadata={
        "name": "exclude_redundant",
        "description": "If the 'exclude_redundant' switch is set to true, this list specifies which assemblies can have the redundancy check performed on.",
    })
    strip_cds: List[bool] = field(default_factory=lambda: [], metadata={
        "name": "strip_cds",
        "description": "List of input assemblies from which the CDS will be stripped. Useful to e.g. remove the CDS from GMAP alignments.",
    })
    source_score: dict = field(default_factory=dict, metadata={
        "name": "source_score",
        "description": "Dictionary linking assemblies with an optional score, for prioritising them during 'pick' (or de-prioritising if the score is negative).",
    })


@dataclass
class PrepareConfiguration:
    exclude_redundant: bool = field(default=False, metadata={
        "name": "exclude_redundant",
        "description": "Boolean. If set to True, fully redundant transcripts across the input files will be removed. Higher scoring transcripts *on the basis of the score associated to a given input file* will be preferentially retained."
    })
    minimum_cdna_length: int = field(default=200, metadata={
        "name": "minimum_cdna_length",
        "description": "Minimum length of a transcript to be retained. Default: 200 bps",
        "validate": validate.Range(min=1),
        "required": True
    })
    max_intron_length: int = field(default=1000000, metadata={
        "name": "max_intron_length",
        "description": "Maximum length of an intron. Transcripts with introns bigger than this value will be split in various sub-transcripts. Default: 1,000,000 bps.",
        "validate": validate.Range(min=20),
        "required": True
    })
    strip_cds: bool = field(default=False, metadata={
        "name": "strip_cds",
        "description": "Boolean flag. If set, the CDS will be stripped from any non-reference input assembly.",
    })
    strip_faulty_cds: bool = field(default=False, metadata={
        "name": "strip_faulty_cds",
        "description": "Boolean flag. If set to false, any transcript - *including reference transcripts* found to have an incorrect CDS will be discarded. If set to to true, these transcripts will be retained but their CDS will be stripped.",
    })
    single: bool = field(default=False, metadata={
        "name": "single",
        "description": "Boolean flag. If set to true, Mikado will run in single-threaded mode, useful for debugging.",
    })
    lenient: bool = field(default=False, metadata={
        "name": "lenient",
        "description": "Boolean flag. If set to true, Mikado will retain transcripts with no canonical junction or with canonical junctions on both strands. If set to false (default), such transcripts will instead be discarded.",
    })
    strand_specific: bool = field(default=False, metadata={
        "name": "strand_specific",
        "description": "Boolean flag. If set to true, all assemblies will be considered as strand-specific. By default, Mikado will consider the strand-specificity of each assembly in isolation, see 'files/strand_specific_assemblies'.",
        "required": True
    })
    canonical: List[List[str]] = field(default_factory=lambda: [["GT", "AG"], ["GC", "AG"], ["AT", "AC"]], metadata={
        "name": "canonical",
        "description": "Accepted canonical splicing junctions for the organism in examination.",
    })
    files: PrepareFilesConfiguration = field(default_factory=PrepareFilesConfiguration, metadata={
        "name": "files",
        "description": "Options related to the input and output files.",
        "required": True
    })
