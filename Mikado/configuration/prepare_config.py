from marshmallow import Schema, fields


class PrepareFilesConfiguration(Schema):
    output_dir: str = fields.Str(missing="")
    out: str = fields.Str(missing="mikado_prepared.gtf")
    out_fasta: str = fields.Str(missing="mikado_prepared.fasta")
    log: str = fields.Str(missing="prepare.log")
    gff: list = fields.List(fields.Str(), missing=list)
    labels: list = fields.List(fields.Str(), missing=list)
    strand_specific_assemblies: list = fields.List(fields.Str(), missing=list)
    reference: list = fields.List(fields.Str(), missing=list)
    exclude_redundant: list = fields.List(fields.Str(), missing=list)
    strip_cds: list = fields.List(fields.Str(), missing=list)
    source_score: dict = fields.Dict(missing=dict)


class PrepareConfiguration(Schema):
    exclude_redundant: bool = fields.Bool(missing=False)
    minimum_cdna_length: int = fields.Int(missing=200)
    max_intron_length: int = fields.Int(missing=1000000)
    strip_cds: bool = fields.Bool(missing=False)
    strip_faulty_cds: bool = fields.Bool(missing=False)
    single: bool = fields.Bool(missing=False)
    lenient: bool = fields.Bool(missing=False)
    strand_specific: bool = fields.Bool(missing=False)
    canonical: list = fields.List(fields.List(fields.Str()), missing=lambda: [["GT", "AG"], ["GC", "AG"], ["AT", "AC"]])
    files: PrepareFilesConfiguration = fields.Nested(PrepareFilesConfiguration)
