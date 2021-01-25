from dataclasses import dataclass, field
from typing import Union

from marshmallow import fields, Schema
import marshmallow_union


class FilesConfiguration(Schema):
    junctions: list = fields.List(fields.Str(), missing=list)
    xml: list = fields.List(fields.Str(), missing=list)
    blast_loading_debug: bool = fields.Bool(missing=False)
    external_scores: str = fields.Str(missing="")
    orfs: list = fields.List(fields.Str(), missing=list)
    transcripts: str = fields.Str(missing="mikado_prepared.fasta")
    log: str = fields.Str(missing="serialise.log")
    blast_targets: list = fields.List(fields.Str(), missing=list)
    output_dir: str = fields.Str(missing="../serializers")


class SerialiseConfiguration(Schema):
    files: FilesConfiguration = fields.Nested(FilesConfiguration)
    substitution_matrix: str = fields.Str(missing="blosum62")
    blast_flavour: str = fields.Str(missing="blastx")
    codon_table: Union[str, int] = marshmallow_union.Union(fields=[fields.Str(), fields.Int()], missing=0)
    max_objects: int = fields.Int(missing=10000000)
    max_regression: float = fields.Float(missing=0.2)
    start_adjustment: bool = fields.Bool(missing=True)
    max_target_seqs: int = fields.Int(missing=100000)
    force: bool = fields.Bool(missing=False)
    single_thread: bool = fields.Bool(missing=False)
