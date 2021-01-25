from dataclasses import dataclass, field
from typing import Union, List, Optional

from marshmallow import validate


@dataclass
class FilesConfiguration:
    junctions: List[str] = field(default_factory=lambda: [""])
    xml: List[str] = field(default_factory=lambda: [""])
    blast_loading_debug: bool = field(default=False)
    external_scores: str = field(default="")
    orfs: List[str] = field(default_factory=lambda: [""])
    transcripts: str = field(default="mikado_prepared.fasta", metadata={"required": True})
    log: str = field(default="serialise.log")
    blast_targets: List[str] = field(default_factory=lambda: [""])
    output_dir: str = field(default="")


@dataclass
class SerialiseConfiguration:
    files: FilesConfiguration = field(default_factory=FilesConfiguration)
    substitution_matrix: Optional[str] = field(default="blosum62",
                                               metadata={"required": True,
                                                         "validate": validate.OneOf(choices=[
                                                             "blosum45", "blosum50", "blosum62",
                                                             "blosum80", "blosum90", "pam250", "pam30", "pam70"])
                                                         })
    blast_flavour: str = field(default="blastx", metadata={"validate": validate.OneOf(choices=["blastx", "blastp"])})
    codon_table: Union[str, int] = field(default=0, metadata={"required": True,
                                                              "validate": validate.OneOf(
                                                                  choices=[0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14,
                                                                           15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                                                           30, 31,
                                                                           "Standard", "SGC0",
                                                                           "Vertebrate Mitochondrial", "SGC1",
                                                                           "Yeast Mitochondrial", "SGC2",
                                                                           "Mold Mitochondrial",
                                                                           "Protozoan Mitochondrial",
                                                                           "Coelenterate Mitochondrial", "Mycoplasma",
                                                                           "Spiroplasma", "SGC3",
                                                                           "Invertebrate Mitochondrial", "SGC4",
                                                                           "Ciliate Nuclear", "Dasycladacean Nuclear",
                                                                           "Hexamita Nuclear", "SGC5",
                                                                           "Echinoderm Mitochondrial",
                                                                           "Flatworm Mitochondrial",
                                                                           "SGC8", "Euplotid Nuclear", "SGC9",
                                                                           "Bacterial", "Archaeal", "Plant Plastid",
                                                                           "Alternative Yeast Nuclear",
                                                                           "Ascidian Mitochondrial",
                                                                           "Alternative Flatworm Mitochondrial",
                                                                           "Blepharisma Macronuclear",
                                                                           "Chlorophycean Mitochondrial",
                                                                           "Trematode Mitochondrial",
                                                                           "Scenedesmus obliquus Mitochondrial",
                                                                           "Thraustochytrium Mitochondrial",
                                                                           "Pterobranchia Mitochondrial",
                                                                           "Candidate Division SR1", "Gracilibacteria",
                                                                           "Pachysolen tannophilus Nuclear",
                                                                           "Karyorelict Nuclear",
                                                                           "Condylostoma Nuclear",
                                                                           "Mesodinium Nuclear", "Peritrich Nuclear",
                                                                           "Blastocrithidia Nuclear"])})

    max_objects: int = field(default=10000000, metadata={"validate": validate.Range(min=1)})
    max_regression: float = field(default=0.2,
                                  metadata={"validate": validate.Range(min=0.0, max=1.0), "required": True})
    start_adjustment: bool = field(default=True)
    max_target_seqs: int = field(default=100000, metadata={"validate": validate.Range(min=1)})
    force: bool = field(default=False)
    single_thread: bool = field(default=False)
