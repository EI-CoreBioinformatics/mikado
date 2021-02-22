from dataclasses import field

from marshmallow import validate
from marshmallow_dataclass import dataclass, List, Optional, Union


@dataclass
class FilesConfiguration:
    junctions: List[str] = field(default_factory=lambda: [], metadata={
                "metadata": {"description": "File of reliable junctions (from e.g. Portcullis), in BED12 format"},
    })
    xml: List[str] = field(default_factory=lambda: [], metadata={
                "metadata": {
                    "description": "BLAST file or folder of files, in XML or tabular format, optionally compressed."},
    })
    blast_loading_debug: bool = field(default=False, metadata={
                "metadata": {
                    "description": "Boolean flag. If True, the loading of BLAST files will happen \
in a single thread and debugging mode will be activated."},
    })
    external_scores: str = field(default="", metadata={
                "metadata": {
                    "description": "File of additional scores related to the input transcripts. \
Please see the documentation."},
    })
    orfs: List[str] = field(default_factory=lambda: [], metadata={
                "metadata": {
                    "description": "File(s) containing the ORF calling for the input transcripts, in GFF3 or BED12 \
*transcriptomic* (ie in cDNA rather than genomic) coordinates."},
    })
    transcripts: str = field(default="mikado_prepared.fasta", metadata={
                "metadata": {"description": "Input transcripts in FASTA format, ie the output of Mikado prepare."},
                "required": True
    })
    log: str = field(default="serialise.log", metadata={"metadata": {"description": "log file."}})
    blast_targets: List[str] = field(default_factory=lambda: [], metadata={
                "metadata": {"description": "FASTA file(s) with the BLAST targets."},
    })
    output_dir: str = field(default=".", metadata={"validate": validate.Length(min=1)})


@dataclass
class SerialiseConfiguration:
    files: Optional[FilesConfiguration] = field(default_factory=FilesConfiguration, metadata={
                "metadata": {"description": "Options related to input files for serialise"},
    })
    substitution_matrix: str = field(default="blosum62", metadata={
                "metadata": {
                    "description": "Substitution matrix used for the BLAST. This value will be derived from the XML \
files, but it must be provided here or on the command line when using BLAST tabular data. \
Default: blosum62, the default for both BLAST and DIAMOND."},
                "required": True,
                "validate": validate.OneOf(choices=["blosum45", "blosum50", "blosum62",
                                                    "blosum80", "blosum90", "pam250", "pam30", "pam70"])
    })
    blast_flavour: str = field(default="blastx", metadata={
                "metadata": {"description": "Type of BLAST used. Either BLASTP or BLASTX. Default: blastx, \
which should be the sane presumption in most instances."},
                "validate": validate.OneOf(choices=["blastx", "blastp"])
    })
    codon_table: Union[int, str] = field(default=0, metadata={
                "metadata": {
                    "description": "codon table to use for verifying/modifying the ORFs. Default: 0, ie the universal \
codon table but enforcing ATG as the only valid start codon."},
                "required": True,
                "validate": validate.OneOf(
                    choices=[0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22,
                             23, 24, 25, 26, 27, 28, 29, 30, 31,
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
                             "Blastocrithidia Nuclear"])
    })

    max_objects: int = field(default=10000000, metadata={
                "metadata": {
                    "description": "Maximum number of objects to keep in memory while loading data into the database. \
Linearly related to memory usage and inversely correlated with runtime."},
                "validate": validate.Range(min=1)
    })
    max_regression: float = field(default=0.2, metadata={
                "metadata": {
                    "description": "if the ORF lacks a valid start site, this percentage indicates how far along \
the sequence Mikado should look for a good start site. Eg. with a value of 0.1, on a 300bp sequence with an open ORF \
Mikado would look for an alternative in-frame start codon in the first 30 bps (10% of the cDNA)."},
            "validate": validate.Range(min=0.0, max=1.0),
            "required": True
    })
    start_adjustment: bool = field(default=True, metadata={
                "metadata": {"description": "Boolean switch. If set to true (default), if an ORF is truncated at \
the 5' Mikado will look for internal start codons. See 'max_regression'."},
    })
    max_target_seqs: int = field(default=100000, metadata={
                "metadata": {
                    "description": "Equivalently to BLAST, it indicates the maximum number of targets \
to keep per blasted sequence."},
                "validate": validate.Range(min=1)
    })
    force: bool = field(default=True, metadata={
                "metadata": {"description": "Whether to drop and reload everything into the DB"},
    })
    single_thread: bool = field(default=False, metadata={})
