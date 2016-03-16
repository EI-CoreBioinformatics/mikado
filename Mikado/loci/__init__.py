# coding: utf-8

"""
    This module defines the objects which rely the information on the transcript
    location on the genome. The most basic construct is the transcript, which holds
    information about a single RNA molecule.
    Transcripts can then be gathered in superloci, subloci, monosubloci or loci;
    all of these are defined as implementations of the blueprint "abstractlocus" class.
    The creation of the loci is delegated to the "Creator" class.
"""

from .clique_methods import *
from .abstractlocus import Abstractlocus
from .transcript import Transcript
from .excluded import Excluded
from .locus import Locus
from .monosublocusholder import MonosublocusHolder
from .superlocus import Superlocus, Sublocus, Monosublocus
from .reference_gene import Gene

# from .picker import Picker
# from . import abstractlocus
# from . import picker
# from . import excluded
# from . import locus
# from . import monosublocus
# from . import sublocus
# from . import superlocus
# from . import transcript
# from . import transcriptchecker

__title__ = "loci"
