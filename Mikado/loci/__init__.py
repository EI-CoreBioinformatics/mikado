# coding: utf-8

"""
    This module defines the objects which rely the information on the transcript
    location on the genome. The most basic construct is the transcript, which holds
    information about a single RNA molecule.
    Transcripts can then be gathered in superloci, subloci, monosubloci or loci;
    all of these are defined as implementations of the blueprint "abstractlocus" class.
    The creation of the loci is delegated to the "Creator" class.
"""


from ..transcripts import Transcript, Gene
from .abstractlocus import Abstractlocus
from .excluded import Excluded
from .locus import Locus
from .monosublocusholder import MonosublocusHolder
from .monosublocus import Monosublocus
from .superlocus import Superlocus, Sublocus


__title__ = "loci"
