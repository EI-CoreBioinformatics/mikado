# coding: utf-8

"""
    This module defines the objects which rely the information on the transcript
    location on the genome. The most basic construct is the transcript, which holds
    information about a single RNA molecule.
    Transcripts can then be gathered in superloci, subloci, monosubloci or loci;
    all of these are defined as implemenations of the blueprint "abstractlocus" class.
    The creation of the loci is delegated to the "Creator" class.
"""

import Mikado.loci_objects.abstractlocus
import Mikado.loci_objects.Creator
import Mikado.loci_objects.excluded
import Mikado.loci_objects.locus
import Mikado.loci_objects.monosublocus
import Mikado.loci_objects.sublocus
import Mikado.loci_objects.superlocus
import Mikado.loci_objects.transcript
import Mikado.loci_objects.transcriptchecker
