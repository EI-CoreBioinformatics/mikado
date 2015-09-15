# coding: utf-8

"""
    This module defines the objects which rely the information on the transcript
    location on the genome. The most basic construct is the transcript, which holds
    information about a single RNA molecule.
    Transcripts can then be gathered in superloci, subloci, monosubloci or loci;
    all of these are defined as implemenations of the blueprint "abstractlocus" class.
    The creation of the loci is delegated to the "Creator" class.
"""

import mikado_lib.loci_objects.abstractlocus
import mikado_lib.loci_objects.Creator
import mikado_lib.loci_objects.excluded
import mikado_lib.loci_objects.locus
import mikado_lib.loci_objects.monosublocus
import mikado_lib.loci_objects.sublocus
import mikado_lib.loci_objects.superlocus
import mikado_lib.loci_objects.transcript
import mikado_lib.loci_objects.transcriptchecker
