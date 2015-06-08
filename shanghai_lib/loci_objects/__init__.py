''' 
    This module defines the objects which rely the information on the transcript location on the genome.
    The most basic construct is the transcript, which holds information about a single RNA molecule.
    Transcripts can then be gathered in superloci, subloci, monosubloci or loci;
    all of these are defined as implemenations of the blueprint "abstractlocus" class.
    The creation of the loci is delegated to the "Creator" class.
'''

import shanghai_lib.loci_objects.abstractlocus
import shanghai_lib.loci_objects.Creator
import shanghai_lib.loci_objects.excluded_locus
import shanghai_lib.loci_objects.locus
import shanghai_lib.loci_objects.monosublocus
import shanghai_lib.loci_objects.sublocus
import shanghai_lib.loci_objects.superlocus
import shanghai_lib.loci_objects.transcript
import shanghai_lib.loci_objects.transcript_checker