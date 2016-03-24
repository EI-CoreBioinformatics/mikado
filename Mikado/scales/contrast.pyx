from .resultstorer import ResultStorer
from ..loci import Transcript
from ..utilities.overlap cimport c_overlap
from .f1 cimport calc_f1
import cython


__author__ = 'Luca Venturini'

@cython.profile(True)
cdef str __assign_monoexonic_ccode(prediction, reference, long nucl_overlap, tuple stats):

    cdef:
        str ccode
        double nucl_recall, nucl_precision, nucl_f1
        double exon_recall, exon_precision, exon_f1
        double junction_recall, junction_precision, junction_f1
        list overlaps
        short p_exon_num, r_exon_num
        long p_start, p_end, r_start, r_end
        long over
        long over_left, over_right
        long p_cdna_length
        long i_length, i_start, i_end
        
    p_exon_num, p_start, p_end = prediction.exon_num, prediction.start, prediction.end
    r_exon_num, r_start, r_end = reference.exon_num, reference.start, reference.end

    ccode = ""
    (nucl_recall, nucl_precision, nucl_f1,
     exon_recall, exon_precision, exon_f1,
     junction_recall, junction_precision, junction_f1) = stats

    if p_exon_num == 1 and r_exon_num > 1:
        if nucl_precision < 1 and nucl_overlap > 0:
            overlaps = []
            for intron in sorted(reference.introns):
                i_start, i_end = intron[0], intron[1]
                over = c_overlap(i_start, i_end, p_start, p_end, flank=0, positive=True)
                if over > 0:
                    i_length = i_end - i_start + 1
                    overlaps.append((over, i_length))
            if len(overlaps) == 0:
                # Completely contained inside
                ccode = "mo"
            elif len(overlaps) == 1:
                over, i_length = overlaps.pop()
                p_cdna_length = prediction.cdna_length
                if over == i_length and over < p_cdna_length:
                    ccode = "mo"
                elif 10 < over < i_length:
                    if (r_start < p_start < p_end < r_end):
                        ccode = "e"
                    else:
                        ccode = "mo"
                else:
                    ccode = "mo"
            elif len(overlaps) > 2:
                ccode = "mo"
            elif len(overlaps) == 2:
                over_left, over_right = overlaps[0][0], overlaps[1][0]
                if over_left < 10 and over_right < 10:
                    ccode = "mo"
                else:
                    ccode = "e"
        elif nucl_overlap > 0:
            ccode = "mo"
        elif (nucl_recall == 0 and r_start < p_start < r_end):
            ccode = "i"  # Monoexonic fragment inside an intron
    elif p_exon_num > 1 and r_exon_num == 1:
        # if nucl_recall == 1:
        #     ccode = "h"  # Extension
        # else:
        ccode = "O"  # Reverse generic overlap
    elif p_exon_num == r_exon_num == 1:
        if nucl_f1 >= 0.95 and reference.strand == prediction.strand:
            reference_exon = reference.exons[0]
            ccode = "_"
        elif nucl_precision == 1:
            ccode = "c"  # contained
        else:
            ccode = "m"  # just a generic exon overlap b/w two monoexonic transcripts

    # stats = (nucl_recall, nucl_precision, nucl_f1,
    #         exon_recall, exon_precision, exon_f1,
    #         junction_recall, junction_precision, junction_f1)

    return ccode


cdef str __assign_multiexonic_ccode(prediction, reference, long nucl_overlap, tuple stats):

    """
    Static method to assign a class code when both transcripts are multiexonic.
    :param prediction: prediction transcript
    :type prediction: Transcript
    :param reference: reference transcript
    :type reference: Transcript
    :param nucl_overlap: overlap between the exonic parts of the two transcripts
    :type nucl_overlap: int
    :param stats: a tuple of 9 statistics (Base-level precision, recall, F1, then exon-level,
    then junction-level)
    :return:
    """

    cdef:
        str ccode
        double nucl_recall, nucl_precision, nucl_f1
        double exon_recall, exon_precision, exon_f1
        double junction_recall, junction_precision, junction_f1
        list overlaps
        int p_exon_num, r_exon_num
        long p_start, p_end, r_start, r_end
        long over
        set r_splices, p_splices
        set r_introns, p_introns
        long intron_start, intron_end
        bint start_in, end_ind
        long min_splice, max_splice
        list corr_exons

    r_splices, p_splices = reference.splices, prediction.splices
    r_introns, p_introns = reference.introns, prediction.introns
    p_exon_num, p_start, p_end = prediction.exon_num, prediction.start, prediction.end
    r_exon_num, r_start, r_end = reference.exon_num, reference.start, reference.end

    (nucl_recall, nucl_precision, nucl_f1,
     exon_recall, exon_precision, exon_f1,
     junction_recall, junction_precision, junction_f1) = stats

    ccode = ""
    if junction_recall == 1 and junction_precision < 1:
        # Check if this is an extension
        new_splices = set(prediction.splices) - set(reference.splices)
        # If one of the new splices is internal to the intron chain, it's a j
        if any(min(reference.splices) <
               splice <
               max(reference.splices) for splice in new_splices):
            ccode = "j"
        else:
            if any(reference.start < splice < reference.end for splice in new_splices):
                # if nucl_recall < 1:
                ccode = "J"
            else:
                ccode = "n"
    elif 0 < junction_recall < 1 and 0 < junction_precision < 1:
        ccode = "j"
    elif junction_precision == 1 and junction_recall < 1:
        if nucl_precision == 1:
            # assert nucl_recall < 1
            ccode = "c"
        else:
            missed_introns = r_introns - p_introns
            min_splice, max_splice = min(prediction.splices), max(prediction.splices)
            start_in = 0
            end_in = 0
            for intron in missed_introns:
                intron_start, intron_end = intron
                if p_start < intron_start and intron_end < min_splice:
                    start_in = 1
                    break
                if p_end > intron_end and intron_start > max_splice:
                    end_in = 1
                    break

            if start_in or end_in:
                ccode = "j"
            else:
                ccode = "C"

    elif junction_recall == 0 and junction_precision == 0:
        if nucl_f1 > 0:
            corr_exons = []
            for pred_index, intron in enumerate(p_introns):
                for ref_index, ref_intron in enumerate(r_introns):
                    if c_overlap(ref_intron[0], ref_intron[1],
                                 intron[0], intron[1],
                                 flank=0, positive=1) > 0:
                        corr_exons.append((ref_index, pred_index))
                if len(corr_exons) >= 1:
                    break
            if len(corr_exons) == 1:
                ccode = "h"
            else:
                ccode = "o"
        else:
            if nucl_overlap == 0:
                # The only explanation for no nucleotide overlap
                # and no junction overlap is that it is inside an intron
                if r_start < p_start < r_end:
                    ccode = "I"
                elif p_start < r_start < p_end:
                    ccode = "rI"  # reverse intron retention

    # stats = (nucl_recall, nucl_precision, nucl_f1,
    #          exon_recall, exon_precision, exon_f1,
    #          junction_recall, junction_precision, junction_f1)

    return ccode


@cython.profile(True)
cpdef tuple compare(prediction, reference):

    """Function to compare two transcripts and determine a ccode.

    :param prediction: the transcript query
    :type prediction: Transcript

    :param reference: the reference transcript against which we desire to
    calculate the ccode and other stats.
    :type reference: Transcript

    :rtype (ResultStorer, (int,int)) | (ResultStorer, None)

    Available ccodes (from Cufflinks documentation):

    - =    Complete intron chain match
    - c    Contained (perfect junction recall and precision, imperfect recall)
    - j    Potentially novel isoform (fragment): at least one splice junction is shared
    with a reference transcript
    - e    Single exon transfrag overlapping a reference exon and at least
    10 bp of a reference intron, indicating a possible pre-mRNA fragment.
    - i    A *monoexonic* transfrag falling entirely within a reference intron
    - o    Generic exonic overlap with a reference transcript
    - p    Possible polymerase run-on fragment (within 2Kbases of a reference transcript)
    - u    Unknown, intergenic transcript
    - x    Exonic overlap with reference on the opposite strand (class codes e, o, m, c, _)
    - X    Overlap on the opposite strand, with some junctions in common (probably a serious mistake,
           unless non-canonical splicing junctions are involved).

    Please note that the description for i is changed from Cufflinks.

    We also provide the following additional classifications:

    - f    gene fusion - in this case, this ccode will be followed by the
    ccodes of the matches for each gene, separated by comma
    - _    Complete match, for monoexonic transcripts
    (nucleotide F1>=80% - i.e. min(precision,recall)>=66.7%
    - m    Exon overlap between two monoexonic transcripts
    - n    Potential extension of the reference - we have added new splice junctions
    *outside* the boundaries of the transcript itself
    - C    Contained transcript with overextensions on either side
    (perfect junction recall, imperfect nucleotide specificity)
    - J    Potentially novel isoform, where all the known junctions
    have been confirmed and we have added others as well *externally*
    - I    *multiexonic* transcript falling completely inside a known transcript
    - h    AS event in which at least a couple of introns overlaps but without any
           junction in common.
    - O    Reverse generic overlap - the reference is monoexonic while the prediction isn't
    - P    Possible polymerase run-on fragment
    - mo   Monoexonic overlap - the prediction is monoexonic and the reference is multiexonic
    (within 2K bases of a reference transcript), on the opposite strand

    This is a class method, and can therefore be used outside of a class instance.
    """

    prediction.finalize()
    reference.finalize()

    cdef:
        long nucl_overlap, distance
        set __pred_exons, __ref_exons
        tuple exon, other_exon, stats
        double p_cdna_length, r_cdna_length  # Cast as doubles to ensure division correctness
        double nucl_recall, nucl_precision, nucl_f1
        double exon_recall, exon_precision, exon_f1
        double junction_recall, junction_precision, junction_f1
        short r_exon_num, p_exon_num
        long p_start, p_end, r_start, r_end


    p_cdna_length, r_cdna_length = prediction.cdna_length, reference.cdna_length
    p_exon_num, r_exon_num = prediction.exon_num, reference.exon_num
    p_start, p_end, r_start, r_end = (prediction.start, prediction.end,
                                      reference.start, reference.end)

    nucl_overlap = 0

    __pred_exons = set()
    __ref_exons = set()

    for exon in prediction.exons:
        exon = (exon[0], exon[1] + 1)
        __pred_exons.add(exon)
        for other_exon in reference.exons:
            other_exon = (other_exon[0], other_exon[1] + 1)
            __ref_exons.add(other_exon)
            nucl_overlap += c_overlap(exon[0], exon[1],
                                      other_exon[0], other_exon[1],
                                      flank=0,
                                      positive=1)

    # Quick verification that the overlap is not too big
    assert nucl_overlap <= min(r_cdna_length, p_cdna_length), \
        (prediction.id, prediction.cdna_length,
         reference.id, reference.cdna_length, nucl_overlap)

    nucl_recall = nucl_overlap / r_cdna_length  # Sensitivity
    nucl_precision = nucl_overlap / p_cdna_length
    nucl_f1 = calc_f1(nucl_recall, nucl_precision)

    # Exon statistics
    recalled_exons = set.intersection(__pred_exons, __ref_exons)
    exon_recall = len(recalled_exons)/len(reference.exons)
    exon_precision = len(recalled_exons)/len(prediction.exons)
    exon_f1 = calc_f1(exon_recall, exon_precision)

    reference_exon = None

    # Both multiexonic
    if p_exon_num > 1 and r_exon_num > 1:
        junction_overlap = len(set.intersection(
            set(prediction.splices),
            set(reference.splices)))
        junction_recall = junction_overlap / len(reference.splices)
        junction_precision = junction_overlap / len(prediction.splices)
        junction_f1 = calc_f1(junction_recall, junction_precision)

    elif p_exon_num == r_exon_num == 1:
        # junction_overlap = junction_f1 = junction_precision = junction_recall = 1
        junction_f1 = junction_precision = junction_recall = 1
    else:
        # junction_overlap = junction_f1 = junction_precision = junction_recall = 0
        junction_f1 = junction_precision = junction_recall = 0

    ccode = ""
    distance = 0
    if junction_f1 == 1 and prediction.exon_num > 1:
        if prediction.strand == reference.strand:
            ccode = "="  # We have recovered all the junctions
        else:
            ccode = "c"  # We will set this to x at the end of the function

    elif junction_f1 == 1 and nucl_f1 >= 0.80:
        reference_exon = reference.exons[0]
        ccode = "_"  # We have recovered all the junctions

    # Outside the transcript - polymerase run-on
    elif prediction.start > reference.end or prediction.end < reference.start:
        if reference.strand == prediction.strand:
            ccode = "p"
        else:
            ccode = "P"
        distance = max(p_start - r_end, r_start - p_end)

    elif nucl_precision == 1:
        if prediction.exon_num == 1 or (prediction.exon_num > 1 and junction_precision == 1):
            ccode = "c"

    if ccode == "":
        stats = (nucl_recall, nucl_precision, nucl_f1,
                 exon_recall, exon_precision, exon_f1,
                 junction_recall, junction_precision, junction_f1)

        if min(p_exon_num, r_exon_num) > 1:
            ccode = __assign_multiexonic_ccode(prediction, reference,
                                               nucl_overlap, stats)
        else:
            ccode = __assign_monoexonic_ccode(prediction, reference,
                                              nucl_overlap, stats)

    if (prediction.strand != reference.strand and
            all([_ is not None for _ in (prediction.strand, reference.strand)])):
        if ccode in ("e", "mo", "c", "m", "_", "C"):
            ccode = "x"  # "x{0}".format(ccode)
        elif ccode not in ("u", "i", "I", "p", "P", "x"):
            ccode = "X"  # "X{0}".format(ccode)

    if prediction.strand != reference.strand:
        reference_exon = None

    result = ResultStorer(reference.id,
                          ",".join(reference.parent),
                          ccode, prediction.id,
                          ",".join(prediction.parent),
                          # Nucleotide stats
                          round(nucl_precision * 100, 2),
                          round(100 * nucl_recall, 2),
                          round(100 * nucl_f1, 2),
                          # Junction stats
                          round(junction_precision * 100, 2),
                          round(100 * junction_recall, 2),
                          round(100 * junction_f1, 2),
                          # Exonic stats
                          round(exon_precision * 100, 2),
                          round(100 * exon_recall, 2),
                          round(100 * exon_f1, 2),
                          distance)

    if ccode == "":
        raise ValueError("Ccode is null;\n{0}".format(repr(result)))

    return (result, reference_exon)
