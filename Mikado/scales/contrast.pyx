from .resultstorer import ResultStorer
from ..loci import Transcript
from ..utilities.overlap cimport c_overlap
from .f1 cimport calc_f1
import cython


__author__ = 'Luca Venturini'

@cython.profile(True)
@cython.boundscheck(False)
cdef str __assign_monoexonic_ccode(prediction, reference, long nucl_overlap, double stats[9]):

    cdef:
        str ccode
        double nucl_recall, nucl_precision, nucl_f1
        double exon_recall, exon_precision, exon_f1
        double junction_recall, junction_precision, junction_f1
        short p_exon_num, r_exon_num
        long p_start, p_end, r_start, r_end
        long over
        long over_left, over_right
        long p_cdna_length
        long i_length, i_start, i_end
        long introns[2]
        long overlaps[3]
        unsigned short int index
        str r_strand, p_strand
        # double stats[9]
        
    p_exon_num, p_start, p_end = prediction.exon_num, prediction.start, prediction.end
    r_exon_num, r_start, r_end = reference.exon_num, reference.start, reference.end
    if reference.strand is None:
        r_strand = ""
    else:
        r_strand = reference.strand

    if prediction.strand is None:
        p_strand = ""
    else:
        p_strand = prediction.strand

    ccode = ""
    nucl_recall = stats[0]
    nucl_precision = stats[1]
    nucl_f1 = stats[2]
    exon_recall = stats[3]
    exon_precision = stats[4]
    exon_f1 = stats[5]
    junction_recall = stats[6]
    junction_precision = stats[7]
    junction_f1 = stats[8]

    if p_exon_num == 1 and r_exon_num > 1:
        if nucl_precision < 1 and nucl_overlap > 0:
            index = -1
            for intron in sorted(reference.introns):
                i_start, i_end = intron[0], intron[1]
                over = c_overlap(i_start, i_end, p_start, p_end, flank=0, positive=True)
                if over > 0:
                    index += 1
                    if index > 1:
                        break
                    i_length = i_end - i_start + 1
                    overlaps[index] = over
                    introns[index] = i_length

            if index == -1:
                # Completely contained inside
                ccode = "mo"
            elif index == 0:
                over = overlaps[0]
                i_length = introns[0]
                # over, i_length = overlaps[0][0], overlaps[0][1]
                p_cdna_length = prediction.cdna_length
                if (10 < over < i_length) and (r_start < p_start < p_end < r_end):
                    ccode = "e"
                else:
                    ccode = "mo"
            elif index > 1:
                ccode = "mo"
            elif index == 1:
                over_left, over_right = overlaps[0], overlaps[1]
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
        if nucl_f1 >= 0.95 and r_strand == p_strand:
            reference_exon = reference.exons[0]
            ccode = "_"
        elif nucl_precision == 1:
            ccode = "c"  # contained
        else:
            ccode = "m"  # just a generic exon overlap b/w two monoexonic transcripts

    return ccode


cdef str __assign_multiexonic_ccode(prediction, reference, long nucl_overlap, double stats[9]):

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
        long min_splice, max_splice, splice
        short unsigned int overlapping_introns
        short unsigned int index

    r_splices, p_splices = reference.splices, prediction.splices
    r_introns, p_introns = reference.introns, prediction.introns
    p_exon_num, p_start, p_end = prediction.exon_num, prediction.start, prediction.end
    r_exon_num, r_start, r_end = reference.exon_num, reference.start, reference.end

    nucl_recall = stats[0]
    nucl_precision = stats[1]
    nucl_f1 = stats[2]
    exon_recall = stats[3]
    exon_precision = stats[4]
    exon_f1 = stats[5]
    junction_recall = stats[6]
    junction_precision = stats[7]
    junction_f1 = stats[8]

    ccode = ""
    if junction_recall == 1 and junction_precision < 1:
        # Check if this is an extension
        # If one of the new splices is internal to the intron chain, it's a j

        ccode = "n"
        max_splice = max(reference.splices)
        min_splice = min(reference.splices)
        for splice in p_splices:
            if splice in r_splices:
                continue
            if min_splice < splice < max_splice:
                ccode = "j"
                break
            elif r_start < splice < r_end:
                ccode = "J"

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
            overlapping_introns = 0
            ccode = "o"
            for intron in p_introns:
                for ref_intron in r_introns:
                    if c_overlap(ref_intron[0], ref_intron[1],
                                 intron[0], intron[1],
                                 0, 1) > 0:
                        # corr_exons.append((ref_index, pred_index))
                        overlapping_introns += 1
                if overlapping_introns >= 1:
                    ccode = "h"
                    break
        else:
            if nucl_overlap == 0:
                # The only explanation for no nucleotide overlap
                # and no junction overlap is that it is inside an intron
                if r_start < p_start < r_end:
                    ccode = "I"
                elif p_start < r_start < p_end:
                    ccode = "rI"  # reverse intron retention

    return ccode


@cython.profile(True)
@cython.cdivision(True)
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
        long exon_a, exon_b, other_exon_a, other_exon_b
        double p_cdna_length, r_cdna_length  # Cast as doubles to ensure division correctness
        double nucl_recall, nucl_precision, nucl_f1
        double exon_recall, exon_precision, exon_f1
        double junction_recall, junction_precision, junction_f1
        long r_exon_num, p_exon_num
        long p_start, p_end, r_start, r_end
        double recalled_exons
        set r_splices, p_splices
        double len_r_splices, len_p_splices
        long junction_overlap
        double stats[9]
        str r_strand, p_strand

    p_cdna_length, r_cdna_length = prediction.cdna_length, reference.cdna_length
    p_exon_num, r_exon_num = prediction.exon_num, reference.exon_num
    p_start, p_end, r_start, r_end = (prediction.start, prediction.end,
                                      reference.start, reference.end)

    if reference.strand is None:
        r_strand = ""
    else:
        r_strand = reference.strand

    if prediction.strand is None:
        p_strand = ""
    else:
        p_strand = prediction.strand

    nucl_overlap = 0

    __pred_exons = set()
    __ref_exons = set()

    for exon in prediction.exons:
        exon_a, exon_b = exon[0], exon[1] + 1
        __pred_exons.add(exon)
        for other_exon in reference.exons:
            other_exon_a, other_exon_b = other_exon[0], other_exon[1] + 1
            __ref_exons.add(other_exon)
            nucl_overlap += c_overlap(exon_a, exon_b,
                                      other_exon_a, other_exon_b,
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
    recalled_exons = len(set.intersection(__pred_exons, __ref_exons))
    exon_recall = recalled_exons / r_exon_num
    exon_precision = recalled_exons / p_exon_num
    exon_f1 = calc_f1(exon_recall, exon_precision)

    reference_exon = None

    # Both multiexonic
    if p_exon_num > 1 and r_exon_num > 1:
        r_splices, p_splices = reference.splices, prediction.splices
        len_rsplices, len_psplices = len(r_splices), len(p_splices)
        junction_overlap = len(set.intersection(p_splices, r_splices))
        junction_recall = junction_overlap / len_rsplices
        junction_precision = junction_overlap / len_psplices
        junction_f1 = calc_f1(junction_recall, junction_precision)

    elif p_exon_num == r_exon_num == 1:
        # junction_overlap = junction_f1 = junction_precision = junction_recall = 1
        junction_f1 = junction_precision = junction_recall = 1
    else:
        # junction_overlap = junction_f1 = junction_precision = junction_recall = 0
        junction_f1 = junction_precision = junction_recall = 0

    ccode = ""
    distance = 0
    if junction_f1 == 1 and p_exon_num > 1:
        if p_strand == r_strand:
            ccode = "="  # We have recovered all the junctions
        else:
            ccode = "c"  # We will set this to x at the end of the function

    elif junction_f1 == 1 and nucl_f1 >= 0.80:
        reference_exon = reference.exons[0]
        ccode = "_"  # We have recovered all the junctions

    # Outside the transcript - polymerase run-on
    elif p_start > r_end or p_end < r_start:
        if r_strand == p_strand:
            ccode = "p"
        else:
            ccode = "P"
        distance = max(p_start - r_end, r_start - p_end)

    elif nucl_precision == 1:
        if p_exon_num == 1 or (p_exon_num > 1 and junction_precision == 1):
            ccode = "c"

    if ccode == "":
        # stats = (nucl_recall, nucl_precision, nucl_f1,
        #         exon_recall, exon_precision, exon_f1,
        #          junction_recall, junction_precision, junction_f1)
        stats[0] = nucl_recall
        stats[1] = nucl_precision
        stats[2] = nucl_f1
        stats[3] = exon_recall
        stats[4] = exon_precision
        stats[5] = exon_f1
        stats[6] = junction_recall
        stats[7] = junction_precision
        stats[8] = junction_f1


        if min(p_exon_num, r_exon_num) > 1:
            ccode = __assign_multiexonic_ccode(prediction, reference,
                                               nucl_overlap, stats)
        else:
            ccode = __assign_monoexonic_ccode(prediction, reference,
                                              nucl_overlap, stats)

    if (p_strand != r_strand and p_strand != "" and r_strand != ""):
        if ccode in ("e", "mo", "c", "m", "_", "C"):
            ccode = "x"  # "x{0}".format(ccode)
        elif ccode not in ("u", "i", "I", "p", "P", "x"):
            ccode = "X"  # "X{0}".format(ccode)

    elif p_strand != r_strand:
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