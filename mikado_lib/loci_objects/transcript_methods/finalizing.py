"""
This module provides the functions needed to check a transcript for consinstency,
e.g. reliability of the CDS/UTR, sanity of borders, etc.
"""

import intervaltree
import operator
import mikado_lib.exceptions

__author__ = 'Luca Venturini'


def __basic_final_checks(transcript):

    """
    Function that verifies minimal criteria of a transcript before finalising.
    :type transcript: mikado_lib.loci_objects.transcript.Transcript
    
    :return:
    """

    if len(transcript.exons) == 0:
        raise mikado_lib.exceptions.InvalidTranscript(
            "No exon defined for the transcript {0}. Aborting".format(transcript.id))

    if not isinstance(transcript.exons[0], intervaltree.Interval):
        _ = [intervaltree.Interval(int(exon[0]), int(exon[1])) for exon in transcript.exons]
        transcript.logger.debug("Converting to interval objects")
        transcript.exons = _

    assert all([isinstance(exon, intervaltree.Interval) for exon in transcript.exons]), transcript.exons
    assert all([(isinstance(exon[0], int) and isinstance(exon[1], int)) for exon in transcript.exons]), transcript.exons

    if len(transcript.exons) > 1 and transcript.strand is None:
        raise mikado_lib.exceptions.InvalidTranscript(
            "Multiexonic transcripts must have a defined strand! Error for {0}".format(
                transcript.id))

    if transcript.combined_utr != [] and transcript.combined_cds == []:
        raise mikado_lib.exceptions.InvalidTranscript(
            "Transcript {tid} has defined UTRs but no CDS feature!".format(
                tid=transcript.id))


def __check_cdna_vs_utr(transcript):

    """
    Verify that cDNA + UTR in the transcript add up.
    :return:
    """

    if transcript.cdna_length > transcript.combined_utr_length + transcript.combined_cds_length:
        if transcript.combined_utr == [] and transcript.combined_cds != []:
            transcript.combined_cds = sorted(transcript.combined_cds,
                                             key=operator.itemgetter(0, 1))
            for exon in transcript.exons:
                assert isinstance(exon, intervaltree.Interval)
                if len(transcript.combined_cds) > 1:
                    assert isinstance(transcript.combined_cds[0], intervaltree.Interval),\
                        type(transcript.combined_cds[0])
                if exon in transcript.combined_cds:
                    continue
                elif (exon[1] < transcript.combined_cds[0][0] or
                      exon[0] > transcript.combined_cds[-1][1]):
                    transcript.combined_utr.append(exon)
                elif (exon[0] < transcript.combined_cds[0][0] and
                      exon[1] == transcript.combined_cds[0][1]):
                    transcript.combined_utr.append(intervaltree.Interval(
                        exon[0], transcript.combined_cds[0][0] - 1))

                elif (exon[1] > transcript.combined_cds[-1][1] and
                      exon[0] == transcript.combined_cds[-1][0]):
                    transcript.combined_utr.append(intervaltree.Interval(
                        transcript.combined_cds[-1][1] + 1, exon[1]))
                else:
                    if len(transcript.combined_cds) == 1:
                        transcript.combined_utr.append(intervaltree.Interval(
                            exon[0], transcript.combined_cds[0][0] - 1))
                        transcript.combined_utr.append(intervaltree.Interval(
                            transcript.combined_cds[-1][1] + 1, exon[1]))
                    else:
                        raise mikado_lib.exceptions.InvalidCDS(
                            "Error while inferring the UTR",
                            exon, transcript.id,
                            transcript.exons, transcript.combined_cds)

            equality_one = (transcript.combined_cds_length == transcript.combined_utr_length == 0)
            equality_two = (transcript.cdna_length ==
                            transcript.combined_utr_length + transcript.combined_cds_length)
            if not (equality_one or equality_two):
                raise mikado_lib.exceptions.InvalidCDS(
                    "Failed to create the UTR",
                    transcript.id, transcript.exons,
                    transcript.combined_cds, transcript.combined_utr)
        else:
            pass


def __calculate_introns(transcript):

    """Private method to create the stores of intron
    and splice sites positions.
    """

    introns = []
    splices = []

    if len(transcript.exons) > 1:
        for index in range(len(transcript.exons) - 1):
            exona, exonb = transcript.exons[index:index + 2]
            if exona[1] >= exonb[0]:
                raise mikado_lib.exceptions.InvalidTranscript(
                    "Overlapping exons found!\n{0} {1}/{2}\n{3}".format(
                        transcript.id, exona, exonb, transcript.exons))
            # Append the splice junction
            introns.append(intervaltree.Interval(exona[1] + 1, exonb[0] - 1))
            # Append the splice locations
            splices.extend([exona[1] + 1, exonb[0] - 1])
    transcript.introns = set(introns)
    transcript.splices = set(splices)


def __check_completeness(transcript):

    """Private method that checks whether a transcript is complete
    or not based solely on the presence of CDS/UTR information."""

    if len(transcript.combined_utr) > 0:
        if transcript.combined_utr[0][0] < transcript.combined_cds[0][0]:
            if transcript.strand == "+":
                transcript.has_start_codon = True
            elif transcript.strand == "-":
                transcript.has_stop_codon = True
        if transcript.combined_utr[-1][1] > transcript.combined_cds[-1][1]:
            if transcript.strand == "+":
                transcript.has_stop_codon = True
            elif transcript.strand == "-":
                transcript.has_start_codon = True


def __verify_boundaries(transcript):

    """
    Method to verify that the start/end of the transcripts are exactly where they should.
    Called from finalise.
    :return:
    """

    try:
        if transcript.exons[0][0] != transcript.start or transcript.exons[-1][1] != transcript.end:
            if (transcript.exons[0][0] > transcript.start and
                    transcript.selected_cds[0][0] == transcript.start):
                transcript.exons[0] = (transcript.start, transcript.exons[0][0])
            if (transcript.exons[-1][1] < transcript.end and
                    transcript.selected_cds[-1][1] == transcript.end):
                transcript.exons[-1] = (transcript.exons[-1][0], transcript.end)

            if (transcript.exons[0][0] != transcript.start or
                    transcript.exons[-1][1] != transcript.end):
                raise mikado_lib.exceptions.InvalidTranscript(
                    """The transcript {id} has coordinates {tstart}:{tend},
                but its first and last exons define it up until {estart}:{eend}!
                Exons: {exons}
                """.format(id=transcript.id,
                           tstart=transcript.start,
                           tend=transcript.end,
                           estart=transcript.exons[0][0],
                           eend=transcript.exons[-1][1],
                           exons=transcript.exons))
    except IndexError as err:
        raise mikado_lib.exceptions.InvalidTranscript(
            err, transcript.id, str(transcript.exons))


def __check_internal_orf(transcript, exons, orf):

    """
    Method that verifies that an internal ORF does not have any internal gap.

    :param exons: list of original exons
    :param orf: internal ORF to check.
    :return:
    """

    orf_segments = sorted([_[1] for _ in orf if _[0] == "CDS"])

    if len(orf_segments) > 0:
        assert isinstance(orf_segments[0], intervaltree.Interval), (orf_segments, type(orf_segments[0]))

    previous_exon_index = None

    for orf_segment in orf_segments:
        exon_found = False
        for exon_position, exon in enumerate(exons):
            if exon[0] <= orf_segment[0] <= orf_segment[1] <= exon[1]:
                if previous_exon_index is not None and previous_exon_index + 1 != exon_position:
                    exc = mikado_lib.exceptions.InvalidTranscript(
                        """Invalid ORF for {0}, invalid index: {1} (for {2}), expected {3}
                        {4} CDS vs. {5} exons""".format(
                            transcript.id,
                            exon_position,
                            orf_segment,
                            previous_exon_index + 1,
                            orf_segments,
                            exons
                        ))
                    transcript.logger.exception(exc)
                    raise exc
                else:
                    previous_exon_index = exon_position
                    exon_found = True
                    break
        if exon_found is False:
            exc = mikado_lib.exceptions.InvalidTranscript(
                "Invalid ORF for {0}, no exon found: {1} CDS vs. {2} exons".format(
                    transcript.id,
                    orf_segments,
                    exons))
            transcript.logger.exception(exc)
            raise exc

    return


def finalize(transcript):
    """Function to calculate the internal introns from the exons.
    In the first step, it will sort the exons by their internal coordinates.

    :param transcript: the Transcript instance to finalize.
    :type transcript: mikado_lib.loci_objects.transcript.Transcript

    """

    if transcript.finalized is True:
        return

    __basic_final_checks(transcript)
    # Sort the exons by start then stop

    transcript.exons = sorted(transcript.exons, key=operator.itemgetter(0, 1))

    __check_cdna_vs_utr(transcript)

    __calculate_introns(transcript)

    transcript.combined_cds = sorted(transcript.combined_cds,
                                     key=operator.itemgetter(0, 1))

    transcript.combined_utr = sorted(transcript.combined_utr,
                                     key=operator.itemgetter(0, 1))
    __check_completeness(transcript)

    # assert self.selected_internal_orf_index > -1
    if len(transcript.internal_orfs) == 0:
        transcript.segments = [("exon", intervaltree.Interval(e[0], e[1])) for e in transcript.exons]
        transcript.segments.extend([("CDS", intervaltree.Interval(c[0], c[1])) for c in transcript.combined_cds])
        transcript.segments.extend([("UTR", intervaltree.Interval(u[0], u[1])) for u in transcript.combined_utr])
        transcript.segments = sorted(transcript.segments, key=operator.itemgetter(1, 0))
        transcript.internal_orfs = [transcript.segments]
    else:
        assert len(transcript.internal_orfs) > 0

    for internal_orf in transcript.internal_orfs:
        __check_internal_orf(transcript, transcript.exons, internal_orf)

    if len(transcript.combined_cds) > 0:
        transcript.selected_internal_orf_index = 0
        # pylint: disable=protected-access
        if len(transcript.phases) > 0:
            transcript._first_phase = sorted(transcript.phases, key=operator.itemgetter(0),
                                             reverse=(transcript.strand == "-"))[0][1]
        else:
            transcript._first_phase = 0
        # pylint: enable=protected-access

    # Necessary to set it to the default value
    _ = transcript.selected_internal_orf

    if len(transcript.combined_cds) > 0:
        transcript.feature = "mRNA"
    else:
        transcript.feature = "transcript"

    __verify_boundaries(transcript)

    if len(transcript.combined_cds) == 0:
        transcript.selected_internal_orf_cds = tuple([])
    else:
        transcript.selected_internal_orf_cds = tuple(
            internal_cds for internal_cds in transcript.internal_orfs[
                transcript.selected_internal_orf_index] if
            internal_cds[0] == "CDS")

    # Create the interval tree
    transcript.cds_tree = intervaltree.IntervalTree([
        intervaltree.Interval(cds[0]-1, cds[1]+1) for cds in transcript.combined_cds])

    transcript.finalized = True
    transcript.logger.debug("Finished finalising %s", transcript.id)

    return
