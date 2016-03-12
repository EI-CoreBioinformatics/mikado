"""
This module provides the functions needed to check a transcript for consinstency,
e.g. reliability of the CDS/UTR, sanity of borders, etc.
"""

import intervaltree
import operator
from ...exceptions import InvalidCDS, InvalidTranscript

__author__ = 'Luca Venturini'


def __basic_final_checks(transcript):

    """
    Function that verifies minimal criteria of a transcript before finalising.
    :type transcript: mikado_lib.loci_objects.transcript.Transcript

    :return:
    """

    if len(transcript.exons) == 0:
        raise InvalidTranscript(
            "No exon defined for the transcript {0}. Aborting".format(transcript.id))

    if not isinstance(transcript.exons[0], intervaltree.Interval):
        _ = [intervaltree.Interval(int(exon[0]), int(exon[1])) for exon in transcript.exons]
        transcript.logger.debug("Converting to interval objects")
        transcript.exons = _

    new_exons = []
    invalid = False
    for exon in transcript.exons:
        if not isinstance(exon, intervaltree.Interval):
            if (isinstance(exon, (tuple, list)) and len(exon) == 2 and
                    isinstance(exon[0], int) and isinstance(exon[1], int)):
                exon = intervaltree.Interval(*exon)
            else:
                raise ValueError("Invalid exon: {0}, type {1}".format(
                    exon, type(exon)))
        if exon.begin < transcript.start or exon.end > transcript.end:
            invalid = True
            break
        new_exons.append(exon)

    transcript.exons = sorted(new_exons)

    if invalid:
        raise InvalidTranscript("""Exons out of bounds of the transcript:
        ({start}, {end})
        Exons: {exons}""".format(start=transcript.start,
                                 end=transcript.end,
                                 exons=transcript.exons))

    if len(transcript.exons) > 1 and transcript.strand is None:
        raise InvalidTranscript(
            "Multiexonic transcripts must have a defined strand! Error for {0}".format(
                transcript.id))

    if transcript.combined_utr != [] and transcript.combined_cds == []:
        raise InvalidTranscript(
            "Transcript {tid} has defined UTRs but no CDS feature!".format(
                tid=transcript.id))


def _check_cdna_vs_utr(transcript):

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
                # If we have more than one CDS segment
                if len(transcript.combined_cds) > 1:
                    assert isinstance(transcript.combined_cds[0], intervaltree.Interval),\
                        type(transcript.combined_cds[0])
                # Ignore, this is a completely CDS ORF
                if exon in transcript.combined_cds:
                    continue
                # The end of the exon is before the first ORF start
                # or the start is after the last ORF segment: UTR segment
                elif (exon[1] < transcript.combined_cds[0][0] or
                      exon[0] > transcript.combined_cds[-1][1]):
                    transcript.combined_utr.append(exon)

                # The last base of the exon is the first ORF base
                elif (exon[0] < transcript.combined_cds[0][0] and
                      exon[1] == transcript.combined_cds[0][1]):
                    transcript.combined_utr.append(intervaltree.Interval(
                        exon[0], transcript.combined_cds[0][0] - 1))
                # The first base of the exon is the first base of the last ORF segment:
                # UTR after
                elif (exon[1] > transcript.combined_cds[-1][1] and
                      exon[0] == transcript.combined_cds[-1][0]):
                    transcript.combined_utr.append(intervaltree.Interval(
                        transcript.combined_cds[-1][1] + 1, exon[1]))
                else:
                    # If the ORF is contained inside a single exon, with UTR
                    # at both sites, then we create the two UTR segments
                    if len(transcript.combined_cds) == 1:
                        transcript.combined_utr.append(intervaltree.Interval(
                            exon[0], transcript.combined_cds[0][0] - 1))
                        transcript.combined_utr.append(intervaltree.Interval(
                            transcript.combined_cds[-1][1] + 1, exon[1]))
                    else:
                        # This means there is an INTERNAL UTR region between
                        # two CDS segments: something is clearly wrong!
                        raise InvalidCDS(
                            "Error while inferring the UTR",
                            exon, transcript.id,
                            transcript.exons, transcript.combined_cds)

            # If no CDS and no UTR are present, all good
            equality_one = (transcript.combined_cds_length == transcript.combined_utr_length == 0)
            # Otherwise, if cDNA length == UTR + CDS, all good
            equality_two = (transcript.cdna_length ==
                            transcript.combined_utr_length + transcript.combined_cds_length)
            if not (equality_one or equality_two):
                # Something fishy going on
                raise InvalidCDS(
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
                raise InvalidTranscript(
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
                raise InvalidTranscript(
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
        raise InvalidTranscript(
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
        assert isinstance(orf_segments[0],
                          intervaltree.Interval), (orf_segments,
                                                   type(orf_segments[0]))

    previous_exon_index = None

    for orf_segment in orf_segments:
        exon_found = False
        for exon_position, exon in enumerate(exons):
            if exon[0] <= orf_segment[0] <= orf_segment[1] <= exon[1]:
                if previous_exon_index is not None and previous_exon_index + 1 != exon_position:
                    exc = InvalidTranscript(
                        """Invalid ORF for {0}, invalid index: {1} (for {2}), expected {3}
                        {4} CDS vs. {5} exons""".format(
                            transcript.id,
                            exon_position,
                            orf_segment,
                            previous_exon_index + 1,
                            orf_segments,
                            exons
                        ))
                    transcript.logger.error(exc)
                    raise exc
                else:
                    previous_exon_index = exon_position
                    exon_found = True
                    break
        if exon_found is False:
            exc = InvalidTranscript(
                "Invalid ORF for {0}, no exon found: {1} CDS vs. {2} exons".format(
                    transcript.id,
                    orf_segments,
                    exons))
            transcript.logger.exception(exc)
            raise exc

    return


def __calc_cds_introns(transcript):

    """
    Function to calculate and memorize the selected CDS and combined CDS introns.
    :param transcript: the transcript instance.
    :type transcript: mikado_lib.loci_objects.transcript.Transcript
    :return:
    """

    if transcript.number_internal_orfs == 0 or \
                len(transcript.selected_cds) < 2 or \
                len(transcript.combined_cds) < 2:
        pass
    else:
        # Start calculating the selected CDS introns
        cintrons = []
        for first, second in zip(transcript.selected_cds[:-1], transcript.selected_cds[1:]):
            cintrons.append(
                intervaltree.Interval(first[1] + 1,
                                      second[0] - 1)
            )

        cintrons = set(cintrons)
        assert len(cintrons) > 0
        transcript._selected_cds_introns = cintrons
        transcript._combined_cds_introns = transcript._selected_cds_introns.copy()
        if transcript.number_internal_orfs > 1:
            cintrons = []
            for position in range(len(transcript.combined_cds) - 1):
                former = transcript.combined_cds[position]
                latter = transcript.combined_cds[position + 1]
                junc = intervaltree.Interval(former[1] + 1, latter[0] - 1)
                if junc in transcript.introns:
                    cintrons.append(junc)
            cintrons = set(cintrons)
            transcript._combined_cds_introns = cintrons
        assert len(transcript._combined_cds_introns) > 0

    assert len(transcript._combined_cds_introns) >= len(transcript._selected_cds_introns)
    return transcript
    

def finalize(transcript):
    """Function to calculate the internal introns from the exons.
    In the first step, it will sort the exons by their internal coordinates.

    :param transcript: the Transcript instance to finalize.
    :type transcript: mikado_lib.loci_objects.transcript.Transcript

    """

    if transcript.finalized is True:
        return

    transcript.exons = sorted(transcript.exons)
    transcript.__cdna_length = None
    __basic_final_checks(transcript)
    # Sort the exons by start then stop

    try:
        _check_cdna_vs_utr(transcript)
    except InvalidCDS:
        if transcript.combined_cds:
            transcript.logger.warning(
                "Possible faulty UTR annotation for %s, trying to recalculate it.",
                transcript.id)
            transcript.combined_utr = []
            try:
                _check_cdna_vs_utr(transcript)
            except InvalidCDS:
                transcript.logger.warning("CDS for %s completely invalid. Removing it.",
                                          transcript.id)
                transcript.combined_cds = []
                transcript.combined_utr = []
                transcript.segments = []
                transcript.internal_orfs = []
                __basic_final_checks(transcript)
                _check_cdna_vs_utr(transcript)

    __calculate_introns(transcript)

    transcript.combined_cds = sorted(transcript.combined_cds,
                                     key=operator.itemgetter(0, 1))

    transcript.combined_utr = sorted(transcript.combined_utr,
                                     key=operator.itemgetter(0, 1))
    __check_completeness(transcript)

    # assert self.selected_internal_orf_index > -1
    if min(len(transcript.segments), len(transcript.internal_orfs)) == 0:
        # Define exons
        transcript.segments = [("exon", intervaltree.Interval(e[0], e[1]))
                               for e in transcript.exons]
        # Define CDS
        transcript.segments.extend([("CDS", intervaltree.Interval(c[0], c[1]))
                                    for c in transcript.combined_cds])
        # Define UTR segments
        transcript.segments.extend([("UTR", intervaltree.Interval(u[0], u[1]))
                                    for u in transcript.combined_utr])
        # Mix and sort
        transcript.segments = sorted(transcript.segments, key=operator.itemgetter(1, 0))
        # Add to the store as a single entity
        transcript.internal_orfs = [transcript.segments]
    else:
        assert len(transcript.internal_orfs) > 0

    assert all([segment[1] in transcript.exons for segment in transcript.segments if
                segment[0] == "exon"]), (transcript.exons, transcript.segments)

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
    try:
        assert transcript.selected_cds_length >= 0
    except AssertionError as _:
        transcript.logger.warning("%s has an invalid CDS; removing it.", transcript.id)
        transcript.strip_cds()

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

    # BUG somewhere ... I am not sorting this properly before (why?)
    transcript.exons = sorted(transcript.exons)
    transcript = __calc_cds_introns(transcript)

    transcript.finalized = True
    transcript.logger.debug("Finished finalising %s", transcript.id)

    return
