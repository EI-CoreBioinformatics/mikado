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

    transcript.logger.debug("Checking the cDNA for %s", transcript.id)
    if transcript.cdna_length > transcript.combined_utr_length + transcript.combined_cds_length:
        if transcript.combined_utr == transcript.combined_cds == []:
            # non-coding transcript
            transcript.logger.debug("%s is non coding, returning", transcript.id)
            return
        assert transcript.combined_cds != []
        transcript.logger.debug("Recalculating the UTR for %s", transcript.id)
        transcript.combined_utr = []  # Reset
        transcript.combined_cds = sorted(transcript.combined_cds,
                                         key=operator.itemgetter(0, 1))
        for exon in transcript.exons:
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


def __calculate_introns(transcript):

    """Private method to create the stores of intron
    and splice sites positions.
    """


    introns = []
    cds_introns = []
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

    if transcript.number_internal_orfs == 0 or \
                len(transcript.selected_cds) < 2 or \
                len(transcript.combined_cds) < 2:
        pass
    else:
        # Start calculating the selected CDS introns

        for first, second in zip(transcript.selected_cds[:-1],
                                 transcript.selected_cds[1:]):
            assert first != second, transcript.selected_cds
            assert first[1] < second[0], (first, second)
            # first, second = sorted([first, second])
            intron = intervaltree.Interval(first[1] + 1, second[0] - 1)
            assert intron in transcript.introns, (intron, first, second)
            cds_introns.append(intron)

        cintrons = set(cds_introns)
        assert len(cintrons) > 0
        transcript._selected_cds_introns = cintrons

        if transcript.number_internal_orfs > 1:
            cds_introns = []
            for position in range(len(transcript.combined_cds) - 1):
                former = transcript.combined_cds[position]
                latter = transcript.combined_cds[position + 1]
                junc = intervaltree.Interval(former[1] + 1, latter[0] - 1)
                if junc in transcript.introns:
                    cds_introns.append(junc)
            cintrons = set(cds_introns)
            transcript._combined_cds_introns = cintrons
            assert len(transcript._combined_cds_introns) > 0

            # for index, orf in enumerate(transcript.internal_orfs):
            #     if index == transcript.selected_internal_orf_index:
            #         continue
            #     cds = sorted([_[1] for _ in orf if _[0] == "CDS"])
            #     for first, second in zip(cds[:-1], cds[1:]):
            #         assert first != second, transcript.selected_cds
            #         assert first[1] < second[0], (first, second)
            #         # first, second = sorted([first, second])
            #         intron = intervaltree.Interval(first[1] + 1, second[0] - 1)
            #         assert intron in transcript.introns, (intron, first, second)
            #         cds_introns.append(intron)
            # cds_introns = set(cds_introns)
            # transcript._combined_cds_introns = cds_introns

        else:
            transcript._combined_cds_introns = transcript._selected_cds_introns.copy()

        assert len(transcript._combined_cds_introns) > 0

    assert len(transcript._combined_cds_introns) >= len(transcript._selected_cds_introns)
    return transcript


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
            transcript.logger.warning(
                """The transcript {id} has coordinates {tstart}:{tend},
            but its first and last exons define it up until {estart}:{eend}!
            Exons: {exons}. Shrinking it""".format(
                    id=transcript.id,
                    tstart=transcript.start,
                    tend=transcript.end,
                    estart=transcript.exons[0][0],
                    eend=transcript.exons[-1][1],
                    exons=transcript.exons))
            transcript.start = transcript.exons[0][0]
            transcript.end = transcript.exons[-1][1]

    except IndexError as err:
        raise InvalidTranscript(
            err, transcript.id, str(transcript.exons))


def __check_internal_orf(transcript, index, phases=None):

    """
    Method that verifies that an internal ORF does not have any internal gap.

    :param transcript: the transcript to analyse
    :type transcript: Mikado.loci.Transcript
    :param index: index of the internal orf to check
    :type index: int
    :param phases: dictionary of the phases derived from the GFF file
    :type phases: dict

    :return: the updated transcript
    :rtype: Mikado.loci.Transcript
    """

    orf, new_orf = transcript.internal_orfs[index], []
    previous_exon_index, exon_found = None, False
    total_cds_length = 0

    exons = sorted(transcript.exons, reverse=(transcript.strand == "-"))

    if phases:
        phases_keys = sorted(phases.keys(), reverse=(transcript.strand == "-"))
        phase_orf =[phases[_] for _ in phases_keys]
        previous = (3 - phases[phases_keys[0]] % 3) % 3
        total_cds_length += previous
    else:
        phase_orf = []
        previous = 0

    transcript.logger.debug("Phases, previous, phase_orf: %s, %d, %s",
                            phases, previous, phase_orf)

    orf = sorted(orf, key=operator.itemgetter(1), reverse=(transcript.strand == "-"))

    first, last = float("inf"), float("-inf")

    for orf_segment in orf:
        first = min(first, orf_segment[1][0])
        last = max(last, orf_segment[1][1])

        if orf_segment[0] != "CDS":
            new_orf.append(orf_segment)
            continue

        for exon_position, exon in enumerate(exons):
            if exon[0] <= orf_segment[1][0] <= orf_segment[1][1] <= exon[1]:
                if previous_exon_index is None:
                    previous_exon_index = exon_position
                else:
                    previous_exon_index += 1

                if previous_exon_index != exon_position:
                    exc = InvalidTranscript(
                        """Invalid ORF for {0}, invalid index: {1} (for {2}), expected {3} ({4})
                        {4} CDS vs. {5} exons""".format(
                            transcript.id,
                            exon_position,
                            orf_segment,
                            previous_exon_index + 1,
                            [_ for _ in enumerate(exons) if _[1][0] == orf_segment[1][0] or
                             _[1][1] == orf_segment[1][1]][0],
                            orf,
                            exons))
                    transcript.logger.error(exc)
                    current = len(transcript.internal_orfs)
                    del transcript.internal_orfs[index]
                    assert len(transcript.internal_orfs) == current - 1
                    raise exc

                exon_found = True
                break

        if exon_found is False:
            exc = InvalidTranscript(
                "Invalid ORF for {0}, no exon found: {1} CDS vs. {2} exons".format(
                    transcript.id, orf, exons))
            transcript.logger.exception(exc)
            raise exc

        exon_found = False

        phase = (3 - (total_cds_length % 3)) % 3
        if phases:
            key = (orf_segment[1][0], orf_segment[1][1])
            if key not in phases:
                transcript.logger.warning("Phase not found for %s key %s, recalculating." % (
                    transcript.id, key))
                phases = None
            elif phases[key] != phase:
                transcript.logger.warning(
                    "Wrong phase in % sfor key %s, %d instead of %d. recalculating." % (
                        transcript.id, key, phases[key], phase))
            else:
                phase = phases[key]

        phase_orf.append(phase)
        new_orf.append((orf_segment[0], orf_segment[1], phase))
        previous = orf_segment[1].length() + 1
        total_cds_length += previous

    transcript.logger.debug("Total CDS length %d", total_cds_length)
    assert first == transcript.start, (transcript.start, first, orf[0], new_orf[0])
    assert last == transcript.end, (transcript.end, last, orf, new_orf)

    if total_cds_length % 3 != 0:
        # The transcript is truncated. Check it makes sense.

        cds, utr = [], []
        for _ in new_orf:
            if _[0] == "CDS":
                cds.append(_[1])
            elif _[0] == "UTR":
                utr.append(_[1])

        total = sum(_.length() + 1 for _ in cds)
        # This is the amount of bases left over
        remainder = total % 3

        orf_start = min(_[0] for _ in cds)
        orf_end = max(_[1] for _ in cds)

        five_utr = [_ for _ in utr if _[1] < orf_start]
        three_utr = [_ for _ in utr if _[0] > orf_end]
        assert (len(five_utr) + len(three_utr)) == len(utr), (utr,
                                                              five_utr,
                                                              three_utr,
                                                              cds)

        if transcript.strand == "-":
            orf_start, orf_end = orf_end, orf_start
            five_utr, three_utr = three_utr, five_utr

        if phase_orf[0] != 0 and five_utr:
            raise InvalidCDS("Truncated ORF at 5' with 5' UTR for {0}".format(
                transcript.id))

        remainder = total % 3 - phase_orf[0]
        transcript.logger.debug("In %s Remainder %d, three_utr %s (UTR %s)",
                                transcript.id, remainder, three_utr, utr)
        if remainder > 0 and three_utr:
            raise InvalidCDS("Truncated ORF at 3' with 3' UTR for {0}".format(
                transcript.id))

    if transcript.strand == "-":
        new_orf = list(reversed(new_orf))

    transcript.internal_orfs[index] = new_orf
    assert all((len(_) == 3) for _ in new_orf if new_orf[0] == "CDS"), new_orf
    return transcript


def __check_phase_correctness(transcript):

    """
    This method verifies that the phases are assigned correctly in the case of a coding transcript.
    :param transcript: the input transcript.
    :type transcript: Mikado.loci.transcript.Transcript
    :return: Mikado.loci.transcript.Transcript
    """

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
        if any(_[0] == "CDS" for _ in transcript.segments):
            transcript.internal_orfs = [transcript.segments]
        else:
            transcript.selected_internal_orf_index = None
    else:
        assert len(transcript.internal_orfs) > 0

    for orf_index in range(len(transcript.internal_orfs)):
        try:
            transcript.logger.debug("ORF #%d: %s", orf_index, transcript.phases)
            transcript = __check_internal_orf(transcript,
                                              orf_index,
                                              phases=transcript.phases)
        except (InvalidTranscript, InvalidCDS) as exc:
            transcript.logger.exception(exc)
            transcript.logger.warning("Stripping the CDS from %s", transcript.id)
            transcript.strip_cds()
            break

    # Necessary to set it to the default value
    if len(transcript.internal_orfs) > 0:
        transcript.selected_internal_orf_index = 0
        _ = transcript.selected_internal_orf


def finalize(transcript):
    """Function to calculate the internal introns from the exons.
    In the first step, it will sort the exons by their internal coordinates.

    :param transcript: the Transcript instance to finalize.
    :type transcript: Mikado.loci.transcript.Transcript

    """

    if transcript.finalized is True:
        return

    __previous = transcript.deepcopy()

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

    transcript.combined_cds = sorted(transcript.combined_cds,
                                     key=operator.itemgetter(0, 1))

    transcript.combined_utr = sorted(transcript.combined_utr,
                                     key=operator.itemgetter(0, 1))

    try:
        __check_completeness(transcript)
        __verify_boundaries(transcript)

        assert all([segment[1] in transcript.exons for segment in transcript.segments if
                    segment[0] == "exon"]), (transcript.exons, transcript.segments)

        __check_phase_correctness(transcript)
        __calculate_introns(transcript)
    except (InvalidTranscript, InvalidCDS):
        del transcript
        transcript = __previous
        return

    if len(transcript.combined_cds) > 0:
        transcript.feature = "mRNA"
    else:
        transcript.feature = "transcript"

    if len(transcript.combined_cds) == 0:
        transcript.selected_internal_orf_cds = tuple([])
    else:
        assert isinstance(transcript.selected_internal_orf_index, int)
        transcript.selected_internal_orf_cds = tuple(
            internal_cds for internal_cds in transcript.internal_orfs[
                transcript.selected_internal_orf_index] if
            internal_cds[0] == "CDS")

    # Create the interval tree
    transcript.cds_tree = intervaltree.IntervalTree([
        intervaltree.Interval(cds[0]-1, cds[1]+1) for cds in transcript.combined_cds])

    # BUG somewhere ... I am not sorting this properly before (why?)
    transcript.exons = sorted(transcript.exons)
    # transcript = __calc_cds_introns(transcript)

    transcript.finalized = True
    transcript.logger.debug("Finished finalising %s", transcript.id)

    return
