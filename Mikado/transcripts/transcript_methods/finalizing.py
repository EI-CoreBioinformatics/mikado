"""
This module provides the functions needed to check a transcript for consinstency,
e.g. reliability of the CDS/UTR, sanity of borders, etc.
"""

from Mikado.utilities.intervaltree import Interval
import operator
from Mikado.exceptions import InvalidCDS, InvalidTranscript

__author__ = 'Luca Venturini'


def __basic_final_checks(transcript):

    """
    Function that verifies minimal criteria of a transcript before finalising.
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :return:
    """

    if len(transcript.exons) == 0:
        if transcript._possibly_without_exons is True:
            transcript.logger.debug("Inferring that %s is a single-exon transcript")
            new_exon = (transcript.start, transcript.end)
            transcript.add_exon(new_exon)

        elif len(transcript.combined_cds) == 0:
            exc=InvalidTranscript(
                "No exon defined for the transcript {0}. Aborting".format(transcript.id))
            transcript.logger.exception(exc)
            raise exc
        else:
            # Let us try to derive exons from CDS ...
            transcript.exons = sorted([tuple([int(exon[0]), int(exon[1])]) for exon in transcript.combined_cds])
            if len(transcript.combined_utr) == 0:
                # Enlarge the terminal exons to include the starts
                if transcript.start is not None:
                    transcript.exons[0] = (transcript.start, transcript.exons[0][1])
                if transcript.end is not None:
                    transcript.exons[-1] = (transcript.exons[-1][0], transcript.end)
            else:
                __utr = sorted([tuple([int(exon[0]), int(exon[1])]) for exon in transcript.combined_utr])
                try:
                    __before = [_ for _ in __utr if _[1] < transcript.exons[0][0]]
                except IndexError:
                    raise IndexError((__utr, transcript.exons))
                if __before[-1][1] == transcript.exons[0][0] - 1:
                    transcript.exons[0] = (__before[-1][0], transcript.exons[0][1])
                    __before.pop()
                __after = [_ for _ in __utr if _[0] > transcript.exons[-1][1]]
                if __after[0][0] == transcript.exons[-1][1] + 1:
                    transcript.exons[-1] = (transcript.exons[-1][0], __after[0][1])
                    __after = __after[1:]
                transcript.exons = __before + transcript.exons + __after

    transcript.logger.debug("Converting to tuples")
    transcript.exons = [tuple([int(exon[0]), int(exon[1])]) for exon in transcript.exons]

    new_exons = []
    invalid = False

    # Set the start and end automatically if none has been explicitly provided
    if transcript.start is None:
        transcript.start = min(_[0] for _ in transcript.exons)
    if transcript.end is None:
        transcript.end = max(_[1] for _ in transcript.exons)

    for exon in transcript.exons:
        if not isinstance(exon, tuple):
            if (isinstance(exon, Interval) or
                    (isinstance(exon, list) and len(exon) == 2 and
                     isinstance(exon[0], int) and isinstance(exon[1], int))):
                exon = tuple([exon])
            else:
                raise ValueError("Invalid exon: {0}, type {1}".format(
                    exon, type(exon)))
        if exon[0] < transcript.start or exon[1] > transcript.end:
            invalid = True
            break
        new_exons.append(exon)

    transcript.exons = sorted(new_exons)

    if invalid:
        exc = InvalidTranscript("""Exons out of bounds of the transcript:
        ({start}, {end})
        Exons: {exons}""".format(start=transcript.start,
                                 end=transcript.end,
                                 exons=transcript.exons))
        transcript.logger.exception(exc)
        raise exc

    if len(transcript.exons) > 1 and transcript.strand is None:

        exc = InvalidTranscript(
            "Multiexonic transcripts must have a defined strand! Error for {0}".format(
                transcript.id))
        transcript.logger.exception(exc)
        raise exc

    if transcript.combined_utr != [] and transcript.combined_cds == []:

        exc = InvalidTranscript(
            "Transcript {tid} has defined UTRs but no CDS feature!".format(
                tid=transcript.id))
        transcript.logger.exception(exc)
        raise exc

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
            assert isinstance(exon, tuple)
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
                transcript.combined_utr.append(tuple([
                    exon[0], transcript.combined_cds[0][0] - 1]))
            # The first base of the exon is the first base of the last ORF segment:
            # UTR after
            elif (exon[1] > transcript.combined_cds[-1][1] and
                  exon[0] == transcript.combined_cds[-1][0]):
                transcript.combined_utr.append(tuple([
                    transcript.combined_cds[-1][1] + 1, exon[1]]))
            else:
                # If the ORF is contained inside a single exon, with UTR
                # at both sites, then we create the two UTR segments
                if len(transcript.combined_cds) == 1:
                    transcript.combined_utr.append(tuple([
                        exon[0], transcript.combined_cds[0][0] - 1]))
                    transcript.combined_utr.append(tuple([
                        transcript.combined_cds[-1][1] + 1, exon[1]]))
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
                """"Failed to create the UTR:
ID: {}
Exons: {}
Combined CDS: {}
Combined UTR: {}
CDS == UTR == 0: {}
CDNA == CDS + UTR: {}
CDNA == {}
CDS == {}
UTR == {}""".format(transcript.id,
                    transcript.exons,
                    transcript.combined_cds,
                    transcript.combined_utr, equality_one, equality_two,
                    transcript.cdna_length, transcript.combined_cds_length, transcript.combined_utr_length))


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
                exc = InvalidTranscript(
                    "Overlapping exons found!\n{0} {1}/{2}\n{3}".format(
                        transcript.id, exona, exonb, transcript.exons))
                transcript.logger.exception(exc)
                raise exc
            # Append the splice junction
            introns.append(tuple([exona[1] + 1, exonb[0] - 1]))
            # Append the splice locations
            splices.extend([exona[1] + 1, exonb[0] - 1])
    transcript.introns = set(introns)
    transcript.splices = set(splices)

    if (transcript.number_internal_orfs == 0 or
            len(transcript.selected_cds) < 2 or
            len(transcript.combined_cds) < 2):
        pass
    else:
        # Start calculating the selected CDS introns
        for first, second in zip(transcript.selected_cds[:-1],
                                 transcript.selected_cds[1:]):
            assert first != second, transcript.selected_cds
            # assert first[1] < second[0], (first, second)
            first, second = sorted([first, second])
            intron = tuple([first[1] + 1, second[0] - 1])
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
                junc = tuple([former[1] + 1, latter[0] - 1])
                if junc in transcript.introns:
                    cds_introns.append(junc)
            cintrons = set(cds_introns)
            assert len(cintrons) > 0
            transcript._combined_cds_introns = cintrons
            # assert len(transcript._combined_cds_introns) > 0

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

        # assert len(transcript._combined_cds_introns) > 0

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


def __calculate_phases(coding, previous):
    """

    :param coding:
    :param previous:
    :return:
    """

    total_cds_length = -previous

    __calculated_phases = dict()
    for cds_segment in coding:
        length = cds_segment[1][1] - cds_segment[1][0] + 1
        phase = (3 - (total_cds_length % 3)) % 3
        __calculated_phases[cds_segment[1]] = phase
        total_cds_length += length

    return total_cds_length, __calculated_phases


def __check_internal_orf(transcript, index):

    """
    Method that verifies that an internal ORF does not have any internal gap.

    :param transcript: the transcript to analyse
    :type transcript: Mikado.loci.Transcript
    :param index: index of the internal orf to check
    :type index: int

    :return: the updated transcript
    :rtype: Mikado.loci.Transcript
    """

    orf, new_orf = transcript.internal_orfs[index], []

    exons = sorted(transcript.exons, reverse=(transcript.strand == "-"))

    coding = sorted([(_[0], _[1]) for _ in orf if _[0] == "CDS"],
                    key=operator.itemgetter(1))

    if not coding:
        raise InvalidCDS("No ORF for {}, index {}!".format(transcript.id, index))
    before = sorted([_ for _ in orf
                     if _[0] == "UTR" and _[1][1] < coding[0][1][0]], key=operator.itemgetter(1))
    after = sorted([_ for _ in orf
                    if _[0] == "UTR" and _[1][0] > coding[-1][1][1]], key=operator.itemgetter(1))

    first = min(coding[0][1][0], float("inf") if not before else before[0][1][0])
    last = max(coding[-1][1][1], float("-inf") if not after else after[-1][1][1])

    if first != transcript.start or last != transcript.end:
        raise InvalidCDS("""Invalid start and stop of the ORF for {}
First: {} Start: {}
Last: {} End {}
Coding: {}
Before: {}
After: {}



dict: {}""".format(transcript.id,
                   first, transcript.start,
                   last, transcript.end,
                   coding,
                   before,
                   after,
                   transcript.__dict__))

    # Check that the number of exons with a coding section is correct and that they are in the correct order.
    coding_exons = [_ for _ in enumerate(exons) if
                    _[1][1] >= coding[0][1][1] and _[1][0] <= coding[-1][1][0]]
    if len(coding_exons) != len(coding) or coding_exons[-1][0] - coding_exons[0][0] + 1 != len(coding):
        raise InvalidCDS(""""Invalid number of coding exons for {} ({} vs {})
Coding: {}
Coding_exons (recalculated): {}""".format(
            transcript.id,
            len(coding), len(coding_exons),
            coding, coding_exons))

    # Now it's time to check the phases
    if transcript.strand == "-":
        coding = list(reversed(coding))
        five_utr = list(reversed(after))
        three_utr = list(reversed(before))
    else:
        five_utr = before
        three_utr = after

    del before, after

    if index == 0 and transcript.phases:
        phases_keys = sorted(transcript.phases.keys(), reverse=(transcript.strand == "-"))
        phase_orf = [transcript.phases[_] for _ in phases_keys]
        # Calculating the complement of the phase so that
        # previous = (3 - phase_orf[0]) % 3
        previous = phase_orf[0]
        # transcript.logger.warning(previous)
    elif index == 0 and transcript._first_phase is not None:
        previous = transcript._first_phase
        phase_orf = []
    else:
        phase_orf = []
        for segment in sorted(orf, key=operator.itemgetter(1), reverse=(transcript.strand == "-")):
            if segment[0] != "CDS":
                continue
            else:
                if len(segment) == 3:
                    phase_orf.append(segment[2])
                else:
                    break
        if phase_orf and len(phase_orf) == len(coding):
            previous = phase_orf[0]
        else:
            previous = 0
            phase_orf = []

    if transcript._trust_orf is True and index == 0 and len(phase_orf) == len(coding):
        total_cds_length = sum([_[1][1] - _[1][0] + 1 for _ in coding])
        __calculated_phases = phase_orf[:]
    else:
        total_cds_length, __calculated_phases = __calculate_phases(coding, previous)
        new_phases_keys = sorted(__calculated_phases.keys(), reverse=(transcript.strand == "-"))
        new_phase_orf = [__calculated_phases[_] for _ in new_phases_keys]

        if len(__calculated_phases) != len(coding):
            # This is a mistake which should crash the program
            raise ValueError("Error in calculating the phases!")

        if phase_orf and new_phase_orf != phase_orf:
            transcript.logger.debug("Wrong phases for %s, using recalculated ones (\n%s\nvs\n%s)",
                                    transcript.id,
                                    phase_orf, __calculated_phases)
        else:
            transcript.logger.debug("Correct phases for %s: %s",
                                    transcript.id, __calculated_phases)
        if total_cds_length % 3 != 0 and three_utr and five_utr:
            # The transcript is truncated.
            raise InvalidCDS(""""Both UTR presents with a truncated ORF in {}
5'UTR: {}
3' UTR: {}""".format(transcript.id, five_utr, three_utr))
        elif total_cds_length % 3 != 0 and three_utr:
            for num in (0, 1, 2):
                total_cds_length, __calculated_phases = __calculate_phases(coding,
                                                                           num)
                if total_cds_length % 3 == 0:
                    break

            if total_cds_length % 3 != 0:
                raise InvalidCDS("Persistently wrong ORF for %s at 5' end", transcript.id)

        # new_phase_orf = [__calculated_phases[_] for _ in phases_keys]
        if __calculated_phases[sorted(__calculated_phases.keys(), reverse=(transcript.strand == "-"))[0]] != 0 and five_utr:
            raise InvalidCDS("5'UTR present with a truncated ORF at 5' end for {}".format(
                             transcript.id))

    transcript.phases = __calculated_phases

    transcript.logger.debug("Total CDS length %d", total_cds_length)

    new_orf = five_utr[:]
    new_orf.extend([(_[0][0], _[0][1], _[1]) for _ in zip(
        coding,
        [__calculated_phases[_]  for _ in sorted(__calculated_phases.keys(), reverse=(transcript.strand == "-"))])])
    new_orf.extend(three_utr)

    new_orf.extend([("exon", _) for _ in transcript.exons])
    new_orf = sorted(new_orf, key=operator.itemgetter(1, 0))

    transcript.internal_orfs[index] = new_orf
    return transcript


def __check_phase_correctness(transcript):

    """
    This method verifies that the phases are assigned correctly in the case of a coding transcript.
    :param transcript: the input transcript.
    :type transcript: Mikado.loci.transcript.Transcript
    :return: Mikado.loci.transcript.Transcript
    """

    if min(len(transcript.segments), len(transcript.internal_orfs)) == 0:
        transcript.logger.debug("Redefining segments for %s", transcript.id)
        # Define exons
        transcript.segments = [("exon", tuple([e[0], e[1]]))
                               for e in transcript.exons]
        # Define CDS
        transcript.segments.extend([("CDS", tuple([c[0], c[1]]))
                                    for c in transcript.combined_cds])
        # Define UTR segments
        transcript.segments.extend([("UTR", tuple([u[0], u[1]]))
                                    for u in transcript.combined_utr])
        # Mix and sort
        transcript.segments = sorted(transcript.segments, key=operator.itemgetter(1, 0))
        # Add to the store as a single entity
        if any(_[0] == "CDS" for _ in transcript.segments):
            transcript.internal_orfs = [transcript.segments]
        else:
            transcript.selected_internal_orf_index = None
    elif len(transcript.internal_orfs) == 0:
        exception = AssertionError("No internal ORF for {}".format(transcript.id))
        transcript.logger.exception(exception)
        raise exception

    transcript.logger.debug("{} has {} internal ORF{}".format(
        transcript.id, len(transcript.internal_orfs),
        "s" if len(transcript.internal_orfs) > 1 else ""))
    for orf_index in range(len(transcript.internal_orfs)):
        transcript.logger.debug("ORF #%d for %s: %s",
                                orf_index, transcript.id, transcript.phases)
        try:
            transcript = __check_internal_orf(transcript,
                                              orf_index)
        except (InvalidTranscript, InvalidCDS) as exc:
            transcript.logger.warning("Stripping the CDS from %s, error: %s",
                                      transcript.id, exc)
            transcript.strip_cds(strand_specific=True)
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

    # __previous = transcript.deepcopy()

    transcript.exons = sorted(transcript.exons)

    # Add the stop codon to the CDS
    if transcript.stop_codon:
        transcript.logger.debug("Adding the stop codon to %s", transcript.id)
        transcript.stop_codon = sorted(transcript.stop_codon)
        transcript.combined_cds = sorted(transcript.combined_cds)
        if transcript.strand == "-":
            transcript.logger.debug("%s: CDS[0]: %s, Stop codon: %s",
                                    transcript.id,
                                    transcript.combined_cds[0],
                                    transcript.stop_codon)
            if transcript.combined_cds[0][0] == transcript.stop_codon[-1][1] + 1:
                transcript.logger.debug("Moving %s last CDS from %d to %d",
                                        transcript.id,
                                        transcript.combined_cds[0][0],
                                        transcript.stop_codon[-1][0]
                                        )
                transcript.combined_cds[0] = (
                    transcript.stop_codon.pop(-1)[0],
                    transcript.combined_cds[0][1])
            transcript.logger.debug("Extend the CDS with: %s", transcript.stop_codon)
            transcript.combined_cds.extend(transcript.stop_codon)
            transcript.logger.debug("Final CDS: %s", transcript.combined_cds)
        else:
            transcript.logger.debug("%s: CDS[-1]: %s, Stop codon: %s",
                                    transcript.id,
                                    transcript.combined_cds[-1],
                                    transcript.stop_codon)
            if transcript.combined_cds[-1][1] == transcript.stop_codon[0][0] - 1:
                transcript.logger.debug("Moving %s last CDS from %d to %d",
                                        transcript.id,
                                        transcript.combined_cds[-1][1],
                                        transcript.stop_codon[0][1]
                                        )
                transcript.combined_cds[-1] = (
                    transcript.combined_cds[-1][0],
                    transcript.stop_codon.pop(0)[1])
            transcript.logger.debug("Extend the CDS with: %s", transcript.stop_codon)
            transcript.combined_cds.extend(transcript.stop_codon)
            transcript.logger.debug("Final CDS: %s", transcript.combined_cds)

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
            except InvalidCDS as exc:
                transcript.logger.warning("CDS for %s completely invalid. Removing it.",
                                          transcript.id)
                transcript.logger.exception(exc)
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
        transcript.logger.debug("Verifying phase correctness for %s", transcript.id)
        __check_phase_correctness(transcript)
        transcript.logger.debug("Calculating intron correctness for %s", transcript.id)
        __calculate_introns(transcript)
    except (InvalidCDS, InvalidTranscript):
        transcript.finalized = True
        transcript.unfinalize()
        return

    if transcript.feature == "transcript":
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

    # Create the internal trees
    _ = transcript.cds_tree
    _ = transcript.cds_introntree
    _ = transcript.segmenttree

    # BUG somewhere ... I am not sorting this properly before (why?)
    transcript.exons = sorted(transcript.exons)
    # transcript = __calc_cds_introns(transcript)

    transcript.finalized = True
    transcript.logger.debug("Finished finalising %s", transcript.id)

    return
