"""
This module provides the functions needed to check a transcript for consinstency,
e.g. reliability of the CDS/UTR, sanity of borders, etc.
"""

from ...utilities.intervaltree import Interval, IntervalTree
from ...utilities.overlap import overlap
from collections import defaultdict
import operator
# from sys import intern
from ...exceptions import InvalidCDS, InvalidTranscript

__author__ = 'Luca Venturini'


def __basic_final_checks(transcript):

    """
    Function that verifies minimal criteria of a transcript before finalising.
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :return:
    """

    _exons = transcript.exons

    if not _exons:
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
            _exons = sorted([tuple([int(exon[0]), int(exon[1])]) for exon in transcript.combined_cds])
            if len(transcript.combined_utr) == 0:
                # Enlarge the terminal exons to include the starts
                if transcript.start is not None:
                    _exons[0] = (transcript.start, _exons[0][1])
                if transcript.end is not None:
                    _exons[-1] = (_exons[-1][0], transcript.end)
            else:
                __utr = sorted([tuple([int(exon[0]), int(exon[1])]) for exon in transcript.combined_utr])
                try:
                    __before = [_ for _ in __utr if _[1] < _exons[0][0]]
                    if __before and __before[-1][1] == _exons[0][0] - 1:
                        _exons[0] = (__before[-1][0], _exons[0][1])
                        __before.pop()
                    __after = [_ for _ in __utr if _[0] > _exons[-1][1]]
                    if __after and __after[0][0] == _exons[-1][1] + 1:
                        _exons[-1] = (_exons[-1][0], __after[0][1])
                        __after = __after[1:]
                    _exons = __before + _exons + __after
                except IndexError:
                    exc = InvalidTranscript("Transcript {} has a mangled CDS/UTR annotation. Please revise it.")
                    transcript.logger.exception(exc)
                    raise exc

    transcript.logger.debug("Converting to tuples")
    _exons = [tuple([int(exon[0]), int(exon[1])]) for exon in _exons]

    new_exons = []
    # invalid = False

    # Set the start and end automatically if none has been explicitly provided
    if transcript.start is None:
        transcript.start = min(_[0] for _ in _exons)
    if transcript.end is None:
        transcript.end = max(_[1] for _ in _exons)

    for exon in _exons:
        if not isinstance(exon, tuple):
            if (isinstance(exon, Interval) or
                    (isinstance(exon, list) and len(exon) == 2 and
                     isinstance(exon[0], int) and isinstance(exon[1], int))):
                exon = tuple([exon])
            else:
                raise ValueError("Invalid exon: {0}, type {1}".format(
                    exon, type(exon)))
        if exon[0] < transcript.start or exon[1] > transcript.end:
            exc = InvalidTranscript("{} for {} is an invalid exon (start {}, end {})".format(
                exon, transcript.id, transcript.start, transcript.end))
            transcript.logger.debug(exc)
            raise exc
        new_exons.append(exon)

    transcript._set_exons(sorted(new_exons))

    if len(transcript.exons) > 1 and transcript.strand is None:
        if transcript._accept_undefined_multi is False:

            exc = InvalidTranscript(
                "Multiexonic transcripts must have a defined strand! Error for {0}".format(
                    transcript.id))
            transcript.logger.exception(exc)
            raise exc
        else:
            transcript.strand = "?"

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
    cdna_length = transcript.cdna_length
    if cdna_length > transcript.combined_utr_length + transcript.combined_cds_length:
        if transcript.combined_utr == transcript.combined_cds == []:
            # non-coding transcript
            transcript.logger.debug("%s is non coding, returning", transcript.id)
            return
        assert transcript.combined_cds != []

        transcript.logger.debug("Recalculating the UTR for %s. Reason: cDNA length %s, UTR %s, CDS %s (total %s)",
                                transcript.id, cdna_length, transcript.combined_utr_length,
                                transcript.combined_cds_length,
                                transcript.combined_utr_length + transcript.combined_cds_length)
        transcript.combined_utr = []  # Reset
        transcript.combined_cds = sorted(transcript.combined_cds,
                                         key=operator.itemgetter(0, 1))

        exons = IntervalTree.from_intervals([Interval(*exon) for exon in transcript.exons])

        mapper = defaultdict(list)
        for cds in transcript.combined_cds:
            fexon = exons.find(cds[0] - 1, cds[1], strict=False)
            if len(fexon) > 1:
                raise InvalidCDS(
                    "{} has a CDS ({}) which straddles {} different exons ({}).".format(
                        transcript.id, cds, len(fexon), fexon
                    )
                )
            elif len(fexon) == 0:
                raise InvalidCDS(
                    "{} has a CDS ({}) which is not mapped to any exon.".format(
                        transcript.id, cds, len(fexon), fexon
                    )
                )
            mapper[fexon[0]].append(cds)

        for exon in transcript.exons:
            if exon not in mapper:
                transcript.combined_utr.append(exon)
                continue
            elif len(mapper[exon]) == 1:
                cds = mapper[exon][0]
                if cds[0] == exon[0] and exon[1] == cds[1]:
                    continue
                else:
                    before = None
                    after = None
                    if cds[0] < exon[0] or cds[1] > exon[1]:
                        raise InvalidCDS(
                            "{} in {} debords its exon {}".format(cds, transcript.id, exon)
                        )
                    if cds[0] > exon[0]:
                        before = (exon[0], max(cds[0] - 1, exon[0]))
                        transcript.combined_utr.append(before)
                    if cds[1] < exon[1]:
                        after = (min(cds[1] + 1, exon[1]), exon[1])
                        transcript.combined_utr.append(after)
                    assert before or after, (exon, cds)
            else:
                transcript.logger.debug("Starting to find the UTRs for %s", exon)
                found = sorted(mapper[exon])
                utrs = []
                for pos, interval in enumerate(found):
                    if pos == len(found) - 1:
                        if exon[1] > interval[1]:
                            utrs.append((min(exon[1], interval[1] + 1), exon[1]))
                        continue
                    if pos == 0 and exon[0] < interval[0]:
                        utrs.append((exon[0], max(exon[0], interval[0] - 1)))
                    next_interval = found[pos + 1]
                    if not (interval[1] + 1 <= next_interval[0] - 1):
                        raise InvalidCDS(
                            "Error while inferring the UTR for a transcript with multiple ORFs: overlapping CDS found.")
                    utrs.append((interval[1] + 1, next_interval[0] - 1))
                assert utrs, found
                utr_sum = sum([_[1] - _[0] + 1 for _ in utrs])
                cds_sum = sum(_[1] - _[0] + 1 for _ in found)
                assert utr_sum + cds_sum == exon[1] - exon[0] + 1, (utr_sum, cds_sum,
                                                                    exon[1] - exon[0] + 1, utrs, found)
                transcript.combined_utr.extend(utrs)

        # If no CDS and no UTR are present, all good
        equality_one = (transcript.combined_cds_length == transcript.combined_utr_length == 0)
        # Otherwise, if cDNA length == UTR + CDS, all good
        equality_two = (cdna_length ==
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
                    "Overlapping exons found for\n{0} {1}/{2}\n{3}".format(
                        transcript.id, exona, exonb, transcript.exons))
                transcript.logger.debug(exc)
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
            assert first != second, (transcript.id, transcript.selected_cds)
            # assert first[1] < second[0], (first, second)
            first, second = sorted([first, second])
            intron = tuple([first[1] + 1, second[0] - 1])
            if intron not in transcript.introns:
                continue
            # assert intron in transcript.introns, (transcript.id, intron, first, second)
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
        else:
            transcript._combined_cds_introns = transcript._selected_cds_introns.copy()

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


def __calculate_phases(coding: list, previous: int) -> (int, dict):
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

    if transcript._trust_orf is True and index == 0:
        if (transcript.is_coding and transcript.phases) or not transcript.is_coding:
            new_orf = []
            for segment in transcript.internal_orfs[index]:
                if segment[0] == "CDS":
                    segment = tuple([segment[0], segment[1], transcript.phases[segment[1]]])
                new_orf.append(segment)
            transcript.internal_orfs[index] = new_orf
            return transcript
        else:
            pass

    orf, new_orf = transcript.internal_orfs[index], []

    exons = sorted(transcript.exons, reverse=(transcript.strand == "-"))

    coding = sorted([_ for _ in orf if _[0] == "CDS"], key=operator.itemgetter(1))
    transcript.logger.debug("ORF for %s: %s", transcript.id, coding)

    if not coding:
        err = "No ORF for {}, index {}!".format(transcript.id, index)
        transcript.logger.debug(err)
        raise InvalidCDS(err)

    before = sorted([_ for _ in orf
                     if _[0] == "UTR" and _[1][1] < coding[0][1][0]], key=operator.itemgetter(1))
    after = sorted([_ for _ in orf
                    if _[0] == "UTR" and _[1][0] > coding[-1][1][1]], key=operator.itemgetter(1))

    first = min(coding[0][1][0], float("inf") if not before else before[0][1][0])
    last = max(coding[-1][1][1], float("-inf") if not after else after[-1][1][1])

    if first != transcript.start or last != transcript.end:
        # Let's try to salvage this.
        err = """Invalid start and stop of the ORF for {}
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
                   {},
                   # transcript.__dict__
                   )
        transcript.logger.debug(err)
        raise InvalidCDS(err)

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

    phase_orf = []
    for _ in coding:
        if len(_) == 3:
            if _[2] not in (None, 0, 1, 2):
                raise ValueError("Invalid phase value for {}".format(transcript.id))
            phase_orf.append(_[2])
        elif len(_) == 2:
            continue
        else:
            raise ValueError("Invalid CDS fragment: {}".format(_))

    if len(phase_orf) != 0 and len(phase_orf) != len(coding):
        transcript.logger.warning("Invalid phases for %s. Resetting.", transcript.id)
        phase_orf = []

    if not phase_orf and transcript.phases:
        phases_keys = sorted(transcript.phases.keys(), reverse=(transcript.strand == "-"))
        phase_orf = [transcript.phases[_] for _ in phases_keys]
        # Calculating the complement of the phase so that
        # previous = (3 - phase_orf[0]) % 3
        previous = phase_orf[0]
        # transcript.logger.warning(previous)
    elif not phase_orf and transcript._first_phase is not None:
        previous = transcript._first_phase
        phase_orf = []
    elif phase_orf:
        previous = phase_orf[0]
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
        raise InvalidCDS(""""Both UTR presents with a truncated ORF (length {}, modulo {}) in {};
5'UTR: {}
3' UTR: {}""".format(total_cds_length, total_cds_length % 3, transcript.id, five_utr, three_utr))
    elif total_cds_length % 3 != 0 and three_utr:
        for num in (0, 1, 2):
            total_cds_length, __calculated_phases = __calculate_phases(coding,
                                                                       num)
            if total_cds_length % 3 == 0:
                break

        if total_cds_length % 3 != 0:
            raise InvalidCDS("Persistently wrong ORF for %s at 5' end", transcript.id)

    # new_phase_orf = [__calculated_phases[_] for _ in phases_keys]
    if ((__calculated_phases[sorted(__calculated_phases.keys(), reverse=(transcript.strand == "-"))[0]] != 0)
            and five_utr):
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

    segments, internal_orfs = transcript.segments, transcript.internal_orfs

    if min(len(segments), len(internal_orfs)) == 0:
        # transcript.logger.debug("Redefining segments for %s", transcript.id)
        # Define exons
        segments = [("exon", tuple([e[0], e[1]])) for e in transcript.exons]
        # Define CDS
        if len(internal_orfs) > 0:
            for orf in internal_orfs:
                for segment in orf:
                    if segment[0] == "exon":
                        continue
                    elif segment[0] == "UTR":
                        segments.append(("UTR", (segment[1][0], segment[1][1])))
                    elif segment[0] == "CDS":
                        segments.append(("CDS", (segment[1][0], segment[1][1])))
        else:
            segments.extend([("CDS", tuple([c[0], c[1]])) for c in transcript.combined_cds])
        # Define UTR segments
        segments.extend([("UTR", tuple([u[0], u[1]])) for u in transcript.combined_utr])
        # Mix and sort
        segments = sorted(segments, key=operator.itemgetter(1, 0))
        # Add to the store as a single entity
        if not internal_orfs and any(_[0] == "CDS" for _ in segments):
            internal_orfs = [segments]
        else:
            transcript.selected_internal_orf_index = None
    elif len(internal_orfs) == 0:
        exception = AssertionError("No internal ORF for {}".format(transcript.id))
        transcript.logger.exception(exception)
        raise exception
    else:
        pass

    transcript.segments, transcript.internal_orfs = segments, internal_orfs

    __orfs_to_remove = []
    for orf_index in range(len(internal_orfs)):
        # transcript.logger.debug("ORF #%d for %s: %s",
        #                         orf_index, transcript.id, transcript.internal_orfs[orf_index])
        # transcript = __check_internal_orf(transcript, orf_index)
        try:
            transcript = __check_internal_orf(transcript, orf_index)
        except (InvalidTranscript, InvalidCDS) as exc:
            transcript.logger.warning("ORF %s of %s is invalid, removing.",
                                      orf_index, transcript.id)
            __orfs_to_remove.append(orf_index)

    __num_orfs = len(internal_orfs)
    if (__num_orfs > 0) and (len(__orfs_to_remove) == __num_orfs):
        transcript.logger.warning("Every ORF of %s is invalid, stripping the CDS", transcript.id)
        transcript.strip_cds(strand_specific=True)
    elif len(__orfs_to_remove):
        transcript.logger.warning("Stripping %s of %s ORFs out of %s",
                                  transcript.id, len(__orfs_to_remove), __num_orfs)
        for orf_index in reversed(sorted(__orfs_to_remove)):
            internal_orfs.pop(orf_index)
        transcript.internal_orfs = internal_orfs
    else:
        pass

    if len(transcript.internal_orfs) > 0:
        transcript.selected_internal_orf_index = 0
        _ = transcript.selected_internal_orf


def _fix_stop_codon(transcript):

    """This private function will fix the CDS and stop codons when the transcript comes from GTF2
    and therefore has, incorrectly, the stop codon outside the CDS."""

    if transcript.strand == "-":
        # We need to check whether the stop codon is actually in the same exon.
        if transcript.stop_codon[-1][1] == transcript.combined_cds[0][0] - 1:
            phase = transcript.phases.pop(transcript.combined_cds[0], None)
            transcript.combined_cds[0] = (transcript.stop_codon.pop(-1)[0],
                                      transcript.combined_cds[0][1])
            transcript.phases[transcript.combined_cds[0]] = phase
        transcript.combined_cds = [tuple(_) for _ in transcript.stop_codon] + transcript.combined_cds
        for pos, utr in enumerate(transcript.combined_utr):
            if utr[0] > transcript.combined_cds[-1][1]:
                continue  # Skip the 5'
            over = overlap(utr, transcript.combined_cds[0])
            if over < 0:
                continue
            elif over > 3:
                raise InvalidCDS("Invalid overlap between UTR and CDS found")
            else:
                if over == utr[1] - utr[0] + 1:  # This is equivalent to a fragment. Remove.
                    transcript.combined_utr[pos] = None
                else:
                    transcript.combined_utr[pos] = (utr[0], max(utr[0], transcript.combined_cds[0][0] - 1))
    else:
        # Expand the last CDS
        if transcript.stop_codon[0][0] == transcript.combined_cds[-1][1] + 1:
            phase = transcript.phases.pop(transcript.combined_cds[-1], None)
            transcript.combined_cds[-1] = (transcript.combined_cds[-1][0],
                                           transcript.stop_codon.pop(0)[1])
            transcript.phases[transcript.combined_cds[-1]] = phase
        transcript.combined_cds.extend([tuple(_) for _ in transcript.stop_codon])
        for pos, utr in enumerate(transcript.combined_utr):
            if utr[1] < transcript.combined_cds[0][0]:
                continue  # Skip the 5'
            over = overlap(utr, transcript.combined_cds[-1])
            if over < 0:
                continue
            elif over > 3:
                raise InvalidCDS("Invalid overlap between UTR and CDS found")
            else:
                if over == utr[1] - utr[0] + 1:  # This is equivalent to a fragment. Remove.
                    transcript.combined_utr[pos] = None
                else:
                    transcript.combined_utr[pos] = (min(utr[1], transcript.combined_cds[-1][1] + 1),
                                                    utr[1])
    transcript.combined_utr = [_ for _ in transcript.combined_utr if _ is not None]  # Remove the deleted UTRs
    return transcript


def finalize(transcript):
    """Function to calculate the internal introns from the exons.
    In the first step, it will sort the exons by their internal coordinates.

    :param transcript: the Transcript instance to finalize.
    :type transcript: Mikado.loci.transcript.Transcript

    """

    if transcript.finalized is True:
        return

    transcript.exons = sorted(transcript.exons)

    # Add the stop codon to the CDS
    if transcript.stop_codon:
        transcript.logger.debug("Adding the stop codon to %s", transcript.id)
        transcript.stop_codon = sorted(transcript.stop_codon)
        transcript.combined_cds = sorted(transcript.combined_cds)
        transcript.combined_utr = sorted(transcript.combined_utr)
        stop_in_cds = True
        if transcript.strand == "-" and transcript.combined_cds[0][0] > transcript.stop_codon[-1][1]:
            stop_in_cds = False
        elif transcript.combined_cds[-1][1] < transcript.stop_codon[0][0]:
            stop_in_cds = False
        if not stop_in_cds:
            # Here comes the complicated part
            try:
                transcript = _fix_stop_codon(transcript)
            except InvalidCDS:
                transcript.strip_cds()

    transcript.__cdna_length = None
    __basic_final_checks(transcript)
    # Sort the exons by start then stop

    try:
        _check_cdna_vs_utr(transcript)
    except InvalidCDS as exc:
        if transcript.combined_cds:
            transcript.logger.debug(
                "Possible faulty UTR annotation for %s, trying to recalculate it.",
                transcript.id)
            transcript.logger.debug(exc)
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
    transcript._calculate_segment_tree()
    transcript._calculate_cds_tree()

    # BUG somewhere ... I am not sorting this properly before (why?)
    # transcript.exons = sorted(transcript.exons)
    transcript.logger.debug("Checking various attributes")
    if transcript.has_start_codon is None:
        if "has_start_codon" in transcript.attributes:
            transcript.logger.debug("%s has start codon attribute (%s)", transcript.id,
                                    transcript.attributes["has_start_codon"])
            transcript.has_start_codon = transcript.attributes["has_start_codon"]
        else:
            transcript.logger.debug("No predetermined has_start_codon attribute for %s. Attributes: %s",
                                    transcript.id, transcript.attributes)
    else:
        transcript.logger.debug("%s has start codon attribute (%s)", transcript.id, transcript.has_start_codon)

    if transcript.has_stop_codon is None:
        if "has_stop_codon" in transcript.attributes:
            transcript.logger.debug("%s has stop codon attribute (%s)", transcript.id,
                                    transcript.attributes["has_stop_codon"])
            transcript.has_stop_codon = transcript.attributes["has_stop_codon"]
        else:
            transcript.logger.debug("No predetermined has_stop_codon attribute for %s. Attributes: %s",
                                    transcript.id, transcript.attributes)
    else:
        transcript.logger.debug("%s has stop codon attribute (%s)", transcript.id, transcript.has_stop_codon)

    for prop in list(transcript.attributes.keys()):
        val = transcript.attributes[prop]
        if hasattr(transcript, prop):
            try:
                setattr(transcript, prop, val)
            except AttributeError:  # Some instance attributes CANNOT be set from the attributes of the GTF
                transcript.attributes.pop(prop)

    _ = transcript.cdna_length
    transcript._set_basic_lengths()
    transcript._set_distances()
    transcript.finalized = True
    transcript.logger.debug("Finished finalising %s", transcript.id)

    return
