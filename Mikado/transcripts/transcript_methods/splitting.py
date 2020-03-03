"""
This module contains the methods used by the Transcript class to split an instance into
multiple transcripts, if the conditions are met (multiple ORFs present and BLAST not
supporting them being part of the same transcript).
"""

from sys import version_info
if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict
import collections
import operator
from ...utilities.intervaltree import IntervalTree, Interval
from ...utilities import overlap
from ...exceptions import InvalidTranscript
from ...parsers.blast_utils import merge
from ...parsers.bed12 import BED12

__author__ = 'Luca Venturini'


def check_split_by_blast(transcript, cds_boundaries):

    """
    This method verifies if a transcript with multiple ORFs has support by BLAST to
    NOT split it into its different components.

    The minimal overlap between ORF and HSP is defined inside the JSON at the key
        ["chimera_split"]["blast_params"]["minimal_hsp_overlap"]
    basically, we consider a HSP a hit only if the overlap is over a certain threshold
    and the HSP evalue under a certain threshold.

    The split by CDS can be executed in three different ways - PERMISSIVE, LENIENT, STRINGENT:

    - PERMISSIVE: split if two CDSs do not have hits in common,
    even when one or both do not have a hit at all.
    - STRINGENT: split only if two CDSs have hits and none
    of those is in common between them.
    - LENIENT: split if *both* lack hits, OR *both* have hits and none
    of those is in common.

    :param transcript: the transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript
    :param cds_boundaries:
    :return: cds_boundaries
    :rtype: dict
    """

    # Establish the minimum overlap between an ORF and a BLAST hit to consider it
    # to establish belongingness

    minimal_overlap = transcript.json_conf[
        "pick"]["chimera_split"]["blast_params"]["minimal_hsp_overlap"]

    cds_hit_dict = SortedDict().fromkeys(cds_boundaries.keys())
    for key in cds_hit_dict:
        cds_hit_dict[key] = collections.defaultdict(list)

    # BUG, this is a hacky fix
    if not hasattr(transcript, "blast_hits"):
        transcript.logger.warning(
            "BLAST hits store lost for %s! Creating a mock one to avoid a crash",

            transcript.id)
        transcript.blast_hits = []

    transcript.logger.debug("%s has %d possible hits", transcript.id, len(transcript.blast_hits))

    # Determine for each CDS which are the hits available
    min_eval = transcript.json_conf["pick"]['chimera_split']['blast_params']['hsp_evalue']
    for hit in transcript.blast_hits:
        for hsp in iter(_hsp for _hsp in hit["hsps"] if
                        _hsp["hsp_evalue"] <= min_eval):
            for cds_run in cds_boundaries:
                # If I have a valid hit b/w the CDS region and the hit,
                # add the name to the set
                overlap_threshold = minimal_overlap * (cds_run[1] + 1 - cds_run[0])
                overl = overlap(cds_run, (hsp['query_hsp_start'], hsp['query_hsp_end']))

                if overl >= overlap_threshold:
                    cds_hit_dict[cds_run][(hit["target"], hit["target_length"])].append(hsp)
                    transcript.logger.debug(
                        "Overlap %s passed for %s between %s CDS and %s HSP (threshold %s)",
                        overlap,
                        transcript.id,
                        cds_run,
                        (hsp['query_hsp_start'], hsp['query_hsp_end']),
                        overlap_threshold)
                else:
                    transcript.logger.debug(
                        "Overlap %s rejected for %s between %s CDS and %s HSP (threshold %s)",
                        overlap,
                        transcript.id,
                        cds_run, (hsp['query_hsp_start'], hsp['query_hsp_end']),
                        overlap_threshold)

    transcript.logger.debug("Final cds_hit_dict for %s: %s", transcript.id, cds_hit_dict)

    final_boundaries = SortedDict()
    for boundary in __get_boundaries_from_blast(transcript, cds_boundaries, cds_hit_dict):
        if len(boundary) == 1:
            assert len(boundary[0]) == 2
            boundary = boundary[0]
            final_boundaries[boundary] = cds_boundaries[boundary]
        else:
            nboun = (boundary[0][0], boundary[-1][1])
            final_boundaries[nboun] = []
            for boun in boundary:
                final_boundaries[nboun].extend(cds_boundaries[boun])
    transcript.logger.debug("Final boundaries for %s: %s",
                            transcript.id, final_boundaries)

    cds_boundaries = final_boundaries.copy()
    return cds_boundaries


def check_common_hits(transcript, cds_hits, old_hits):
    """
    This private method verifies whether we have to split a transcript
    if there are hits for both ORFs and some of them refer to the same target.
    To do so, we check whether the two CDS runs actually share at least one HSPs
    (in which case we do NOT want to split); if not, we verify whether the HSPs
    cover a large fraction of the target length. If this is the case, we decide to
    break down the transcript because we are probably in the presence of a tandem
    duplication.

    :param transcript: the transcript instance to analyse
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param cds_hits:
    :param old_hits:
    :return:
    """

    in_common = set.intersection(set(cds_hits.keys()),
                                 set(old_hits.keys()))
    # We do not have any hit in common
    min_overlap_duplication = transcript.json_conf[
        "pick"]['chimera_split']['blast_params']['min_overlap_duplication']
    to_break = True
    for common_hit in in_common:
        old_hsps = old_hits[common_hit]
        cds_hsps = cds_hits[common_hit]
        # First check ... do we have HSPs in common?
        # if len(set.intersection(old_hsps, cds_hsps)) > 0:
        #     to_break = False
        #     break
        old_query_boundaries = IntervalTree.from_tuples(
            [(h["query_hsp_start"], h["query_hsp_end"]) for h in old_hsps])

        # old_query_boundaries = IntervalTree([Interval(h["query_hsp_start"],
        #                                               h["query_hsp_end"])
        #                                      for h in old_hsps])
        # Look for HSPs that span the two ORFs
        if any([len(
                old_query_boundaries.search(new_cds["query_hsp_start"],
                                            new_cds["query_hsp_end"])) > 0]
               for new_cds in cds_hsps):
            to_break = False

        old_target_boundaries = IntervalTree.from_tuples([
            (h["target_hsp_start"], h["target_hsp_end"]) for h in old_hsps])
        for cds_hsp in cds_hsps:
            boundary = (cds_hsp["target_hsp_start"], cds_hsp["target_hsp_end"])
            for target_hit in old_target_boundaries.search(*boundary):
                overlap_fraction = transcript.overlap(boundary,
                                                      target_hit)/common_hit[1]
                transcript.logger.debug("Checking overlap duplication for %s; OF %s, minimum %s",
                                        transcript.id, overlap_fraction, min_overlap_duplication)
                if overlap_fraction >= min_overlap_duplication:
                    to_break = True and to_break
                else:
                    to_break = False
    return to_break


def __get_boundaries_from_blast(transcript, cds_boundaries, cds_hit_dict):

    """
    Private method that calculates the CDS boundaries to keep
    given the blast hits. Called by check_split_by_blast
    :param cds_boundaries:
    :return:
    """
    new_boundaries = []
    leniency = transcript.json_conf["pick"]['chimera_split']['blast_params']['leniency']
    for cds_boundary in cds_boundaries:
        if not new_boundaries:
            new_boundaries.append([cds_boundary])
        else:
            old_boundary = new_boundaries[-1][-1]
            cds_hits = cds_hit_dict[cds_boundary]
            old_hits = cds_hit_dict[old_boundary]
            if len(cds_hits) == len(old_hits) == 0:  # No hit found for either CDS
                # If we are stringent, we DO NOT split
                if leniency == "STRINGENT":
                    new_boundaries[-1].append(cds_boundary)
                else:  # Otherwise, we do split
                    new_boundaries.append([cds_boundary])
            elif min(len(cds_hits), len(old_hits)) == 0:  # We have hits for only one
                # If we are permissive, we split
                if leniency == "PERMISSIVE":
                    new_boundaries.append([cds_boundary])
                else:
                    new_boundaries[-1].append(cds_boundary)
            else:
                if check_common_hits(transcript, cds_hits, old_hits) is True:
                    new_boundaries.append([cds_boundary])
                # We have hits in common
                else:
                    new_boundaries[-1].append(cds_boundary)
        # } # Finish BLAST check
    return new_boundaries


def __split_complex_exon(transcript, exon, texon, sentinel, boundary, invert=False):

    """
    Private method used to split an exon when it is only partially coding,
    :param exon: Exon to be analysed
    :param texon: Transcriptomic coordinates of the exon
    :param sentinel: tuple of boolean flags, it indicates whether there are transcripts on the left
    (first member) and/or on the right (second member) of the current instance.

    :param boundary: Transcriptomic coordinates of the ORF boundary.
    :return:
    """

    left, right = sentinel

    to_discard = None
    new_exon = list(exon)

    if texon[1] == boundary[0]:
        # In this case we have that the exon ends exactly at the end of the
        # UTR, so we have to keep a one-base exon
        if left is False:
            transcript.logger.debug("Appending mixed UTR/CDS 5' exon %s", exon)
        else:
            if transcript.strand == "+":
                # Keep only the LAST base
                to_discard = (exon[0], exon[1]-1)
                new_exon = (exon[1]-1, exon[1])
                texon = (texon[1]-1, texon[1])
                transcript.logger.debug("Appending monobase CDS exon %s (Texon %s)",
                                        new_exon,
                                        texon)
            else:
                # Keep only the FIRST base
                to_discard = (exon[0]+1, exon[1])
                new_exon = (exon[0], exon[0]+1)
                texon = (texon[1]-1, texon[1])
                transcript.logger.debug(
                    "Appending monobase CDS exon %s (Texon %s)",
                    new_exon,
                    texon)
    elif texon[0] == boundary[1]:
        # In this case we have that the exon ends exactly at the end of the
        # CDS, so we have to keep a one-base exon
        if right is False:
            transcript.logger.debug(
                "Appending mixed UTR/CDS right exon %s",
                exon)
        else:
            if transcript.strand == "+":
                # In this case we have to keep only the FIRST base
                to_discard = (exon[0]+1, exon[1])
                new_exon = (exon[0], exon[0]+1)
                texon = (texon[0], texon[0]+1)
                transcript.logger.debug(
                    "Appending monobase CDS exon %s (Texon %s)",
                    new_exon,
                    texon)
            else:
                # In this case we have to keep only the LAST base
                to_discard = (exon[0], exon[1]-1)
                new_exon = (exon[1]-1, exon[1])
                texon = (texon[0], texon[0]+1)
                transcript.logger.debug(
                    "Appending monobase CDS exon %s (Texon %s)",
                    new_exon,
                    texon)

    elif texon[0] <= boundary[0] <= boundary[1] <= texon[1]:
        # Monoexonic
        transcript.logger.debug("Exon %s, case 3.1", exon)

        # if transcript.monoexonic is False:
        if invert is True:
            if left is True:
                new_exon[1] = exon[0] + (texon[1] - boundary[0])
                transcript.logger.debug(
                    "Case 3.1: Negative strand, another transcript on the left, new exon: %d, %d",
                    new_exon[0], new_exon[1])
            if right is True:
                new_exon[0] = exon[1] - (boundary[1] - texon[0])
                transcript.logger.debug(
                    "Case 3.1: Negative strand, another transcript on the right, new exon: %d, %d",
                    new_exon[0], new_exon[1])
        else:
            if left is True:
                new_exon[0] = exon[1] - (texon[1] - boundary[0])
                transcript.logger.debug(
                    "Case 3.1: Positive strand, another transcript on the left, new exon: %d, %d",
                    new_exon[0], new_exon[1])
            if right is True:
                new_exon[1] = exon[0] + (boundary[1] - texon[0])
                transcript.logger.debug(
                    "Case 3.1: Positive strand, another transcript on the right, new exon: %d, %d",
                    new_exon[0], new_exon[1])

        transcript.logger.debug(
            "[Monoexonic] Tstart shifted for %s, %d to %d",
            transcript.id, texon[0], boundary[0])
        transcript.logger.debug(
            "[Monoexonic] GStart shifted for %s, %d to %d",
            transcript.id, exon[0], new_exon[1])
        transcript.logger.debug(
            "[Monoexonic] Tend shifted for %s, %d to %d",
            transcript.id, texon[1], boundary[1])
        transcript.logger.debug(
            "[Monoexonic] Gend shifted for %s, %d to %d",
            transcript.id, exon[1], new_exon[1])

        if left is True:
            texon[0] = boundary[0]
        if right is True:
            texon[1] = boundary[1]

    elif texon[0] <= boundary[0] <= texon[1] <= boundary[1]:
        # In this case we have that exon is sitting halfway
        # i.e. there is a partial 5'UTR
        if left is True:
            if transcript.strand == "-":
                new_exon[1] = exon[0] + (texon[1] - boundary[0])
            else:
                new_exon[0] = exon[1] - (texon[1] - boundary[0])
            transcript.logger.debug(
                "Tstart shifted for %s, %d to %d", transcript.id, texon[0], boundary[0])
            transcript.logger.debug(
                "GStart shifted for %s, %d to %d", transcript.id, exon[0], new_exon[1])
            texon[0] = boundary[0]

    elif texon[1] >= boundary[1] >= texon[0] >= boundary[0]:
        # In this case we have that exon is sitting halfway
        # i.e. there is a partial 3'UTR
        if right is True:
            if transcript.strand == "-":
                new_exon[0] = exon[1] - (boundary[1] - texon[0])
            else:
                new_exon[1] = exon[0] + (boundary[1] - texon[0])
            transcript.logger.debug(
                "Tend shifted for %s, %d to %d",
                transcript.id, texon[1], boundary[1])
            transcript.logger.debug(
                "Gend shifted for %s, %d to %d",
                transcript.id, exon[1], new_exon[1])
            texon[1] = boundary[1]
        else:
            transcript.logger.debug("New exon: %s", new_exon)
            transcript.logger.debug("New texon: %s", texon)

    # Prevent monobase exons
    # if new_exon[0] == new_exon[1]:
    #     new_exon[1] += 1
    new_exon = tuple([new_exon[0], new_exon[1]])

    return new_exon, texon, to_discard


def __create_splitted_exons(transcript, boundary, left, right, orf_strand):

    """
    Given a boundary in transcriptomic coordinates, this method will extract the
    exons retained in the splitted part of the model.

    :param boundary: the *transcriptomic* coordinates of start/end of the ORF(s)
    to be included in the new transcript
    :type boundary: (int,int)

    :param left: boolean flag indicating whether there is another sub-transcript
    to the left of the one we mean to create, irrespective of *genomic* strand
    :type left: bool

    :param left: boolean flag indicating whether there is another sub-transcript
    to the right of the one we mean to create, irrespective of *genomic* strand
    :type right: bool


    :return: my_exons (final exons), discarded_exons (eventual discarded exons),
    tstart (new transcript start), tend (new transcript end)
    :rtype: (list(int,int),list(int,int),int,int)
    """

    my_exons = []

    discarded_exons = []
    tlength = 0
    tstart = float("Inf")
    tend = float("-Inf")

    if transcript.strand == "-":
        reversal = True
    else:
        reversal = False

    transcript.logger.debug("Starting analysis on %s, boundaries %s, left: %s, \
right: %s, reversal: %s",
                            transcript.id, boundary, left, right, reversal)

    for exon in sorted(transcript.exons, key=operator.itemgetter(0), reverse=reversal):
        # Translate into transcript coordinates
        elength = exon[1] - exon[0] + 1
        texon = [tlength + 1, tlength + elength]
        tlength += elength
        transcript.logger.debug("Analysing exon %s [%s] for %s",
                                exon, texon, transcript.id)

        # SIMPLE CASES
        # Exon completely contained in the ORF
        if boundary[0] <= texon[0] < texon[1] <= boundary[1]:
            transcript.logger.debug("Appending CDS exon %s", exon)
            my_exons.append(exon)
        # Exon on the left of the CDS
        elif texon[1] < boundary[0]:
            if left is False:
                transcript.logger.debug("Appending 5'UTR exon %s", exon)
                my_exons.append(exon)
            else:
                transcript.logger.debug("Discarding 5'UTR exon %s", exon)
                discarded_exons.append(exon)
                continue
        elif texon[0] > boundary[1]:
            if right is False:
                transcript.logger.debug("Appending 3'UTR exon %s", exon)
                my_exons.append(exon)
            else:
                transcript.logger.debug("Discarding 3'UTR exon %s", exon)
                discarded_exons.append(exon)
                continue
        else:
            # exon with partial UTR, go to the relative function
            # to handle these complex cases
            assert transcript.strand is not None
            exon, texon, to_discard = __split_complex_exon(
                transcript, exon, texon, (left, right), boundary,
                invert=(transcript.strand != orf_strand))

            my_exons.append(exon)

            if to_discard is not None:
                discarded_exons.append(to_discard)

        tstart = min(tstart, texon[0])
        tend = max(tend, texon[1])

    return my_exons, discarded_exons, tstart, tend


def __create_splitted_transcripts(transcript, cds_boundaries):

    """
    Private method called by split_by_cds to create the various (N>1) transcripts
    that are its output.
    :param cds_boundaries: a list of int tuples, containing the boundaries
     of the new transcripts.
    :return:
    """

    spans = []
    new_transcripts = []

    for counter, (boundary, bed12_objects) in enumerate(
            sorted(cds_boundaries.items(), key=operator.itemgetter(0))):
        new_transcript = transcript.__class__()
        new_transcript.feature = "mRNA"
        for attribute in ["chrom", "source", "score", "strand", "attributes"]:
            setattr(new_transcript, attribute, getattr(transcript, attribute))
        # Determine which ORFs I have on my right and left
        new_transcript.parent = transcript.parent
        left = True
        right = True
        if counter == 0:  # leftmost
            left = False
        if 1 + counter == len(cds_boundaries):  # rightmost
            right = False
        counter += 1  # Otherwise they start from 0
        new_transcript.id = "{0}.split{1}".format(transcript.id, counter)
        new_transcript.external_scores.update(transcript.external_scores.items())
        new_transcript.logger = transcript.logger
        bed12_strand = set(_.strand for _ in bed12_objects)
        assert len(bed12_strand) == 1
        bed12_strand = bed12_strand.pop()

        transcript.logger.debug("Splitting exons for %s", new_transcript.id)
        my_exons, discarded_exons, tstart, tend = __create_splitted_exons(
            transcript, boundary, left, right, bed12_strand)

        transcript.logger.debug("""TID %s counter %d, boundary %s, left %s right %s""",
                                transcript.id,
                                counter,
                                boundary,
                                left,
                                right)

        if right is True:
            transcript.logger.debug("TID %s TEND %d Boun[1] %s",
                                    transcript.id, tend, boundary[1])
        if left is True:
            transcript.logger.debug("TID %s TSTART %d Boun[0] %s",
                                    transcript.id, tstart, boundary[0])

        assert len(my_exons) > 0, (discarded_exons, boundary)

        new_transcript.exons = my_exons

        new_transcript.start = min(exon[0] for exon in new_transcript.exons)
        new_transcript.end = max(exon[1] for exon in new_transcript.exons)
        assert new_transcript.end <= transcript.end
        assert new_transcript.start >= transcript.start
        assert new_transcript.is_coding is False
        new_transcript.json_conf = transcript.json_conf
        # Now we have to modify the BED12s to reflect
        # the fact that we are starting/ending earlier
        new_transcript.finalize()
        # if transcript.monoexonic is True:
        # if new_transcript.monoexonic is True:
        #     new_transcript.strand = None

        transcript.logger.debug("Relocating %d ORFs into the new transcript (%d, %d), \
tcoordinates (%d, %d)",
                                len(bed12_objects),
                                new_transcript.start, new_transcript.end,
                                tstart, tend
                                )
        new_bed12s = __relocate_orfs(transcript, bed12_objects, tstart, tend)
        assert len([_ for _ in new_bed12s if _.strand == "+"]) > 0
        transcript.logger.debug("Loading %d ORFs into the new transcript (%d, %d): %s",
                                len(new_bed12s),
                                new_transcript.start, new_transcript.end,
                                "\n\t"+"\n".join([str(_) for _ in new_bed12s]))
        new_transcript.logger = transcript.logger
        new_transcript.load_orfs(new_bed12s)

        if new_transcript.selected_cds_length <= 0:
            err_message = "No CDS information retained for {0} split {1}\n".format(
                transcript.id, counter)
            err_message += "BED: {0}".format("\n\t".join([str(x) for x in new_bed12s]))
            raise InvalidTranscript(err_message)

        # Load the blast hits
        __load_blast_hits(new_transcript, boundary, transcript)
        new_transcript.finalize()
        new_transcripts.append(new_transcript)
        nspan = (new_transcript.start, new_transcript.end)
        transcript.logger.debug(
            "Transcript {0} split {1}, discarded exons: {2}".format(
                transcript.id, counter, discarded_exons))
        __check_collisions(transcript, nspan, spans)
        spans.append([new_transcript.start, new_transcript.end])

    return new_transcripts


def __load_blast_hits(new_transcript, boundary, transcript):

    """
    Function to load the BLAST hits into the new splitted transcript.
    :param new_transcript: the splitted transcript
    :type new_transcript: Mikado.loci_objects.Transcript
    :param boundary: tuple(start, end) of the boundary of the new transcript
    :type boundary: tuple(int, int)
    :param transcript:  the original transcript
    :type transcript: Mikado.loci_objects.Transcript
    :return:
    """

    for hit in transcript.blast_hits:
        if overlap((hit["query_start"], hit["query_end"]), boundary) > 0:

            minimal_overlap = transcript.json_conf[
                "pick"]["chimera_split"]["blast_params"]["minimal_hsp_overlap"]
            new_hit = __recalculate_hit(hit, boundary, minimal_overlap)
            if new_hit is not None:
                transcript.logger.debug("""Hit %s,
                                        previous id/query_al_length/t_al_length %f/%f/%f,
                                        novel %f/%f/%f""",
                                        new_hit["target"],
                                        hit["global_identity"],
                                        hit["query_aligned_length"],
                                        hit["target_aligned_length"],
                                        new_hit["global_identity"],
                                        new_hit["query_aligned_length"],
                                        new_hit["target_aligned_length"])

                new_transcript.blast_hits.append(new_hit)
            else:
                transcript.logger.debug("Hit %s did not pass overlap checks for %s",
                                        hit["target"], new_transcript.id)
        else:
            transcript.logger.debug("Ignoring hit %s as it is not intersecting", hit)
            continue


def __check_collisions(transcript, nspan, spans):

    """
    This method checks whether a new transcript collides with a previously
    defined transcript.
    :param nspan:
    :param spans:
    :return:
    """

    if len(spans) == 0:
        return
    for span in spans:
        overl = overlap(span, nspan)

        transcript.logger.debug(
            "Comparing start-ends for split of %s. SpanA: %s SpanB: %s Overlap: %d",
            transcript.id, span,
            nspan, overl)

        if overl > 0:
            err_message = "Invalid overlap for {0}! T1: {1}. T2: {2}".format(
                transcript.id, span, nspan)
            transcript.logger.error(err_message)
            raise InvalidTranscript(err_message)


def __recalculate_hit(hit, boundary, minimal_overlap):
    """Static method to recalculate coverage/identity for new hits."""

    __valid_matches = set([chr(x) for x in range(65, 91)] + [chr(x) for x in range(97, 123)] +
                          ["|"])

    hit_dict = dict()
    for key in iter(k for k in hit.keys() if k not in ("hsps",)):
        hit_dict[key] = hit[key]

    hsp_dict_list = []
    # hit_dict["global_identity"] = []
    q_intervals = []
    t_intervals = []

    identical_positions, positives = set(), set()

    best_hsp = (float("inf"), float("-inf"))

    for hsp in hit["hsps"]:
        _ = overlap((hsp["query_hsp_start"], hsp["query_hsp_end"]), boundary)
        if _ >= minimal_overlap * (boundary[1] + 1 - boundary[0]):
            hsp_dict_list.append(hsp)
            if hsp["hsp_evalue"] < best_hsp[0]:
                best_hsp = (hsp["hsp_evalue"], hsp["hsp_bits"])

            q_intervals.append((hsp["query_hsp_start"], hsp["query_hsp_end"]))
            t_intervals.append((hsp["target_hsp_start"], hsp["target_hsp_end"]))

            query_pos = hsp["query_hsp_start"] - 1

            for amino in hsp["match"]:
                if amino in __valid_matches or amino == "+":
                    query_pos += 1
                    positives.add(query_pos)
                    if amino != "+":
                        identical_positions.add(query_pos)
                elif amino == "_":  # Gap in the target sequence
                    query_pos += 1

    if len(hsp_dict_list) == 0:
        return None

    q_merged_intervals, q_aligned = merge(q_intervals)
    hit_dict["query_aligned_length"] = q_aligned
    hit_dict["query_start"] = q_merged_intervals[0][0]
    hit_dict["query_end"] = q_merged_intervals[-1][1]

    t_merged_intervals, t_aligned = merge(t_intervals)
    hit_dict["target_aligned_length"] = t_aligned
    hit_dict["target_start"] = t_merged_intervals[0][0]
    hit_dict["target_end"] = t_merged_intervals[-1][1]
    hit_dict["global_identity"] = len(identical_positions) * 100 / q_aligned
    hit_dict["global_positives"] = len(positives) * 100 / q_aligned
    hit_dict["hsps"] = hsp_dict_list
    hit_dict["bits"] = max(x["hsp_bits"] for x in hit_dict["hsps"])
    hit_dict["evalue"] = min(x["hsp_evalue"] for x in hit_dict["hsps"])

    return hit_dict


def __relocate_orfs(transcript, bed12_objects, tstart, tend):
    """
    Function to recalculate the coordinates of BED12 objects based on
    the new transcriptomic start/end
    :param bed12_objects: list of the BED12 ORFs to relocate
    :param tstart: New transcriptomic start
    :param tend: New transcriptomic end
    :return:
    """
    new_bed12s = []
    tstart, tend = sorted([tstart, tend])

    for obj in bed12_objects:
        # import copy
        # obj = copy.deepcopy(obj)
        new = BED12(table=transcript.codon_table)
        new.transcriptomic = True
        # Phase is necessary for truncated models
        for attr in ["chrom", "start", "end", "strand", "thick_start", "thick_end",
                     "block_starts", "block_sizes", "phase"]:
            setattr(new, attr, getattr(obj, attr))
        # obj.transcriptomic = True
        if new.strand == "-":
            thick_start = obj.end - obj.thick_end + 1
            thick_end = obj.end - obj.thick_start + 1
            old_start, old_end = tstart, tend
            local_tstart = obj.end - old_end + 1
            local_tend = obj.end - old_start + 1
            assert (old_end - old_start) == (local_tend - local_tstart), (
                (old_start, old_end), (local_tstart, local_tend))
            assert (thick_end - thick_start) == (obj.thick_end - obj.thick_start)
            new.strand = "+"
            new.start = 1
            new.end = local_tend - local_tstart + 1
            new.fasta_length = new.end
            assert new.end > 0 and new.end > new.start
            new.thick_start = max(thick_start - local_tstart + 1, 1)
            assert new.thick_start > 0
            new.thick_end = thick_end - local_tstart + 1
            assert new.thick_end > new.thick_start > 0
            new.block_sizes = [new.end]
            new.block_starts = [new.block_starts]
            transcript.logger.debug("Inverting negative ORF in %s",
                                    transcript.id)
        else:
            new.start = 1
            new.end = min(obj.end, tend) - tstart + 1
            new.fasta_length = new.end
            new.thick_start = min(new.thick_start, tend) - tstart + 1
            new.thick_end = min(obj.thick_end, tend) - tstart + 1
            new.block_sizes = [new.end]
            new.block_starts = new.block_starts[:]

        assert new.thick_start > 0, new.thick_start
        assert new.thick_end > 0, new.thick_end
        assert new.invalid is False, (len(new), new.cds_len, new.fasta_length,
                                      new.invalid_reason,
                                      str(new))

        new_bed12s.append(new)
    assert len(new_bed12s) > 0
    return new_bed12s


def split_by_cds(transcript):
    """This method is used for transcripts that have multiple ORFs.
    It will split them according to the CDS information into multiple transcripts.
    UTR information will be retained only if no ORF is down/upstream.

    :param transcript: the transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript
    """

    transcript.finalize()

    # List of the transcript that will be retained

    if transcript.number_internal_orfs < 2:
        new_transcripts = [transcript]  # If we only have one ORF this is easy
    elif (transcript.json_conf and
              transcript.json_conf.get("pick", None) and
              transcript.json_conf["pick"].get("chimera_split", None) and
              transcript.json_conf["pick"]["chimera_split"].get("skip", None) and
              transcript.source in transcript.json_conf["pick"]["chimera_split"]["skip"]):
        # Disable splitting for transcripts with certain tags
        transcript.logger.warning("%s (label %s) to be skipped for splitting",
                                  transcript.id,
                                  transcript.id.split("_")[0])
        new_transcripts = [transcript]
    else:
        cds_boundaries = SortedDict()
        for orf in sorted(transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        # Check whether we have to split or not based on BLAST data
        if transcript.json_conf is not None:
            if transcript.json_conf["pick"]["chimera_split"]["blast_check"] is True:
                cds_boundaries = check_split_by_blast(transcript, cds_boundaries)

        if len(cds_boundaries) == 1:
            # Recheck how many boundaries we have - after the BLAST check
            # we might have determined that the transcript has not to be split
            new_transcripts = [transcript]
        else:
            try:
                new_transcripts = __create_splitted_transcripts(transcript, cds_boundaries)
            except InvalidTranscript as err:
                exc = InvalidTranscript(err)
                transcript.logger.error("Error in splitting %s by ORF",
                                        transcript.id)
                transcript.logger.exception(exc)
                transcript.logger.error("ORFs: %s",
                                        "\n".join([str(_) for _ in transcript.internal_orfs]))
                transcript.logger.error("BED12: %s",
                                        "\n".join([str(_) for _ in transcript.loaded_bed12]))
                transcript.logger.error("Stripping %s of its CDS.",
                                        transcript.id)
                transcript.strip_cds()
                new_transcripts = [transcript]

    assert len(new_transcripts) > 0, str(transcript)
    __original = set()
    for internal in transcript.internal_orfs:
        __original.add(tuple([_[1] for _ in internal if _[0] == "CDS"]))

    for new_transc in new_transcripts:
        for internal in new_transc.internal_orfs:
            internal = tuple([_[1] for _ in internal if _[0] == "CDS"])
            assert internal in __original, (transcript.id, __original, internal)

        new_transc.verified_introns = set.intersection(set(new_transc.introns),
                                                       transcript.verified_introns)
        new_transc.attributes["split"] = True
        yield new_transc

    return
