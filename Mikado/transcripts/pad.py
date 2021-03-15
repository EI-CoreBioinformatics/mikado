from .transcript import Transcript
import pysam
from ..exceptions import InvalidTranscript
from .transcriptchecker import TranscriptChecker


def pad_transcript(transcript: Transcript,
                   backup: Transcript,
                   start_transcript: [Transcript, bool],
                   end_transcript: [Transcript, bool],
                   fai: pysam.libcfaidx.FastaFile,
                   logger):

    """This method will enlarge the coordinates and exon structure of a transcript, given:
    :param transcript: the transcript to modify.
    :type transcript: Transcript
    :param backup: a copy of the transcript to be modified.
    :type backup: Transcript
    :param start_transcript: the template transcript for the 5' end.
    :param end_transcript: the template transcript for the 3' end.
    :param fai: the indexed genomic sequence.
    :param logger: the logger to be used in the function.
    """

    # If there is nothing to do, just get out
    assert transcript == backup
    transcript.finalize()
    if start_transcript not in (False, None):
        start_transcript.finalize()
    if end_transcript not in (False, None):
        end_transcript.finalize()

    if start_transcript in (False, None) and end_transcript in (False, None):
        logger.debug("%s does not need to be expanded, exiting", transcript.id)
        return transcript

    if transcript.strand == "-":
        start_transcript, end_transcript = end_transcript, start_transcript

    # Make a backup copy of the transcript
    # First get the ORFs
    # Remove the CDS and unfinalize
    logger.debug("Starting expansion of %s", transcript.id)
    strand = transcript.strand
    transcript.strip_cds()
    transcript.unfinalize()
    assert strand == transcript.strand

    upstream, up_exons, new_first_exon, up_remove = _enlarge_start(transcript, backup, start_transcript)
    downstream, up_exons, down_exons, down_remove = _enlarge_end(transcript,
                                                                 backup, end_transcript, up_exons, new_first_exon)

    first_exon, last_exon = transcript.exons[0], transcript.exons[-1]

    assert upstream >= 0 and downstream >= 0

    if up_remove is True:
        # Remove the first exon
        transcript.remove_exon(first_exon)
    if down_remove is True:
        if not (up_remove is True and first_exon == last_exon):
            transcript.remove_exon(last_exon)

    new_exons = up_exons + down_exons
    if not new_exons:
        logger.debug("%s does not need to be expanded, exiting", transcript.id)
        return backup

    transcript.add_exons(new_exons)
    transcript.start, transcript.end = None, None
    transcript.finalize()

    if transcript.strand == "-":
        downstream, upstream = upstream, downstream

    if backup.is_coding:
        seq = check_expanded(transcript, backup, start_transcript, end_transcript,
                             fai, upstream, downstream, logger)
        transcript = enlarge_orfs(transcript, backup, seq, upstream, downstream, logger)
        transcript.finalize()

    logger.debug("%s: start (before %s, now %s, %s), end (before %s, now %s, %s)",
                 transcript.id,
                 backup.start, transcript.start, transcript.start < backup.start,
                 backup.end, transcript.end, transcript.end > backup.end)
    if transcript.start < backup.start or transcript.end > backup.end:
        transcript.attributes["padded"] = True

    # Now check that we have a valid expansion
    if backup.is_coding and not transcript.is_coding:
        # Something has gone wrong. Just return the original transcript.
        assert new_exons
        logger.info("Padding %s would lead to an invalid CDS (up exons: %s). Aborting.",
                    transcript.id, up_exons)
        return backup
    elif backup.is_coding:
        abort = False
        if backup.strand == "-" and backup.combined_cds_end < transcript.combined_cds_end:
            abort = True
        elif backup.strand != "-" and backup.combined_cds_end > transcript.combined_cds_end:
            abort = True
        if abort is True:
            msg = "Padding {} (strand: {}) would lead to an in-frame stop codon ({} to {}, \
vs original {} to {}. Aborting.".format(
                transcript.id, backup.strand, transcript.combined_cds_start, transcript.combined_cds_end,
                backup.combined_cds_start, backup.combined_cds_end)
            logger.info(msg)
            return backup

    return transcript


def _enlarge_start(transcript: Transcript,
                   backup: Transcript,
                   start_transcript: Transcript) -> (int, list, [None, tuple], bool):

    """This method will enlarge the transcript at the 5' end, using another transcript as the template.
    :param transcript: the original transcript to modify.
    :param backup: a copy of the transcript. As we are modifying the original one, we do need a hard copy.
    :param start_transcript: the template transcript.

    The function returns the following:
    :returns: the upstream modification, the list of upstream exons to add, the new first exon (if any),
              a boolean flag indicating whether the first exon of the transcript should be removed.
    """

    upstream = 0
    up_exons = []
    new_first_exon = None
    to_remove = False
    if start_transcript:
        transcript.start = start_transcript.start
        upstream_exons = sorted(
            [_ for _ in start_transcript.find_upstream(transcript.exons[0][0], transcript.exons[0][1])
             if _.value == "exon"])
        intersecting_upstream = sorted(start_transcript.search(
            transcript.exons[0][0], transcript.exons[0][1]))

        if not intersecting_upstream:
            raise KeyError("No exon or intron found to be intersecting with %s vs %s, this is a mistake",
                           transcript.id, start_transcript.id)

        if intersecting_upstream[0].value == "exon":
            new_first_exon = (min(intersecting_upstream[0][0], backup.start),
                              transcript.exons[0][1])
            if new_first_exon != transcript.exons[0]:
                upstream += backup.start - new_first_exon[0]
                up_exons.append(new_first_exon)
                to_remove = True
            else:
                new_first_exon = None
            if intersecting_upstream[0] in upstream_exons:
                upstream_exons.remove(intersecting_upstream[0])
            upstream += sum(_[1] - _[0] + 1 for _ in upstream_exons)
            up_exons.extend([(_[0], _[1]) for _ in upstream_exons])
        elif intersecting_upstream[0].value == "intron":
            # Check whether the first exon of the model *ends* within an *intron* of the template
            # If that is the case, we have to keep the first exon in place and
            # just expand it until the end
            # Now we have to expand until the first exon in the upstream_exons
            if intersecting_upstream[0][1] == transcript.exons[0][0] - 1:
                assert upstream_exons
                to_remove = False
            elif upstream_exons:
                to_remove = True
                upstream_exon = upstream_exons[-1]
                new_first_exon = (upstream_exon[0], transcript.exons[0][1])
                upstream_exons.remove(upstream_exon)
                upstream += backup.start - new_first_exon[0]
                up_exons.append(new_first_exon)
            else:
                # Something fishy going on here. Let us double check everything.
                if start_transcript.exons[0][0] == transcript.start:
                    raise ValueError(
                        "Something has gone wrong. The template transcript should have returned upstream exons."
                    )
                elif start_transcript.exons[0][0] < transcript.start:
                    raise ValueError(
                        "Something has gone wrong. We should have found the correct exons."
                    )
                else:
                    pass

            upstream += sum(_[1] - _[0] + 1 for _ in upstream_exons)
            up_exons.extend([(_[0], _[1]) for _ in upstream_exons])

    return upstream, up_exons, new_first_exon, to_remove


def _enlarge_end(transcript: Transcript,
                 backup: Transcript,
                 end_transcript: Transcript,
                 up_exons: list,
                 new_first_exon: [None, tuple]) -> [int, list, list, bool]:

    """
    This method will enlarge the transcript at the 5' end, using another transcript as the template.
    :param transcript: the original transcript to modify.
    :param backup: a copy of the transcript. As we are modifying the original one, we do need a hard copy.
    :param end_transcript: the template transcript.
    :param up_exons: the list of exons added at the 5' end.
    :param new_first_exon: the new coordinates of what used to be the first exon of the transcript.
                           This is necessary because if the transcript is monoexonic, we might need to re-modify it.

    The function returns the following:
    :returns: the downstream modification, the (potentially modified) list of upstream exons to add,
              the list of downstream exons to add, a boolean flag indicating whether the last exon of the transcript
              should be removed.
    """

    downstream = 0
    down_exons = []
    to_remove = False

    if end_transcript:
        transcript.end = end_transcript.end
        downstream_exons = sorted(
            [_ for _ in end_transcript.find_downstream(transcript.exons[-1][0], transcript.exons[-1][1])
             if _.value == "exon"])
        intersecting_downstream = sorted(end_transcript.search(
            transcript.exons[-1][0], transcript.exons[-1][1]))
        if not intersecting_downstream:
            raise KeyError("No exon or intron found to be intersecting with %s vs %s, this is a mistake",
                           transcript.id, end_transcript.id)
        # We are taking the right-most intersecting element.
        if intersecting_downstream[-1].value == "exon":
            if transcript.monoexonic and new_first_exon is not None:
                new_exon = (new_first_exon[0], max(intersecting_downstream[-1][1], new_first_exon[1]))
                if new_exon != new_first_exon:
                    up_exons.remove(new_first_exon)
                    downstream += new_exon[1] - backup.end
                    down_exons.append(new_exon)
                    to_remove = True
            else:
                new_exon = (transcript.exons[-1][0],
                            max(intersecting_downstream[-1][1], transcript.exons[-1][1]))
                if new_exon != transcript.exons[-1]:
                    downstream += new_exon[1] - backup.end
                    down_exons.append(new_exon)
                    to_remove = True

            if intersecting_downstream[-1] in downstream_exons:
                downstream_exons.remove(intersecting_downstream[-1])
            downstream += sum(_[1] - _[0] + 1 for _ in downstream_exons)
            down_exons.extend([(_[0], _[1]) for _ in downstream_exons])
        elif intersecting_downstream[-1].value == "intron":
            # Now we have to expand until the first exon in the upstream_exons
            if intersecting_downstream[-1][0] == transcript.exons[-1][1] + 1:
                assert downstream_exons
                to_remove = False
            elif downstream_exons:
                downstream_exon = downstream_exons[0]
                assert downstream_exon[1] > backup.end
                assert downstream_exon[0] > backup.end
                if transcript.monoexonic and new_first_exon is not None:
                    new_exon = (new_first_exon[0], downstream_exon[1])
                    up_exons.remove(new_first_exon)
                    to_remove = True
                else:
                    new_exon = (transcript.exons[-1][0], downstream_exon[1])
                    to_remove = True
                downstream_exons.remove(downstream_exon)
                downstream += new_exon[1] - backup.end
                down_exons.append(new_exon)
            else:
                # Something fishy going on here. Let us double check everything.
                if end_transcript.exons[-1][1] == transcript.end:
                    raise ValueError(
                        "Something has gone wrong. The template transcript should have returned upstream exons."
                    )
                elif end_transcript.exons[-1][1] > transcript.end:
                    raise ValueError(
                        "Something has gone wrong. We should have found the correct exons."
                    )
            downstream += sum(_[1] - _[0] + 1 for _ in downstream_exons)
            down_exons.extend([(_[0], _[1]) for _ in downstream_exons])

    return downstream, up_exons, down_exons, to_remove


def check_expanded(transcript, backup, start_transcript, end_transcript, fai, upstream, downstream, logger) -> str:

    """
    This function checks that the expanded transcript is valid, and it also calculates and returns its cDNA sequence.
    :param transcript: the modified transcript.
    :param backup: The original transcript, before expansion.
    :param start_transcript: the transcript used as template at the 5' end.
    :param end_transcript: the transcript used as template at the 3' end.
    :param fai: The pysam.libcfaidx.FastaFile object indexing the genome.
    :param upstream: the amount of transcriptomic base-pairs added to the transcript at its 5' end.
    :param downstream: the amount of transcriptomic base-pairs added to the transcript at its 3' end.
    :param logger: the logger to use.
    :returns: the cDNA of the modified transcript, as a standard Python string.
    """

    assert transcript.exons != backup.exons
    assert transcript.end <= fai.get_reference_length(transcript.chrom), (
        transcript.end, fai.get_reference_length(transcript.chrom))
    genome_seq = fai.fetch(transcript.chrom, transcript.start - 1, transcript.end)

    if not (transcript.exons[-1][1] - transcript.start + 1 == len(genome_seq)):
        error = "{} should have a sequence of length {} ({} start, {} end), but one of length {} has been given"
        error = error.format(transcript.id, transcript.exons[-1][1] - transcript.start + 1,
                             transcript.start, transcript.end, len(genome_seq))
        logger.error(error)
        raise InvalidTranscript(error)
    seq = TranscriptChecker(transcript, genome_seq, is_reference=True).cdna
    assert len(seq) == transcript.cdna_length, (len(seq), transcript.cdna_length, transcript.exons)
    if not len(seq) == backup.cdna_length + upstream + downstream:
        error = [len(seq), backup.cdna_length + upstream + downstream,
                 backup.cdna_length, upstream, downstream,
                 (transcript.start, transcript.end), (backup.id, backup.start, backup.end),
                 (None if not start_transcript else (start_transcript.id, (start_transcript.start,
                                                                           start_transcript.end))),
                 (None if not end_transcript else (end_transcript.id, (end_transcript.start,
                                                                       end_transcript.end))),
                 (backup.id, backup.exons),
                 None if not start_transcript else (start_transcript.id, start_transcript.exons),
                 None if not end_transcript else (end_transcript.id, end_transcript.exons),
                 (transcript.id + "_expanded", transcript.exons),
                 set.difference(set(transcript.exons), set(backup.exons)),
                 set.difference(set(backup.exons), set(transcript.exons))
                 ]
        error = "\n".join([str(_) for _ in error])
        raise AssertionError(error)
    return seq


def enlarge_orfs(transcript: Transcript,
                 backup: Transcript,
                 seq: str,
                 upstream: int,
                 downstream: int,
                 logger) -> Transcript:

    """
    This method will take an expanded transcript and recalculate its ORF(s). As a consequence of the expansion,
    truncated transcripts might become whole.
    :param transcript: the expanded transcript.
    :param backup: the original transcript. Used to extract the original ORF(s).
    :param seq: the new cDNA sequence of the expanded transcript.
    :param upstream: the amount of expansion that happened at the 5'.
    :param downstream: the amount of expansion that happened at the 3'.
    :param logger: the logger.
    :returns: the modified transcript with the ORF(s) recalculated.
    """

    if backup.combined_cds_length > 0:
        try:
            internal_orfs = list(backup.get_internal_orf_beds())
        except (ValueError, TypeError, AssertionError):
            logger.error("Something went wrong with the CDS extraction for %s. Stripping it.",
                         backup.id)
            internal_orfs = []
    else:
        internal_orfs = []

    if not internal_orfs:
        return transcript

    new_orfs = []
    for orf in internal_orfs:
        logger.debug("Old ORF: %s", str(orf))
        try:
            logger.debug("Sequence for %s: %s[..]%s (upstream %s, downstream %s)",
                         transcript.id, seq[:10], seq[-10:], upstream, downstream)
            orf.expand(seq, upstream, downstream, expand_orf=True, logger=logger)
        except AssertionError as err:
            logger.error(err)
            logger.error("%s, %s, %s, %s",
                         upstream,
                         downstream,
                         transcript.exons,
                         transcript.cdna_length)
            raise AssertionError(err)
        logger.debug("New ORF: %s", str(orf))
        if orf.coding is False:
            raise ValueError(orf)
        elif orf.invalid:
            raise InvalidTranscript(orf.invalid_reason)

        new_orfs.append(orf)

    transcript.load_orfs(new_orfs)
    transcript.finalize()
    if backup.is_coding and not transcript.is_coding:
        raise InvalidTranscript(new_orfs)
    return transcript
