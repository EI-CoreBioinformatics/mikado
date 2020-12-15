from ...exceptions import InvalidTranscript
from ...transcripts import Transcript
import functools
from .transmission import transmit_transcript


def parse_prediction_gtf(args, queue, queue_logger):
    """Method to parse GTF files."""
    invalids = set()
    done = 0
    lastdone = 1
    transcript = None
    __found_with_orf = set()
    constructor = functools.partial(Transcript, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
    rows = []

    for row in args.prediction:
        if row.header is True:
            continue
        if row.is_transcript is True:
            if transcript is not None:
                rows, done, lastdone, __found_with_orf = transmit_transcript(
                    transcript=transcript, done=done, lastdone=lastdone,
                    rows=rows, __found_with_orf=__found_with_orf,
                    queue=queue, queue_logger=queue_logger)
            try:
                transcript = constructor(row)
            except (InvalidTranscript, AssertionError, TypeError, ValueError):
                queue_logger.warning("Row %s is invalid, skipping.", row)
                transcript = None
                invalids.add(row.id)
                continue
        elif row.is_exon is True:
            # Case 1: we are talking about cDNA_match and GFF
            if any(_ in invalids for _ in row.parent):
                # Skip children of invalid things
                continue
            elif transcript is None or (transcript is not None and transcript.id != row.transcript):
                rows, done, lastdone, __found_with_orf = transmit_transcript(
                    transcript=transcript, done=done, lastdone=lastdone,
                    rows=rows, __found_with_orf=__found_with_orf,
                    queue=queue, queue_logger=queue_logger)
                queue_logger.debug("New transcript: %s", row.transcript)
                transcript = constructor(row)
                transcript.add_exon(row)
            elif transcript.id == row.transcript:
                transcript.add_exon(row)
            else:
                raise TypeError("Unmatched exon: {}".format(row))
        else:
            queue_logger.debug("Skipped row: {}".format(row))

    rows, done, lastdone, __found_with_orf = transmit_transcript(
        transcript=transcript, done=done, lastdone=lastdone,
        rows=rows, __found_with_orf=__found_with_orf,
        queue=queue, queue_logger=queue_logger, send_all=True)
    return done, lastdone
