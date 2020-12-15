from ...exceptions import InvalidTranscript
from ...transcripts import Transcript
import functools
from .transmission import transmit_transcript, send_transcripts


def parse_prediction_bed12(args, queue, queue_logger):
    """"""

    constructor = functools.partial(Transcript, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
    transcript = None
    done = 0
    lastdone = 1
    invalids = set()
    __found_with_orf = set()
    rows = []
    for row in args.prediction:
        if row.header is True:
            continue
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
        transcript.parent = transcript.id

    rows, done, lastdone, __found_with_orf = transmit_transcript(
        transcript=transcript, done=done, lastdone=lastdone,
        rows=rows, __found_with_orf=__found_with_orf,
        queue=queue, queue_logger=queue_logger, send_all=True)

    return done, lastdone
