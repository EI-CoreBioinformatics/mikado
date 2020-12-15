from ...exceptions import InvalidTranscript
from ...transcripts import Transcript
import functools


def parse_prediction_bed12(args, queue_logger):
    """"""

    constructor = functools.partial(Transcript, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
    transcript = None
    invalids = set()
    for row in args.prediction:
        if row.header is True:
            continue
        yield transcript
        try:
            transcript = constructor(row)
        except (InvalidTranscript, AssertionError, TypeError, ValueError):
            queue_logger.warning("Row %s is invalid, skipping.", row)
            transcript = None
            invalids.add(row.id)
            continue
        transcript.parent = transcript.id

    yield transcript
