import collections
from ...parsers.bam_parser import BamParser
from ...exceptions import InvalidTranscript
from ...transcripts import Transcript
import functools
from .transmission import transmit_transcript


def parse_prediction_bam(args, queue, queue_logger):
    constructor = functools.partial(Transcript, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)

    transcript = None
    done = 0
    lastdone = 1
    __found_with_orf = set()
    name_counter = collections.Counter()  # This is needed for BAMs
    invalids = set()
    rows = []
    if args.prediction.__annot_type__ == BamParser.__annot_type__:
        for row in args.prediction:
            if row.is_unmapped is True:
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
            if name_counter.get(row.query_name):
                name = "{}_{}".format(row.query_name, name_counter.get(row.query_name))
            else:
                name = row.query_name
            transcript.id = transcript.name = transcript.alias = name
            transcript.parent = transcript.attributes["gene_id"] = "{0}.gene".format(name)
    rows, done, lastdone, __found_with_orf = transmit_transcript(
        transcript=transcript, done=done, lastdone=lastdone,
        rows=rows, __found_with_orf=__found_with_orf,
        queue=queue, queue_logger=queue_logger, send_all=True)

    return done, lastdone
