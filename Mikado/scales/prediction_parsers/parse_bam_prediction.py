import collections
from ...parsers.bam_parser import BamParser
from ...exceptions import InvalidTranscript
from ...transcripts import Transcript
import functools


def parse_prediction_bam(args, queue_logger):
    constructor = functools.partial(Transcript, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)

    transcript = None
    name_counter = collections.Counter()  # This is needed for BAMs
    invalids = set()
    if args.prediction.__annot_type__ == BamParser.__annot_type__:
        for row in args.prediction:
            if row.is_unmapped is True:
                continue
            yield transcript
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
    yield transcript
