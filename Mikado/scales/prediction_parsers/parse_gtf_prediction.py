from ...exceptions import InvalidTranscript


def parse_prediction_gtf(args, queue_logger, transmit_wrapper, constructor):
    """Method to parse GTF files."""
    invalids = set()
    done = 0
    lastdone = 1
    transcript = None
    __found_with_orf = set()

    for row in args.prediction:
        if row.header is True:
            continue
        if row.is_transcript is True:
            if transcript is not None:
                done, lastdone, __found_with_orf = transmit_wrapper(
                    transcript=transcript,
                    __found_with_orf=__found_with_orf,
                    done=done, lastdone=lastdone)
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
                done, lastdone, __found_with_orf = transmit_wrapper(
                    transcript=transcript,
                    __found_with_orf=__found_with_orf,
                    done=done,
                    lastdone=lastdone)
                queue_logger.debug("New transcript: %s", row.transcript)
                transcript = constructor(row)
                transcript.add_exon(row)
            elif transcript.id == row.transcript:
                transcript.add_exon(row)
            else:
                raise TypeError("Unmatched exon: {}".format(row))
        else:
            queue_logger.debug("Skipped row: {}".format(row))
    done, lastdone, __found_with_orf = transmit_wrapper(
        transcript=transcript,
        done=done,
        lastdone=lastdone,
        __found_with_orf=__found_with_orf)
    return done, lastdone
