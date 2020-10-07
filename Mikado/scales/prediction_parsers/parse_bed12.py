from ...exceptions import InvalidTranscript


def parse_prediction_bed12(args, queue_logger, transmit_wrapper, constructor):
    """"""

    transcript = None
    done = 0
    lastdone = 1
    invalids = set()
    __found_with_orf = set()
    coord_list = None
    for row in args.prediction:
        if row.header is True:
            continue
        done, lastdone, coord_list, __found_with_orf = transmit_wrapper(transcript=transcript,
                                                            done=done,
                                                            lastdone=lastdone,
                                                                        coord_list=coord_list,
                                                            __found_with_orf=__found_with_orf)
        try:
            transcript = constructor(row)
        except (InvalidTranscript, AssertionError, TypeError, ValueError):
            queue_logger.warning("Row %s is invalid, skipping.", row)
            transcript = None
            invalids.add(row.id)
            continue
        transcript.parent = transcript.id

    done, lastdone, coord_list, __found_with_orf = transmit_wrapper(transcript=transcript,
                                                        done=done,
                                                        lastdone=lastdone,
                                                                    coord_list=coord_list,
                                                        __found_with_orf=__found_with_orf)
    return done, lastdone, coord_list
