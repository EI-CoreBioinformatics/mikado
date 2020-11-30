from ...transcripts import Transcript
import msgpack
import zlib
from ..assignment.assigner import Assigner
import re


def transmit_transcript(transcript: Transcript, queue):
    transcript.finalize()
    d = transcript.as_dict(remove_attributes=False)
    to_write = msgpack.dumps(d)
    queue.put(to_write)


def get_best_result(transcript, assigner_instance: Assigner):
    transcript.finalize()
    assigner_instance.get_best(transcript)


orf_pattern = re.compile(r"\.orf[0-9]+$", re.IGNORECASE)


def _transmit_transcript(transcript, done, lastdone,
                         transmitter, queue_logger,
                         __found_with_orf):

    if transcript is not None:
        if orf_pattern.search(transcript.id):
            __name = orf_pattern.sub("", transcript.id)
            if __name not in __found_with_orf:
                __found_with_orf.add(__name)
                done += 1
                if done and done % 10000 == 0:
                    queue_logger.info("Parsed %s transcripts", done)
                    lastdone = done
                transmitter(transcript)
            else:
                pass
        else:
            done += 1
            if done and done % 10000 == 0:
                queue_logger.info("Parsed %s transcripts", done)
                lastdone = done
            try:
                transmitter(transcript)
            except AssertionError:
                raise AssertionError((transcript.id, ))

    return done, lastdone, __found_with_orf
