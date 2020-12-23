from ...transcripts import Transcript
import msgpack
import zlib
from ..assignment.assigner import Assigner
import re
from queue import Queue
from logging import Logger
from typing import Union


def send_transcripts(transcript: Union[None,Transcript], rows: list, queue: Queue, maxsize=1000):
    if transcript is not None:
        transcript.finalize()
        d = transcript.as_dict(remove_attributes=False)
        rows.append(d)
        send_all = False
    else:
        send_all = True
    if len(rows) >= maxsize or send_all is True:
        to_write = msgpack.dumps(rows)
        queue.put(to_write)
        rows = []
    return rows


def get_best_result(transcript, assigner_instance: Assigner):
    transcript.finalize()
    assigner_instance.get_best(transcript)


orf_pattern = re.compile(r"\.orf[0-9]+$", re.IGNORECASE)


def transmit_transcript(transcript: Union[None,Transcript], done: int, lastdone: int,
                         rows: list, queue: Queue,
                         queue_logger: Logger,
                         __found_with_orf: set, send_all=False):

    if transcript is not None:
        if orf_pattern.search(transcript.id):
            __name = orf_pattern.sub("", transcript.id)
            if __name not in __found_with_orf:
                __found_with_orf.add(__name)
                done += 1
                if done and done % 10000 == 0:
                    queue_logger.info("Parsed %s transcripts", done)
                    lastdone = done
                rows = send_transcripts(transcript, rows, queue)
            else:
                pass
        else:
            done += 1
            if done and done % 10000 == 0:
                queue_logger.info("Parsed %s transcripts", done)
                lastdone = done
            try:
                rows = send_transcripts(transcript, rows, queue)
            except AssertionError:
                raise AssertionError((transcript.id, ))

    if send_all is True:
        send_transcripts(None, rows, queue)

    return rows, done, lastdone, __found_with_orf
