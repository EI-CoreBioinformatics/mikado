import os
from ...transcripts import Transcript
import msgpack
from ..assignment.assigner import msgpack_default, Assigner
import re


def transmit_transcript(transcript: Transcript, connection):
    transcript.finalize()
    d = transcript.as_dict(remove_attributes=False)
    start = connection.tell()
    to_write = msgpack.dumps(d)
    connection.write(to_write)
    connection.flush()
    assert len(to_write) == 0 or os.stat(connection.name).st_size > 0
    end = connection.tell()
    return tuple([start, end])
    # connection.execute("INSERT INTO dump VALUES (?, ?)", (index, json.dumps(d)))


def get_best_result(transcript, assigner_instance: Assigner):
    transcript.finalize()
    assigner_instance.get_best(transcript)


orf_pattern = re.compile(r"\.orf[0-9]+$", re.IGNORECASE)


def _transmit_transcript(transcript, done, lastdone,
                         assigner_instance, transmitter, queue_logger,
                         queue, dump_db, __found_with_orf,
                         coord_list=None):

    if coord_list is None:
        coord_list = []

    if transcript is not None:
        if orf_pattern.search(transcript.id):
            __name = orf_pattern.sub("", transcript.id)
            if __name not in __found_with_orf:
                __found_with_orf.add(__name)
                done += 1
                if done and done % 10000 == 0:
                    queue_logger.info("Parsed %s transcripts", done)
                    if assigner_instance is None:
                        dump_db.flush()
                        [queue.put(coords) for coords in coord_list]
                        coord_list = []
                    lastdone = done
                coords = transmitter(transcript)
                coord_list.append(coords)
            else:
                pass
        else:
            done += 1
            if done and done % 10000 == 0:
                queue_logger.info("Parsed %s transcripts", done)
                if assigner_instance is None:
                    dump_db.flush()
                    [queue.put(coords) for coords in coord_list]
                    coord_list = []
                lastdone = done
            try:
                coords = transmitter(transcript)
                coord_list.append(coords)
            except AssertionError:
                raise AssertionError((transcript.id, ))

    return done, lastdone, coord_list, __found_with_orf
