from . import Query, Target, Hsp, Hit, prepare_hit, InvalidHit
import sqlalchemy.exc


def load_into_db(self, hits, hsps, force=False, lock=None):
    """
    :param hits:
    :param hsps:

    :param force: boolean flag. If set, data will be loaded no matter what.
    To be used at the end of the serialisation to load the final batch of data.
    :type force: bool

    :return:
    """

    self.logger.debug("Checking whether to load %d hits and %d hsps",
                      len(hits), len(hsps))

    tot_objects = len(hits) + len(hsps)
    if len(hits) == 0:
        self.logger.debug("No hits to serialise. Exiting")
        return hits, hsps

    if tot_objects >= self.maxobjects or force:
        # Bulk load
        self.logger.debug("Loading %d BLAST objects into database", tot_objects)

        if lock is not None:
            lock.acquire()
        try:
            # pylint: disable=no-member
            self.session.begin(subtransactions=True)
            self.engine.execute(Hit.__table__.insert(), hits)
            self.engine.execute(Hsp.__table__.insert(), hsps)
            # pylint: enable=no-member
            self.session.commit()
        except sqlalchemy.exc.IntegrityError as err:
            self.logger.critical("Failed to serialise BLAST!")
            self.logger.exception(err)
            raise err
        finally:
            if lock is not None:
                lock.release()
        self.logger.debug("Loaded %d BLAST objects into database", tot_objects)
        hits, hsps = [], []
    return hits, hsps
