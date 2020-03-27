from . import Hsp, Hit
import sqlalchemy.exc


def load_into_db(self, hits, hsps, force=False, raw=False):
    """
    :param hits:
    :param hsps:

    :param force: boolean flag. If set, data will be loaded no matter what.
    To be used at the end of the serialisation to load the final batch of data.
    :type force: bool

    :return:
    """

    # self.logger.debug("Checking whether to load %d hits and %d hsps", len(hits), len(hsps))

    tot_objects = len(hits) + len(hsps)
    if len(hits) == 0:
        self.logger.debug("No hits to serialise. Exiting")
        return hits, hsps

    if tot_objects >= self.maxobjects or force:
        # Bulk load
        self.logger.debug("Loading %d BLAST objects into database", tot_objects)

        if not hasattr(self, "hit_i_string"):
            self.hit_i_string = str(Hit.__table__.insert(bind=self.engine).compile())
            self.hsp_i_string = str(Hsp.__table__.insert(bind=self.engine).compile())

        try:
            # pylint: disable=no-member
            # self.session.begin(subtransactions=True)
            if hasattr(self, "lock") and self.lock is not None:
                self.lock.acquire()
            if raw is True:
                self.engine.execute(self.hit_i_string, hits)
                self.engine.execute(self.hsp_i_string, hsps)
            else:
                self.engine.execute(Hit.__table__.insert(), hits)
                self.engine.execute(Hsp.__table__.insert(), hsps)
            if hasattr(self, "lock") and self.lock is not None:
                self.lock.release()
            # pylint: enable=no-member
        except sqlalchemy.exc.IntegrityError as err:
            self.logger.critical("Failed to serialise BLAST!")
            self.logger.exception(err)
            raise err
        self.logger.debug("Loaded %d BLAST objects into database", tot_objects)
        hits, hsps = [], []

    if force is True: # Final push
        if hasattr(self, "lock") and self.lock is not None:
            self.lock.acquire()
        self.session.commit()
        if hasattr(self, "lock") and self.lock is not None:
            self.lock.release()

    return hits, hsps
