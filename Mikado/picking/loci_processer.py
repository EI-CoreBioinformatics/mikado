from multiprocessing import Process
from multiprocessing.managers import AutoProxy
import logging
from itertools import product
import logging.handlers as logging_handlers
import functools
from ..utilities import dbutils
from ..scales.assigner import Assigner
from ..loci.superlocus import Superlocus
import os
import collections
import sys
import pickle
import sqlite3
from ..transcripts import Transcript
from ..exceptions import InvalidTranscript
from ..parsers.GTF import GtfLine
import msgpack
try:
    import rapidjson as json
except ImportError:
    import json
from copy import deepcopy

__author__ = 'Luca Venturini'


def default_for_serialisation(obj):
    if isinstance(obj, set):
        return tuple(obj)


def print_locus(stranded_locus,
                gene_counter,
                handles,
                logger=None,
                json_conf=None):
    """
    Method that handles a single superlocus for printing.
    It also detects and flags/discard fragmentary loci.
    :param stranded_locus: the stranded locus to analyse
    :return:
    """

    locus_metrics, locus_scores, locus_out = handles[0]
    sub_metrics, sub_scores, sub_out = handles[1]
    mono_metrics, mono_scores, mono_out = handles[2]

    if json_conf is None:
        from ..configuration.configurator import to_json
        json_conf = to_json(None)

    stranded_locus.logger = logger
    if sub_out is not None:  # Skip this section if no sub_out is defined
        sub_lines = stranded_locus.__str__(
            level="subloci",
            print_cds=not json_conf["pick"]["run_options"]["exclude_cds"])
        if sub_lines != '':
            print(sub_lines, file=sub_out)
        sub_metrics_rows = [_ for _ in stranded_locus.print_subloci_metrics()
                            if _ != {} and "tid" in _]
        sub_scores_rows = [_ for _ in stranded_locus.print_subloci_scores()
                           if _ != {} and "tid" in _]
        for row in sub_metrics_rows:
            [row.__delitem__(key) for key in list(row.keys()) if key not in sub_metrics.fieldnames]
            sub_metrics.writerow(row)
        for row in sub_scores_rows:
            assert "alias" in sub_scores.fieldnames
            [row.__delitem__(key) for key in list(row.keys()) if key not in sub_scores.fieldnames]
            sub_scores.writerow(row)
    if mono_out is not None:
        mono_lines = stranded_locus.__str__(
            level="monosubloci",
            print_cds=not json_conf["pick"]["run_options"]["exclude_cds"])
        if mono_lines != '':
            print(mono_lines, file=mono_out)
        mono_metrics_rows = [_ for _ in stranded_locus.print_monoholder_metrics()
                             if _ != {} and "tid" in _]
        mono_scores_rows = [_ for _ in stranded_locus.print_monoholder_scores()
                            if _ != {} and "tid" in _]
        for row in mono_metrics_rows:
            [row.__delitem__(key) for key in list(row.keys()) if key not in mono_metrics.fieldnames]
            mono_metrics.writerow(row)
        for row in mono_scores_rows:
            [row.__delitem__(key) for key in list(row.keys()) if key not in mono_scores.fieldnames]
            mono_scores.writerow(row)

    for locus in stranded_locus.loci:
        gene_counter += 1
        fragment_test = (
            json_conf["pick"]["fragments"]["remove"]
            is True and stranded_locus.loci[locus].is_fragment is True)

        if fragment_test is True:
            continue
        gene_counter += 1
        new_id = "{0}.{1}G{2}".format(
            json_conf["pick"]["output_format"]["id_prefix"],
            stranded_locus.chrom, gene_counter)
        stranded_locus.loci[locus].logger = logger
        stranded_locus.loci[locus].id = new_id

    locus_lines = stranded_locus.__str__(
        print_cds=not json_conf["pick"]["run_options"]["exclude_cds"],
        level="loci")

    locus_metrics_rows = [x for x in stranded_locus.print_loci_metrics()]
    locus_scores_rows = [x for x in stranded_locus.print_loci_scores()]

    if locus_lines:
        assert len(locus_metrics_rows) > 0
        print(locus_lines, file=locus_out)

    # assert len(locus_metrics_rows) == len(locus_scores_rows)

    for row in locus_metrics_rows:
        [row.__delitem__(key) for key in list(row.keys()) if key not in locus_metrics.fieldnames]
        locus_metrics.writerow(row)
    for row in locus_scores_rows:
        [row.__delitem__(key) for key in list(row.keys()) if key not in locus_scores.fieldnames]
        locus_scores.writerow(row)
    # Necessary to flush out all the files
    [_.flush() for _ in handles if hasattr(_, "close")]
    return gene_counter


def serialise_locus(stranded_loci: [Superlocus],
                    conn: sqlite3.Connection,
                    counter):

    loci = msgpack.dumps([json.dumps(stranded_locus.as_dict(), default=default_for_serialisation,
                                     number_mode=json.NM_NATIVE)
                          for stranded_locus in stranded_loci])
    conn.execute("INSERT INTO loci VALUES (?, ?)", (counter, loci))
    conn.commit()
    return


def merge_loci(num_temp, out_handles,
               json_conf=None,
               logger=None,
               tempdir="mikado_pick_tmp"):

    """ Function to merge the temporary loci files into single output files,
      renaming the genes according to the preferred style.
    :param num_temp: number of temporary files.
    :param out_handles: The names of the output loci files.
    :param json_conf: configuration store
    :param logger: logger to use
    :param tempdir: Temporary directory where the temporary files are located.
    :return:
    """

    dumps = [os.path.join(tempdir, "output-{}.db".format(num)) for num in range(1, num_temp + 1)]
    locus_printer = functools.partial(print_locus,
                                      handles=out_handles,
                                      logger=logger,
                                      json_conf=json_conf)

    conns = [sqlite3.connect(dump) for dump in dumps]
    cursors = [conn.cursor() for conn in conns]
    curr_chrom = None
    common_index = dict()

    # cursor.execute("CREATE TABLE loci (counter integer, json blob)")
    temp_conn = sqlite3.connect(os.path.join(tempdir, "temp_store.db"))
    max_counter = temp_conn.execute("SELECT MAX(counter) FROM transcripts").fetchone()[0]

    decoder = json.Decoder(number_mode=json.NM_NATIVE)

    gene_counter = 0
    for dbindex, cursor in enumerate(cursors):
        d = dict((index[0], dbindex) for index in cursor.execute("SELECT counter FROM loci").fetchall())
        common_index.update(d)

    print_subloci = (out_handles[1][0] is not None)
    print_monoloci = (out_handles[2][0] is not None)
    if max_counter != max(common_index.keys()):
        raise KeyError("I am missing some loci! {} vs {}".format(
            max_counter, max(common_index.keys())))

    factory = Superlocus(None)

    for index in range(1, max(common_index.keys()) + 1):
        cursor = cursors[common_index[index]]
        try:
            stranded_loci = cursor.execute("SELECT json FROM loci WHERE counter=?", (str(index),)).fetchone()
        except ValueError:
            raise ValueError((index, type(index)))
        for stranded_locus_json in msgpack.loads(stranded_loci[0], raw=False):
            stranded_locus = factory
            stranded_locus.load_dict(decoder(stranded_locus_json),
                                     load_transcripts=False,
                                     print_monoloci=print_monoloci, print_subloci=print_subloci
                                     )
            if stranded_locus.chrom != curr_chrom:
                curr_chrom = stranded_locus.chrom
                gene_counter = 0
            gene_counter = locus_printer(stranded_locus, gene_counter)

    return


def remove_fragments(stranded_loci, json_conf, logger):

    """This method checks which loci are possible fragments, according to the
    parameters provided in the configuration file, and tags/remove them according
    to the configuration specification.

    :param stranded_loci: a list of the loci to consider for fragment removal
    :type stranded_loci: list[Superlocus]

    :param json_conf: the configuration dictionary
    :type json_conf: dict

    :param logger: the logger
    :type logger: logging.Logger

    """

    loci_to_check = {True: set(), False: set()}
    # mcdl = json_conf["pick"]["run_options"]["fragments_maximal_cds"]
    # mexons = json_conf["pick"]["run_options"]["fragments_maximal_exons"]
    # mcdna = json_conf["pick"]["run_options"]["fragments_maximal_cdna"]

    total = 0

    stranded_loci_dict = dict()
    loci_to_superloci = dict()

    start, end = float("inf"), float("-inf")
    chrom = None

    for stranded_locus in stranded_loci:
        stranded_loci_dict[stranded_locus.id] = stranded_locus
        start, end = min(stranded_locus.start, start), max(stranded_locus.end, end)
        chrom = chrom or stranded_locus.chrom
        for _, locus_instance in stranded_locus.loci.items():
            loci_to_superloci[locus_instance.id] = stranded_locus.id
            logger.debug("Assessing whether %s could be a fragment", _)
            total += 1
            is_fragment = locus_instance.is_putative_fragment()
            logger.debug("%s is a putative fragment: %s", _, is_fragment)
            loci_to_check[is_fragment].add(locus_instance)

    if len(loci_to_check[True]) == total:
        loci_to_check[False] = loci_to_check.pop(True)
        loci_to_check[True] = set()

    comparisons = collections.defaultdict(list)
    # Produce a list of duples

    for locus_to_check, gene in product(loci_to_check[True], loci_to_check[False]):
        is_to_be_filtered, comparison = gene.other_is_fragment(locus_to_check)
        if is_to_be_filtered is True:
            comparisons[locus_to_check.id].append(comparison)

    for locus in comparisons:
        if json_conf["pick"]["fragments"]["remove"] is True:
            # A bit convoluted: use the locus ID to find the correct superlocus, then delete the ID inside the SL.
            if locus not in stranded_loci_dict[loci_to_superloci[locus]].loci:
                logger.error("Locus %s has been lost from superlocus %s!",
                             locus,
                             loci_to_superloci[locus])
                continue
            del stranded_loci_dict[loci_to_superloci[locus]].loci[locus]
            # stranded_loci_dict[loci_to_superloci[locus]].loci.pop(locus, None)
        else:
            best_comparison = sorted(comparisons[locus], reverse=True, key=Assigner.get_f1)[0]
            stranded_loci_dict[loci_to_superloci[locus]].loci[locus].is_fragment = True
            stranded_loci_dict[loci_to_superloci[locus]].loci[locus].attributes[
                "fragment_of"] = best_comparison.ref_id[0]
            stranded_loci_dict[loci_to_superloci[locus]].loci[locus].attributes[
                "fragment_class_code"] = best_comparison.ccode[0]
            if best_comparison.distance[0] > 0:
                stranded_loci_dict[loci_to_superloci[locus]].loci[locus].attributes[
                    "distance"] = best_comparison.distance[0]

    for stranded_locus in stranded_loci:
        yield stranded_locus


def analyse_locus(slocus: Superlocus,
                  counter: int,
                  json_conf: dict,
                  printer_queue: [AutoProxy, None],
                  logging_queue: AutoProxy,
                  engine=None,
                  data_dict=None) -> [Superlocus]:

    """
    :param slocus: a superlocus instance
    :type slocus: Mikado.loci_objects.superlocus.Superlocus

    :param counter: an integer which is used to create the proper name for the locus.
    :type counter: int

    :param json_conf: the configuration dictionary
    :type json_conf: dict

    :param logging_queue: the logging queue
    :type logging_queue: multiprocessing.managers.AutoProxy

    :param printer_queue: the printing queue
    :type printer_queue: multiprocessing.managers.AutoProxy

    :param engine: an optional engine to connect to the database.
    :type data_dict: sqlalchemy.engine.engine

    :param data_dict: a dictionary of preloaded data
    :type data_dict: (None|dict)

    This function takes as input a "superlocus" instance and the pipeline configuration.
    It also accepts as optional keywords a dictionary with the CDS information
    (derived from a Bed12Parser) and a "lock" used for avoiding writing collisions
    during multithreading.
    The function splits the superlocus into its strand components and calls the relevant methods
    to define the loci.
    When it is finished, it transmits the superloci to the printer function.
    """

    # Define the logger
    if slocus is None:
        # printer_dict[counter] = []
        if printer_queue:
            while printer_queue.qsize() >= json_conf["threads"] * 10:
                continue
            # printer_queue.put_nowait(([], counter))
            return
        else:
            return []

    handler = logging_handlers.QueueHandler(logging_queue)
    logger = logging.getLogger("{0}:{1}-{2}".format(
        slocus.chrom, slocus.start, slocus.end))
    logger.addHandler(handler)

    # We need to set this to the lowest possible level,
    # otherwise we overwrite the global configuration
    logger.setLevel(json_conf["log_settings"]["log_level"])
    logger.propagate = False
    logger.debug("Started with %s, counter %d",
                 slocus.id, counter)
    if slocus.stranded is True:
        logger.warning("%s is stranded already! Resetting", slocus.id)
        slocus.stranded = False

    slocus.logger = logger
    slocus.source = json_conf["pick"]["output_format"]["source"]

    # if engine is None and data_dict is None:
    #     raise ValueError("Error for locus {id}".format(id=slocus.id))

    try:
        slocus.load_all_transcript_data(engine=engine,
                                        data_dict=data_dict)
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        logger.error("Error while loading data for %s", slocus.id)
        logger.exception(exc)
    logger.debug("Loading transcript data for %s", slocus.id)

    # Load the CDS information if necessary
    if slocus.initialized is False:
        # This happens when all transcripts have been removed from the locus,
        # due to errors that have been hopefully logged
        logger.warning(
            "%s had all transcripts failing checks, ignoring it",
            slocus.id)
        # printer_dict[counter] = []
        if printer_queue:
            while printer_queue.qsize >= json_conf["threads"] * 10:
                continue
            # printer_queue.put_nowait(([], counter))
            return
        else:
            return []

    # Split the superlocus in the stranded components
    logger.debug("Splitting by strand")
    stranded_loci = sorted([_ for _ in slocus.split_strands()])
    # Define the loci
    logger.debug("Divided into %d loci", len(stranded_loci))

    for stranded_locus in stranded_loci:
        stranded_locus.logger = logger
        try:
            stranded_locus.define_loci()
        except KeyboardInterrupt:
            raise
        except OSError:
            raise
        except Exception as exc:
            logger.exception(exc)
            logger.error("Removing failed locus %s", stranded_locus.name)
            stranded_loci.remove(stranded_locus)
        logger.debug("Defined loci for %s:%f-%f, strand: %s",
                     stranded_locus.chrom,
                     stranded_locus.start,
                     stranded_locus.end,
                     stranded_locus.strand)

    # Check if any locus is a fragment, if so, tag/remove it
    stranded_loci = sorted(list(remove_fragments(stranded_loci, json_conf, logger)))
    try:
        logger.debug("Size of the loci to send: {0}, for {1} loci".format(
            sys.getsizeof(stranded_loci),
            len(stranded_loci)))
    except Exception as err:
        logger.error(err)
        pass
    if printer_queue:
        while printer_queue.qsize() >= json_conf["threads"] * 10:
            continue
        logger.debug("Finished with %s, counter %d", slocus.id, counter)
        logger.removeHandler(handler)
        handler.close()
        return
    else:
        logger.debug("Finished with %s, counter %d", slocus.id, counter)
        logger.removeHandler(handler)
        handler.close()
        return stranded_loci


class LociProcesser(Process):

    """This process class takes care of getting from the queue the loci,
    analyse them, and print them to the output files."""

    def __init__(self,
                 json_conf,
                 locus_queue,
                 logging_queue,
                 identifier,
                 tempdir="mikado_pick_tmp"
                 ):

        # current_counter, gene_counter, current_chrom = shared_values
        super(LociProcesser, self).__init__()
        json_conf = msgpack.loads(json_conf, raw=False)
        self.logging_queue = logging_queue
        self.__identifier = identifier  # Property directly unsettable
        self.name = "LociProcesser-{0}".format(self.identifier)
        self.json_conf = json_conf
        self.engine = None
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.json_conf["log_settings"]["log_level"])

        self.logger.propagate = False
        self._tempdir = tempdir
        self.locus_queue = locus_queue
        self.regressor = None
        self.dump_db, self.dump_conn, self.dump_cursor = self._create_temporary_store(self._tempdir, self.identifier)

        if self.json_conf["pick"]["scoring_file"].endswith((".pickle", ".model")):
            with open(self.json_conf["pick"]["scoring_file"], "rb") as forest:
                self.regressor = pickle.load(forest)
            from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
            if not isinstance(self.regressor["scoring"], (RandomForestRegressor, RandomForestClassifier)):
                exc = TypeError("Invalid regressor provided, type: %s", type(self.regressor))
                self.logger.critical(exc)
                self.exitcode = 9
                self.join()
        self.logger.debug("Starting Process %s", self.name)

        self.logger.debug("Starting the pool for {0}".format(self.name))
        try:
            self.engine = dbutils.connect(self.json_conf, self.logger)
        except KeyboardInterrupt:
            raise
        except EOFError:
            raise
        except Exception as exc:
            self.logger.exception(exc)
            return

        self.analyse_locus = functools.partial(analyse_locus,
                                               printer_queue=None,
                                               json_conf=self.json_conf,
                                               engine=self.engine,
                                               logging_queue=self.logging_queue)

    @property
    def identifier(self):
        return self.__identifier

    @staticmethod
    def _create_temporary_store(tempdirectory, identifier):

        db = os.path.join(tempdirectory, "output-{}.db".format(identifier))
        conn = sqlite3.connect(db, isolation_level=None,
                               check_same_thread=False)
        cursor = conn.cursor()

        cursor.execute("CREATE TABLE IF NOT EXISTS loci (counter integer, json blob)")
        cursor.execute("CREATE INDEX IF NOT EXISTS loci_idx on loci(counter)")

        return db, conn, cursor

    def __getstate__(self):

        state = self.__dict__.copy()
        state["engine"] = None
        state["analyse_locus"] = None
        state["logger"].removeHandler(self.handler)
        del state["handler"]
        if "logger" in state:
            del state["logger"]
        if "dump_conn" in state:
            del state["dump_conn"]
        if "dump_cursor" in state:
            del state["dump_cursor"]
        return state

    def terminate(self):
        self.__close_handles()
        super().terminate()

    def __close_handles(self):
        """Private method to flush and close all handles."""
        if self.engine is not None:
            self.engine.dispose()
        if hasattr(self, "dump_conn"):
            self.dump_conn.commit()
            self.dump_conn.close()

    def close(self):
        self.__close_handles()
        self.join()

    def __setstate__(self, state):

        self.__dict__.update(state)
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.logger.propagate = False
        self.engine = dbutils.connect(self.json_conf, self.logger)
        self.analyse_locus = functools.partial(analyse_locus,
                                               printer_queue=None,
                                               json_conf=self.json_conf,
                                               engine=self.engine,
                                               logging_queue=self.logging_queue)
        self.dump_db, self.dump_conn, self.dump_cursor = self._create_temporary_store(self._tempdir, self.identifier)
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger.addHandler(self.handler)

    def join(self, timeout=None):
        self.__close_handles()
        # self.terminate()
        super().join(timeout=timeout)

    def run(self):
        """Start polling the queue, analyse the loci, and send them to the printer process."""
        self.logger.debug("Starting to parse data for {0}".format(self.name))
        # Read-only connection

        conn = sqlite3.connect("file:{}?mode=ro".format(os.path.join(self._tempdir, "temp_store.db")),
                               uri=True,  # Necessary to use the Read-only mode from file string
                               isolation_level="DEFERRED",
                               timeout=60,
                               check_same_thread=False  # Necessary for SQLite3 to function in multiprocessing
                               )
        cursor = conn.cursor()

        while True:
            counter = self.locus_queue.get()[0]
            if counter == "EXIT":
                self.logger.debug("EXIT received for %s", self.name)
                self.locus_queue.task_done()
                self.locus_queue.put((counter, ))
                # self.__close_handles()
                break
                # self.join()
            else:
                try:
                    transcripts = cursor.execute(
                        "SELECT json FROM transcripts WHERE counter=?", (str(counter),)).fetchone()
                except sqlite3.ProgrammingError as exc:
                    self.logger.exception(sqlite3.ProgrammingError((exc, counter, str(counter), (str(counter),))))
                    # self.__close_handles()
                    break
                    
                if transcripts is None:
                    raise KeyError("Nothing found in the database for %s", counter)

                transcripts = msgpack.loads(transcripts[0], raw=False)
                if len(transcripts) == 0:
                    stranded_loci = []
                else:
                    tobjects = []
                    for tjson in transcripts:
                        definition = GtfLine(tjson["definition"]).as_dict()
                        transcript = Transcript(logger=self.logger,
                                                source=definition["source"],
                                                intron_range=self.json_conf["pick"]["run_options"]["intron_range"])
                        transcript.chrom, transcript.start, transcript.end = (definition["chrom"],
                                                                              definition["start"], definition["end"])
                        try:
                            transcript.id = definition["transcript"]
                        except KeyError:
                            raise KeyError(definition)
                        transcript.strand, transcript.feature = definition["strand"], definition["feature"]
                        transcript.attributes = definition["attributes"]
                        try:
                            for exon in tjson["exon_lines"]:
                                start, end, feature, phase = exon
                                transcript.add_exon((start, end), feature=feature, phase=phase)
                            transcript.finalize()
                            tobjects.append(transcript)
                        except InvalidTranscript as exc:
                            self.logger.exception("Transcript %s is invalid. Ignoring. Error: %s",
                                                  transcript.id, exc)

                    slocus = Superlocus(tobjects.pop(),
                                        stranded=False,
                                        json_conf=self.json_conf,
                                        source=self.json_conf["pick"]["output_format"]["source"])
                    while len(tobjects) > 0:
                        slocus.add_transcript_to_locus(tobjects.pop(),
                                                       check_in_locus=False)

                    if self.regressor is not None:
                        slocus.regressor = self.regressor
                    stranded_loci = self.analyse_locus(slocus, counter)

                serialise_locus(stranded_loci, self.dump_conn, counter)

                self.locus_queue.task_done()

        return
