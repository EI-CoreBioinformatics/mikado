import zlib
from multiprocessing import Process
from typing import Union
from multiprocessing.managers import AutoProxy
import logging
from itertools import product
import logging.handlers as logging_handlers
import functools
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session, sessionmaker
from ..configuration import MikadoConfiguration, DaijinConfiguration
from ..utilities import dbutils
from ..scales.assignment.assigner import Assigner
from ..loci.superlocus import Superlocus
from ._merge_loci_utils import __create_gene_counters, manage_index
import collections
import sys
from ..transcripts import Transcript
from ..exceptions import InvalidTranscript, InvalidConfiguration, InvalidCDS, NotInLocusError
from ..parsers.GTF import GtfLine
from ..configuration.configurator import load_and_validate_config
import msgpack
from ._loci_serialiser import serialise_locus


__author__ = 'Luca Venturini'


def merge_loci(mapper,
               out_handles,
               total,
               configuration,
               logger=None,
               source="Mikado"):

    """ Function to merge the temporary loci files into single output files,
      renaming the genes according to the preferred style.
    :param num_temp: number of temporary files.
    :param out_handles: The names of the output loci files.
    :param logger: logger to use
    :return:
    """

    locus_metrics, locus_scores, locus_out = out_handles[0]
    sub_metrics, sub_scores, sub_out = out_handles[1]
    mono_metrics, mono_scores, mono_out = out_handles[2]

    if len(mapper["results"]) != total:
        raise KeyError("I am missing some loci! {} vs {}".format(len(mapper["results"]), total))

    # Start iterating the output dictionaries ("cursors")
    print_subloci = (out_handles[1][0] is not None)
    print_monoloci = (out_handles[2][0] is not None)
    __valid = set(range(1, total + 1))

    if set(mapper["results"].keys()) != __valid:
        missing = set.difference(__valid, set(mapper["results"].keys()))
        raise AssertionError("Missing the following loci: {}".format(missing))

    new_common, total_genes = __create_gene_counters(mapper["results"])

    done = 0
    logger.info("We have a total of %d genes", total_genes)
    tot_loci = set()
    for index in sorted(mapper["results"]):
        loci, batch = manage_index((index, new_common[index]), mapper["results"], source=source,
                                   configuration=configuration)
        done += 1
        if loci and set(loci).issubset(tot_loci):
            raise ValueError("Duplicated loci! {}".format(loci))

        tot_loci.update(set(loci))
        for minibatch in batch:
            if minibatch[0] is not None:
                locus_lines, locus_metrics_rows, locus_scores_rows = minibatch[0]
                if locus_lines:
                    assert len(locus_metrics_rows) > 0
                    print(locus_lines, file=locus_out)
                for row in locus_metrics_rows:
                    print(*[row[key] for key in locus_metrics.fieldnames],
                          sep="\t", file=locus_metrics.handle)
                for row in locus_scores_rows:
                    print(*[row[key] for key in locus_scores.fieldnames],
                          sep="\t", file=locus_scores.handle)

            if print_subloci and minibatch[1]:
                sub_lines, sub_metrics_rows, sub_scores_rows = minibatch[1]
                if sub_lines != '':
                    print(sub_lines, file=sub_out)

                for row in sub_metrics_rows:
                    try:
                        print(*[row[key] for key in sub_metrics.fieldnames],
                              sep="\t", file=sub_metrics.handle)
                    except KeyError:
                        raise KeyError("\n".join([str((key, str(row.get(key, "MISSING"))))
                                                  for key in sub_metrics.fieldnames]))
                for row in sub_scores_rows:
                    print(*[row[key] for key in sub_scores.fieldnames],
                              sep="\t", file=sub_scores.handle)
            if print_monoloci and minibatch[2]:
                mono_lines, mono_metrics_rows, mono_scores_rows = minibatch[2]
                if mono_lines != '':
                    print(mono_lines, file=mono_out)
                for row in mono_metrics_rows:
                    print(*[row[key] for key in mono_metrics.fieldnames],
                          sep="\t", file=mono_metrics.handle)
                for row in mono_scores_rows:
                    print(*[row[key] for key in mono_scores.fieldnames],
                          sep="\t", file=mono_scores.handle)

    if done != total:
        raise KeyError("Something has been lost")

    return


def remove_fragments(stranded_loci, configuration, logger):

    """This method checks which loci are possible fragments, according to the
    parameters provided in the configuration file, and tags/remove them according
    to the configuration specification.

    :param stranded_loci: a list of the loci to consider for fragment removal
    :type stranded_loci: list[Superlocus]

    :param configuration: the configuration dictionary
    :type configuration: (MikadoConfiguration|DaijinConfiguration)

    :param logger: the logger
    :type logger: logging.Logger

    """

    loci_to_check = {True: list(), False: list()}

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
            loci_to_check[is_fragment].append(locus_instance)

    if len(loci_to_check[True]) == total:
        loci_to_check[False] = loci_to_check.pop(True)
        loci_to_check[True] = list()

    comparisons = collections.defaultdict(list)
    # Produce a list of duples

    for locus_to_check, gene in product(loci_to_check[True], loci_to_check[False]):
        is_to_be_filtered, comparison = gene.other_is_fragment(locus_to_check)
        if is_to_be_filtered is True:
            comparisons[locus_to_check.id].append(comparison)

    for locus in comparisons:
        if configuration.pick.fragments.remove is True:
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
                  configuration: Union[MikadoConfiguration,DaijinConfiguration],
                  logging_queue: AutoProxy,
                  engine=None,
                  session=None) -> [Superlocus]:

    """
    :param slocus: a superlocus instance
    :type slocus: Mikado.loci_objects.superlocus.Superlocus

    :param counter: an integer which is used to create the proper name for the locus.
    :type counter: int

    :param configuration: the configuration dictionary
    :type configuration: (MikadoConfiguration|DaijinConfiguration)

    :param logging_queue: the logging queue
    :type logging_queue: multiprocessing.managers.AutoProxy

    :param engine: an optional engine to connect to the database.
    :type data_dict: sqlalchemy.engine.engine

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
        return []

    handler = logging_handlers.QueueHandler(logging_queue)
    logger = logging.getLogger("{0}:{1}-{2}".format(
        slocus.chrom, slocus.start, slocus.end))
    logger.addHandler(handler)

    # We need to set this to the lowest possible level,
    # otherwise we overwrite the global configuration
    logger.setLevel(configuration.log_settings.log_level)
    logger.propagate = False
    logger.debug("Started with %s, counter %d", slocus.id, counter)
    if slocus.stranded is True:
        logger.warning("%s is stranded already! Resetting", slocus.id)
        slocus.stranded = False

    slocus.logger = logger
    slocus.source = configuration.pick.output_format.source
    slocus.load_all_transcript_data(engine=engine, session=session)
    logger.debug("Loading transcript data for %s", slocus.id)

    # Load the CDS information if necessary
    if slocus.initialized is False:
        # This happens when all transcripts have been removed from the locus,
        # due to errors that have been hopefully logged
        logger.warning(
            "%s had all transcripts failing checks, ignoring it",
            slocus.id)
        return []

    # Split the superlocus in the stranded components
    logger.debug("Splitting by strand")
    stranded_loci = sorted([_ for _ in slocus.split_strands()])
    # Define the loci
    logger.debug("Divided into %d loci", len(stranded_loci))

    failed = 0
    original_size = len(stranded_loci)
    for stranded_locus in stranded_loci:
        stranded_locus.logger = logger
        try:
            stranded_locus.define_loci()
        except (InvalidTranscript, InvalidConfiguration, InvalidCDS, NotInLocusError) as exc:
            logger.exception(exc)
            logger.error("Removing failed locus %s", stranded_locus.name)
            stranded_loci.remove(stranded_locus)
            assert stranded_locus not in stranded_loci
            failed += 1
        logger.debug("Defined loci for %s:%f-%f, strand: %s",
                     stranded_locus.chrom,
                     stranded_locus.start,
                     stranded_locus.end,
                     stranded_locus.strand)

    assert len(stranded_loci) + failed == original_size, (len(stranded_loci), failed, original_size)
    # Check if any locus is a fragment, if so, tag/remove it
    if len(stranded_loci) > 0:
        stranded_loci = sorted(list(remove_fragments(stranded_loci, configuration, logger)))
    logger.debug("Size of the loci to send: {0}, for {1} loci".format(sys.getsizeof(stranded_loci), len(stranded_loci)))
    logger.debug("Finished with %s, counter %d", slocus.id, counter)
    logger.removeHandler(handler)
    handler.close()
    return stranded_loci


class LociProcesser(Process):

    """This process class takes care of getting from the queue the loci,
    analyse them, and print them to the output files."""

    def __init__(self,
                 configuration,
                 locus_queue,
                 logging_queue,
                 status_queue,
                 identifier,
                 ):

        super(LociProcesser, self).__init__()
        configuration = load_and_validate_config(msgpack.loads(configuration, raw=False))
        self.logging_queue = logging_queue
        self.status_queue = status_queue
        self.__identifier = identifier  # Property directly unsettable
        self.name = "LociProcesser-{0}".format(self.identifier)
        self.configuration = configuration
        for section in (self.configuration.scoring.requirements, self.configuration.scoring.as_requirements,
                        self.configuration.scoring.cds_requirements, self.configuration.scoring.not_fragmentary):
            # Compile the expression
            _ = section.compiled

        self.engine = None
        self.session = None
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.configuration.log_settings.log_level)

        self.logger.propagate = False
        self.locus_queue = locus_queue
        self.logger.debug("Starting Process %s", self.name)

        self.logger.debug("Starting the pool for {0}".format(self.name))

    @property
    def identifier(self):
        return self.__identifier

    def __getstate__(self):

        state = self.__dict__.copy()
        state["engine"] = None
        state["session"] = None
        state["analyse_locus"] = None
        state["logger"].removeHandler(self.handler)
        del state["handler"]
        if "logger" in state:
            del state["logger"]
        return state

    def terminate(self):
        self.__close_handles()
        super().terminate()

    def __close_handles(self):
        """Private method to flush and close all handles."""
        if self.engine is not None:
            self.engine.dispose()
        if self.session is not None:
            self.session.close_all()

    def close(self):
        self.__close_handles()
        self.join()

    def __setstate__(self, state):

        self.__dict__.update(state)
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.configuration.log_settings.log_level)
        self.logger.propagate = False
        self.engine = dbutils.connect(self.configuration, self.logger)
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger.addHandler(self.handler)
        self.connect_to_db(engine=None, session=None)

    def join(self, timeout=None):
        self.__close_handles()
        # self.terminate()
        super().join(timeout=timeout)

    def connect_to_db(self, engine: Union[Engine, None], session: Union[Session, None]):

        """
        :param engine: the connection pool
        :type engine: Engine

        :param session: a preformed session
        :type session: Session

        This method will connect to the database using the information
        contained in the JSON configuration.
        """

        if isinstance(session, Session):
            self.session = session
            self.sessionmaker = sessionmaker()
            self.sessionmaker.configure(bind=self.session.bind)
            self.engine = self.session.bind

        if engine is None:
            self.engine = dbutils.connect(self.configuration)
        else:
            self.engine = engine

        self.sessionmaker = sessionmaker()
        self.sessionmaker.configure(bind=self.engine)
        self.session = self.sessionmaker()

    def run(self):
        """Start polling the queue, analyse the loci, and send them to the printer process."""
        self.connect_to_db(None, None)
        self.logger.debug("Starting to parse data for {0}".format(self.name))
        self.analyse_locus = functools.partial(analyse_locus,
                                               configuration=self.configuration,
                                               engine=self.engine,
                                               session=self.session,
                                               logging_queue=self.logging_queue)

        print_cds = (not self.configuration.pick.run_options.exclude_cds)

        print_monoloci = (self.configuration.pick.files.monoloci_out is not None and
                          len(self.configuration.pick.files.monoloci_out) > 0)
        print_subloci = (self.configuration.pick.files.subloci_out is not None and
                         len(self.configuration.pick.files.subloci_out) > 0)

        accumulator = []
        while True:
            vals = self.locus_queue.get()
            try:
                counter, transcripts = vals
            except ValueError:
                raise ValueError(vals)

            if counter == "EXIT":
                self.logger.debug("EXIT received for %s", self.name)
                accumulator = zlib.compress(msgpack.dumps(accumulator))
                self.status_queue.put(accumulator)
                self.logger.debug("Put all the remaining results into the queue for %s (size %s)", self.name,
                                  sys.getsizeof(accumulator))
                self.locus_queue.task_done()
                self.locus_queue.put((counter, None))
                break
            else:
                try:
                    transcripts = msgpack.loads(transcripts, raw=False)
                except TypeError as err:
                    raise TypeError("{}, {}".format(err, transcripts))
                if len(transcripts) == 0:
                    stranded_loci = []
                    self.logger.warning("No transcript found for index %d", counter)
                else:
                    tobjects = []
                    chroms = set()
                    for tjson in transcripts:
                        definition = GtfLine(tjson["definition"]).as_dict()
                        is_reference = definition["source"] in self.configuration.prepare.files.reference
                        transcript = Transcript(logger=self.logger,
                                                source=definition["source"],
                                                intron_range=self.configuration.pick.run_options.intron_range,
                                                is_reference=is_reference)
                        transcript.chrom, transcript.start, transcript.end = (definition["chrom"],
                                                                              definition["start"], definition["end"])
                        chroms.add(transcript.chrom)
                        assert len(chroms) == 1, chroms
                        transcript.id = definition["transcript"]
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
                                        configuration=self.configuration,
                                        source=self.configuration.pick.output_format.source)
                    while len(tobjects) > 0:
                        slocus.add_transcript_to_locus(tobjects.pop(),
                                                       check_in_locus=False)

                    stranded_loci = self.analyse_locus(slocus, counter)

                if len(stranded_loci) == 0:
                    self.logger.warning("No loci left for index %d", counter)
                else:
                    if len(accumulator) >= 100:
                        accumulator = zlib.compress(msgpack.dumps(accumulator))
                        self.status_queue.put(accumulator)
                        accumulator = []

                    accumulator = serialise_locus(
                        stranded_loci, accumulator, counter, print_cds=print_cds,
                        print_monosubloci=print_monoloci, print_subloci=print_subloci)
                self.locus_queue.task_done()

        return
