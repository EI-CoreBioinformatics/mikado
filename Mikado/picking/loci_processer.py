from multiprocessing import Process
from multiprocessing.managers import AutoProxy
import logging
import logging.handlers as logging_handlers
import functools
from ..utilities import dbutils
from ..loci.superlocus import Superlocus
from ..parsers.GFF import GffLine
import os
import collections
import csv
import re
import sys
import pickle
from itertools import zip_longest

__author__ = 'Luca Venturini'


def print_gene(current_gene, gene_counter, handle, prefix):

    """
    This function takes a gene and reformats it using the new derived name.
    :param current_gene: a dictionary with the data for the current gene.
    :type current_gene: dict
    :param gene_counter: A counter to be used inside the genes
    :type gene_counter: int
    :param handle: The handle to the file to print to
    :type handle: io.TextIOWrapper
    :param prefix: the prefix to add to the name of the genes
    :type prefix: str
    :return: a dictionary with the correspondeces for the transcripts
    :rtype: dict
    """

    # Reset the name
    current_gene["gene"].name = current_gene["gene"].id
    print(current_gene["gene"], file=handle)
    tid_corrs = dict()
    chrom = current_gene["gene"].chrom

    transcripts = [_ for _ in current_gene["transcripts"]]

    transcripts = sorted(transcripts,
                         key=lambda _:
                         (current_gene["transcripts"][_]["transcript"].start,
                          current_gene["transcripts"][_]["transcript"].end))
    # transcript_counter = 1

    for transcript in transcripts:
        current_transcript = current_gene["transcripts"][transcript]["transcript"]
        # Get the original transcript counter
        try:
            # first = re.sub("{0}\.".format(current_transcript.parent[0]), "", other)
            # transcript_counter = int(re.sub("\.orf[0-9]+", "", first))
            transcript_counter = int(re.sub("\.orf[0-9]+", "",
                                            current_transcript.id).split(".")[-1])
        except ValueError:
            assert isinstance(current_transcript.parent, list),\
                type(current_transcript.parent)
            raise
        assert transcript_counter >= 1

        tid = "{0}.{1}G{2}.{3}".format(prefix,
                                       chrom,
                                       gene_counter,
                                       transcript_counter)
        if ".orf" in transcript:
            tid = "{0}{1}".format(tid,
                                  transcript[transcript.find(".orf"):])

        current_transcript.parent = current_gene["gene"].id
        current_exons = current_gene["transcripts"][transcript]["exons"]
        current_transcript.attributes["Alias"] = current_transcript.id[:]
        # name = re.sub("\.orf[0-9]+", "", tid)
        # tid_corrs[re.sub("\.orf[0-9]+", "", current_transcript.id)] = name
        tid_corrs[current_transcript.id] = tid
        current_transcript.id = tid
        current_transcript.name = tid
        print(current_transcript, file=handle)
        for exon in current_exons:
            exon.parent = tid
            exon.id = re.sub(current_transcript.attributes["Alias"],
                             tid, exon.id)
            exon.name = re.sub(current_transcript.attributes["Alias"],
                               tid, exon.id)
            print(exon, file=handle)
    return tid_corrs


def merge_loci_gff(gff_filenames, gff_handle, prefix=""):

    """
    This function will merge different partial GFF files into a single loci file,
    while changing the names to reflect the ordering.
    :param gff_filenames:
    :param gff_handle:
    :param prefix:
    :return:
    """

    current_lines = dict()
    gffs = [open(_) for _ in gff_filenames]

    for lines in zip_longest(*gffs):
        for num, line in enumerate(lines):
            if line is None:
                continue
            _ = line.split("/")
            index = int(_[0])
            if index not in current_lines:
                current_lines[index] = {"filenum": num,
                                        "lines": []}
            else:
                assert current_lines[index]["filenum"] == num, (num,
                                                                current_lines[index])
            current_lines[index]["lines"].append("/".join(_[1:]))

    [_.close() for _ in gffs]

    gid_to_new = dict()
    tid_to_new = dict()

    with open(gff_handle, "a") as gff_handle:
        current_gene = dict()
        current_chrom = None
        gene_counter = 0
        for index in sorted(current_lines.keys()):
            file_index = current_lines[index]["filenum"]
            lines = [GffLine(_) for _ in current_lines[index]["lines"]]
            for line in lines:
                if line.header is True:
                    continue
                if current_chrom is not None and current_chrom != line.chrom:
                    gene_counter = 0
                    current_chrom = line.chrom
                elif current_chrom is None:
                    current_chrom = line.chrom
                # Start the printing process
                if line.is_gene:
                    if current_gene != dict():
                        tid_corrs = print_gene(current_gene,
                                               gene_counter,
                                               gff_handle,
                                               prefix)
                        for tid in tid_corrs:
                            assert (file_index, tid) not in tid_to_new, (file_index, tid)
                            tid_to_new[(file_index, tid)] = tid_corrs[tid]
                        current_gene = dict()
                    current_gene["transcripts"] = dict()

                    # Create the correspondence for the new gene
                    gene_counter += 1
                    new_id = "{0}.{1}G{2}".format(prefix, line.chrom, gene_counter)
                    assert (file_index, line.id) not in gid_to_new, ((file_index, line.id),
                                                                     gid_to_new)
                    gid_to_new[(file_index, line.id)] = new_id
                    line.id = new_id
                    current_gene["gene"] = line
                elif line.is_transcript:
                    assert current_gene != dict()
                    current_gene["transcripts"][line.id] = dict()
                    current_gene["transcripts"][line.id]["transcript"] = line
                    current_gene["transcripts"][line.id]["exons"] = []
                    if line.attributes["primary"].lower() in ("true", "false"):
                        if line.attributes["primary"].lower() == "true":
                            primary = True
                        else:
                            primary = False
                        current_gene["transcripts"][line.id]["primary"] = primary
                    else:
                        raise ValueError("Invalid value for \"primary\" field: {0}".format(
                            line.attributes["primary"]))
                elif line.is_exon:
                    for parent in line.parent:
                        assert parent in current_gene["transcripts"]
                        current_gene["transcripts"][parent]["exons"].append(line)
                else:
                    if current_gene != dict():
                        tid_corrs = print_gene(current_gene,
                                               gene_counter,
                                               gff_handle,
                                               prefix)
                        for tid in tid_corrs:
                            assert (file_index, tid) not in tid_to_new, (file_index, tid)
                            tid_to_new[(file_index, tid)] = tid_corrs[tid]
                        current_gene = dict()
                        print("###", file=gff_handle)

                    print(line, file=gff_handle)
                    continue

            if current_gene != dict():
                tid_corrs = print_gene(current_gene, gene_counter, gff_handle, prefix)
                for tid in tid_corrs:
                    assert (file_index, tid) not in tid_to_new, (file_index, tid)
                    tid_to_new[(file_index, tid)] = tid_corrs[tid]
                current_gene = dict()
                print("###", file=gff_handle)
            del current_lines[index]

    [os.remove(_) for _ in gff_filenames]
    return gid_to_new, tid_to_new


def merge_loci(num_temp, out_handles, prefix="", tempdir="mikado_pick_tmp"):

    """ Function to merge the temporary loci files into single output files,
      renaming the genes according to the preferred style.
    :param num_temp: number of temporary files.
    :param out_handles: The names of the output loci files.
    :param prefix: Prefix to use for the gene names.
    :param tempdir: Temporary directory where the temporary files are located.
    :return:
    """

    metrics_handle, scores_handle, gff_handle = out_handles

    gff_filenames = [os.path.join(tempdir,
                                  "{0}-{1}".format(os.path.basename(gff_handle),
                                                   _))
                     for _ in range(1, num_temp + 1)]

    gid_to_new, tid_to_new = merge_loci_gff(gff_filenames, gff_handle, prefix)

    for handle in metrics_handle, scores_handle:
        filenames = [os.path.join(tempdir,
                                  "{0}-{1}".format(os.path.basename(handle), _))
                     for _ in range(1, num_temp + 1)]
        handle = open(handle, "a")
        current_lines = collections.defaultdict(list)
        filenames = [open(_) for _ in filenames]
        finished = set()
        while len(finished) < len(filenames):
            for num, _ in enumerate(filenames):
                if _.name in finished:
                    continue
                else:
                    try:
                        line = next(_)
                        fields = line.split("/")
                        current_lines[int(fields[0])].append((num, "/".join(fields[1:])))
                    except StopIteration:
                        _.close()
                        finished.add(_.name)

        while len(current_lines) > 0:
            current = min(current_lines.keys())
            for index, line in current_lines[current]:
                fields = line.split("\t")
                tid, gid = fields[:2]
                assert (index, gid) in gid_to_new, ((index, gid),
                                                    "\n".join(str(_) for _ in gid_to_new.items()))
                assert (index, tid) in tid_to_new, (index, tid, tid_to_new)

                fields[0] = tid_to_new[(index, tid)]
                fields[1] = gid_to_new[(index, gid)]
                line = "\t".join(fields)
                print(line, file=handle, end="")
            del current_lines[current]
        [os.remove(_) for _ in finished]
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
    mcdl = json_conf["pick"]["run_options"]["fragments_maximal_cds"]
    mexons = json_conf["pick"]["run_options"]["fragments_maximal_exons"]

    total = 0
    for stranded_locus in stranded_loci:
        for _, locus_instance in stranded_locus.loci.items():
            locus_instance.logger = logger
            total += 1
            to_check = (locus_instance.primary_transcript.combined_cds_length < mcdl)
            loci_to_check[to_check].add(locus_instance)

    if len(loci_to_check[True]) == total:
        loci_to_check[False] = loci_to_check[True].copy()

    bool_remove_fragments = json_conf["pick"]["run_options"]["remove_overlapping_fragments"]
    for stranded_locus in stranded_loci:
        to_remove = set()
        for locus_id, locus_instance in stranded_locus.loci.items():
            if locus_instance in loci_to_check[True]:
                logger.debug("Checking if %s is a fragment", locus_instance.id)

                for other_locus in iter(
                        olocus for olocus in loci_to_check[False]
                        if olocus.primary_transcript_id != locus_instance.primary_transcript_id):
                    if other_locus.other_is_fragment(
                            locus_instance,
                            minimal_cds_length=mcdl,
                            minimal_exons=mexons) is True:
                        if bool_remove_fragments is False:
                            # Just mark it as a fragment
                            stranded_locus.loci[locus_id].is_fragment = True
                        else:
                            to_remove.add(locus_id)
                            # del stranded_locus.loci[locus_id]
                        break
        for locus_id in to_remove:
            del stranded_locus.loci[locus_id]
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
            while printer_queue.qsize() >= json_conf["pick"]["run_options"]["procs"] * 10:
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
        logger.warn("%s is stranded already! Resetting",
                    slocus.id)
        slocus.stranded = False

    slocus.logger = logger
    slocus.source = json_conf["pick"]["output_format"]["source"]

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
            while printer_queue.qsize >= json_conf["pick"]["run_options"]["procs"] * 10:
                continue
            # printer_queue.put_nowait(([], counter))
            return
        else:
            return []

    # Split the superlocus in the stranded components
    logger.debug("Splitting by strand")
    stranded_loci = sorted(list(slocus.split_strands()))
    # Define the loci
    logger.debug("Divided into %d loci", len(stranded_loci))

    for stranded_locus in stranded_loci:
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

    # Remove overlapping fragments.
    loci_to_check = {True: set(), False: set()}
    for stranded_locus in stranded_loci:
        for _, locus_instance in stranded_locus.loci.items():
            locus_instance.logger = logger
            loci_to_check[locus_instance.monoexonic].add(locus_instance)

    # Check if any locus is a fragment, if so, tag/remove it
    stranded_loci = sorted(list(remove_fragments(stranded_loci, json_conf, logger)))
    try:
        logger.debug("Size of the loci to send: {0}, for {1} loci".format(
            sys.getsizeof(stranded_loci),
            len(stranded_loci)))
    except Exception as err:
        logger.error(err)
        pass
    # printer_dict[counter] = stranded_loci
    if printer_queue:
        while printer_queue.qsize() >= json_conf["pick"]["run_options"]["procs"] * 10:
            continue
        # printer_queue.put_nowait((stranded_loci, counter))
        # printer_queue.put((stranded_loci, counter))
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
                 data_dict,
                 output_files,
                 locus_queue,
                 logging_queue,
                 identifier,
                 tempdir="mikado_pick_tmp"
                 ):

        # current_counter, gene_counter, current_chrom = shared_values
        super(LociProcesser, self).__init__()
        self.logging_queue = logging_queue
        self.__identifier = identifier  # Property directly unsettable
        self.name = "LociProcesser-{0}".format(self.identifier)
        self.json_conf = json_conf
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.logger.propagate = False
        self._tempdir = tempdir

        self.__data_dict = data_dict
        self.locus_queue = locus_queue
        # self.lock = lock
        self.__output_files = output_files
        self.locus_metrics, self.locus_scores, self.locus_out = [None] * 3
        self.sub_metrics, self.sub_scores, self.sub_out = [None] * 3
        self.mono_metrics, self.mono_scores, self.mono_out = [None] * 3
        self._handles = []
        self.regressor = None

        if self.json_conf["pick"]["scoring_file"].endswith((".pickle", ".model")):
            with open(self.json_conf["pick"]["scoring_file"], "rb") as forest:
                self.regressor = pickle.load(forest)
            from sklearn.ensemble import RandomForestRegressor
            if not isinstance(self.regressor, RandomForestRegressor):
                exc = TypeError("Invalid regressor provided, type: %s", type(self.regressor))
                self.logger.critical(exc)
                self.exitcode = 9
                self.join()

        self._create_handles(self.__output_files)
        self.__gene_counter = 0
        assert self.locus_out is not None
        self.logger.debug("Starting Process %s", self.name)

        self.logger.debug("Starting the pool for {0}".format(self.name))
        try:
            if self.json_conf["pick"]["run_options"]["preload"] is False:
                self.engine = dbutils.connect(self.json_conf, self.logger)
            else:
                self.engine = None
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
                                               data_dict=self.__data_dict,
                                               engine=self.engine,
                                               logging_queue=self.logging_queue)

    @property
    def identifier(self):
        return self.__identifier

    def __getstate__(self):

        state = self.__dict__.copy()
        for h in state["_handles"]:
            h.close()

        state["_handles"] = []

        for name in ["locus_metrics", "locus_scores", "locus_out",
                     "sub_metrics", "sub_scores", "sub_out",
                     "mono_metrics", "mono_scores", "mono_out"]:
            state[name] = None
        state["engine"] = None
        state["analyse_locus"] = None
        del state["handler"]
        del state["logger"]
        return state

    def __setstate__(self, state):

        self.__dict__.update(state)
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.logger.propagate = False

        self._create_handles(self.__output_files)
        if self.json_conf["pick"]["run_options"]["preload"] is False:
            self.engine = dbutils.connect(self.json_conf, self.logger)
        else:
            self.engine = None
        self.analyse_locus = functools.partial(analyse_locus,
                                               printer_queue=None,
                                               json_conf=self.json_conf,
                                               data_dict=self.__data_dict,
                                               engine=self.engine,
                                               logging_queue=self.logging_queue)

    def _create_handles(self, handles):

        (locus_metrics_file,
         locus_scores_file,
         locus_out_file) = [os.path.join(self._tempdir,
                                         "{0}-{1}".format(os.path.basename(_),
                                                          self.identifier))
                            for _ in handles[0]]
        locus_metrics_file = open(locus_metrics_file, "w")
        locus_scores_file = open(locus_scores_file, "w")

        if self.regressor is None:
            score_keys = sorted(list(self.json_conf["scoring"].keys()))
        else:
            score_keys = self.regressor.metrics
        # Define mandatory output files

        score_keys = ["tid", "parent", "score"] + sorted(score_keys + ["source_score"])
        self.locus_metrics = csv.DictWriter(
            locus_metrics_file,
            Superlocus.available_metrics,
            delimiter="\t")

        self.locus_scores = csv.DictWriter(locus_scores_file, score_keys, delimiter="\t")

        self.locus_metrics.handle = locus_metrics_file
        self.locus_metrics.flush = self.locus_metrics.handle.flush
        self.locus_metrics.close = self.locus_metrics.handle.close
        self.locus_scores.handle = locus_scores_file
        self.locus_scores.flush = self.locus_scores.handle.flush
        self.locus_scores.close = self.locus_scores.handle.close

        self.locus_out = open(locus_out_file, 'w')
        self._handles.extend((self.locus_out, self.locus_metrics, self.locus_out))

        if handles[1][0]:
            (sub_metrics_file,
             sub_scores_file,
             sub_out_file) = [os.path.join(self._tempdir,
                                           "{0}-{1}".format(os.path.basename(_),
                                                            self.identifier))
                              for _ in handles[1]]
            sub_metrics_file = open(sub_metrics_file, "w")
            sub_scores_file = open(sub_scores_file, "w")
            self.sub_metrics = csv.DictWriter(
                sub_metrics_file,
                Superlocus.available_metrics,
                delimiter="\t")
            self.sub_metrics.handle = sub_metrics_file
            self.sub_metrics.flush = self.sub_metrics.handle.flush
            self.sub_metrics.close = self.sub_metrics.handle.close
            self.sub_scores = csv.DictWriter(
                sub_scores_file, score_keys, delimiter="\t")
            self.sub_scores.handle = sub_scores_file
            self.sub_scores.flush = self.sub_scores.handle.flush
            self.sub_scores.close = self.sub_scores.handle.close
            self.sub_out = open(sub_out_file, "w")
            self._handles.extend([self.sub_metrics, self.sub_scores, self.sub_out])

        if handles[2][0]:
            (mono_metrics_file,
             mono_scores_file,
             mono_out_file) = [os.path.join(self._tempdir,
                                            "{0}-{1}".format(os.path.basename(_),
                                                             self.identifier))
                               for _ in handles[2]]
            mono_metrics_file = open(mono_metrics_file, "w")
            mono_scores_file = open(mono_scores_file, "w")
            self.mono_metrics = csv.DictWriter(
                mono_metrics_file,
                Superlocus.available_metrics,
                delimiter="\t")
            self.mono_metrics.handle = mono_metrics_file
            self.mono_metrics.flush = self.mono_metrics.handle.flush
            self.mono_metrics.close = self.mono_metrics.handle.close
            self.mono_scores = csv.DictWriter(
                mono_scores_file, score_keys, delimiter="\t")
            self.mono_scores.handle = mono_scores_file
            self.mono_scores.flush = self.mono_scores.handle.flush
            self.mono_scores.close = self.mono_scores.handle.close
            self.mono_out = open(mono_out_file, "w")
            self._handles.extend([self.mono_metrics, self.mono_scores, self.mono_out])            

        return

    def run(self):
        """Start polling the queue, analyse the loci, and send them to the printer process."""
        self.logger.debug("Starting to parse data for {0}".format(self.name))
        current_chrom = None
        while True:
            slocus, counter = self.locus_queue.get()
            if slocus == "EXIT":
                self.logger.debug("EXIT received for %s", self.name)
                self.locus_queue.put((slocus, counter))
                if self.engine is not None:
                    self.engine.dispose()
                self.locus_out.close()
                self.locus_scores.close()
                self.locus_metrics.close()
                if self.sub_metrics is not None:
                    self.sub_metrics.close()
                    self.sub_scores.close()
                    self.sub_out.close()
                if self.mono_out is not None:
                    self.mono_out.close()

                return
            else:
                if slocus is not None:
                    if current_chrom != slocus.chrom:
                        self.__gene_counter = 0
                        current_chrom = slocus.chrom
                    if self.regressor is not None:
                        slocus.regressor = self.regressor
                    stranded_loci = self.analyse_locus(slocus, counter)
                else:
                    stranded_loci = []
                for stranded_locus in stranded_loci:
                    self._print_locus(stranded_locus, counter)

    def _print_locus(self, stranded_locus, counter):

        """
        Private method that handles a single superlocus for printing.
        It also detects and flags/discard fragmentary loci.
        :param stranded_locus: the stranded locus to analyse
        :return:
        """

        if self.sub_out is not None:  # Skip this section if no sub_out is defined
            sub_lines = stranded_locus.__str__(
                level="subloci",
                print_cds=not self.json_conf["pick"]["run_options"]["exclude_cds"])
            if sub_lines != '':
                sub_lines = "\n".join(
                    ["{0}/{1}".format(counter, line) for line in sub_lines.split("\n")])
                print(sub_lines, file=self.sub_out)
            sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()
                                if x != {} and "tid" in x]
            sub_scores_rows = [x for x in stranded_locus.print_subloci_scores()
                               if x != {} and "tid" in x]
            for row in sub_metrics_rows:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
                self.sub_metrics.writerow(row)
            for row in sub_scores_rows:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
                self.sub_scores.writerow(row)
        if self.mono_out is not None:
            mono_lines = stranded_locus.__str__(
                level="monosubloci",
                print_cds=not self.json_conf["pick"]["run_options"]["exclude_cds"])
            if mono_lines != '':
                mono_lines = "\n".join(
                    ["{0}/{1}".format(counter, line) for line in mono_lines.split("\n")])
                print(mono_lines, file=self.mono_out)
            mono_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()
                                 if x != {} and "tid" in x]
            mono_scores_rows = [x for x in stranded_locus.print_subloci_scores()
                                if x != {} and "tid" in x]
            for row in mono_metrics_rows:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
                self.mono_metrics.writerow(row)
            for row in mono_scores_rows:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
                self.mono_scores.writerow(row)

        for locus in stranded_locus.loci:
            fragment_test = (
                self.json_conf["pick"]["run_options"]["remove_overlapping_fragments"]
                is True and stranded_locus.loci[locus].is_fragment is True)

            if fragment_test is True:
                continue
            self.__gene_counter += 1
            new_id = "{0}.{1}G{2}".format(
                self.json_conf["pick"]["output_format"]["id_prefix"],
                stranded_locus.chrom, self.__gene_counter)
            stranded_locus.loci[locus].id = new_id

        locus_lines = stranded_locus.__str__(
            print_cds=not self.json_conf["pick"]["run_options"]["exclude_cds"],
            level="loci")

        locus_metrics_rows = [x for x in stranded_locus.print_loci_metrics()]
        locus_scores_rows = [x for x in stranded_locus.print_loci_scores()]

        # assert len(locus_metrics_rows) == len(locus_scores_rows)

        for row in locus_metrics_rows:
            row["tid"] = "{0}/{1}".format(counter, row["tid"])
            self.locus_metrics.writerow(row)
        for row in locus_scores_rows:
            row["tid"] = "{0}/{1}".format(counter, row["tid"])
            self.locus_scores.writerow(row)

        if locus_lines != '':
            locus_lines = "\n".join(
                    ["{0}/{1}".format(counter, line) for line in locus_lines.split("\n")])
            print(locus_lines, file=self.locus_out)
