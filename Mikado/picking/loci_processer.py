from multiprocessing import Process
from multiprocessing.managers import AutoProxy
import logging
from itertools import product
import logging.handlers as logging_handlers
import functools
from ..utilities import dbutils
from ..scales.assigner import Assigner
from ..loci.superlocus import Superlocus
from ..parsers.GFF import GffLine
from ..serializers.external import ExternalSource
import os
import collections
import csv
import re
import sys
import pickle
from itertools import zip_longest
from sqlalchemy.engine import create_engine  # SQLAlchemy/DB imports
import sqlalchemy.orm.session

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

    foo_counter = ["failed", 0]

    for transcript in transcripts:
        current_transcript = current_gene["transcripts"][transcript]["transcript"]
        # Get the original transcript counter
        try:
            # first = re.sub("{0}\.".format(current_transcript.parent[0]), "", other)
            # transcript_counter = int(re.sub("\.orf[0-9]+", "", first))
            transcript_counter = int(re.sub("\.orf[0-9]+", "",
                                            current_transcript.id).split(".")[-1])
        except ValueError as exc:
            # Patch. For
            if isinstance(current_transcript.parent, list):
                foo_counter[1] += 1
                transcript_counter = "_".join(foo_counter)
            else:
                raise ValueError((exc, str(current_transcript)))
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


def print_locus(stranded_locus,
                gene_counter,
                handles,
                counter=None,
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
            if counter is not None:
                sub_lines = "\n".join(
                    ["{0}/{1}".format(counter, line) for line in sub_lines.split("\n")])
            print(sub_lines, file=sub_out)
        sub_metrics_rows = [_ for _ in stranded_locus.print_subloci_metrics()
                            if _ != {} and "tid" in _]
        sub_scores_rows = [_ for _ in stranded_locus.print_subloci_scores()
                           if _ != {} and "tid" in _]
        for row in sub_metrics_rows:
            if counter is not None:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
            sub_metrics.writerow(row)
        for row in sub_scores_rows:
            if counter is not None:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
            sub_scores.writerow(row)
    if mono_out is not None:
        mono_lines = stranded_locus.__str__(
            level="monosubloci",
            print_cds=not json_conf["pick"]["run_options"]["exclude_cds"])
        if mono_lines != '':
            mono_lines = "\n".join(
                ["{0}/{1}".format(counter, line) for line in mono_lines.split("\n")])
            print(mono_lines, file=mono_out)
        mono_metrics_rows = [_ for _ in stranded_locus.print_monoholder_metrics()
                             if _ != {} and "tid" in _]
        mono_scores_rows = [_ for _ in stranded_locus.print_monoholder_scores()
                            if _ != {} and "tid" in _]
        for row in mono_metrics_rows:
            if counter is not None:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
            mono_metrics.writerow(row)
        for row in mono_scores_rows:
            if counter is not None:
                row["tid"] = "{0}/{1}".format(counter, row["tid"])
            mono_scores.writerow(row)

    for locus in stranded_locus.loci:
        gene_counter += 1
        fragment_test = (
            json_conf["pick"]["clustering"]["remove_overlapping_fragments"]
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
        if counter is not None:
            locus_lines = "\n".join(
                ["{0}/{1}".format(counter, line) for line in locus_lines.split("\n")])
        print(locus_lines, file=locus_out)

    # assert len(locus_metrics_rows) == len(locus_scores_rows)

    for row in locus_metrics_rows:
        if counter is not None:
            row["tid"] = "{0}/{1}".format(counter, row["tid"])
        locus_metrics.writerow(row)
    for row in locus_scores_rows:
        if counter is not None:
            row["tid"] = "{0}/{1}".format(counter, row["tid"])
        locus_scores.writerow(row)
    # Necessary to flush out all the files
    [_.flush() for _ in handles if hasattr(_, "close")]
    return gene_counter


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

        # Parsing scores and metrics
        for current_index in sorted(current_lines.keys()):
            # current = min(current_lines.keys())
            for index, line in current_lines[current_index]:
                fields = line.split("\t")
                tid, gid = fields[:2]
                if (index, gid) not in gid_to_new:
                    raise KeyError("GID {} not found in {}!".format(
                        (index, gid), handle.name))
                if (index, tid) not in tid_to_new:
                    raise KeyError("TID {} not found in {}!".format(
                        (index, tid), handle.name))

                fields[0] = tid_to_new[(index, tid)]
                fields[1] = gid_to_new[(index, gid)]
                line = "\t".join(fields)
                print(line, file=handle, end="")
            # del current_lines[current]
        [_.close() for _ in filenames]
        [os.remove(_) for _ in finished]
        handle.close()
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

    for stranded_locus in stranded_loci:
        stranded_loci_dict[stranded_locus.id] = stranded_locus
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
        if json_conf["pick"]["clustering"]["remove_overlapping_fragments"] is True:
            # A bit convoluted: use the locus ID to find the correct superlocus, then delete the ID inside the SL.
            del stranded_loci_dict[loci_to_superloci[locus]].loci[locus]
        else:
            best_comparison = sorted(comparisons[locus], reverse=True, key=Assigner.get_f1)[0]
            stranded_loci_dict[loci_to_superloci[locus]].loci[locus].is_fragment = True
            stranded_loci_dict[loci_to_superloci[locus]].loci[locus].attributes["fragment_of"] = best_comparison.ref_id[0]
            stranded_loci_dict[loci_to_superloci[locus]].loci[locus].attributes["fragment_class_code"] = best_comparison.ccode[0]
            if best_comparison.distance[0] > 0:
                stranded_loci_dict[loci_to_superloci[locus]].loci[locus].attributes["distance"] = best_comparison.distance[0]

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
    stranded_loci = sorted([_ for _ in slocus.split_strands()])
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
            from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
            if not isinstance(self.regressor["scoring"], (RandomForestRegressor, RandomForestClassifier)):
                exc = TypeError("Invalid regressor provided, type: %s", type(self.regressor))
                self.logger.critical(exc)
                self.exitcode = 9
                self.join()

        self._create_handles(self.__output_files)
        self.__gene_counter = 0
        assert len(self._handles) > 0

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

        for name in ["locus_metrics", "locus_scores", "locus_out",
                     "sub_metrics", "sub_scores", "sub_out",
                     "mono_metrics", "mono_scores", "mono_out"]:
            state[name] = None
        state["engine"] = None
        state["analyse_locus"] = None
        del state["handler"]
        del state["logger"]
        return state

    def terminate(self):
        self.__close_handles()
        super().terminate()

    def __close_handles(self):
        """Private method to flush and close all handles."""

        for group in self._handles:
            [_.flush() for _ in group if hasattr(_, "flush") and _.closed is False]
            [_.close() for _ in group if hasattr(_, "close") and _.closed is False]
        if self.engine is not None:
            self.engine.dispose()

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

    def __create_step_handles(self, handles, metrics, score_keys):

        """Private method to create the handles for a given step (eg Locus).

        :param handles: the list with the filename prefixes
        :type handles: [list|tuple]

        :param metrics: list of metrics name, to be used as header for the metrics file
        :type metrics: list

        :param score_keys: list of metrics names used for scoring, to be used as header
        for the score file
        :type score_keys: list

        :returns: a list of handles to be used for writing
        :rtype: list
        """

        (locus_metrics_file,
         locus_scores_file,
         locus_out_file) = [os.path.join(self._tempdir,
                                         "{0}-{1}".format(os.path.basename(_),
                                                          self.identifier))
                            for _ in handles]
        locus_metrics_handle = open(locus_metrics_file, "a")
        locus_scores_handle = open(locus_scores_file, "a")
        locus_metrics = csv.DictWriter(
            locus_metrics_handle,
            metrics,
            delimiter="\t")
        locus_metrics.handle = locus_metrics_handle
        locus_metrics.close = locus_metrics.handle.close
        locus_metrics.closed = locus_metrics.handle.closed
        locus_metrics.flush = locus_metrics.handle.flush

        locus_scores = csv.DictWriter(locus_scores_handle, score_keys, delimiter="\t")
        locus_scores.handle = locus_scores_handle
        locus_scores.close = locus_scores.handle.close
        locus_scores.closed = locus_scores.handle.closed
        locus_scores.flush = locus_scores.handle.flush

        locus_out = open(locus_out_file, 'w')

        return [locus_metrics, locus_scores, locus_out]

    def _create_handles(self, handles):

        if self.regressor is None:
            score_keys = sorted(list(self.json_conf["scoring"].keys()))
        else:
            score_keys = self.regressor["scoring"].metrics
        # Define mandatory output files

        db_connection = functools.partial(
            dbutils.create_connector,
            self.json_conf,
            self.logger)

        engine = create_engine("{0}://".format(self.json_conf["db_settings"]["dbtype"]),
                               creator=db_connection)
        session = sqlalchemy.orm.sessionmaker(bind=engine)()

        score_keys = ["tid", "parent", "score"] + sorted(score_keys + ["source_score"])
        metrics = Superlocus.available_metrics[3:]
        metrics.extend(["external.{}".format(_.source) for _ in session.query(ExternalSource.source).all()])
        metrics = Superlocus.available_metrics[:3] + sorted(metrics)

        self._handles.append(self.__create_step_handles(handles[0],
                                                        metrics, score_keys))

        # Subloci
        if handles[1][0]:
            self._handles.append(self.__create_step_handles(handles[1],
                                                        metrics, score_keys))
        else:
            self._handles.append([None, None, None])

        # Monoloci
        if handles[2][0]:
            self._handles.append(self.__create_step_handles(handles[2],
                                                        metrics, score_keys))
        else:
            self._handles.append([None, None, None])

        return

    def join(self, timeout=None):
        self.__close_handles()
        # self.terminate()
        super().join(timeout=timeout)

    def run(self):
        """Start polling the queue, analyse the loci, and send them to the printer process."""
        self.logger.debug("Starting to parse data for {0}".format(self.name))
        current_chrom = None
        while True:
            slocus, counter = self.locus_queue.get()
            if slocus == "EXIT":
                self.logger.debug("EXIT received for %s", self.name)
                self.locus_queue.put((slocus, counter))
                self.__close_handles()
                break
                # self.join()
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
                    self.__gene_counter = print_locus(
                        stranded_locus, self.__gene_counter, self._handles,
                        counter=counter, logger=self.logger, json_conf=self.json_conf)

        return

