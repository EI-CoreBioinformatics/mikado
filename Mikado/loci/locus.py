# coding: utf-8

"""
This module defines the last object to be created during the picking,
i.e. the locus.
"""

import collections
import itertools
import operator
from collections import defaultdict
import pysam
from ..transcripts.transcript import Transcript
from ..transcripts.transcriptchecker import TranscriptChecker
from .abstractlocus import Abstractlocus, rgetattr
from ..parsers.GFF import GffLine
from ..scales.assignment.assigner import Assigner
from ..exceptions import InvalidTranscript
import networkx as nx
import io
from pkg_resources import resource_stream
import rapidjson as json
import jsonschema
import random


with io.TextIOWrapper(resource_stream("Mikado.configuration",
                                      "configuration_blueprint.json")) as blue:
        blue_print = json.loads(blue.read())


_valid_ccodes = jsonschema.Draft7Validator(
            blue_print["properties"]["pick"]["properties"]["alternative_splicing"]["properties"]["valid_ccodes"]
        )


_valid_redundant = jsonschema.Draft7Validator(
            blue_print["properties"]["pick"]["properties"]["alternative_splicing"]["properties"]["redundant_ccodes"]
        )


class Locus(Abstractlocus):
    """Class that defines the final loci.
    It is a child of monosublocus, but it also has the possibility of adding
    additional transcripts if they are valid splicing isoforms.
    """

    def __init__(self, transcript=None, logger=None, json_conf=None,
                 pad_transcripts=None, **kwargs):
        """
        Constructor class. Like all loci, also Locus is defined starting from a transcript.

        :param transcript: the transcript which is used to initialize the Locus
        :type transcript: [Transcript|None]

        :param logger: the logger instance.
        :type logger: None | logging.Logger
        """

        self.counter = 0  # simple tag to avoid collisions
        Abstractlocus.__init__(self, logger=logger, json_conf=json_conf, **kwargs)
        if transcript is not None:
            transcript.attributes["primary"] = True
            if transcript.is_coding:
                transcript.feature = "mRNA"
            else:
                transcript.feature = "ncRNA"
            self.monoexonic = transcript.monoexonic
            Abstractlocus.add_transcript_to_locus(self, transcript)
            self.locus_verified_introns = transcript.verified_introns
            self.tid = transcript.id
            self.logger.debug("Created Locus object with {0}".format(transcript.id))
            self.primary_transcript_id = transcript.id

        # this must be defined straight away
        # A set of the transcript we will ignore during printing
        # because they are duplications of the original instance. Done solely to
        # get the metrics right.
        self.__orf_doubles = collections.defaultdict(set)
        self.excluded = None
        self.parent = None

        self.attributes = dict()
        self.attributes["is_fragment"] = False
        self.metric_lines_store = []
        self.__id = None
        self.fai = None
        self.__finalized = False

        self._reference_sources = set(source for source, is_reference in
                                      zip(self.json_conf["prepare"]["files"]["labels"],
                                          self.json_conf["prepare"]["files"]["reference"]) if is_reference is True)

        if pad_transcripts in (False, True):
            # self._set_padding(pad_transcripts)
            self.json_conf["pick"]["alternative_splicing"]["pad"] = pad_transcripts
            assert self.perform_padding == pad_transcripts

        if self.perform_padding is True and self.reference_update is True:
            self._add_to_alternative_splicing_codes("=")
            self._add_to_alternative_splicing_codes("_")
            self._add_to_alternative_splicing_codes("n")

    def __str__(self, print_cds=True) -> str:

        self.feature = self.__name__
        assert self.feature != "Monosublocus"
        # Hacky fix to make sure that the primary transcript has the attribute
        # Set to True in any case.
        self.primary_transcript.attributes["primary"] = True
        # BF, really just a hack.
        for transcript in self.transcripts:
            if transcript == self.primary_transcript_id:
                continue
            self.transcripts[transcript].attributes["primary"] = False

        lines = []

        self_line = GffLine('')
        for attr in ["chrom", 'feature', 'source', 'start', 'end', 'strand']:
            setattr(self_line, attr, getattr(self, attr))
        self_line.phase, self_line.score = None, self.score
        self_line.id = self.id
        self_line.name = self.name
        self_line.attributes["superlocus"] = self.parent
        self_line.attributes.update(self.attributes)
        if "is_fragment" in self.attributes and self.attributes["is_fragment"] is False:
            del self_line.attributes["is_fragment"]
        self_line.attributes["multiexonic"] = (not self.monoexonic)
        lines.append(str(self_line))

        for tid in self.transcripts:
            transcript_instance = self.transcripts[tid]
            transcript_instance.source = self.source
            transcript_instance.parent = self_line.id
            self.logger.debug(self.attributes)
            if self.transcripts[tid].is_coding:
                self.transcripts[tid].feature = "mRNA"
            else:
                self.transcripts[tid].feature = "ncRNA"
            for attribute in self.attributes:
                if attribute not in transcript_instance.attributes:
                    if attribute == "is_fragment" and self.attributes[attribute] is False:
                        continue
                    transcript_instance.attributes[attribute] = self.attributes[attribute]

            lines.append(transcript_instance.format(
                "gff", with_cds=print_cds,
                all_orfs=self.json_conf["pick"]["output_format"]["report_all_orfs"]
            ).rstrip())

        return "\n".join(lines)

    def __setstate__(self, state):
        super(Locus, self).__setstate__(state)
        self._not_passing = set(self._not_passing)

    def finalize_alternative_splicing(self, _scores=None):

        """"This method ensures that all the transcripts retained in the locus
        are within the score threshold. This is due to the fact that the score
        changes depending on the transcript considered together; so that a transcript
        that might have scored relatively well on its own will score pretty badly when
        brought inside the locus.

        Please note that, **when doing a reference update**, reference transcripts that have passed the previous checks
        will *always* be retained and will not count for the purposes of determining the maximum number of isoforms.

        :param _scores:only for testing. Load directly a score mock-up.
        """

        if self._finalized is True:
            self.logger.debug("Locus already finalised, returning")
            return
        else:
            self.logger.debug("Starting to finalise AS for %s", self.id)

        if _scores is not None:
            self._load_scores(_scores)
        else:
            self.metrics_calculated = False
            self.scores_calculated = False
            self.logger.debug("Re-calculating metrics and scores for %s", self.id)
            self.filter_and_calculate_scores()
            self.logger.debug("Re-calculated metrics and scores for %s", self.id)

        max_isoforms = self.json_conf["pick"]["alternative_splicing"]["max_isoforms"]
        original = dict((tid, self.transcripts[tid].copy()) for tid in self.transcripts)

        # *Never* lose the primary transcript
        reference_sources = set(source for source, is_reference in zip(
            self.json_conf["prepare"]["files"]["labels"],
            self.json_conf["prepare"]["files"]["reference"]
        ))
        reference_transcripts = dict(
            (tid, self.transcripts[tid].is_reference or self.transcripts[tid].original_source in reference_sources)
            for tid in self.transcripts)

        order = sorted([[tid, self.transcripts[tid].score] for tid in self.transcripts
                       if tid != self.primary_transcript_id and
                        (reference_transcripts[tid] is False or
                         self.json_conf["pick"]["run_options"]["reference_update"] is False)],
                       key=operator.itemgetter(1), reverse=True)
        threshold = self.json_conf["pick"]["alternative_splicing"]["min_score_perc"] * self.primary_transcript.score

        score_passing = [_ for _ in order if _[1] >= threshold]
        removed = set([_[0] for _ in order if _[1] < threshold])

        remainder = score_passing[max_isoforms:]
        score_passing = score_passing[:max_isoforms]

        [self.remove_transcript_from_locus(tid) for tid in removed]
        [self.remove_transcript_from_locus(tup[0]) for tup in remainder]

        self.logger.debug("Kept %s transcripts out of %s (threshold: %s)", len(score_passing), len(order), threshold)

        iteration = 0
        while True:
            iteration += 1
            self.logger.debug("Starting iteration %s", iteration)
            if _scores is None:
                    self.metrics_calculated = False
                    self.scores_calculated = False
                    self.filter_and_calculate_scores()
            else:
                pass
            with_retained = self._remove_retained_introns()
            if len(with_retained) > 0:
                self.logger.debug("Transcripts with retained introns: %s", ", ".join(list(with_retained)))
            else:
                self.logger.debug("No transcripts with retained introns found.")
            if self.perform_padding is True and len(self.transcripts) > 1:
                self.logger.debug("Starting padding procedure for %s", self.id)
                failed = self.__launch_padding()
                if failed:
                    # Restart the padding procedure
                    continue
                else:
                    removed.update(self.__remove_redundant_after_padding())
                    self.primary_transcript.attributes.pop("ccode", None)
                    self.logger.debug("Updated removal set: %s", ", ".join(removed))
            missing = len(self.transcripts) - max_isoforms
            if missing > 0 and remainder:
                __to_add = remainder[:missing]
                for tid, score, start, obj in __to_add:
                    self.add_transcript_to_locus(obj, check_in_locus=False)
                remainder = remainder[missing:]
            else:
                break

        for tid in self.transcripts:
            assert tid in original
            backup = original[tid]
            transcript = self.transcripts[tid]
            if (transcript.start, transcript.end) != (backup.start, backup.end):
                self.transcripts[tid].attributes["padded"] = True
                message = "{transcript.id} is now padded and has now start {transcript.start}, end {transcript.end}"
                if (backup.is_coding and ((backup.combined_cds_end != transcript.combined_cds_end) or
                                          (backup.combined_cds_start != transcript.combined_cds_start))):
                    transcript.attributes["cds_padded"] = True
                    message += "; CDS moved to {transcript.combined_cds_start}, end {transcript.combined_cds_end}"
                elif backup.is_coding:
                    transcript.attributes["cds_padded"] = False
                message += "."
                self.logger.info(message.format(**locals()))

        # Now that we have added the padding ... time to remove redundant alternative splicing events.
        self.logger.debug("%s has %d transcripts (%s)", self.id, len(self.transcripts),
                          ", ".join(list(self.transcripts.keys())))
        self._finalized = True
        return

    def __launch_padding(self):

        self.logger.debug("Launched padding for %s", self.id)
        failed = False
        backup = dict()
        for tid in self.transcripts:
            backup[tid] = self.transcripts[tid].deepcopy()

        # The "templates" are the transcripts that we used to expand the others.
        templates = self.pad_transcripts()
        # First off, let us update the transcripts.
        tid_keys = list(self.transcripts.keys())
        for tid in tid_keys:
            self.logger.debug("Swapping %s", tid)
            self._swap_transcript(backup[tid], self.transcripts[tid])

        self.logger.debug("Done padding for %s", self.id)
        self.metrics_calculated = False
        self.scores_calculated = False
        self.filter_and_calculate_scores()
        self.logger.debug("Recalculated metrics after padding in %s", self.id)

        self._not_passing = set()
        self._check_requirements()
        if self.primary_transcript_id in self._not_passing or set.intersection(templates, self._not_passing):
            self.logger.debug(
                "Either the primary or some template transcript has not passed the muster. Removing, restarting.")
            if self._not_passing == {self.primary_transcript_id}:
                self.logger.info("The primary transcript %s is invalidated by the other transcripts in the locus.\
                Leaving only the main transcript in %s.", self.primary_transcript_id, self.id)
                self._not_passing = set(self.transcripts.keys()) - {self.primary_transcript_id}
            # Templates are clearly wrong. Remove them
            for tid in self._not_passing - {self.primary_transcript_id}:
                self.remove_transcript_from_locus(tid)
            for tid in set(backup.keys()) - (self._not_passing - {self.primary_transcript_id}):
                self.logger.debug("Swapping the old transcript for %s", tid)
                self._swap_transcript(self.transcripts[tid], backup[tid])
            self.metrics_calculated = False
            self.scores_calculated = False
            self.filter_and_calculate_scores()
            failed = True
            return failed

        order = sorted([(tid, self.transcripts[tid].score) for tid in self.transcripts
                        if tid != self.primary_transcript_id],
                       key=operator.itemgetter(1), reverse=True)

        self.logger.debug("Transcripts potentially kept in the locus: %s", ",".join([_[0] for _ in order]))

        # Now that we are sure that we have not ruined the primary transcript, let us see whether
        # we should discard any other transcript.
        [self._add_to_alternative_splicing_codes(ccode) for ccode in ("=", "_")]
        removed = set()
        added = set()

        for tid, score in order:
            is_valid, ccode, _ = self.is_alternative_splicing(self.transcripts[tid], others=added)
            if is_valid:
                self.logger.debug("Keeping %s in the locus, ccode: %s", tid, ccode)
                self.transcripts[tid].attributes["ccode"] = ccode
                added.add(tid)
            else:
                self.logger.debug("Removing %s from the locus after padding, ccode: %s", tid, ccode)
                self.remove_transcript_from_locus(tid)
                removed.add(tid)

        removed.update(self._remove_retained_introns())
        # Now let us check whether we have removed any template transcript.
        # If we have, remove the offending ones and restart

        if len(set.intersection(set.union(set(templates), {self.primary_transcript_id}), removed)) > 0:
            self.logger.debug("Removed: %s; Templates: %s; Primary: %s", ",".join(removed), ",".join(templates),
                              self.primary_transcript_id)
            self.transcripts = backup
            [self.remove_transcript_from_locus(tid) for tid in set.intersection(templates, removed)]
            self.metrics_calculated = False
            self.scores_calculated = False
            self.filter_and_calculate_scores()
            failed = True
            self.logger.debug("Padding failed for %s (removed: %s), restarting", self.id, removed)
            return failed

        return failed

    def _remove_retained_introns(self):
        self.logger.debug("Now checking the retained introns for %s", self.id)
        removed = set()
        while True:
            cds_disrupted = set()
            retained_introns = set()
            to_remove = set()
            for tid, transcript in self.transcripts.items():
                if tid == self.primary_transcript_id:
                    continue
                self.find_retained_introns(transcript)
                self.transcripts[tid].attributes["retained_intron"] = (transcript.retained_intron_num > 0)
                if transcript.retained_intron_num > 0:
                    retained_introns.add(tid)
                if transcript.cds_disrupted_by_ri is True:
                    assert transcript.retained_intron_num > 0
                    cds_disrupted.add(tid)
            if max(len(retained_introns), len(cds_disrupted)) == 0:
                break
            if self.json_conf["pick"]["alternative_splicing"]["keep_cds_disrupted_by_ri"] is False:
                self.logger.debug("Removing {} because their CDS is disrupted by retained introns".format(
                    ", ".join(list(cds_disrupted))))
                to_remove.update(cds_disrupted)
                retained_introns -= cds_disrupted
            if self.json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] is False:
                self.logger.debug("Removing {} because they contain retained introns".format(
                    ", ".join(list(retained_introns))))
                to_remove.update(retained_introns)
            if len(to_remove) > 0:
                removed.update(to_remove)
                for tid in to_remove:
                    self.remove_transcript_from_locus(tid)
                self.__segmenttree = self._calculate_segment_tree(self.exons, self.introns)
                self.metrics_calculated = False
                self.scores_calculated = False
                self.filter_and_calculate_scores()
            else:
                break

        return removed

    def __remove_redundant_after_padding(self):

        # First thing: calculate the class codes
        to_remove = set()
        if len(self.transcripts) == 1:
            return to_remove

        self.logger.debug("Starting to remove redundant transcripts from %s", self.id)

        class_codes = dict()
        ichains = collections.defaultdict(list)

        for tid in self.transcripts:
            if self.transcripts[tid].introns:
                key = tuple(sorted(self.transcripts[tid].introns))
            else:
                key = None
            self.logger.debug("Intron key for %s: %s", tid, key)
            ichains[key].append(tid)

        for ichain in ichains:
            for t1, t2 in itertools.combinations(ichains[ichain], 2):
                class_codes[(t1, t2)] = Assigner.compare(self[t1], self[t2])[0]
                self.logger.debug("Comparing ichain %s, %s vs %s: nF1 %s", ichain, t1, t2, class_codes[(t1, t2)].n_f1)

        for couple, comparison in class_codes.items():
            # if self[couple[0]].is_reference and self[couple[1]].is_reference:
            #     self.logger.debug("Both %s and %s are references, continuing", couple[0], couple[1])
            #     continue

            removal = None
            if comparison.n_f1[0] == 100:
                try:
                    if couple[0] == self.primary_transcript_id:
                        removal = couple[1]
                    elif couple[1] == self.primary_transcript_id:
                        removal = couple[0]
                    elif self[couple[0]].is_reference and not self[couple[1]].is_reference:
                        removal = couple[1]
                    elif self[couple[1]].is_reference and not self[couple[0]].is_reference:
                        removal = couple[0]
                    elif self[couple[0]].score > self[couple[1]].score:
                        removal = couple[1]
                    elif self[couple[1]].score > self[couple[0]].score:
                        removal = couple[0]
                    else:
                        # removal = np.random.choice(sorted(couple))
                        removal = random.choice(sorted(couple))
                except (TypeError, ValueError):
                    raise ValueError((couple, self[couple[0]].score, self[couple[1]].score))
                finally:
                    if removal:
                        to_remove.add(removal)
                        self.logger.debug("Removing %s from locus %s because after padding it is redundant with %s",
                                          removal, self.id, (set(couple) - {removal}).pop())

        if to_remove:
            self.logger.debug("Removing from %s: %s", self.id, ", ".join(to_remove))
            for tid in to_remove:
                self.remove_transcript_from_locus(tid)

        [self._add_to_redundant_splicing_codes(_) for _ in ("=", "_", "n", "c")]
        [self._remove_from_alternative_splicing_codes(_) for _ in ("=", "_", "n", "c")]

        order = sorted([(tid, self[tid].score) for tid in self if tid != self.primary_transcript_id],
                       key=operator.itemgetter(1))

        others = set(self.transcripts.keys())

        for tid, score in order:
            if self[tid].is_reference:
                continue
            is_valid, main_ccode, main_result = self.is_alternative_splicing(self[tid], others=others)
            if is_valid is False:
                self.logger.debug("Removing %s from %s as it is a redundant splicing isoform after padding.",
                                  tid, self.id)
                to_remove.add(tid)
                others.remove(tid)

        if to_remove:
            self.logger.debug("Removing from %s: %s", self.id, ", ".join(to_remove))
            for tid in to_remove:
                self.remove_transcript_from_locus(tid)

        return to_remove

    def as_dict(self):
        # self.calculate_scores()
        return super().as_dict()

    def remove_transcript_from_locus(self, tid: str):

        """Overloading of the AbstractLocus class, in order to ensure that the primary transcript will *not*
        be removed."""

        if tid == self.primary_transcript_id:
            raise KeyError("%s is the primary transcript of %s!" % (tid, self.id))

        super().remove_transcript_from_locus(tid)
        self._finalized = False

    def add_transcript_to_locus(self, transcript: Transcript, check_in_locus=True,
                                **kwargs):
        """Implementation of the add_transcript_to_locus method.
        Before a transcript is added, the class checks that it is a valid splicing isoform
        and that we have not exceeded already the maximum number of isoforms for the Locus.

        The checks performed are, in order:

        #. whether the locus already has the maximum number of acceptable
           isoforms ("max_isoforms")
        #. (optional) whether all the introns *specific to the transcript when
           compared with the primary transcript* are confirmed by external
           validation tools (eg Portcullis)
        #. Whether the score of the proposed AS event has a score over the
           minimum percentage of the primary transcript score (eg if the minimum
           percentage is 0.6 and the primary is scored 20, a model with a score
           of 11 would be rejected and one with a score of 12 would be accepted)
        #. Whether the strand of the candidate is the same as the one of the locus
        #. Whether the AS event is classified (ie has a class code) which is acceptable as valid AS
        #. Whether the transcript shares enough cDNA with the primary transcript ("min_cdna_overlap")
        #. Whether the proposed model has too much UTR
        #. (optional) Whether the proposed model has a retained intron compared
           to the primary, ie part of its 3' non-coding regions overlaps
           one intron of the primary model
        #. Whether the proposed model shares enough CDS with the primary model (min_cds_overlap)

        :param transcript: the candidate transcript
        :type transcript: Transcript

        :param check_in_locus: boolean flag to bypass *all* checks. Use with caution!
        :type check_in_locus: bool

        :param kwargs: optional keyword arguments are ignored.
        """

        _ = kwargs
        if check_in_locus is False:
            Abstractlocus.add_transcript_to_locus(self, transcript, check_in_locus=False)
            self.locus_verified_introns.update(transcript.verified_introns)
            return

        to_be_added = True

        if to_be_added and transcript.strand != self.strand:
            self.logger.debug("%s not added because it has a different strand from %s (%s vs. %s)",
                              transcript.id, self.id, transcript.strand, self.strand)
            to_be_added = False

        reference_pass = False
        if self.reference_update is True and \
                (transcript.is_reference is True or transcript.original_source in self._reference_sources):
            reference_pass = True

        if self.json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] is True and reference_pass is False:
            to_check = (transcript.introns - transcript.verified_introns) - self.primary_transcript.introns
            to_be_added = len(to_check) == 0
            if not to_be_added:
                self.logger.debug(
                    "%s not added because it has %d non-confirmed intron%s",
                    transcript.id,
                    len(to_check),
                    "s" * min(1, len(to_check) - 1))

        # Add a check similar to what we do for the minimum requirements and the fragments
        if to_be_added and "as_requirements" in self.json_conf:
            to_be_added = self.__check_as_requirements(transcript, is_reference=reference_pass)

        if to_be_added is True:
            is_alternative, ccode, _ = self.is_alternative_splicing(transcript)
            if is_alternative is False:
                self.logger.debug("%s not added because it is not a \
                valid splicing isoform. Ccode: %s",
                                  transcript.id, ccode)
                to_be_added = False
            else:
                transcript.attributes["ccode"] = ccode
                self.logger.debug("%s is a valid splicing isoform; Ccode: %s", transcript.id, ccode)

        if to_be_added is False:
            return

        self.logger.debug("Keeping %s as a valid alternative isoform for %s",
                          transcript.id, self.id)
        transcript.attributes["primary"] = False

        if transcript.is_coding:
            transcript.feature = "mRNA"
        else:
            transcript.feature = "ncRNA"
        Abstractlocus.add_transcript_to_locus(self, transcript)

        self.locus_verified_introns.update(transcript.verified_introns)

    def __check_as_requirements(self, transcript: Transcript, is_reference=False) -> bool:

        to_be_added = True
        if is_reference is True and self.json_conf["pick"]["run_options"]["check_references"] is False:
            return True
        if ("compiled" not in self.json_conf["as_requirements"] or
                self.json_conf["as_requirements"]["compiled"] is None):
            self.json_conf["as_requirements"]["compiled"] = compile(
                self.json_conf["as_requirements"]["expression"], "<json>",
                "eval")
        evaluated = dict()
        for key in self.json_conf["as_requirements"]["parameters"]:
            value = rgetattr(transcript,
                             self.json_conf["as_requirements"]["parameters"][key]["name"])
            if "external" in key:
                value = value[0]

            evaluated[key] = self.evaluate(
                value,
                self.json_conf["as_requirements"]["parameters"][key])
            # pylint: disable=eval-used
        if eval(self.json_conf["as_requirements"]["compiled"]) is False:
            self.logger.debug("%s fails the minimum requirements for AS events", transcript.id)
            to_be_added = False
        return to_be_added

    def is_intersecting(self, *args):
        """Not implemented: this function makes no sense for a single-transcript container.
        :param args: any argument to this nethod will be ignored.
        """
        raise NotImplementedError("""Loci do not use this method, but rather
        assess whether a transcript is a splicing isoform or not.""")

    def is_putative_fragment(self):

        """This method will use the expression in the "not_fragmentary" section
        of the configuration to determine whether it is itself a putative fragment."""

        if not self.json_conf["pick"]["run_options"]["check_references"] and \
                any(self.transcripts[tid].is_reference is True for tid in self.transcripts):
            return False

        self.json_conf["not_fragmentary"]["compiled"] = compile(
            self.json_conf["not_fragmentary"]["expression"], "<json>",
            "eval")

        current_id = self.id[:]

        evaluated = dict()
        for key in self.json_conf["not_fragmentary"]["parameters"]:
            value = rgetattr(self.primary_transcript,
                             self.json_conf["not_fragmentary"]["parameters"][key]["name"])
            if "external" in key:
                value = value[0]
            try:
                evaluated[key] = self.evaluate(
                    value,
                    self.json_conf["not_fragmentary"]["parameters"][key])
            except Exception as err:
                self.logger.error(
                    """Exception while calculating putative fragments. Key: {}, \
                    Transcript value: {} (type {}) \
                    configuration value: {} (type {}).""".format(
                        key, value, type(value), self.json_conf["not_fragmentary"]["parameters"][key],
                        type(self.json_conf["not_fragmentary"]["parameters"][key])
                    ))
                self.logger.exception(err)
                raise err
        if eval(self.json_conf["not_fragmentary"]["compiled"]) is True:
            self.logger.debug("%s cannot be a fragment according to the definitions, keeping it",
                              self.id)
            fragment = False
        else:
            self.logger.debug(
                "%s could be a fragment according to the definitions, tagging it for analysis",
                self.id)
            fragment = True

        self.id = current_id
        assert self.id == current_id
        return fragment

    def other_is_fragment(self,
                          other):
        """
        :param other: another Locus to compare against
        :type other: Locus

        This function checks whether another *monoexonic* Locus
        *on the opposite strand* is a fragment,by checking its classification
        according to Assigner.compare.
        """

        if not isinstance(self, type(other)):
            raise TypeError("I can compare only loci.")

        if other.primary_transcript_id == self.primary_transcript_id:
            self.logger.debug("Self-comparisons are not allowed!")
            return False, None

        if any(other.transcripts[tid].is_reference is True for tid in other.transcripts.keys()):
            self.logger.debug("Locus %s has a reference transcript, hence it will not be discarded", other.id)
            return False, None

        self.logger.debug("Comparing %s with %s",
                          self.primary_transcript_id,
                          other.primary_transcript_id)

        result, _ = Assigner.compare(other.primary_transcript,
                                     self.primary_transcript,
                                     strict_strandedness=True)

        # We should get rid of the "max_distance" value here.
        max_distance = max(0,
                           min(self.json_conf["pick"].get("fragments", {}).get("max_distance", float("inf")),
                               self.json_conf["pick"]["clustering"]["flank"]))

        self.logger.debug("Comparison between {0} (strand {3}) and {1}: class code \"{2}\"".format(
            self.primary_transcript.id,
            other.primary_transcript.id,
            result.ccode[0],
            other.strand))
        if (result.ccode[0] in self.json_conf["pick"]["fragments"]["valid_class_codes"] and
                    result.distance[0] <= max_distance):
            self.logger.debug("{0} is a fragment (ccode {1})".format(
                other.primary_transcript.id, result.ccode[0]))
            return True, result

        return False, None

    def set_json_conf(self, jconf: dict):
        """
        Setter for the configuration dictionary.
        :param jconf:
        :type jconf: dict
        """
        if not isinstance(jconf, dict):
            raise TypeError("Invalid configuration of type {0}".format(type(jconf)))
        self.json_conf = jconf

    def calculate_metrics(self, tid: str):
        """
        :param tid: the name of the transcript to be analysed
        :type tid: str

        This function will calculate the metrics for a transcript which are relative in nature
        i.e. that depend on the other transcripts in the sublocus. Examples include the fraction
        of introns or exons in the sublocus, or the number/fraction of retained introns.
        """

        self.logger.debug("Calculating metrics for %s", tid)
        self.transcripts[tid].finalize()
        if (self.transcripts[tid].number_internal_orfs <= 1 or
                    self.json_conf["pick"]["output_format"]["report_all_orfs"] is False):
            super().calculate_metrics(tid)
        else:
            transcript = self.transcripts[tid]
            selected = transcript.selected_internal_orf
            new_transcript = transcript.copy()
            new_transcript.id = "{0}.orf1".format(new_transcript.id)
            self.transcripts[new_transcript.id] = new_transcript
            super().calculate_metrics(new_transcript.id)
            self.__orf_doubles[tid].add(new_transcript.id)

            for num, orf in enumerate([_ for _ in transcript.internal_orfs if
                                       _ != selected]):
                new_transcript = transcript.copy()
                assert isinstance(new_transcript, Transcript)
                new_transcript.internal_orfs = [orf]
                new_transcript.internal_orfs.extend([_ for _ in transcript.internal_orfs if
                                                     orf != _])
                new_transcript.id = "{0}.orf{1}".format(new_transcript.id, num + 2)
                self.transcripts[new_transcript.id] = new_transcript
                super().calculate_metrics(new_transcript.id)
                self.__orf_doubles[tid].add(new_transcript.id)

        self.logger.debug("Calculated metrics for {0}".format(tid))

    def filter_and_calculate_scores(self):
        """
        Function to calculate a score for each transcript, given the metrics derived
        with the calculate_metrics method and the scoring scheme provided in the JSON configuration.
        If any requirements have been specified, all transcripts which do not pass them
        will be assigned a score of 0 and subsequently ignored.
        Scores are rounded to the nearest integer.
        """

        super().filter_and_calculate_scores()

        self.logger.debug("Calculated scores for %s, now checking for double IDs", self.id)
        for index, item in enumerate(reversed(self.metric_lines_store)):
            if item["tid"] in self.__orf_doubles:
                del self.metric_lines_store[index]
            else:
                continue

        for doubled in self.__orf_doubles:
            for partial in self.__orf_doubles[doubled]:
                if partial in self.transcripts:
                    del self.transcripts[partial]
                if partial in self.scores:
                    del self.scores[partial]

    def print_scores(self):
        """This method yields dictionary rows that are given to a csv.DictWriter class."""
        self.filter_and_calculate_scores()
        if self.regressor is None:
            score_keys = sorted(list(self.json_conf["scoring"].keys()) + ["source_score"])
        else:
            score_keys = sorted(self.regressor.metrics + ["source_score"])
        keys = ["tid", "alias", "parent", "score"] + score_keys

        for tid in self.scores:
            row = dict().fromkeys(keys)
            row["tid"] = tid
            row["parent"] = self.id
            row["alias"] = self.transcripts[tid].alias
            if tid in self._not_passing:
                row["score"] = 0
            else:
                row["score"] = round(self.scores[tid]["score"], 2)
            calculate_total = (self.regressor is None)
            for key in score_keys:
                if calculate_total:
                    assert self.scores[tid][key] != "NA" and self.scores[tid][key] is not None
                    row[key] = round(self.scores[tid][key], 2)

            if calculate_total is True:
                score_sum = sum(row[key] for key in score_keys)
                if tid not in self._not_passing and self.scores[tid]["score"] > 0:
                    assert round(score_sum, 2) == round(self.scores[tid]["score"], 2), (
                        score_sum,
                        self.transcripts[tid].score,
                        tid)
                else:
                    assert self.scores[tid]["score"] == 0

            yield row

    def is_alternative_splicing(self, other, others=None):

        """This function defines whether another transcript could be a
        putative alternative splice variant of the primary Locus
        transcript.
        To do so, it compares the candidate against all transcripts in the Locus, and calculates
        the class code using scales.Assigner.compare.
        If all the matches are "n" or "j", the transcript is considered as an AS event.

        :param other: another transcript to compare against
        :type other: Transcript

        :param others: a set of AS transcripts used to check whether the new transcript should be included.
        This allows to avoid creating a whole new locus for testing the padding.
        :type others: (None|set)
        """

        is_valid = True

        valid_ccodes = self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"]
        redundant_ccodes = self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"]
        cds_only = self.json_conf["pick"]["alternative_splicing"]["cds_only"]

        if other.is_coding and not self.primary_transcript.is_coding:
            reason = "{} is coding, and cannot be added to a non-coding locus.".format(other.id)
            enough_overlap, overlap_reason = False, reason
            main_ccode = "NA"
            main_result = None
        elif (other.is_coding is False) and (self.primary_transcript.is_coding is True):
            reason = "{} is non-coding, and cannot be added to a coding locus.".format(other.id)
            enough_overlap, overlap_reason = False, reason
            main_ccode = "NA"
            main_result = None
        else:
            if cds_only is True:
                self.logger.debug("Checking whether the CDS of %s and %s are overlapping enough",
                                  other.id, self.primary_transcript_id)
                main_result, _ = Assigner.compare(other._selected_orf_transcript,
                                                  self.primary_transcript._selected_orf_transcript)

                enough_overlap, overlap_reason = self._evaluate_transcript_overlap(
                    other._selected_orf_transcript,
                    self.primary_transcript._selected_orf_transcript,
                    min_cdna_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"],
                    min_cds_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"],
                    comparison=main_result,
                    check_references=self.json_conf["pick"]["run_options"]["check_references"],
                    fixed_perspective=True)
            else:
                main_result, _ = Assigner.compare(other,
                                                  self.primary_transcript)
                enough_overlap, overlap_reason = self._evaluate_transcript_overlap(
                    other,
                    self.primary_transcript,
                    min_cdna_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"],
                    min_cds_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"],
                    comparison=main_result,
                    check_references=self.json_conf["pick"]["run_options"]["check_references"],
                    fixed_perspective=True)

            self.logger.debug(overlap_reason)

            main_ccode = main_result.ccode[0]
            if main_ccode not in valid_ccodes:
                self.logger.debug("%s is not a valid splicing isoform. Ccode: %s",
                                  other.id,
                                  main_result.ccode[0])
                is_valid = False

        if not enough_overlap:
            self.logger.debug("%s is not a valid splicing isoform. Reason: %s",
                              other.id, overlap_reason)
            is_valid = False

        if (is_valid and self.json_conf["pick"]["run_options"]["check_references"] is False and
                other.is_reference is True):
            pass
        elif is_valid:
            if others is None:
                others = self.transcripts.keys()
            else:
                diff = set.difference(set(others), set(self.transcripts.keys()))
                assert diff == set(), diff

            for tid in iter(tid for tid in others if tid not in (self.primary_transcript_id, other.id)):
                candidate = self.transcripts[tid]
                if cds_only is True:
                    result, _ = Assigner.compare(other._selected_orf_transcript, candidate._selected_orf_transcript)
                else:
                    result, _ = Assigner.compare(other, candidate)
                if result.ccode[0] in redundant_ccodes:
                    self.logger.debug("%s is a redundant isoform of %s (ccode %s)",
                                      other.id, candidate.id, result.ccode[0])
                    is_valid = False
                    break

        if is_valid:
            self.logger.debug("%s is a valid splicing isoform of %s (class code %s, overlap: %s)",
                              other.id, self.id, main_ccode, overlap_reason)

        return is_valid, main_ccode, main_result

    def pad_transcripts(self) -> set:

        """
        """

        try:
            if isinstance(self.json_conf["reference"]["genome"], pysam.FastaFile):
                self.fai = self.json_conf["reference"]["genome"]
            else:
                self.fai = pysam.FastaFile(self.json_conf["reference"]["genome"])
        except KeyError:
            raise KeyError(self.json_conf.keys())

        five_graph = self.define_graph(objects=self.transcripts, inters=self._share_extreme, three_prime=False)
        three_graph = self.define_graph(objects=self.transcripts, inters=self._share_extreme, three_prime=True)

        self.logger.debug("5' graph: %s", five_graph.edges)
        self.logger.debug("3' graph: %s", three_graph.edges)
        # TODO: Tie breaks!

        __to_modify = self._find_communities_boundaries(five_graph, three_graph)

        self.logger.debug("To modify: %s", __to_modify)
        templates = set()

        # Now we can do the proper modification
        for tid in sorted(__to_modify.keys()):
            if __to_modify[tid][0]:
                templates.add(__to_modify[tid][0].id)
            if __to_modify[tid][1]:
                templates.add(__to_modify[tid][1].id)

            self.logger.debug("Expanding %s to have start %s (from %s) and end %s (from %s)",
                              tid, __to_modify[tid][0] if not __to_modify[tid][0] else __to_modify[tid][0].start,
                              self[tid].start,
                              __to_modify[tid][1] if not __to_modify[tid][1] else __to_modify[tid][1].end,
                              self[tid].end)
            try:
                new_transcript = expand_transcript(self[tid].deepcopy(),
                                                   __to_modify[tid][0],
                                                   __to_modify[tid][1],
                                                   self.fai,
                                                   self.logger)
            except KeyboardInterrupt:
                raise
            except Exception as exc:
                self.logger.exception(exc)
                raise
            if (new_transcript.start == self.transcripts[tid].end) and (new_transcript.end == self.transcripts[tid].end):
                self.logger.debug("No expansion took place for %s!", tid)
            else:
                self.logger.debug("Expansion took place for %s!", tid)
            self.transcripts[tid] = new_transcript

        self.exons = set()
        for tid in self:
            self.exons.update(self[tid].exons)
        return templates

    def _find_communities_boundaries(self, five_graph, three_graph):

        five_found = set()

        __to_modify = dict()

        while len(five_graph) > 0:
            # Find the sinks
            sinks = {node for node in five_graph.nodes() if node not in
                     {edge[0] for edge in five_graph.edges()}}
            __putative = defaultdict(list)
            for sink in sinks:
                for ancestor in nx.ancestors(five_graph, sink):
                    __putative[ancestor].append((sink, self[sink].score))

            for ancestor in __putative:
                best = sorted(__putative[ancestor], key=operator.itemgetter(1), reverse=True)[0][0]
                self.logger.debug("Putative 5' for %s: %s. Best: %s", ancestor, __putative[ancestor], best)
                __to_modify[ancestor] = [self[best], False]
                five_found.add(ancestor)

            five_graph.remove_nodes_from(set.union(sinks, __putative.keys()))

        three_found = set()
        while len(three_graph) > 0:
            sinks = {node for node in three_graph.nodes() if node not in
                     {edge[0] for edge in three_graph.edges()}}
            __putative = defaultdict(list)
            for sink in sinks:
                for ancestor in nx.ancestors(three_graph, sink):
                    __putative[ancestor].append((sink, self[sink].score))

            for ancestor in __putative:
                best = sorted(__putative[ancestor], key=operator.itemgetter(1), reverse=True)[0][0]
                if ancestor in __to_modify:
                    __to_modify[ancestor][1] = self[best]
                else:
                    __to_modify[ancestor] = [False, self[best]]
                three_found.add(ancestor)
                self.logger.debug("Putative 3' for %s: %s. Best: %s", ancestor, __putative[ancestor], best)
            three_graph.remove_nodes_from(set.union(sinks, __putative.keys()))

        self.logger.debug("Communities for modifications: %s", __to_modify)

        return __to_modify

    def define_graph(self, objects: dict, inters=None, three_prime=False):

        graph = nx.DiGraph()
        graph.add_nodes_from(objects.keys())
        if inters is None:
            inters = self._share_extreme

        if len(objects) >= 2:
            if (three_prime is True and self.strand != "-") or (three_prime is False and self.strand == "-"):
                reverse = True
            else:
                reverse = False
            order = sorted([(objects[tid].start, objects[tid].end, tid) for tid in objects], reverse=reverse)

            for pos in range(len(order) - 1):
                obj = order[pos]
                for other_obj in order[pos + 1:]:
                    if obj == other_obj:
                        continue
                    elif self.overlap(obj[:2], obj[:2], positive=False, flank=0) == 0:
                        break
                    else:
                        self.logger.debug("Comparing %s to %s (%s')", obj[2], other_obj[2],
                                          "5" if not three_prime else "3")
                        edge = inters(objects[obj[2]], objects[other_obj[2]], three_prime=three_prime)
                        if edge:
                            assert edge[0].id in self
                            assert edge[1].id in self
                            # assert edge[1].id in self.scores
                            graph.add_edge(edge[0].id, edge[1].id)
        else:

            self.logger.debug("No comparison to be made (objects: %s)", objects)

        return graph

    def _share_extreme(self, first: Transcript, second: Transcript, three_prime=False):

        """
        This function will determine whether two transcripts "overlap" at the 3' or 5' end.
        The things to be considered are:
        - whether their extremes are within the maximum distance
        - whether unifying them would bridge too many splice sites within the locus.
        :param first:
        :param second:
        :return:
        """

        if self.strand == "-":  # Remember we have to invert on the negative strand!
            three_prime = not three_prime

        if three_prime:
            return self._share_three_prime(first, second)
        else:
            return self._share_five_prime(first, second)

    def _share_five_prime(self, first: Transcript, second: Transcript):

        reason = None
        ts_splices = 0
        ts_distance = 0

        if second.start == first.start:
            self.logger.debug("%s and %s start at the same coordinate. No expanding.", first.id, second.id)
            return False

        decision = False
        first, second = sorted([first, second], key=operator.attrgetter("start"))
        if self.overlap((first.start, first.end), (second.start, second.end)) <= 0:
            decision = False
            reason = "{first.id} and {second.id} are not intersecting each other."
            self.logger.debug(reason)
            return decision            

        # Now let us check whether the second falls within an intron
        matched = first.segmenttree.find(second.exons[0][0], second.exons[0][1])
        self.logger.debug("{second.id} last exon {second.exons[0]} intersects in {first.id}: {matched}".format(
            **locals()))
        if len(matched) > 0 and (matched[0].value == "intron" or second.exons[0][0] < matched[0].start):
            decision = False
            reason = "{second.id} first exon ends within an intron of {first.id}".format(**locals())
        elif len(matched) > 0:
            upstream = [_ for _ in first.find_upstream(second.exons[0][0], second.exons[0][1])
                        if _.value == "exon" and _ not in matched]
            if matched[0][0] < second.start:
                if upstream:
                    ts_splices += 1
                ts_distance += second.start - matched[0][0] + 1
            for up in upstream:
                if up.start == first.start:
                    ts_splices += 1
                else:
                    ts_splices += 2
                ts_distance += up.end - up.start - 1

        if reason is None:
            decision = (ts_distance <= self.ts_distance) and (ts_splices <= self.ts_max_splices)
            if decision:
                decision = (second, first)
            reason = "{first.id} {doesit} overlap {second.id} (distance {ts_distance} max {self.ts_distance}, splices \
{ts_splices} max {self.ts_max_splices})".format(doesit="does" if decision else "does not", **locals())
        self.logger.debug(reason)
        return decision

    def _share_three_prime(self, first: Transcript, second: Transcript):

        if second.end == first.end:
            self.logger.debug("%s and %s end at the same coordinate. No expanding.", first.id, second.id)
            return False

        reason = None
        ts_splices = 0
        ts_distance = 0
        decision = False
        first, second = sorted([first, second], key=operator.attrgetter("end"), reverse=False)
        # Now let us check whether the second falls within an intron
        if self.overlap((first.start, first.end), (second.start, second.end)) <= 0:
            decision = False
            reason = "{first.id} and {second.id} are not intersecting each other."
            self.logger.debug(reason)
            return decision            
        
        matched = second.segmenttree.find(first.exons[-1][0], first.exons[-1][1])
          
        if len(matched) > 0 and (matched[-1].value == "intron" or first.exons[-1][1] > matched[-1].end):
            decision = False
            reason = "{first.id} last exon ends within an intron of {second.id}".format(**locals())
        elif len(matched) > 0:
            downstream = [_ for _ in second.find_downstream(first.exons[-1][0], first.exons[-1][1])
                          if _.value == "exon" and _ not in matched]

            if matched[-1][1] > first.end:
                if downstream:
                    ts_splices += 1
                ts_distance += matched[-1][1] - first.end + 1

            for down in downstream:
                if down.end == second.end:
                    ts_splices += 1
                else:
                    ts_splices += 2
                ts_distance += down.end - down.start - 1

        if reason is None:
            decision = (ts_distance <= self.ts_distance) and (ts_splices <= self.ts_max_splices)
            if decision:
                decision = (first, second)
            reason = "{second.id} {doesit} overlap {first.id} (distance {ts_distance} max \
{self.ts_distance}, splices {ts_splices} max {self.ts_max_splices})".format(
                doesit="does" if decision else "does not", **locals())
        self.logger.debug(reason)
        return decision

    @property
    def __name__(self):
        if len(self.transcripts) == 0:
            return "locus"
        elif any(transcript.selected_cds_length > 0 for
                 transcript in self.transcripts.values()):
            return "gene"
        else:
            return "ncRNA_gene"

    # pylint: disable=invalid-name
    @property
    def id(self):
        """
        Override of the abstractlocus method.
        :rtype str
        """
        if self.__id is None:
            myid = Abstractlocus.id.fget(self)  # @UndefinedVariable

            if self.counter > 0:
                myid = "{0}.{1}".format(myid, self.counter)
            # self.__set_id(myid)
            self.__id = myid

        return self.__id

    # pylint: disable=arguments-differ
    @id.setter
    def id(self, string):
        """
        Override of the original method from AbstractLocus. This override allows to
        create proper IDs for the final annotation to be output by Mikado.
        :param string:
        :return:
        """

        self.logger.debug("Setting new ID for %s to %s", self.id, string)
        self.__set_id(string)
        self.__id = string

    def __set_id(self, string):

        if string == self.__id:
            return
        primary_id = "{0}.1".format(string)
        old_primary = self.primary_transcript.id
        self.primary_transcript.attributes["alias"] = self.primary_transcript.id
        self.primary_transcript.id = primary_id
        self.transcripts[primary_id] = self.primary_transcript
        self.primary_transcript_id = primary_id
        del self.transcripts[old_primary]

        order = sorted([k for k in self.transcripts.keys() if k != primary_id],
                       key=lambda xtid: self.transcripts[xtid])

        mapper = {old_primary: primary_id}

        for counter, tid in enumerate(order):
            counter += 2
            self.transcripts[tid].attributes["alias"] = tid
            new_id = "{0}.{1}".format(string, counter)
            self.transcripts[tid].id = new_id
            self.transcripts[new_id] = self.transcripts.pop(tid)
            mapper[tid] = new_id

        if self.scores_calculated is True:
            for tid in mapper:
                self.scores[mapper[tid]] = self.scores.pop(tid)
                self._metrics[mapper[tid]] = self._metrics.pop(tid)
        if self.metrics_calculated is True:
            for index in range(len(self.metric_lines_store)):
                self.metric_lines_store[index]["tid"] = mapper[self.metric_lines_store[index]["tid"]]
                self.metric_lines_store[index]["parent"] = self.id


    # pylint: enable=invalid-name,arguments-differ

    @property
    def is_fragment(self):
        """
        :rtype : bool
        Flag. It returns the value of self.attributes["is_fragment"]

        """
        return self.attributes["is_fragment"]

    @is_fragment.setter
    def is_fragment(self, val: bool):
        """
        Setter for is_fragment. Only boolean values are accepted.

        :param val: flag
        :type val: bool
        """
        if not isinstance(val, bool):
            raise ValueError(val)
        self.attributes["is_fragment"] = val

    @property
    def primary_transcript(self):
        """
        This property returns the primary transcript of the Locus
        (i.e. the one which has been used for creation and which has the highest score).
        :rtype : Transcript
        """

        return self.transcripts[self.primary_transcript_id]

    @property
    def purge(self):

        """Overloading of the base property. Loci should never purge."""

        return False

    @property
    def _finalized(self):
        return self.__finalized

    @_finalized.setter
    def _finalized(self, value):
        if not isinstance(value, bool):
            raise ValueError(value)
        self.__finalized = value

    @property
    def ts_distance(self):
        return self.json_conf["pick"]["alternative_splicing"]["ts_distance"]

    @property
    def ts_max_splices(self):
        return self.json_conf["pick"]["alternative_splicing"]["ts_max_splices"]

    @property
    def has_reference_transcript(self):
        return any(self.transcripts[transcript].is_reference for transcript in self)

    @property
    def _redundant_ccodes(self):
        return self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"]

    def _get_alternative_splicing_codes(self):
        """Method to retrieve the currently valid alternative splicing event codes"""
        return self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"]

    def _add_to_alternative_splicing_codes(self, code):
        """Method to retrieve the currently valid alternative splicing event codes"""
        self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"] = list(set(
            list(self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"]) + [code]
        ))
        _valid_ccodes.validate(self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"])
        self._remove_from_redundant_splicing_codes(code)

    def _add_to_redundant_splicing_codes(self, code):
        """Method to retrieve the currently valid alternative splicing event codes"""
        self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"] = list(set(
            list(self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"]) + [code]
        ))
        _valid_redundant.validate(self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"])
        self._remove_from_alternative_splicing_codes(code)

    def _remove_from_alternative_splicing_codes(self, *ccodes):
        sub = self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"]
        for ccode in ccodes:
            if ccode in sub:
                sub.remove(ccode)
        self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"] = sub

    def _remove_from_redundant_splicing_codes(self, *ccodes):
        self.logger.debug("Removing from redundant ccodes: %s. Current: %s", ccodes,
                          self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"])
        sub = self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"]
        sub = [_ for _ in sub if _ not in ccodes]
        self.logger.debug("New redundant ccodes: %s", sub)
        self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"] = sub

        
def expand_transcript(transcript: Transcript,
                      start_transcript: [Transcript, bool],
                      end_transcript: [Transcript, bool],
                      fai: pysam.libcfaidx.FastaFile,
                      logger):

    """This method will enlarge the coordinates and exon structure of a transcript, given:
    :param transcript: the transcript to modify.
    :type transcript: Transcript
    :param start_transcript: the template transcript for the 5' end.
    :param end_transcript: the template transcript for the 3' end.
    :param fai: the indexed genomic sequence.
    :param logger: the logger to be used in the function.
    """

    # If there is nothing to do, just get out
    transcript.finalize()
    if start_transcript not in (False, None):
        start_transcript.finalize()
    if end_transcript not in (False, None):
        end_transcript.finalize()

    if start_transcript in (False, None) and end_transcript in (False, None):
        logger.debug("%s does not need to be expanded, exiting", transcript.id)
        return transcript

    if transcript.strand == "-":
        start_transcript, end_transcript = end_transcript, start_transcript

    # Make a backup copy of the transcript
    backup = transcript.deepcopy()
    # First get the ORFs
    # Remove the CDS and unfinalize
    logger.debug("Starting expansion of %s", transcript.id)
    strand = transcript.strand
    transcript.strip_cds()
    transcript.unfinalize()
    assert strand == transcript.strand
    downstream = 0
    down_exons = []

    upstream, up_exons, new_first_exon, up_remove = _enlarge_start(transcript, backup, start_transcript)
    downstream, up_exons, down_exons, down_remove = _enlarge_end(transcript,
                                                                 backup, end_transcript, up_exons, new_first_exon)

    first_exon, last_exon = transcript.exons[0], transcript.exons[-1]

    assert upstream >= 0 and downstream >= 0

    if up_remove is True:
        # Remove the first exon
        transcript.remove_exon(first_exon)
    if down_remove is True:
        if not (up_remove is True and first_exon == last_exon):
            transcript.remove_exon(last_exon)

    new_exons = up_exons + down_exons
    if not new_exons:
        logger.debug("%s does not need to be expanded, exiting", transcript.id)
        return backup

    transcript.add_exons(new_exons)
    transcript.start, transcript.end = None, None
    transcript.finalize()

    if transcript.strand == "-":
        downstream, upstream = upstream, downstream

    if (up_exons or down_exons):
        if backup.is_coding:
            seq = check_expanded(transcript, backup, start_transcript, end_transcript,
                                 fai, upstream, downstream, logger)
            transcript = enlarge_orfs(transcript, backup, seq, upstream, downstream, logger)
            transcript.finalize()
    else:
        return backup

    # Now finalize again
    logger.debug("%s: start (before %s, now %s, %s), end (before %s, now %s, %s)",
                 transcript.id,
                 backup.start, transcript.start, transcript.start < backup.start,
                 backup.end, transcript.end, transcript.end > backup.end)
    if transcript.start < backup.start or transcript.end > backup.end:
        transcript.attributes["padded"] = True

    # Now check that we have a valid expansion
    if backup.is_coding and not transcript.is_coding:
        # Something has gone wrong. Just return the original transcript.
        assert new_exons
        logger.info("Padding %s would lead to an invalid CDS (up exons: %s). Aborting.",
                    transcript.id, up_exons)
        return backup
    elif backup.is_coding:
        abort = False
        if backup.strand == "-" and backup.combined_cds_end < transcript.combined_cds_end:
            abort = True
        elif backup.strand != "-" and backup.combined_cds_end > transcript.combined_cds_end:
            abort = True
        if abort is True:
            msg = "Padding {} (strand: {}) would lead to an in-frame stop codon ({} to {}, \
vs original {} to {}. Aborting.".format(
                transcript.id, backup.strand, transcript.combined_cds_start, transcript.combined_cds_end,
                backup.combined_cds_start, backup.combined_cds_end)
            logger.info(msg)
            return backup

    return transcript


def _enlarge_start(transcript: Transcript,
                   backup: Transcript,
                   start_transcript: Transcript) -> (int, list, [None, tuple], bool):

    """This method will enlarge the transcript at the 5' end, using another transcript as the template.
    :param transcript: the original transcript to modify.
    :param backup: a copy of the transcript. As we are modifying the original one, we do need a hard copy.
    :param start_transcript: the template transcript.

    The function returns the following:
    :returns: the upstream modification, the list of upstream exons to add, the new first exon (if any),
              a boolean flag indicating whether the first exon of the transcript should be removed.
    """

    upstream = 0
    up_exons = []
    new_first_exon = None
    to_remove = False
    if start_transcript:
        transcript.start = start_transcript.start
        upstream_exons = sorted([_ for _ in
                                 start_transcript.find_upstream(transcript.exons[0][0], transcript.exons[0][1])
                                      if _.value == "exon"])
        intersecting_upstream = sorted(start_transcript.search(
            transcript.exons[0][0], transcript.exons[0][1]))

        if not intersecting_upstream:
            raise KeyError("No exon or intron found to be intersecting with %s vs %s, this is a mistake",
                           transcript.id, start_transcript.id)

        if intersecting_upstream[0].value == "exon":
            new_first_exon = (min(intersecting_upstream[0][0], backup.start),
                              transcript.exons[0][1])
            if new_first_exon != transcript.exons[0]:
                upstream += backup.start - new_first_exon[0]
                up_exons.append(new_first_exon)
                to_remove = True
            else:
                new_first_exon = None
            if intersecting_upstream[0] in upstream_exons:
                upstream_exons.remove(intersecting_upstream[0])
            upstream += sum(_[1] - _[0] + 1 for _ in upstream_exons)
            up_exons.extend([(_[0], _[1]) for _ in upstream_exons])
        elif intersecting_upstream[0].value == "intron":
            # Check whether the first exon of the model *ends* within an *intron* of the template
            # If that is the case, we have to keep the first exon in place and
            # just expand it until the end
            # Now we have to expand until the first exon in the upstream_exons
            if intersecting_upstream[0][1] == transcript.exons[0][0] - 1:
                assert upstream_exons
                to_remove = False
            elif upstream_exons:
                to_remove = True
                upstream_exon = upstream_exons[-1]
                new_first_exon = (upstream_exon[0], transcript.exons[0][1])
                upstream_exons.remove(upstream_exon)
                upstream += backup.start - new_first_exon[0]
                up_exons.append(new_first_exon)
            else:
                # Something fishy going on here. Let us double check everything.
                if start_transcript.exons[0][0] == transcript.start:
                    raise ValueError(
                        "Something has gone wrong. The template transcript should have returned upstream exons."
                    )
                elif start_transcript.exons[0][0] < transcript.start:
                    raise ValueError(
                        "Something has gone wrong. We should have found the correct exons."
                    )
                else:
                    pass

            upstream += sum(_[1] - _[0] + 1 for _ in upstream_exons)
            up_exons.extend([(_[0], _[1]) for _ in upstream_exons])

    return upstream, up_exons, new_first_exon, to_remove


def _enlarge_end(transcript: Transcript,
                 backup: Transcript,
                 end_transcript: Transcript,
                 up_exons: list,
                 new_first_exon: [None, tuple]) -> [int, list, list, bool]:

    """
    This method will enlarge the transcript at the 5' end, using another transcript as the template.
    :param transcript: the original transcript to modify.
    :param backup: a copy of the transcript. As we are modifying the original one, we do need a hard copy.
    :param end_transcript: the template transcript.
    :param up_exons: the list of exons added at the 5' end.
    :param new_first_exon: the new coordinates of what used to be the first exon of the transcript.
                           This is necessary because if the transcript is monoexonic, we might need to re-modify it.

    The function returns the following:
    :returns: the downstream modification, the (potentially modified) list of upstream exons to add,
              the list of downstream exons to add, a boolean flag indicating whether the last exon of the transcript
              should be removed.
    """

    downstream = 0
    down_exons = []
    to_remove = False

    if end_transcript:
        transcript.end = end_transcript.end
        downstream_exons = sorted([_ for _ in
                                    end_transcript.find_downstream(transcript.exons[-1][0], transcript.exons[-1][1])
                                    if _.value == "exon"])
        intersecting_downstream = sorted(end_transcript.search(
            transcript.exons[-1][0], transcript.exons[-1][1]))
        if not intersecting_downstream:
            raise KeyError("No exon or intron found to be intersecting with %s vs %s, this is a mistake",
                           transcript.id, end_transcript.id)
        # We are taking the right-most intersecting element.
        if intersecting_downstream[-1].value == "exon":
            if transcript.monoexonic and new_first_exon is not None:
                new_exon = (new_first_exon[0], max(intersecting_downstream[-1][1], new_first_exon[1]))
                if new_exon != new_first_exon:
                    up_exons.remove(new_first_exon)
                    downstream += new_exon[1] - backup.end
                    down_exons.append(new_exon)
                    to_remove = True
            else:
                new_exon = (transcript.exons[-1][0],
                            max(intersecting_downstream[-1][1], transcript.exons[-1][1]))
                if new_exon != transcript.exons[-1]:
                    downstream += new_exon[1] - backup.end
                    down_exons.append(new_exon)
                    to_remove = True

            if intersecting_downstream[-1] in downstream_exons:
                downstream_exons.remove(intersecting_downstream[-1])
            downstream += sum(_[1] - _[0] + 1 for _ in downstream_exons)
            down_exons.extend([(_[0], _[1]) for _ in downstream_exons])
        elif intersecting_downstream[-1].value == "intron":
            # Now we have to expand until the first exon in the upstream_exons
            if intersecting_downstream[-1][0] == transcript.exons[-1][1] + 1:
                assert downstream_exons
                to_remove = False
            elif downstream_exons:
                downstream_exon = downstream_exons[0]
                assert downstream_exon[1] > backup.end
                assert downstream_exon[0] > backup.end
                if transcript.monoexonic and new_first_exon is not None:
                    new_exon = (new_first_exon[0], downstream_exon[1])
                    up_exons.remove(new_first_exon)
                    to_remove = True
                else:
                    new_exon = (transcript.exons[-1][0], downstream_exon[1])
                    to_remove = True
                downstream_exons.remove(downstream_exon)
                downstream += new_exon[1] - backup.end
                down_exons.append(new_exon)
            else:
                # Something fishy going on here. Let us double check everything.
                if end_transcript.exons[-1][1] == transcript.end:
                    raise ValueError(
                        "Something has gone wrong. The template transcript should have returned upstream exons."
                    )
                elif end_transcript.exons[-1][1] > transcript.end:
                    raise ValueError(
                        "Something has gone wrong. We should have found the correct exons."
                    )
            downstream += sum(_[1] - _[0] + 1 for _ in downstream_exons)
            down_exons.extend([(_[0], _[1]) for _ in downstream_exons])

    return downstream, up_exons, down_exons, to_remove


def check_expanded(transcript, backup, start_transcript, end_transcript, fai, upstream, downstream, logger) -> str:

    """
    This function checks that the expanded transcript is valid, and it also calculates and returns its cDNA sequence.
    :param transcript: the modified transcript.
    :param backup: The original transcript, before expansion.
    :param start_transcript: the transcript used as template at the 5' end.
    :param end_transcript: the transcript used as template at the 3' end.
    :param fai: The pysam.libcfaidx.FastaFile object indexing the genome.
    :param upstream: the amount of transcriptomic base-pairs added to the transcript at its 5' end.
    :param downstream: the amount of transcriptomic base-pairs added to the transcript at its 3' end.
    :param logger: the logger to use.
    :returns: the cDNA of the modified transcript, as a standard Python string.
    """

    assert transcript.exons != backup.exons
    assert transcript.end <= fai.get_reference_length(transcript.chrom), (
        transcript.end, fai.get_reference_length(transcript.chrom))
    genome_seq = fai.fetch(transcript.chrom, transcript.start - 1, transcript.end)

    if not (transcript.exons[-1][1] - transcript.start + 1 == len(genome_seq)):
        error = "{} should have a sequence of length {} ({} start, {} end), but one of length {} has been given"
        error = error.format(transcript.id, transcript.exons[-1][1] - transcript.start + 1,
                             transcript.start, transcript.end, len(genome_seq))
        logger.error(error)
        raise InvalidTranscript(error)
    seq = TranscriptChecker(transcript, genome_seq, is_reference=True).cdna
    assert len(seq) == transcript.cdna_length, (len(seq), transcript.cdna_length, transcript.exons)
    if not len(seq) == backup.cdna_length + upstream + downstream:
        error = [len(seq), backup.cdna_length + upstream + downstream,
                 backup.cdna_length, upstream, downstream,
                 (transcript.start, transcript.end), (backup.id, backup.start, backup.end),
                 (None if not start_transcript else (start_transcript.id, (start_transcript.start,
                                                                           start_transcript.end))),
                 (None if not end_transcript else (end_transcript.id, (end_transcript.start,
                                                                       end_transcript.end))),
                 (backup.id, backup.exons),
                 None if not start_transcript else (start_transcript.id, start_transcript.exons),
                 None if not end_transcript else (end_transcript.id, end_transcript.exons),
                 (transcript.id + "_expanded", transcript.exons),
                 set.difference(set(transcript.exons), set(backup.exons)),
                 set.difference(set(backup.exons), set(transcript.exons))
                 ]
        error = "\n".join([str(_) for _ in error])
        raise AssertionError(error)
    return seq


def enlarge_orfs(transcript: Transcript,
                 backup: Transcript,
                 seq: str,
                 upstream: int,
                 downstream: int,
                 logger) -> Transcript:

    """
    This method will take an expanded transcript and recalculate its ORF(s). As a consequence of the expansion,
    truncated transcripts might become whole.
    :param transcript: the expanded transcript.
    :param backup: the original transcript. Used to extract the original ORF(s).
    :param seq: the new cDNA sequence of the expanded transcript.
    :param upstream: the amount of expansion that happened at the 5'.
    :param downstream: the amount of expansion that happened at the 3'.
    :param logger: the logger.
    :returns: the modified transcript with the ORF(s) recalculated.
    """

    if backup.combined_cds_length > 0:
        try:
            internal_orfs = list(backup.get_internal_orf_beds())
        except (ValueError, TypeError, AssertionError):
            logger.error("Something went wrong with the CDS extraction for %s. Stripping it.",
                         backup.id)
            internal_orfs = []
    else:
        internal_orfs = []

    if not internal_orfs:
        return transcript

    new_orfs = []
    for orf in internal_orfs:
        logger.debug("Old ORF: %s", str(orf))
        try:
            logger.debug("Sequence for %s: %s[..]%s (upstream %s, downstream %s)",
                         transcript.id, seq[:10], seq[-10:], upstream, downstream)
            orf.expand(seq, upstream, downstream, expand_orf=True, logger=logger)
        except AssertionError as err:
            logger.error(err)
            logger.error("%s, %s, %s, %s",
                         upstream,
                         downstream,
                         transcript.exons,
                         transcript.cdna_length)
            raise AssertionError(err)
        logger.debug("New ORF: %s", str(orf))
        if orf.coding is False:
            raise ValueError(orf)
        elif orf.invalid:
            raise InvalidTranscript(orf.invalid_reason)

        new_orfs.append(orf)

    transcript.load_orfs(new_orfs)
    transcript.finalize()
    if backup.is_coding and not transcript.is_coding:
        raise InvalidTranscript(new_orfs)
    return transcript
