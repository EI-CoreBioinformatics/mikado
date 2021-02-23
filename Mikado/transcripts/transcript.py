from .._transcripts.transcript_base import TranscriptBase, Metric
from ..configuration.configuration import MikadoConfiguration
from ..configuration.daijin_configuration import DaijinConfiguration
from .transcript_methods import splitting, retrieval
from typing import List
from ..serializers.blast_serializer import Hit
from sqlalchemy import and_
from sqlalchemy import bindparam
from sqlalchemy.ext import baked
from sqlalchemy.sql.expression import desc, asc  # SQLAlchemy imports
import pysam
import functools
import inspect
default_config = MikadoConfiguration()


class Transcript(TranscriptBase):
    # Query baking to minimize overhead
    bakery = baked.bakery()
    blast_baked = bakery(lambda session: session.query(Hit))
    blast_baked += lambda q: q.filter(and_(Hit.query == bindparam("query"),
                                           Hit.evalue <= bindparam("evalue")), )

    blast_baked += lambda q: q.order_by(asc(Hit.evalue))
    # blast_baked += lambda q: q.limit(bindparam("max_target_seqs"))

    def __init__(self, *args, configuration=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.__configuration = None
        self.configuration = configuration

    @property
    def configuration(self):
        """
        Configuration dictionary. It can be None.
        :return:
        """
        if self.__configuration is None:
            self.__configuration = default_config.copy()

        return self.__configuration

    @configuration.setter
    def configuration(self, configuration):

        """
        Setter for the configuration dictionary.
        :param configuration: None or a dictionary
        :type configuration: (None | MikadoConfiguration | DaijinConfiguration)
        :return:
        """

        if configuration is None:
            configuration = default_config.copy()

        assert isinstance(configuration, (MikadoConfiguration, DaijinConfiguration))
        self.__configuration = configuration

    def __getstate__(self):
        state = super().__getstate__()
        if hasattr(self, "configuration") and self.configuration is not None:
            state["configuration"] = self.configuration.copy()
            assert isinstance(state["configuration"], (MikadoConfiguration, DaijinConfiguration)), type(
                self.configuration)
            if isinstance(state["configuration"].reference.genome, pysam.FastaFile):
                state["configuration"]["reference"]["genome"] = state["configuration"].reference.genome.filename
        return state

    def __setstate__(self, state):
        self.configuration = state.pop("configuration", None)
        self.__dict__.update(state)
        self._calculate_cds_tree()
        self._calculate_segment_tree()
        self.logger = None

    def split_by_cds(self) -> List:
        """This method is used for transcripts that have multiple ORFs.
        It will split them according to the CDS information into multiple transcripts.
        UTR information will be retained only if no ORF is down/upstream.
        """

        for new_transcript in splitting.split_by_cds(self):
            yield new_transcript

        return

    def load_information_from_db(self, configuration, introns=None, data_dict=None):
        """This method will load information regarding the transcript from the provided database.

        :param configuration: Necessary configuration file
        :type configuration: (MikadoConfiguration|DaijinConfiguration)

        :param introns: the verified introns in the Locus
        :type introns: None,set

        :param data_dict: a dictionary containing the information directly
        :type data_dict: dict

        Verified introns can be provided from outside using the keyword.
        Otherwise, they will be extracted from the database directly.
        """

        retrieval.load_information_from_db(self,
                                           configuration,
                                           introns=introns,
                                           data_dict=data_dict)

    def load_orfs(self, candidate_orfs):

        """
        Thin layer over the load_orfs method from the retrieval module.
        :param candidate_orfs: list of candidate ORFs in BED12 format.
        :return:
        """

        retrieval.load_orfs(self, candidate_orfs)

    def find_overlapping_cds(self, candidate_orfs):

        """
        Thin wrapper for the homonym function in retrieval
        :param candidate_orfs: List of candidate ORFs
        :return:
        """

        return retrieval.find_overlapping_cds(self, candidate_orfs)

    # We need to overload this because otherwise we won't get the metrics from the base class.
    @classmethod
    @functools.lru_cache(maxsize=None, typed=True)
    def get_available_metrics(cls) -> list:
        """This function retrieves all metrics available for the class."""

        metrics = TranscriptBase.get_available_metrics()
        for member in inspect.getmembers(cls):
            if not member[0].startswith("__") and member[0] in cls.__dict__ and isinstance(
                    cls.__dict__[member[0]], Metric):
                metrics.append(member[0])

        _metrics = sorted(set([metric for metric in metrics]))
        final_metrics = ["tid", "alias", "parent", "original_source", "score"] + _metrics
        return final_metrics

    # We need to overload this because otherwise we won't get the metrics from the base class.
    @classmethod
    @functools.lru_cache(maxsize=None, typed=True)
    def get_modifiable_metrics(cls) -> set:

        metrics = TranscriptBase.get_modifiable_metrics()
        for member in inspect.getmembers(cls):
            not_private = (not member[0].startswith("_" + cls.__name__ + "__") and not member[0].startswith("__"))
            in_dict = (member[0] in cls.__dict__)
            if in_dict:
                is_metric = isinstance(cls.__dict__[member[0]], Metric)
                has_fset = (getattr(cls.__dict__[member[0]], "fset", None) is not None)
            else:
                is_metric = None
                has_fset = None
            if all([not_private, in_dict, is_metric, has_fset]):
                metrics.append(member[0])
        return set(metrics)
