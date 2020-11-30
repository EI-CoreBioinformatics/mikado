"""
This module contains the methods used by the Transcript class to retrieve information
from the database/dictionary provided during the pick operation.
"""

import operator
from itertools import groupby

from sqlalchemy import and_
from sqlalchemy.orm.session import sessionmaker

from ...serializers.junction import Junction
from ..clique_methods import define_graph, find_cliques, find_communities
from ...utilities import dbutils

__author__ = 'Luca Venturini'


def load_orfs(transcript, candidate_orfs):

    """
    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param candidate_orfs: The ORFs to be inspected for loading.
    :type candidate_orfs: list[Mikado.serializers.orf.Orf|Mikado.parsers.bed12.BED12]

    This method replicates what is done internally by the
    "cdna_alignment_orf_to_genome_orf.pl"
    utility in the TransDecoder suite. It takes as argument "candidate_orfs"
    i.e. a list of BED12 serialised objects.
    The method expects as argument a dictionary containing BED entries,
    and indexed by the transcript name. The indexed name *must* equal the
    "id" property, otherwise the method returns immediately.
    If no entry is found for the transcript, the method exits immediately.
    Otherwise, any CDS information present in the original GFF/GTF
    file is completely overwritten.
    Briefly, it follows this logic:
    - Finalise the transcript
    - Retrieve from the dictionary (input) the CDS object
    - Sort in decreasing order the CDSs on the basis of:
        - Presence of start/stop codon
        - CDS length (useful for monoexonic transcripts where we might want to set the strand)
    - For each CDS:
        - If the ORF is on the + strand:
            - all good
        - If the ORF is on the - strand:
            - if the transcript is monoexonic: invert its genomic strand
            - if the transcript is multiexonic: skip
        - Start looking at the exons
    """

    # Prepare the transcript
    transcript.finalize()

    transcript.logger.debug("Finding the candidate ORFS out of %d .. ",
                            len(candidate_orfs))
    # This will also exclude invalid ORFs
    for orf in candidate_orfs:
        orf.logger = transcript.logger

    candidate_orfs = find_overlapping_cds(transcript, candidate_orfs)
    transcript.logger.debug("Retained %d candidate ORFS",
                            len(candidate_orfs))

    if candidate_orfs is None or len(candidate_orfs) == 0:
        transcript.logger.debug("No ORF for %s; %s", transcript.id, candidate_orfs)
        return

    transcript.combined_utr = []
    transcript.combined_cds = []
    transcript.internal_orfs = []
    transcript.finalized = False
    # Token to be set to False after the first CDS is exhausted
    primary_orf = True
    primary_strand = None
    # This will keep in memory the original BED12 objects
    transcript.loaded_bed12 = []

    primary_phase = None
    for orf in candidate_orfs:
        # Minimal check
        transcript.logger.debug("ORF for %s: start %s, end %s (%s), phase %s",
                                transcript.id, orf.thick_start, orf.thick_end, orf.end, orf.phase)
        if primary_orf is False and orf.strand != primary_strand:
            continue

        check_sanity = (orf.thick_start >= 1 and orf.thick_end <= transcript.cdna_length)
        if len(orf) != transcript.cdna_length or not check_sanity:
            message = "Wrong ORF for {0}: ".format(orf.id)
            message += "cDNA length: {0}; ".format(transcript.cdna_length)
            message += "orf length: {0}; ".format(len(orf))
            message += "CDS: {0}-{1}".format(orf.thick_start, orf.thick_end)
            transcript.logger.warning(message)
            continue

        if transcript.strand is None:
            transcript.strand = orf.strand

        transcript.loaded_bed12.append(orf)
        cds_exons = __create_internal_orf(transcript, orf)
        cds_exons = sorted(cds_exons,
                           key=operator.itemgetter(1))

        transcript.internal_orfs.append(cds_exons)

        if primary_orf is True:
            transcript.logger.debug("%s has start codon: %s, %s",
                                    transcript.id,
                                    orf.has_start_codon,
                                    orf.has_stop_codon)
            (transcript.has_start_codon, transcript.has_stop_codon) = (orf.has_start_codon,
                                                                       orf.has_stop_codon)
            primary_phase = orf.phase
            primary_orf = False
            primary_strand = orf.strand
            transcript.selected_internal_orf_index = transcript.internal_orfs.index(cds_exons)
        # transcript.phases.append(cds_exons[0][0], 0)

    # Now verify the loaded content
    check_loaded_orfs(transcript, primary_phase=primary_phase)


def check_loaded_orfs(transcript, primary_phase=0):

    """
    This function verifies the ORF status after
    loading from an external data structure/database.

    :param transcript: the transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param primary_phase: the starting phase of the primary ORF.
    :type primary_phase: int

    :return:
    """

    if len(transcript.internal_orf_lengths) == 0:
        transcript.logger.warning("No candidate ORF retained for %s",
                                  transcript.id)

    if len(transcript.internal_orfs) == 1:
        transcript.logger.debug("Found 1 ORF for %s", transcript.id)
        transcript.combined_cds = sorted(
            [a[1] for a in iter(_ for _ in transcript.internal_orfs[0]
                                if _[0] == "CDS")],
            key=operator.itemgetter(0, 1))
        transcript.combined_utr = sorted(
            [a[1] for a in iter(_ for _ in transcript.internal_orfs[0]
                                if _[0] == "UTR")],
            key=operator.itemgetter(0, 1)
        )

    elif len(transcript.internal_orfs) > 1:
        transcript.logger.debug("Found %d ORFs for %s",
                                len(transcript.internal_orfs),
                                transcript.id)
        cds_spans = []
        candidates = []
        for internal_cds in transcript.internal_orfs:
            candidates.extend(
                [a[1] for a in iter(tup for tup in internal_cds
                                    if tup[0] == "CDS")])

        for comm in transcript.find_communities(candidates):
            span = tuple([min(t[0] for t in comm), max(t[1] for t in comm)])
            cds_spans.append(span)

        transcript.combined_cds = sorted(cds_spans, key=operator.itemgetter(0, 1))

        # This method is probably OBSCENELY inefficient,
        # but I cannot think of a better one for the moment.

        # curr_utr_segment = None  # temporary token with the UTR positions

        # Calculate the positions inside the UTR
        utr_pos = sorted(list(set.difference(
            set.union(*[set(range(exon[0], exon[1] + 1)) for exon in transcript.exons]),
            set.union(*[set(range(cds[0], cds[1] + 1)) for cds in transcript.combined_cds])
        )))

        # This code snippet using groupby is massively more efficient
        # than the previous naive implementation
        getter = operator.itemgetter(1)
        for _, group in groupby(enumerate(utr_pos), lambda items: items[0] - items[1]):
            group = [getter(element) for element in group]
            # group = list(map(operator.itemgetter(1), group))
            if len(group) > 1:
                transcript.combined_utr.append(tuple([group[0], group[-1]]))
            else:
                transcript.combined_utr.append(tuple([group[0], group[0]]))

        # Check everything is alright
        equality = (transcript.cdna_length ==
                    transcript.combined_cds_length + transcript.combined_utr_length)

        assert equality, (transcript.cdna_length, transcript.combined_cds, transcript.combined_utr)

    if transcript.internal_orfs:
        transcript.feature = "mRNA"

    transcript.phases = dict()
    transcript._first_phase = primary_phase
    transcript._trust_orf = False
    transcript.finalize()


def __load_blast(transcript, data_dict=None, reverse=False):

    """This method looks into the DB for hits corresponding to the desired requirements.
    Hits will be loaded into the "blast_hits" list;
    we will not store the SQLAlchemy query object,
    but rather its representation as a dictionary
    (using the Hit.as_dict() method).

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript
    """

    # if self.query_id is None:
    #     return

    # if self.json_conf["pick"]["chimera_split"]["blast_check"] is False:
    #     return

    max_target_seqs = transcript.json_conf[
        "pick"]["chimera_split"]["blast_params"]["max_target_seqs"]
    maximum_evalue = transcript.json_conf["pick"]["chimera_split"]["blast_params"]["evalue"]

    if data_dict is None:
        blast_hits_query = [_.as_dict() for _ in transcript.blast_baked(transcript.session).params(
            query=transcript.id,
            evalue=maximum_evalue)]
    else:
        blast_hits_query = data_dict.get("hits", dict()).get(transcript.id, [])

    transcript.logger.debug("Starting to load BLAST data for %s",
                              transcript.id)
    # variables which take care of defining
    # the maximum evalue and target seqs
    # We do not trust the limit in the sqlite because
    # we might lose legitimate hits with the same evalue as the
    # best ones. So we collapse all hits with the same evalue.
    previous_evalue = -1
    counter = 0
    for hit in blast_hits_query:

        if counter > max_target_seqs and previous_evalue < hit["evalue"]:
            break
        elif previous_evalue < hit["evalue"]:
            previous_evalue = hit["evalue"]

        query_frames = [_["query_frame"] for _ in hit["hsps"]]
        transcript.logger.debug("Query frames for %s: %s", transcript.id, query_frames)
        if reverse is True:
            query_frames = [_ * -1 for _ in query_frames]
            transcript.logger.debug("Query frames for %s after reversal: %s", transcript.id, query_frames)

        if any(_ < 0 for _ in query_frames):
            transcript.logger.debug("Hit %s skipped for %s as it is on opposite strand",
                                    hit["target"], transcript.id)
            continue

        counter += 1

        transcript.blast_hits.append(hit)

    transcript.logger.debug("Loaded %d BLAST hits for %s",
                            counter, transcript.id)


def _connect_to_db(transcript):

    """This method will connect to the database using the information
    contained in the JSON configuration.
    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    """

    transcript.engine = dbutils.connect(
        transcript.json_conf, transcript.logger)

    transcript.sessionmaker = sessionmaker()
    transcript.sessionmaker.configure(bind=transcript.engine)
    transcript.session = transcript.sessionmaker()


def load_information_from_db(transcript, json_conf, introns=None, session=None,
                             data_dict=None):
    """This method will invoke the check for:

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param json_conf: Necessary configuration file
    :type json_conf: dict

    :param introns: the verified introns in the Locus
    :type introns: None,set

    :param session: an SQLAlchemy session
    :type session: sqlalchemy.orm.session

    :param data_dict: a dictionary containing the information directly
    :type data_dict: dict

    Verified introns can be provided from outside using the keyword.
    Otherwise, they will be extracted from the database directly.
    """

    transcript.logger.debug("Loading {0}".format(transcript.id))
    transcript.json_conf = json_conf

    __load_verified_introns(transcript, data_dict, introns)
    if data_dict is not None:
        retrieve_from_dict(transcript, data_dict)
    else:
        if session is None:
            _connect_to_db(transcript)
        else:
            transcript.session = session
        candidate_orfs = []
        if transcript.is_reference is False:
            transcript.logger.debug("Retrieving the ORFs for %s", transcript.id)
            ext_results = transcript.external_baked(transcript.session).params(query=transcript.id).all()
            for row in ext_results:
                if row.rtype == "int":
                    score = int(row.score)
                elif row.rtype == "bool":
                    score = bool(int(row.score))
                elif row.rtype == "complex":
                    score = complex(row.score)
                elif row.rtype == "float":
                    score = float(row.score)
                else:
                    raise ValueError("Invalid rtype: {}".format(row.rtype))

                transcript.external_scores[row.source] = (score, row.valid_raw)

            for row in transcript.external_sources(transcript.session).all():
                if row.source not in transcript.external_scores:
                    if row.rtype in ("int", "complex", "float"):
                        score = 0
                    elif row.rtype == "bool":
                        score = False
                    else:
                        raise ValueError("Invalid rtype: {}".format(row.rtype))
                    transcript.external_scores[row.source] = (score, row.valid_raw)

            for orf in retrieve_orfs(transcript):
                candidate_orfs.append(orf)
            transcript.logger.debug("Retrieved the ORFs for %s", transcript.id)

            transcript.logger.debug("Loading the ORFs for %s", transcript.id)
            old_strand = transcript.strand
            load_orfs(transcript, candidate_orfs)

            if transcript.monoexonic is False:
                is_reversed = False
            elif old_strand != transcript.strand and ((old_strand is None and transcript.strand == "-")
                or (old_strand is not None)):
                is_reversed = True
            else:
                is_reversed = False
            transcript.logger.debug("Loaded the ORFs for %s", transcript.id)

        else:
            transcript.logger.debug("Skipping ORF loading for reference %s", transcript.id)
            is_reversed = False

        transcript.logger.debug("Loading the BLAST data for %s", transcript.id)
        __load_blast(transcript, reverse=is_reversed)
        transcript.logger.debug("Loaded the BLAST data for %s", transcript.id)
    # Finally load introns, separately
    transcript.logger.debug("Loaded data for %s", transcript.id)


def retrieve_from_dict(transcript, data_dict):
    """
    Method to retrieve transcript data directly from a dictionary.

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param data_dict: the dictionary with loaded data from DB
    :type data_dict: (None | dict)
    """

    transcript.logger.debug(
        "Retrieving information from DB dictionary for %s",
        transcript.id)
    # Intron data
    transcript.logger.debug("Checking introns for %s, candidates %s",
                            transcript.id,
                            sorted(transcript.introns))

    # ORF data
    trust_strand = transcript.json_conf["pick"]["orf_loading"]["strand_specific"]
    min_cds_len = transcript.json_conf["pick"]["orf_loading"]["minimal_orf_length"]

    transcript.logger.debug("Retrieving ORF information from DB dictionary for %s",
                            transcript.id)

    if transcript.id in data_dict.get("external", dict()):
        ext_score = data_dict["external"][transcript.id]
        transcript.external_scores.update(ext_score)
        assert set(transcript.external_scores.keys()).issubset(set(ext_score.keys()))
    elif data_dict.get("external", dict()):
        raise KeyError("{} not found in external scores".format(transcript.id))

    if transcript.is_reference is False:
        if transcript.id in data_dict.get("orfs", dict()):
            candidate_orfs = list(orf for orf in data_dict["orfs"][transcript.id] if
                                  orf.cds_len >= min_cds_len)
        else:
            candidate_orfs = []

        # They must already be as ORFs
        if (transcript.monoexonic is False) or (transcript.monoexonic is True and trust_strand is True and
                                                transcript.strand is not None):
            # Remove negative strand ORFs for multiexonic transcripts,
            # or monoexonic strand-specific transcripts
            candidate_orfs = list(orf for orf in candidate_orfs if orf.strand != "-")

        old_strand = transcript.strand
        load_orfs(transcript, candidate_orfs)

        if transcript.monoexonic is False:
            is_reversed = False
        elif old_strand != transcript.strand and ((old_strand is None and transcript.strand == "-")
                                                  or (old_strand is not None)):
            is_reversed = True
        else:
            is_reversed = False
    else:
        transcript.logger.debug("Skipping ORF loading for reference %s", transcript.id)
        is_reversed = False

    # if transcript.json_conf["pick"]["chimera_split"]["blast_check"] is True:
    transcript.logger.debug("Retrieving BLAST hits for %s",
                            transcript.id)

    __load_blast(transcript, data_dict=data_dict, reverse=is_reversed)

    transcript.logger.debug("Retrieved information from DB dictionary for %s",
                            transcript.id)


def find_overlapping_cds(transcript, candidates: list) -> list:
    """
    Wrapper for the Abstractlocus method, used for finding overlapping ORFs.
    It will pass to the function the class's "is_overlapping_cds" method
    (which would be otherwise be inaccessible from the Abstractlocus class method).
    As we are interested only in the communities, not the cliques,
    this wrapper discards the cliques
    (first element of the Abstractlocus.find_communities results)

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param candidates: candidate ORFs to analyse
    :type candidates: list[Mikado.serializers.orf.Orf]

    """

    # If we are looking at a multiexonic transcript
    if not (transcript.monoexonic is True and transcript.strand is None):
        candidates = list(corf for corf in candidates if corf.strand == "+")

    # Prepare the minimal secondary length parameter
    if transcript.json_conf is not None:
        minimal_secondary_orf_length = \
            transcript.json_conf["pick"]["orf_loading"]["minimal_secondary_orf_length"]
    else:
        minimal_secondary_orf_length = 0
    transcript.logger.debug("Minimal orf loading: %d", minimal_secondary_orf_length)

    transcript.logger.debug("{0} input ORFs for {1}".format(len(candidates), transcript.id))
    if any(corf.transcriptomic is False for corf in candidates):
        transcript.logger.debug("%d non-transcriptomic ORFs in the candidates",
                                len([corf.transcriptomic is False for corf in candidates]))
    if any(corf.invalid is True for corf in candidates):
        for corf in candidates:
            if not corf.invalid is True:
                continue
            else:
                transcript.logger.debug("Invalid ORF, reason: %s", corf.invalid_reason)

    candidates = list(corf for corf in candidates if (
        corf.invalid is False and corf.transcriptomic is True))

    ids = set(_.name for _ in candidates)
    if len(ids) < len(candidates):
        transcript.logger.debug("Colliding IDs found for the ORFs. Redefining them.")
        for pos in range(1, len(candidates) + 1):
            candidates[pos - 1].name = "{transcript.id}.orf{pos}".format(**locals())

    transcript.logger.debug("{0} filtered ORFs for {1}".format(len(candidates), transcript.id))
    if len(candidates) == 0:
        return []

    orf_dictionary = dict((x.name, x) for x in candidates)

    # First define the graph
    graph = define_graph(orf_dictionary, inters=transcript.is_overlapping_cds)
    candidate_orfs = find_candidate_orfs(transcript, graph, orf_dictionary)

    transcript.logger.debug("{0} candidate retained ORFs for {1}: {2}".format(
        len(candidate_orfs),
        transcript.id,
        [x.name for x in candidate_orfs]))
    final_orfs = [candidate_orfs[0]]
    if len(candidate_orfs) > 1:
        others = list(corf for corf in candidate_orfs[1:] if
                      corf.cds_len >= minimal_secondary_orf_length)
        transcript.logger.debug("Found {0} secondary ORFs for {1} of length >= {2}".format(
            len(others), transcript.id,
            minimal_secondary_orf_length
        ))
        final_orfs.extend(others)

    transcript.logger.debug("Retained %d ORFs for %s: %s",
                            len(final_orfs),
                            transcript.id,
                            [orf.name for orf in final_orfs])
    return final_orfs


def __create_internal_orf(transcript, orf):

    """
    Private method that calculates the assignment of the exons given the
    coordinates of the transcriptomic ORF.

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param orf: candidate ORF to transform into an internal ORF
    :type orf: Mikado.serializers.orf.Orf

    """

    cds_exons = []
    current_start, current_end = 0, 0

    transcript.logger.debug("Initial phase for %s: %s", transcript.id, orf.phase)
    if orf.strand == "-":
        phase = orf.phase
        # We might decide to remove this check
        assert transcript.monoexonic is True
        current_end = transcript.start + (orf.thick_start - 1)
        current_start = current_end + orf.cds_len - 1

        assert transcript.start <= current_end <= transcript.end
        assert transcript.start <= current_end < current_start <= transcript.end, (
            transcript.start, current_end, current_start, transcript.end
        )
        # This is not true in the case of truncated ORFs!
        # assert (current_start - current_end + 1) % 3 == 0, (
        #     current_end, current_start,
        #     current_start - current_end + 1
        # )
        cds_exons.append(("exon", transcript.exons[0]))
        if current_end > transcript.start:
            cds_exons.append(("UTR", tuple([transcript.start, current_end - 1])))
        cds_exons.append(("CDS", tuple([current_end, current_start]), phase))
        if current_start < transcript.end:
            cds_exons.append(("UTR", tuple([current_start + 1, transcript.end])))
        transcript.strand = "-"
    else:
        previous = -orf.phase
        for exon in sorted(transcript.exons, key=operator.itemgetter(0, 1),
                           reverse=(transcript.strand == "-")):
            cds_exons.append(("exon", tuple([exon[0], exon[1]])))
            current_start += 1
            current_end += exon[1] - exon[0] + 1
            # Whole UTR
            if current_end < orf.thick_start or current_start > orf.thick_end:
                cds_exons.append(("UTR", tuple([exon[0], exon[1]])))
            else:
                if transcript.strand == "+":
                    c_start = exon[0] + max(0, orf.thick_start - current_start)
                    c_end = exon[1] - max(0, current_end - orf.thick_end)
                else:
                    c_start = exon[0] + max(0, current_end - orf.thick_end)
                    c_end = exon[1] - max(0, orf.thick_start - current_start)

                if c_start > exon[0]:
                    u_end = c_start - 1
                    cds_exons.append(("UTR", tuple([exon[0], u_end])))
                if c_start <= c_end:
                    phase = (3 - (previous % 3)) % 3
                    previous += c_end - c_start + 1
                    cds_exons.append(("CDS", tuple([c_start, c_end]), phase))
                if c_end < exon[1]:
                    cds_exons.append(("UTR", tuple([c_end + 1, exon[1]])))
            current_start = current_end
        if orf.phase != 0:
            transcript.logger.debug("Non-0 phase (%d) for %s [orf: %s]",
                                    orf.phase,
                                    transcript.id,
                                    cds_exons)

    return cds_exons


def __load_verified_introns(transcript, data_dict=None, introns=None):

    """This method will load verified junctions from the external
    (usually the superlocus class).

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param data_dict: the dictionary with data to load
    :type data_dict: (dict | None)

    :param introns: verified introns
    :type introns: (set | None)
    """

    transcript.logger.debug("Checking introns; candidates %s", transcript.introns)
    if data_dict is None:
        transcript.logger.debug("Checking introns using the database for %s",
                                transcript.id)

        for intron in transcript.introns:
            # Disable checks as the hybridproperties confuse
            # both pycharm and pylint
            # noinspection PyCallByClass,PyTypeChecker
            # pylint: disable=no-value-for-parameter
            # chrom_id = self.session.query(Chrom).filter(Chrom.name == self.chrom).one().chrom_id
            # import sqlalchemy
            for ver_intron in transcript.session.query(Junction).filter(and_(
                    Junction.chrom == transcript.chrom,
                    intron[0] == Junction.junction_start,
                    intron[1] == Junction.junction_end)):
                if ver_intron.strand in (transcript.strand, None):
                    transcript.logger.debug("Verified intron %s%s:%d-%d for %s",
                                            transcript.chrom, ver_intron.strand,
                                            intron[0], intron[1], transcript.id)
                    transcript.verified_introns.add(intron)

    else:
        transcript.logger.debug("Checking introns using data structure for %s; introns: %s",
                                transcript.id, introns)
        for intron in transcript.introns:
            transcript.logger.debug("Checking intron %s%s:%d-%d for %s",
                                    transcript.chrom, transcript.strand,
                                    intron[0], intron[1], transcript.id)
            if (intron[0], intron[1], transcript.strand) in introns:
                transcript.logger.debug("Verified intron %s%s:%d-%d for %s",
                                        transcript.chrom, transcript.strand,
                                        intron[0], intron[1], transcript.id)
                transcript.verified_introns.add(intron)
            elif (intron[0], intron[1], None) in introns:
                transcript.logger.debug("Verified intron %s%s:%d-%d for %s",
                                        transcript.chrom, None,
                                        intron[0], intron[1], transcript.id)
                transcript.verified_introns.add(intron)

    transcript.logger.debug("Found these introns for %s: %s",
                            transcript.id, transcript.verified_introns)
    return


def retrieve_orfs(transcript):

    """This method will look up the ORFs loaded inside the database.
    During the selection, the function will also remove overlapping ORFs.

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    """

    # if self.query_id is None:
    #     return []

    trust_strand = transcript.json_conf["pick"]["orf_loading"]["strand_specific"]
    min_cds_len = transcript.json_conf["pick"]["orf_loading"]["minimal_orf_length"]

    orf_results = transcript.orf_baked(transcript.session).params(query=transcript.id,
                                                                  cds_len=min_cds_len)
    transcript.logger.debug("Retrieving ORFs from database for %s",
                            transcript.id)

    assert orf_results is not None

    if (transcript.monoexonic is False) or (transcript.monoexonic is True and trust_strand is True and
                                            transcript.strand is not None):
        # Remove negative strand ORFs for multiexonic transcripts,
        # or monoexonic strand-specific transcripts
        assert orf_results is not None
        candidate_orfs = list(orf for orf in orf_results if orf.strand != "-")
    else:
        candidate_orfs = orf_results.all()

    transcript.logger.debug("Found %d ORFs for %s",
                            len(candidate_orfs), transcript.id)
    assert isinstance(candidate_orfs, list)

    if len(candidate_orfs) == 0:
        return []
    else:
        result = [orf.as_bed12() for orf in candidate_orfs]
        for orf in result:
            assert orf.chrom == transcript.id, (orf.chrom, transcript.id)
        return result


def orf_sorter(orf):
    """Sorting function for the ORFs."

    :param orf: an ORF to sort
    :type orf: Mikado.serializers.orf.Orf
    """
    return (orf.cds_len,
            (orf.has_start_codon and orf.has_stop_codon),
            (orf.has_start_codon or orf.has_stop_codon),
            )


def find_candidate_orfs(transcript, graph, orf_dictionary) -> list:

    """
    Function that returns the best non-overlapping ORFs

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param graph: The NetworkX graph to be analysed
    :type graph: networkx.Graph

    :param orf_dictionary: a dictionary which contains the orf indexed by name
    :type orf_dictionary: dict

    :return:
    """

    candidate_orfs = []

    while len(graph) > 0:
        cliques = find_cliques(graph, logger=transcript.logger)
        communities = find_communities(graph, logger=transcript.logger)
        clique_str = []
        for clique in cliques:
            clique_str.append(str([(orf_dictionary[x].thick_start,
                                    orf_dictionary[x].thick_end) for x in clique]))
        comm_str = []
        for comm in communities:
            comm_str.append(str([(orf_dictionary[x].thick_start,
                                  orf_dictionary[x].thick_end) for x in comm]))
        transcript.logger.debug("{0} communities for {1}:\n\t{2}".format(
            len(communities),
            transcript.id,
            "\n\t".join(comm_str)))
        transcript.logger.debug("{0} cliques for {1}:\n\t{2}".format(
            len(cliques),
            transcript.id,
            "\n\t".join(clique_str)))

        to_remove = set()
        for comm in communities:
            comm = [orf_dictionary[x] for x in comm]
            best_orf = sorted(comm, key=orf_sorter, reverse=True)[0]
            candidate_orfs.append(best_orf)
            for clique in iter(cl for cl in cliques if best_orf.name in cl):
                to_remove.update(clique)
        graph.remove_nodes_from(to_remove)

    candidate_orfs = sorted(candidate_orfs, key=orf_sorter, reverse=True)
    return candidate_orfs
