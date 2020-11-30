from ...parsers.bed12 import Bed12Parser
import collections
from ...transcripts import Transcript, Gene


def prepare_reference(reference, queue_logger, ref_gff=False,
                      exclude_utr=False, protein_coding=False) -> (dict, collections.defaultdict):

    """
    Method to prepare the data structures that hold the reference
    information for the parsing.
    :param reference:
    :param queue_logger:
    :param ref_gff:
    :return: genes, positions
    :param exclude_utr:
    :param protein_coding:
    """

    # from ..parsers.GFF import GFF3
    # from ..parsers.GTF import GTF
    # from ..parsers.bam_parser import BamParser

    genes = dict()
    positions = collections.defaultdict(dict)
    transcript2gene = dict()

    for row in reference:
        # Assume we are going to use GTF for the moment
        if row.header is True:
            continue
        elif reference.__annot_type__ == Bed12Parser.__annot_type__:
            transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
            if transcript.parent:
                gid = transcript.parent[0]
            else:
                transcript.parent = gid = row.id
            transcript2gene[row.id] = gid
            if gid not in genes:
                genes[gid] = Gene(transcript, gid=gid, logger=queue_logger)
            genes[gid].add(transcript)
        elif row.is_gene is True:
            gid = row.id
            genes[gid] = Gene(row, gid=gid, logger=queue_logger)
        elif row.is_transcript is True or (ref_gff is True and row.feature == "match"):
            queue_logger.debug("Transcript\n%s", str(row))
            transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
            if row.feature == "match":
                gid = row.id
            else:
                gid = row.gene

            if gid is None:
                queue_logger.warning("No gene ID found for %s, creating a mock one.", row.id)
                row.parent = f"{row.id}.gene"
                gid = row.parent[0]

            transcript2gene[row.id] = gid
            assert gid is not None
            if gid not in genes:
                genes[gid] = Gene(transcript, gid=gid, logger=queue_logger)
            genes[gid].add(transcript)
            assert transcript.id in genes[gid].transcripts
        elif row.is_exon is True:
            if ref_gff is True:
                if "cDNA_match" in row.feature:
                    row.parent = row.id
                    # row.gene = row.id
                    if row.id not in transcript2gene:
                        genes[row.id] = Gene(None, gid=row.id, logger=queue_logger)
                        transcript2gene[row.id] = row.id
                        transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
                        genes[row.id].add(transcript)
                found = False
                for transcript in row.transcript:
                    if transcript in transcript2gene:
                        # We have to perform the check because there are some GFFs
                        # e.g. TAIR
                        # where CDSs are defined within a spurious "Protein" feature
                        found = True
                        gid = transcript2gene[transcript]
                        genes[gid][transcript].add_exon(row)
                if found is False:
                    for transcript in row.transcript:
                        if transcript in genes:  # for pseudogenes and the like
                            found = True
                            genes[transcript].add_exon(row)
                if found is False:
                    queue_logger.warn("This feature has no corresponding transcript! %s",
                                      str(row))
            else:
                if row.gene in genes and row.transcript in genes[row.gene].transcripts:
                    genes[row.gene][row.transcript].add_exon(row)
                else:
                    if row.gene not in genes:
                        genes[row.gene] = Gene(None, gid=row.gene, logger=queue_logger)
                    if row.transcript not in genes[row.gene]:
                        transcript = Transcript(row, logger=queue_logger, trust_orf=True, accept_undefined_multi=True)
                        transcript2gene[row.id] = row.gene
                        genes[row.gene].add(transcript)
                    try:
                        genes[row.gene][row.transcript].add_exon(row)
                    except TypeError as exc:
                        queue_logger.critical("Failed with %s", exc)
                        queue_logger.critical("%s, %s", row.transcript, row)
                        raise

    genes, positions = finalize_reference(genes, positions, queue_logger, exclude_utr=exclude_utr,
                                          protein_coding=protein_coding)

    if len(genes) == 0:
        raise KeyError("No genes remained for the reference!")
    return genes, positions


def finalize_reference(genes, positions, queue_logger, exclude_utr=False, protein_coding=False) \
        -> (dict, collections.defaultdict):

    """
:param genes:
:param positions:
:param queue_logger:
:param exclude_utr:
:param protein_coding:
:return:
"""

    non_coding_to_remove = set()
    genes_to_remove = set()
    for gid in genes:
        genes[gid].logger = queue_logger
        genes[gid].finalize(exclude_utr=exclude_utr)
        if len(genes[gid]) == 0:
            genes_to_remove.add(gid)
            continue
        if protein_coding is True:
            to_remove = []
            for tid in genes[gid].transcripts:
                if genes[gid].transcripts[tid].combined_cds_length == 0:
                    to_remove.append(tid)
                    queue_logger.debug("No CDS for %s", tid)
            if len(to_remove) == len(genes[gid].transcripts):
                non_coding_to_remove.add(gid)
                queue_logger.debug("Noncoding gene: %s", gid)
                continue
            elif len(to_remove) > 0:
                for tid in to_remove:
                    genes[gid].remove(tid)
        key = (genes[gid].start, genes[gid].end)
        if key not in positions[genes[gid].chrom]:
            positions[genes[gid].chrom][key] = []
        positions[genes[gid].chrom][key].append(gid)

    for gid in genes_to_remove:
        queue_logger.warn("Removed from reference: %s; error: %s",
                          gid, genes[gid].exception_message)
        del genes[gid]
    for gid in non_coding_to_remove:
        del genes[gid]
    return genes, positions
