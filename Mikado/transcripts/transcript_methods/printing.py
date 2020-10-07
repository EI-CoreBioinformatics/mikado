"""
This module contains the methods related to creating the proper output lines for printing
GFFs/GTFs starting from the transcript class.
"""


from itertools import zip_longest
import functools
from ...parsers.GTF import GtfLine
from ...parsers.GFF import GffLine
from ...parsers.bed12 import BED12
import numpy as np


__author__ = 'Luca Venturini'


gff_constructor = GffLine.string_from_dict
gtf_constructor = GtfLine.string_from_dict


def __create_cds_lines(transcript,
                       cds_run,
                       tid,
                       to_gtf=False,
                       with_introns=False):

    """
    Private method to create the exon/UTR/CDS lines for printing
    out in GTF/GFF format.
    :param transcript: the transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript
    :param cds_run: the internal orf run we are preparing
    :param tid: name of the transcript
    :param to_gtf: boolean, indicates whether the lines should be GTF or GFF
    :param with_introns: boolean, if True introns will be added to the output
    we want GTF or GFF output
    :return:
    """

    exon_lines = []
    cds_begin = False
    counter = dict()

    line_creator = functools.partial(__create_exon_line,
                                     transcript,
                                     **{"to_gtf": to_gtf,
                                        "tid": tid})

    if with_introns is True:
        cds_run = cds_run[:]
        for intron in transcript.introns:
            cds_run.append(("intron", intron))

    cds_run = sorted(cds_run, key=lambda segment: (segment[1][0],
                                                   segment[0].lower()))

    for segment in cds_run:
        try:
            exon_line, counter, cds_begin = line_creator(segment,
                                                         counter,
                                                         cds_begin)
        except IndexError:
            raise IndexError(cds_run)
        exon_lines.append(exon_line)

    # assert not any(True for x in exon_lines if x.feature == "CDS" and x.phase is None), [str(_) for _ in exon_lines]

    return [str(line) for line in exon_lines]


# pylint: disable=too-many-arguments
def __create_exon_line(transcript, segment, counter, cds_begin,
                       tid="", to_gtf=False):
    """
    Private method that creates an exon line for printing.
    :param transcript: the transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param segment: a segment of the form (feature, start, end)
    :type segment: list(str, tuple)

    :param counter: a dict object that keeps track of how many exons,
    CDS, UTR segments we have already seen
    :type counter: dict

    :param cds_begin: boolean flag that indicates whether the CDS has already begun
    :type cds_begin: bool

    :param tid: name of the transcript
    :param to_gtf: boolean flag

    :return: exon_line, counter, cds_begin
    :rtype: str, dict, bool
    """

    if to_gtf is False:
        constructor = gff_constructor
        utr3_feature = "three_prime_UTR"
        utr5_feature = "five_prime_UTR"
    else:
        constructor = gtf_constructor
        utr3_feature = "3UTR"
        utr5_feature = "5UTR"

    assert segment[0] in ("UTR", "CDS", "exon", "intron"), segment

    if hasattr(transcript, "chrom"):
        chrom, source, strand = transcript.chrom, transcript.source, transcript.strand
        parent = transcript.parent
    else:
        chrom, source, strand = transcript["chrom"], transcript["source"], transcript["strand"]
        parent = transcript["parent"]

    if not source:
        source = "Mikado"
    if not strand:
        strand = "."

    phase = None
    if segment[0] == "UTR":
        if (cds_begin is True and strand == "-") or \
                (strand == "+" and cds_begin is False):
            feature = utr5_feature
            counter["five"] = counter.get("five", 0) + 1
            index = counter["five"]
        else:
            feature = utr3_feature
            counter["three"] = counter.get("three", 0) + 1
            index = counter["three"]
    elif segment[0] == "CDS":
        cds_begin = True
        counter["CDS"] = counter.get("CDS", 0) + 1
        index = counter["CDS"]
        feature = "CDS"
        try:
            phase = segment[2]
        except IndexError:
            raise IndexError(segment)
    else:
        counter[segment[0]] = counter.get(segment[0], 0) + 1
        index = counter[segment[0]]
        feature = segment[0]

    if phase is None:
        phase = "."

    data = {
        "chrom": chrom,
        "source": source if source else "Mikado",
        "feature": feature,
        "start": segment[1][0],
        "end": segment[1][1],
        "strand": strand if strand else ".",
        "phase": phase,
        "score": ".",
        "attributes": dict()
    }

    if to_gtf is True:
        # noinspection PyPropertyAccess
        # assert transcript.parent
        assert tid
        exon_line = constructor(data, gene=parent, transcript=tid)
    else:
        exon_line = constructor(data,
                                parent=tid,
                                mid="{0}.{1}{2}".format(tid, feature, index),
                                name=None,
                                attribute_order=None)

    return exon_line, counter, cds_begin
# pylint: enable=too-many-arguments


def create_lines_cds(transcript,
                     to_gtf=False,
                     with_introns=False,
                     all_orfs=False,
                     transcriptomic=False):

    """
    Method to create the GTF/GFF lines for printing in the presence of CDS information.
    WARNING: at the moment, the phase support is disabled.
    :param transcript: the transcript instance
    :type transcript: Mikado.loci.transcript.Transcript

    :param to_gtf: boolean, it indicates whether the output is GTF (True) or GFF3 (False)

    :param with_introns: boolean, if set to True, introns will be printed as well.
    :return:
    """

    if to_gtf is False:
        constructor = gff_constructor
    else:
        constructor = gtf_constructor

    lines = []
    transcript_counter = 0

    if transcript.is_coding is False:
        lines = create_lines_no_cds(transcript, to_gtf=to_gtf, with_introns=with_introns)
    else:
        if all_orfs is True:
            iterable = transcript.internal_orfs
        else:
            iterable = [transcript.selected_internal_orf]

        for index, cds_run in enumerate(iterable):
            transcript.logger.debug("CDS run for %s: %s", transcript.id, cds_run)
            if transcript.number_internal_orfs > 1 and all_orfs is True:
                transcript_counter += 1
                tid = "{0}.orf{1}".format(transcript.id, transcript_counter)

                if index == transcript.selected_internal_orf_index:
                    transcript.attributes["maximal"] = True
                else:
                    transcript.attributes["maximal"] = False
            else:
                tid = transcript.id
            cds_run = transcript.internal_orfs[index]

            if transcriptomic is False:
                parent_line = dict(
                    (attr, getattr(transcript, attr))
                    for attr in ("chrom", "source", "feature", "start", "end",
                                 "score", "strand", "attributes")
                )
                if parent_line["score"] is None:
                    parent_line["score"] = "."

                parent_line["phase"] = '.'
                if parent_line["source"] is None:
                    parent_line["source"] = "Mikado"
                if parent_line["strand"] is None:
                    parent_line["strand"] = "."

                if to_gtf is True:
                    parent_line["attributes"]["gene_id"] = transcript.parent
                    parent_line["attributes"]["transcript_id"] = tid
                else:
                    parent_line["attributes"]["parent"] = transcript.parent
                    parent_line["attributes"]["ID"] = tid

                parent_line["attributes"]["name"] = transcript.id

                exon_lines = __create_cds_lines(transcript,
                                                cds_run,
                                                tid,
                                                to_gtf=to_gtf,
                                                with_introns=with_introns)

                lines.append(constructor(parent_line))
                lines.extend(exon_lines)
            else:

                parent_line = dict(
                    (attr, getattr(transcript, attr))
                    for attr in ["source", "feature", "score", "attributes"]
                )

                if parent_line["score"] is None:
                    parent_line["score"] = "."

                parent_line["chrom"] = transcript.id
                parent_line["start"] = 1
                parent_line["end"] = transcript.cdna_length
                parent_line["strand"] = "+"
                parent_line["phase"] = "."
                # data["id = tid
                # parent_line.parent = "{}_gene".format(tid)

                if to_gtf is True:
                    lines.append(
                        constructor(parent_line, gene=transcript.parent, transcript=tid)
                    )
                else:
                    lines.append(
                        constructor(parent_line, parent=transcript.parent, mid=tid,
                                    name=transcript.name)
                    )

                new_cds_run = []

                cds = sorted([_ for _ in cds_run if _[0] == "CDS"])

                if transcript.strand == "+":
                    cds_start = cds[0][1][0]
                    phase = cds[0][2]
                    five_utr = [_[1] for _ in cds_run if _[0] == "UTR" and _[1][1] < cds_start]
                else:
                    cds_start = cds[-1][1][1]
                    phase = cds[-1][2]
                    five_utr = [_[1] for _ in cds_run if _[0] == "UTR" and _[1][0] > cds_start]
                if five_utr:
                    five_utr = sum([_[1] - _[0] + 1 for _ in five_utr])
                else:
                    five_utr = 0
                if five_utr:
                    cds_start = five_utr + 1
                else:
                    cds_start = 1
                cds_end = sum([_[1][1] - _[1][0] + 1 for _ in cds]) + cds_start - 1

                if five_utr:
                    new_cds_run.append(("UTR", (1, five_utr)))
                new_cds_run.append(("CDS", (cds_start, cds_end), phase))
                if cds_end < transcript.cdna_length:
                    new_cds_run.append(("UTR", (cds_end + 1, transcript.cdna_length)))
                exon_lines = __create_cds_lines(parent_line,
                                                new_cds_run,
                                                tid,
                                                to_gtf=to_gtf,
                                                with_introns=False)
                lines.extend(exon_lines)

    return lines


def as_bed12(transcript, transcriptomic=False):
    """
    Method to create a BED12 object for printing
    :param transcript: Mikado.loci.transcript.Transcript
    :return:
    """

    transcript.finalize()
    bed12 = BED12(table=transcript.codon_table)
    bed12.transcriptomic = False
    bed12.header = False
    bed12.chrom = transcript.chrom
    bed12.start = transcript.start
    bed12.end = transcript.end

    if transcript.is_coding is True:
        if transcript.strand != "-":
            try:
                phase = transcript.phases[transcript.selected_cds[0]]
            except KeyError:
                raise KeyError((transcript.selected_cds[0], transcript.phases))
        else:
            try:
                phase = transcript.phases[transcript.selected_cds[-1]]
            except KeyError:
                raise KeyError((transcript.selected_cds[-1], transcript.phases))

        name = "ID={ID};coding={coding};phase={phase}".format(
            ID=transcript.id,
            coding=transcript.is_coding,
            # Now we have to get the phase of the first CDS exon ..
            phase=phase)
    else:
        name = "ID={ID};coding={coding}".format(
            ID=transcript.id,
            coding=transcript.is_coding,
            # Now we have to get the phase of the first CDS exon ..
            )

    if transcript.alias is not None and transcript.alias != transcript.id:
        name += ";alias={}".format(transcript.alias)

    bed12.name = name
    bed12.score = transcript.score if transcript.score else 0
    bed12.strand = transcript.strand
    if transcript.is_coding:
        bed12.coding = True
        first_exon = [_ for _ in transcript.selected_cds if transcript.selected_cds_start in _]
        assert len(first_exon) == 1
        bed12.phase = transcript.phases[first_exon.pop()]
        bed12.thick_start = transcript.selected_cds[0][0]
        bed12.thick_end = transcript.selected_cds[-1][1]
    else:
        bed12.thick_start = bed12.thick_end = bed12.start
    bed12.block_count = transcript.exon_num
    bed12.block_sizes = [exon[1] - exon[0] + 1 for exon in transcript.exons]
    _introns = np.concatenate([np.array([intron[1] - intron[0] + 1 for intron in sorted(transcript.introns)],
                                        dtype=np.int64),
                               np.zeros(1, dtype=np.int64)])
    bed12.block_starts = np.concatenate([np.zeros(1, dtype=np.int64),
                                         (bed12.block_sizes + _introns).cumsum()[:-1]], axis=0)
    assert bed12.block_starts[0] == 0, bed12.block_starts
    if transcriptomic:
        bed12 = bed12.to_transcriptomic(alias=transcript.alias, start_adjustment=False,
                                        coding=transcript.is_coding)
        bed12.chrom = transcript.id
    return bed12


def create_lines_bed(transcript, transcriptomic=False):

    """
    Method to return the BED12 format of the transcript.
    :param transcript:
    :return:
    """

    return str(as_bed12(transcript, transcriptomic=transcriptomic))


def create_lines_no_cds(transcript,
                        to_gtf=False,
                        with_introns=False,
                        transcriptomic=False):

    """
    Method to create the GTF/GFF lines for printing in the absence of CDS information.

    :param transcript: the Transcript instance
    :type transcript: Mikado.loci_objects.transcript.Transcript

    :param to_gtf: boolean, it indicates whether the output is GTF (True) or GFF3 (False)
    :type to_gtf: bool
    """

    if to_gtf is True:
        constructor = gtf_constructor
    else:
        constructor = gff_constructor

    if transcriptomic is False:
        parent_line = {
            "chrom": transcript.chrom,
            "source": transcript.source if transcript.source else "Mikado",
            "feature": transcript.feature,
            "start": transcript.start,
            "end": transcript.end,
            "score": transcript.score if transcript.score is not None else ".",
            "strand": transcript.strand if transcript.strand else ".",
            "phase": ".",
            "attributes": transcript.attributes
        }

        if to_gtf is True:
            # This prevents a strange bug when converting TAIR10 files
            parent_line["attributes"]["transcript_id"] = transcript.id
            parent_line["attributes"]["gene_id"] = transcript.parent
        else:
            parent_line["attributes"]["ID"] = transcript.id
            parent_line["attributes"]["Parent"] = transcript.parent

        parent_line["attributes"]["name"] = transcript.name

        lines = [constructor(parent_line)]
        exon_lines = []

        if with_introns is False:
            intron_list = [None] * len(transcript.exons)
        else:
            intron_list = sorted(transcript.introns)
        counter = dict()

        line_creator = functools.partial(__create_exon_line,
                                         transcript,
                                         **{"to_gtf": to_gtf,
                                            "tid": transcript.id,
                                            "cds_begin": False})

        for exon, intron in zip_longest(sorted(transcript.exons),
                                        intron_list):
            exon_line, counter, _ = line_creator(("exon", exon), counter)

            exon_lines.append(str(exon_line))
            if intron is not None:
                intron_line, counter, _ = line_creator(("intron", intron), counter)
                exon_lines.append(str(intron_line))

        lines.extend(exon_lines)
    else:
        parent_line = {
            "chrom": transcript.id,
            "source": transcript.source if transcript.source else "Mikado",
            "feature": transcript.feature,
            "start": 1,
            "end": transcript.cdna_length,
            "score": transcript.score if transcript.score else ".",
            "strand": "+",
            "phase": ".",
            "attributes": transcript.attributes
        }


        if to_gtf is True:
            parent_line["attributes"]["transcript_id"] = transcript.id
            parent_line["attributes"]["gene_id"] = transcript.parent
        else:
            parent_line["attributes"]["ID"] = transcript.id
            parent_line["attributes"]["Parent"] = transcript.parent

        parent_line["attributes"]["Name"] = transcript.name

        lines = [constructor(parent_line)]

        line_creator = functools.partial(__create_exon_line,
                                         transcript,
                                         **{"to_gtf": to_gtf,
                                            "tid": transcript.id,
                                            "cds_begin": False})

        exon_line, _, _ = line_creator(("exon", (1, transcript.cdna_length)), 0)
        lines.append(exon_line)

    return lines
