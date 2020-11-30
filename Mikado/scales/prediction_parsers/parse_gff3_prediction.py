import functools
from ...transcripts import Gene


def parse_prediction_gff3(args, queue_logger, transmit_wrapper, constructor):
    """Method to parse GFF files. This will use the Gene, rather than Transcript, class."""

    gene = None
    done = 0
    lastdone = 1
    invalids = set()
    __found_with_orf = set()

    gconstructor = functools.partial(Gene,
                                     only_coding=args.protein_coding,
                                     logger=queue_logger,
                                     use_computer=False)
    for row in args.prediction:
        if row.header is True:
            queue_logger.debug("Skipping row %s", row)
            continue
        elif row.is_gene is True:
            if gene is not None:
                gene.finalize(exclude_utr=args.exclude_utr)
                for transcript in gene:
                    transcript.finalize()
                    done, lastdone, __found_with_orf = transmit_wrapper(transcript=transcript,
                        done=done, lastdone=lastdone,
                        __found_with_orf=__found_with_orf)
            gene = gconstructor(row)
            queue_logger.debug("Creating gene %s", gene.id)
        elif row.is_transcript is True or row.feature == "match":
            queue_logger.debug("Analysing %s", row)
            transcript = constructor(row)
            if gene is None or gene.id not in transcript.parent:  # Orphan transcript
                queue_logger.debug(
                    "%s is an orphan transcript (Parent: %s, GeneID: %s)",
                    transcript.id, ",".join(transcript.parent), None if gene is None else gene.id)
                if gene is not None:
                    gene.finalize(exclude_utr=args.exclude_utr)
                    queue_logger.debug("Sending transcripts of %s", gene.id)
                    for gtranscript in gene:
                        done, lastdone, __found_with_orf = transmit_wrapper(
                            transcript=gtranscript,
                            done=done, lastdone=lastdone,
                            __found_with_orf=__found_with_orf)
                gene = gconstructor(transcript)
                queue_logger.debug("Creating gene %s", gene.id)
            else:
                gene.add(transcript)
        elif row.is_exon is True:
            if any(_ in invalids for _ in row.parent):
                # Skip children of invalid things
                continue
            if gene is None:
                gene = gconstructor(row)
            elif any(parent in gene.transcripts for parent in row.parent):
                gene.add_exon(row)
                assert gene.transcripts
            elif any(parent == gene.id for parent in row.parent) or gene.id == row.id:
                gene.add_exon(row)
                assert gene.transcripts
            else:
                queue_logger.debug("Creating gene from %s", row.id)
                if gene is not None:
                    gene.finalize(exclude_utr=args.exclude_utr)
                    assert gene.transcripts, (gene.id, str(row))
                    for gtranscript in gene:
                        queue_logger.debug("Sending %s", gtranscript.id)
                        done, lastdone, __found_with_orf = transmit_wrapper(
                            transcript=gtranscript,
                            done=done, lastdone=lastdone,
                            __found_with_orf=__found_with_orf)
                gene = gconstructor(row)
        else:
            queue_logger.warning("Skipped row: {}".format(row))

    if gene is not None:
        gene.finalize(exclude_utr=args.exclude_utr)
        for transcript in gene:
            transcript.finalize()
            done, lastdone, __found_with_orf = transmit_wrapper(
                transcript=transcript,
                done=done, lastdone=lastdone,
                __found_with_orf=__found_with_orf)
    return done, lastdone
