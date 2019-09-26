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

    :param gene_counter: integer used to keep track of the gene count.
    :type gene_counter: int

    :param handles: list of lists of the handles used for printing.
    :type handles: [list]

    :param logger: logger

    :param json_conf: configuration dictionary
    :type json_conf: [dict|None]
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
            print(*[row[key] for key in sub_metrics.fieldnames],
                  sep="\t", file=sub_metrics.handle)
        for row in sub_scores_rows:
            print(*[row[key] for key in sub_scores.fieldnames],
                  sep="\t", file=sub_scores.handle)
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
            print(*[row[key] for key in mono_metrics.fieldnames],
                  sep="\t", file=mono_metrics.handle)
        for row in mono_scores_rows:
            print(*[row[key] for key in mono_scores.fieldnames],
                  sep="\t", file=mono_scores.handle)

    for locus in stranded_locus.loci:
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

    for row in locus_metrics_rows:
        print(*[row[key] for key in locus_metrics.fieldnames],
              sep="\t", file=locus_metrics.handle)
    for row in locus_scores_rows:
        print(*[row[key] for key in locus_scores.fieldnames],
              sep="\t", file=locus_scores.handle)
    # Necessary to flush out all the files
    [_.flush() for _ in handles if hasattr(_, "close")]
    return gene_counter
