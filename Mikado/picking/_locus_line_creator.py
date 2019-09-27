from ..loci import Superlocus
import sys


def _create_locus_lines(stranded_locus: Superlocus,
                        gene_counter: int):
    json_conf = stranded_locus.json_conf

    if not stranded_locus.loci:
        return None, gene_counter

    for locus in stranded_locus.loci:
        gene_counter += 1
        new_id = "{0}.{1}G{2}".format(
            json_conf["pick"]["output_format"]["id_prefix"],
            stranded_locus.chrom, gene_counter)
        stranded_locus.loci[locus].id = new_id

    if len(stranded_locus.loci) > 0:
        assert stranded_locus.start != sys.maxsize
        assert stranded_locus.end != -sys.maxsize

    if stranded_locus.start != sys.maxsize:
        assert not stranded_locus.id.endswith("{0}--{0}".format(sys.maxsize))

    locus_lines = stranded_locus.__str__(
        print_cds=not json_conf["pick"]["run_options"]["exclude_cds"],
        level="loci")
    locus_metrics_rows = [x for x in stranded_locus.print_loci_metrics()]
    locus_scores_rows = [x for x in stranded_locus.print_loci_scores()]
    batch = [locus_lines, locus_metrics_rows, locus_scores_rows]
    return batch, gene_counter


