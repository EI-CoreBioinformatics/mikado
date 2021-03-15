from ..loci import Superlocus
from functools import partial
from ..utilities import default_for_serialisation
import rapidjson as json
dumper = partial(json.dumps, default=default_for_serialisation)
import msgpack


def serialise_locus(stranded_loci: [Superlocus],
                    accumulator,
                    counter,
                    print_subloci=True,
                    print_cds=True,
                    print_monosubloci=True):

    loci = []
    subloci = []
    monoloci = []
    for stranded_locus in sorted(stranded_loci):
        loci.append([json.dumps(locus.as_dict(), default=default_for_serialisation)
                     for locus in stranded_locus.loci.values()])

        if print_subloci is True:
            batch = []
            batch.append(stranded_locus.__str__(level="subloci", print_cds=print_cds))
            for sublocus in stranded_locus.subloci:
                assert sublocus.scores_calculated is True, (sublocus.id, sublocus.scores_calculated)

            sub_metrics_rows = [_ for _ in stranded_locus.print_subloci_metrics() if _ != {} and "tid" in _]
            sub_scores_rows = [_ for _ in stranded_locus.print_subloci_scores() if _ != {} and "tid" in _]
            batch.append(sub_metrics_rows)
            batch.append(sub_scores_rows)
            subloci.append(batch)

        if print_monosubloci is True:
            batch = []
            batch.append(stranded_locus.__str__(level="monosubloci", print_cds=print_cds))
            mono_metrics_rows = [_ for _ in stranded_locus.print_monoholder_metrics() if _ is not None
                                and _ != {} and "tid" in _]
            mono_scores_rows = [_ for _ in stranded_locus.print_monoholder_scores() if _ is not None
                                and _ != {} and "tid" in _]
            batch.append(mono_metrics_rows)
            batch.append(mono_scores_rows)
            monoloci.append(batch)

    loci = msgpack.dumps(loci)
    subloci = msgpack.dumps(subloci)
    monoloci = msgpack.dumps(monoloci)
    if not stranded_loci:
        chrom = ""
        num_genes = 0
    else:
        chrom = stranded_loci[0].chrom
        num_genes = sum(len(slid.loci) for slid in stranded_loci)
    accumulator.append((counter, chrom, num_genes, loci, subloci, monoloci))
    return accumulator
