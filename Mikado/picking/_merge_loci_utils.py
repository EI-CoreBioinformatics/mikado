import rapidjson as json
import msgpack
from ..loci import Superlocus
from ..loci import Locus
import sys
import collections
import itertools
import numpy as np
from ._locus_line_creator import _create_locus_lines


decoder = json.Decoder()


def manage_index(index, data, source):
    index, (chrom, gene_counter, gene_max) = index[0], index[1]
    orig_gene_counter = gene_counter
    batch = []
    chrom, num_genes, stranded_loci, sublocus_dump, monolocus_dump = data[index]

    loci = []
    sublocus_dump = decoder(msgpack.loads(sublocus_dump, raw=False))

    sub_length = len(sublocus_dump)

    monolocus_dump = decoder(msgpack.loads(monolocus_dump, raw=False))

    mono_length = len(monolocus_dump)

    for pos, stranded_locus_json in enumerate(msgpack.loads(stranded_loci, raw=False)):
        stranded_locus = Superlocus(None)
        for locus_string in stranded_locus_json:
            locus_dict = decoder(locus_string)
            locus = Locus(None)
            locus.load_dict(locus_dict)
            if locus is not None:
                stranded_locus.add_locus(locus)
        stranded_locus.source = source

        if not stranded_locus.id.endswith(str(sys.maxsize)):
            loci.append(stranded_locus.id)

        minibatch, gene_counter = _create_locus_lines(stranded_locus,
                                                      gene_counter)
        minibatch = [minibatch]

        if pos < sub_length:
            assert len(sublocus_dump[pos]) == 3, sublocus_dump
            minibatch.append(sublocus_dump[pos])
        else:
            minibatch.append([])

        if pos < mono_length:
            minibatch.append(monolocus_dump[pos])
        else:
            minibatch.append([])

        batch.append(minibatch)

    assert (gene_counter - orig_gene_counter) == gene_max, (orig_gene_counter, gene_counter, gene_max)
    if len(set(loci)) != len(loci):
        seen = set()
        duplicated = []
        for lid in loci:
            if lid in seen:
                duplicated.append(lid)
            else:
                seen.add(lid)
        raise ValueError("Duplicated loci in counter {}! {}".format(index, duplicated))
    batch = [loci, batch]
    # batch = msgpack.dumps(batch)
    return batch


def __create_gene_counters(common_index: dict) -> (dict, int):
    """Function to assign to each counter in the database the correct base and maximum number of genes.
    This allows to parallelise the printing.
    The common index has the following structure:

    d[counter] = (database index, chrom, number of genes in locus)
    """

    chroms = []
    num_genes = []

    for index in range(1, max(common_index.keys()) + 1):
        chrom, n_genes = common_index[index][:2]
        chroms.append(chrom)
        num_genes.append(n_genes)

    chroms = np.array(chroms)
    num_genes = np.array(num_genes)

    gene_counters = dict()
    total_genes = sum(num_genes)

    chrom_tots = collections.defaultdict(list)
    for chrom in np.unique(chroms):
        index = np.where(chroms == chrom)
        totals = num_genes[index]
        cumu = totals.cumsum()
        for counter, former, num in zip(index[0], itertools.chain([0], cumu[:-1]), totals):
            gene_counters[counter + 1] = (former, num)
            if chrom:
                chrom_tots[chrom].extend(list(range(former + 1, former + num + 1)))

    tot_found = 0
    for chrom in chrom_tots:
        if len(set(chrom_tots[chrom])) != len(chrom_tots[chrom]):
            seen = set()
            duplicated = set()
            for num in chrom_tots[chrom]:
                if num in seen:
                    duplicated.add(num)
                else:
                    seen.add(num)
            raise AssertionError((chrom,
                                  len(set(chrom_tots[chrom])),
                                  len(chrom_tots[chrom]), max(chrom_tots[chrom]),
                                  duplicated,
                                  chrom_tots[chrom]))
        if len(chrom_tots[chrom]) > 0:
            assert len(list(range(1, chrom_tots[chrom][-1] + 1))) == len(chrom_tots[chrom])
            tot_found += chrom_tots[chrom][-1]

    assert tot_found == total_genes, (tot_found, total_genes)

    assert min(common_index) == 1

    new_common = dict()
    for key in common_index:
        # chrom, former id, new id
        new_common[key] = (common_index[key][0], gene_counters[key][0], gene_counters[key][1])
    return new_common, total_genes
