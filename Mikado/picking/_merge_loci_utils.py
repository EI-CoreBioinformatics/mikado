import rapidjson as json
import sqlite3
import msgpack
from ..loci import Superlocus
from ..loci import Locus
import sys
import collections
from ._locus_line_creator import _create_locus_lines


decoder = json.Decoder(number_mode=json.NM_NATIVE)
backup_decoder = json.Decoder()


def manage_index(data, dumps, source):
    index, data = data
    dump_index, gene_counter, gene_max = data
    orig_gene_counter = gene_counter
    conn = sqlite3.connect(dumps[dump_index])
    cursor = conn.cursor()
    batch = []
    try:
        stranded_loci = cursor.execute("SELECT json FROM loci WHERE counter=?", (str(index),)).fetchone()
    except ValueError:
        raise ValueError((index, type(index)))
    loci = []
    sublocus_dump = decoder(msgpack.loads(
        cursor.execute("SELECT json FROM subloci WHERE counter=?", (str(index),)).fetchone()[0],
        raw=False))

    sub_length = len(sublocus_dump)

    monolocus_dump = decoder(msgpack.loads(
        cursor.execute("SELECT json FROM subloci WHERE counter=?", (str(index),)).fetchone()[0],
        raw=False))

    mono_length = len(monolocus_dump)

    for pos, stranded_locus_json in enumerate(msgpack.loads(stranded_loci[0], raw=False)):
        try:
            stranded_locus = Superlocus(None)
            for locus_string in stranded_locus_json:
                locus_dict = decoder(locus_string)
                locus = Locus(None)
                locus.load_dict(locus_dict)
                stranded_locus.add_locus(locus)
            stranded_locus.source = source
        except (ValueError, json.JSONDecodeError):
            stranded_locus = Superlocus(None)
            for locus_string in stranded_locus_json:
                locus_dict = backup_decoder(locus_string)
                locus = Locus(None)
                locus = locus.load_dict(locus_dict)
                stranded_locus.add_locus(locus)

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
    batch = msgpack.dumps(batch)
    return batch


def __create_gene_counters(common_index: dict) -> (dict, int):
    """Function to assign to each counter in the database the correct base and maximum number of genes.
    This allows to parallelise the printing.
    """

    chroms, nums = list(zip(*[common_index[index][1:3] for index in range(1, max(common_index.keys()) + 1)]))
    total_genes = sum(nums)
    gene_counters = dict()
    chrom_tots = collections.defaultdict(list)
    assert len(chroms) == len(common_index), (len(chroms), len(common_index))
    for pos in range(len(chroms)):
        key = pos + 1
        chrom, num = chroms[pos], nums[pos]
        if chrom == '' and pos > 0:
            assert num == 0
            former = gene_counters[pos][0]
        elif pos == 0 or chrom != chroms[pos - 1]:
            if chroms[pos - 1] != "":
                former = 0
            else:  # The previous one is wrong ..
                prev_pos = pos - 1
                prev_chrom = chroms[prev_pos]
                while prev_chrom == "":
                    prev_pos -= 1
                    if prev_pos < 0:
                        break
                    prev_chrom = chroms[prev_pos]
                if prev_chrom == "" or prev_chrom != chrom:
                    former = 0
                else:
                    former = gene_counters[pos][0] + gene_counters[pos][1]
        else:
            former = gene_counters[pos][0] + gene_counters[pos][1]
        gene_counters[key] = (former, num)
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
    new_common = dict()
    assert min(common_index) == 1

    for key in common_index:
        new_common[key] = (common_index[key][0], gene_counters[key][0], gene_counters[key][1])
    return new_common, total_genes
