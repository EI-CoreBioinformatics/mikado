#!/usr/bin/env python3

from Mikado.parsers.GFF import GFF3
from Mikado.loci import Gene, Transcript
import sys
import argparse
from collections import defaultdict, deque
from Mikado.utilities.log_utils import create_default_logger
from collections import OrderedDict


__doc__ = """Script that tries to emulate GenomeTools sort"""



def tree():
    return defaultdict(tree)


def tset(t, path, value=None):
  for node in path[:-1]:
    t = t[node]
  if value is not None:
    t[path[-1]] = value
  else:
    t = t[path[-1]]


def tget(t, path):
    for node in path[:-1]:
        t = t[node]
    if path[-1] in t:
        return t[path[-1]]
    else:
        return None


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("gff3", type=GFF3,
                        help="Input GFF3 file")
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?",
                        default=sys.stdout)
    args = parser.parse_args()

    heads = dict()
    child2parent = dict()
    positions = dict()

    logger = create_default_logger("sorter", level="ERROR")

    for row in args.gff3:
        if row.header is True:
            continue
        elif row.is_parent is True:
            if row.id in heads:
                exc=KeyError("Duplicated parent ID: {}".format(row.id))
                logger.error(exc)
                raise exc
            heads[row.id] = tree()
            if row.chrom not in positions:
                positions[row.chrom] = defaultdict(list)
            positions[row.chrom][(row.start, row.end, row.strand)].append(row.id)
            # Use bytes to decrease the memory footprint
            heads[row.id][b"parent_line"] = row
            heads[row.id][b"children"] = dict()
        else:
            for parent in row.parent:
                if row.id is not None:
                    child2parent[row.id] = parent
                parent_order = deque()
                parent_order.appendleft(parent)
                while parent not in heads:
                    if parent not in child2parent:
                        exc=KeyError("Parent of {} not found!".format(row))
                        logger.error(exc)
                        raise exc
                    parent = child2parent[parent]
                    parent_order.appendleft(parent)
                node = None
                while parent_order:
                    if node is None:
                        node = heads[parent_order.popleft()][b"children"]
                    else:
                        node = node[parent_order.popleft()][b"children"]
                if row.id is not None:
                    assert row.id not in node
                    node[row.id] = dict()
                    node[row.id][b"parent_line"] = row
                    node[row.id][b"children"] = dict()
                else:
                    if b"lines" not in node:
                        node[b"lines"] = []
                    node[b"lines"].append(row)

    print("#gff-version\t3", file=args.out)
    chrom_keys = sorted(positions.keys())
    for chrom in chrom_keys:
        # positions[chrom] = sorted(positions[chrom])
        pos_keys = sorted(positions[chrom])
        new_dict = OrderedDict()
        for key in pos_keys:
            new_dict[key] = positions[chrom][key]
        positions[chrom] = new_dict
        start = pos_keys[0][0]
        end = pos_keys[-1][1]
        print("##sequence-region  {} {} {}".format(chrom, start, end), file=args.out)

    for chrom in chrom_keys:
        for position in positions[chrom]:
            # Now start to iterate over the parents ..
            for pid in positions[chrom][position]:
                pdict = heads[pid]
                print_children(pdict, args.out)
                print("###", file=args.out)

    args.out.close()
    args.gff3.close()
    return

main()