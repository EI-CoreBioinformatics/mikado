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


def tdel(t, path):
    for node in path:
        t = t[node]
    del t


def tget(t, path):
    for node in path[:-1]:
        t = t[node]
    if path[-1] in t:
        return t[path[-1]]
    else:
        return None


def print_tree(t, out):

    print(t[b"lines"][0], file=out)
    if len(t[b"lines"]) > 1:
        print(sorted(t[b"lines"][1:]), file=out)

    for node in t:
        if node == b"lines":
            continue
        newt = t[node]
        assert isinstance(t, defaultdict), (newt, node)
        print(newt)
        print_tree(newt, out)


def tappend(t, path, value):

    try:
        while len(path) > 1:
            node = path.popleft()
            t = t[node]

        node = path.popleft()

        if node in t:
            t[node].append(value)
        else:
            t[node] = [value]
    except TypeError as exc:
        raise TypeError("{} {} {}\n{}".format(exc, type(path), path, t))


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("gff3", type=GFF3,
                        help="Input GFF3 file")
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?",
                        default=sys.stdout)
    args = parser.parse_args()

    index = 0

    heads = tree()
    child2parent = dict()
    positions = tree()
    lines = dict()

    logger = create_default_logger("sorter", level="ERROR")

    for row in args.gff3:
        if row.header is True:
            continue
        else:
            index += 1
            if row.is_parent is True:
                key = (row.start, row.end, row.strand)
                tappend(positions, deque([row.chrom, key]), row.id)
                tappend(heads, deque([row.id, b"lines"]), row)
            else:
                if isinstance(row.parent, (str, bytes)):
                    parents = [row.parent]
                else:
                    parents = row.parent

                for parent in parents:
                    parent_path = deque()
                    parent_path.append(b"lines")
                    parent_path.append(parent)
                    if row.id is not None:
                        child2parent[row.id] = parent
                    while parent not in heads:
                        assert parent in child2parent, (parent, child2parent)
                        parent = child2parent[parent]
                        parent_path.appendleft(parent)
                    tappend(heads, parent_path, row)

    for chrom in sorted(positions.keys()):
        for key in sorted(positions[chrom]):
            for gid in positions[chrom][key]:
                print_tree(heads[gid], args.out)

                print("###")

    args.out.close()
    args.gff3.close()
    return

main()