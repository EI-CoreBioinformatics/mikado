#!/usr/bin/env python3

import csv
import re
import argparse
from collections import defaultdict, Counter

__author__ = 'Luca Venturini'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metrics", type=argparse.FileType("r"))
    parser.add_argument("refmap", type=argparse.FileType("r"))
    args = parser.parse_args()

    refmap = set()

    for row in csv.DictReader(args.refmap, delimiter="\t"):
        refmap.add(row["tid"])

    csvver = csv.DictReader(args.metrics, delimiter="\t")
    fieldn = csvver.fieldnames

    loci = defaultdict(dict)
    for row in csvver:
        for metric in fieldn[3:]:
            # Kids do NOT do this at home ...
            row[metric] = eval(row[metric])
        tid = re.sub(r".orf[0-9]$", "", row["tid"])
        assert ".orf" not in tid
        loci[row["parent"]][tid] = row

    results = defaultdict(list)

    # For the moment I am going to assume that
    # there is only one possible winner per locus.
    for locus in loci:
        bests = set.intersection(set(loci[locus].keys()),
                                 refmap)
        if len(bests) == 0:
            continue

        for metric in fieldn[3:]:
            # Get all possible metrics
            metric_vals = [loci[locus][tid][metric] for tid in loci[locus]]
            for best in bests:
                if isinstance(loci[locus][best][metric], bool):
                    results[metric].append(("target", loci[locus][best][metric]))
                else:
                    maximum = max(metric_vals)
                    assert isinstance(maximum, (int,float)), (maximum, type(maximum))
                    minimum = min(metric_vals)
                    assert isinstance(minimum, (int,float))

                    if loci[locus][best][metric] == max(metric_vals) and max(metric_vals) > 0:
                        results[metric].append(("max", loci[locus][best][metric]))
                    elif loci[locus][best][metric] == min(metric_vals):
                        results[metric].append(("min", loci[locus][best][metric]))
                    else:
                        results[metric].append(("target", loci[locus][best][metric]))

    print()
    print("Metric", "Rescaler", "Target", "Vals", sep="\t")
    for metric in sorted(results):
        instances = Counter([x[0] for x in results[metric]])
        most_common = instances.most_common(1)[0][0]
        counting = Counter([x[1] for x in results[metric] if x[0]==most_common])
        if most_common == "target":
            target = Counter([x[1] for x in results[metric] if x[0]=="target"]).most_common(1)[0][0]
        else:
            target = None
        print(metric, instances.most_common(1)[0][0], target, counting, sep="\t")


if __name__ == "__main__":
    main()