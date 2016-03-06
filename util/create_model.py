#!/usr/bin/env python3

import re
import argparse
import csv
import pickle
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from Mikado.scales.resultstorer import ResultStorer
from Mikado.loci_objects import Transcript
import operator
# import collections
from scipy.stats import hmean


class MetricEntry:

    """
    Basic class that defines the metrics loaded from the Mikado file.
    """

    metrics = [_ for _ in Transcript.get_available_metrics() if
               _ not in ["tid", "parent", "score", "best_bits", "snowy_blast_score"]]

    def __init__(self, row):

        self.__tid = row["tid"]
        self.__locus = row["parent"]

        for key in self.metrics:
            if row[key].lower() == "true":
                row[key] = 1.0
            elif row[key].lower() == "false":
                row[key] = 0.0
            else:
                try:
                    row[key] = float(row[key])
                except ValueError as exc:
                    raise ValueError("Invalid value for key {0}: {1}.\n{2}".format(
                        key, row[key],
                        exc))
            setattr(self, key, row[key])

    @property
    def matrix_row(self):
        return [getattr(self, key) for key in self.metrics]

    @property
    def features(self):
        return self.metrics

    @property
    def tid(self):
        return self.__tid

    @property
    def locus(self):
        return self.__locus


def load_tmap(tmap_file) -> dict:

    """
    Function to serialise the TMAP file into the original ResultStorer objects.
    :param tmap_file:
    :return:
    """

    results = dict()
    with open(tmap_file) as refmap:
        for row in csv.DictReader(refmap, delimiter="\t"):
            vals = []
            if re.search(r"\.orf[0-9]*$", row["tid"]):
                row["tid"] = re.sub(r"\.orf[0-9]*$", "", row["tid"])

            for key in ResultStorer.__slots__:
                if key in ("n_prec", "n_recall", "n_f1",
                           "j_prec", "j_recall", "j_f1",
                           "e_prec", "e_recall", "e_f1"):
                    row[key] = tuple(float(_) for _ in row[key].split(","))
                elif key == "distance":
                    if row[key] == "-":
                        row[key] = 5001
                    else:
                        row[key] = int(row[key])
                elif key in ("ccode", "tid", "gid",):
                    row[key] = tuple(row[key].split(","))
                vals.append(row[key])
            result = ResultStorer(*vals)
            results[result.tid[0]] = result

    return results


def load_metrics(metrics_file) -> [MetricEntry]:

    metrics = []
    with open(metrics_file) as metrics_handle:
        reader = csv.DictReader(metrics_handle, delimiter="\t")
        for row in reader:
            row = MetricEntry(row)
            metrics.append(row)

    return metrics


def main():

    """
    Main script function.
    :return:
    """

    parser = argparse.ArgumentParser("Script to build a model for Mikado starting from TMAP and metric files.")
    parser.add_argument("-t", "--tmap", help="The TMAP file with the comparison results.", required=True)
    parser.add_argument("-m", "--metrics", help="The metrics file.", required=True)
    parser.add_argument("-o", "--out", help="Output file.", default="forest.model")
    args = parser.parse_args()

    # X should contain a matrix of features derived from the portcullis tab file
    # y should contain the labels (0 not a valid junction, 1 a valid junction).
    # Confirmed with the reference.

    # Load tab file and produce matrix
    # bed, tab = loadtab(args.input)
    tmap_results = load_tmap(args.tmap)

    print("# TMAP results: " + str(len(tmap_results)))

    # Load reference and add labels
    # ref = bed12.loadbed(args.reference, False, False)
    metrics = load_metrics(args.metrics)
    print("# metered transcripts:", len(metrics))

    X = np.zeros((len(metrics), len(MetricEntry.metrics)))
    y = []

    for index, entry in enumerate(metrics):
        X[index] = entry.matrix_row
        score = hmean([tmap_results[entry.tid].j_f1,
                       tmap_results[entry.tid].e_f1,
                       tmap_results[entry.tid].n_f1])
        y.append(score)

    clf = RandomForestRegressor(n_estimators=int(len(MetricEntry.metrics)/3),
                                max_depth=None,
                                n_jobs=10,
                                random_state=0)

    clf.fit(X, y)
    clf.metrics = MetricEntry.metrics
    importances = clf.feature_importances_
    # std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
    # indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print("Feature ranking:")

    ordered = sorted([(MetricEntry.metrics[_], importances[_]) for _ in range(len(importances))],
                     key=operator.itemgetter(1), reverse=True)

    for rank, couple in enumerate(ordered, start=1):
        print(rank, "feature", couple[0], couple[1] * 100)

    print("Total contribution", 100 * sum([_[1] for _ in ordered]))

    # for f in range(X.shape[1]):
    #     print("{0}, feature {1} ({2})".format(
    #         f + 1,
    #         MetricEntry.metrics[f],
    #         importances[f]
    #     ))

    with open(args.out, "wb") as forest:
        pickle.dump(clf, forest)
    # scores = cross_val_score(clf, X, y)
    # scores.mean()

    # with open("Xy.pickle", "wb") as xy:
    #     pickle.dump([X, y], xy)

    # print("#" * 30)
    #
    # loci = collections.defaultdict(list)
    # for metric in metrics:
    #     loci[metric.locus].append(metric)
    #
    # for locus, locus_metrics in loci.items():
    #     scores = clf.predict([_.matrix_row for _ in locus_metrics])
    #
    #     scores = zip([_.tid for _ in locus_metrics], scores)
    #     print(locus)
    #     for score in sorted(scores, key=operator.itemgetter(1), reverse=True):
    #         print("\t", *score)
    # Plot the feature importances of the forest
    # plt.figure()
    # plt.title("Feature importances")
    # plt.bar(range(X.shape[1]), importances[indices],
    #        color="r", yerr=std[indices], align="center")
    # plt.xticks(range(X.shape[1]), TabEntry.sortedFeatures(indices))
    # locs, labels = plt.xticks()
    # plt.setp(labels, rotation=90)
    # plt.xlim([-1, X.shape[1]])
    # plt.tight_layout()
    #
    # plt.savefig(args.output + ".png")

if __name__ == "__main__":
    main()
