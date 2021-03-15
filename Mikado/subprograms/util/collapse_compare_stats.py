#!/usr/bin/env python3

import argparse
import collections
import sys
import tabulate
import pandas as pd
import re

__doc__ = "Script to collapse various stat files into one."


# 1 Command line:
# 2 /home/lucve/miniconda3/envs/mikado2/bin/mikado compare -r reference.gff3 -p Daijin/5-mikado/pick/permissive/mikado-permissive.loci.gff3 -o compare -l compare.log
# 3 18 reference RNAs in 12 genes
# 4 22 predicted RNAs in  15 genes
# 5 --------------------------------- |   Sn |   Pr |   F1 |
# 6                         Base level: 94.90  83.22  88.68
# 7             Exon level (stringent): 80.56  71.60  75.82
# 8               Exon level (lenient): 91.18  76.54  83.22
# 9                  Splice site level: 95.19  81.15  87.61
# 10                       Intron level: 96.84  88.19  92.31
# 11                  Intron level (NR): 94.34  79.37  86.21
# 12                 Intron chain level: 69.23  50.00  58.06
# 13            Intron chain level (NR): 69.23  50.00  58.06
# 14       Transcript level (stringent): 55.56  45.45  50.00
# 15   Transcript level (>=95% base F1): 72.22  59.09  65.00
# 16   Transcript level (>=80% base F1): 72.22  59.09  65.00
# 17          Gene level (100% base F1): 75.00  60.00  66.67
# 18         Gene level (>=95% base F1): 83.33  66.67  74.07
# 19         Gene level (>=80% base F1): 83.33  66.67  74.07
# 20
# 21 #   Matching: in prediction; matched: in reference.
# 22
# 23             Matching intron chains: 9
# 24              Matched intron chains: 9
# 25    Matching monoexonic transcripts: 4
# 26     Matched monoexonic transcripts: 4
# 27         Total matching transcripts: 13
# 28          Total matched transcripts: 13
# 29
# 30           Missed exons (stringent): 14/72  (19.44%)
# 31            Novel exons (stringent): 23/81  (28.40%)
# 32             Missed exons (lenient): 6/68  (8.82%)
# 33              Novel exons (lenient): 19/81  (23.46%)
# 34                     Missed introns: 3/53  (5.66%)
# 35                      Novel introns: 13/63  (20.63%)
# 36
# 37        Missed transcripts (0% nF1): 0/18  (0.00%)
# 38         Novel transcripts (0% nF1): 3/22  (13.64%)
# 39              Missed genes (0% nF1): 0/12  (0.00%)
# 40               Novel genes (0% nF1): 3/15  (20.00%)

digit = r"[0-9]*\.?[0-9]+"
accuracy_pat = re.compile(r"^\s+(.*):\s+({digit})\s+({digit})\s+({digit})\s+$".format(digit=digit))
match_pat = re.compile(r"\s+(.*): (\d+)")
missed_novel_pat = re.compile(r"\s+(.*): (\d+)/(\d+)\s+\(({digit})%\)".format(digit=digit))


accuracy_stats = {
    6: "Base level", 7: "Exon level (stringent)", 8: "Exon level (lenient)", 9: "Splice site level",
    10: "Intron level", 11: "Intron level (NR)", 12: "Intron chain level", 13: "Intron chain level (NR)",
    14: "Transcript level (stringent)", 15: "Transcript level (>=95% base F1)",
    16: "Transcript level (>=80% base F1)", 17: "Gene level (100% base F1)", 18: "Gene level (>=95% base F1)",
    19: "Gene level (>=80% base F1)"
}


match_stats = {
    23: "Matching intron chains", 24: "Matched intron chains",
    25: "Matching monoexonic transcripts", 26: "Matched monoexonic transcripts",
    27: "Total matching transcripts", 28: "Total matched transcripts"
}


novel_missed_stats = {
    30: "Missed exons (stringent)", 31: "Novel exons (stringent)", 32: "Missed exons (lenient)",
    33: "Novel exons (lenient)", 34: "Missed introns", 35: "Novel introns",
    37: "Missed transcripts (0% nF1)",  38: "Novel transcripts (0% nF1)",
    39: "Missed genes (0% nF1)", 40: "Novel genes (0% nF1)"
}


def accuracy_retrieval(line: str, requested: str) -> (float, float, float):
    level, sn, pr, f1 = accuracy_pat.search(line).groups()
    assert level == requested, (line, requested)
    sn, pr, f1 = [float(val) for val in (sn, pr, f1)]
    assert 0 <= min([sn, pr, f1]) <= max([sn, pr, f1]) <= 100
    return sn, pr, f1


def match_retrieval(line: str, requested: str) -> int:
    level, match = match_pat.search(line).groups()
    assert requested == level, (line, requested)
    match = int(match)
    assert match >= 0
    return match


def missed_novel_retrieval(line: str, requested: str) -> (int, int, float):

    level, found, maximum, proportion = missed_novel_pat.search(line).groups()
    assert level == requested, (line, requested)
    found, maximum = int(found), int(maximum)
    proportion = float(proportion)
    assert 0 <= proportion <= 100
    proportion = f"{proportion}%"
    return found, maximum, proportion


def launch(args):
    if args.avf is True:
        print("Available formats for this script:")
        print(*tabulate.tabulate_formats, sep="\n")
        sys.exit(0)

    data = {
        "Sn": dict((key, collections.defaultdict(dict)) for counter, key in sorted(accuracy_stats.items())),
        "Pr": dict((key, collections.defaultdict(dict)) for counter, key in sorted(accuracy_stats.items())),
        "F1": dict((key, collections.defaultdict(dict)) for counter, key in sorted(accuracy_stats.items())),
        "Matches": dict((key, collections.defaultdict(dict)) for counter, key in sorted(match_stats.items())),
        "MNR":{
            "Found": dict((key, collections.defaultdict(dict)) for counter, key
                          in sorted(novel_missed_stats.items())),
            "Maximum": dict((key, collections.defaultdict(dict)) for counter, key
                              in sorted(novel_missed_stats.items())),
            "Proportion": dict((key, collections.defaultdict(dict)) for counter, key in
                               sorted(novel_missed_stats.items()))
        }
    }

    for stat in args.stat:
        for line_counter, line in enumerate(stat, start=1):
            if line_counter in accuracy_stats:
                level = accuracy_stats[line_counter]
                sn, pr, f1 = accuracy_retrieval(line, level)
                data["Sn"][level][stat.name] = sn
                data["Pr"][level][stat.name] = pr
                data["F1"][level][stat.name] = f1
            elif line_counter in match_stats:
                level = match_stats[line_counter]
                match = match_retrieval(line, level)
                data["Matches"][level][stat.name] = match
                continue
            elif line_counter in novel_missed_stats:
                level = novel_missed_stats[line_counter]
                found, maximum, proportion = missed_novel_retrieval(line, level)
                if stat.name not in data["MNR"]:
                    data["MNR"][stat.name] = collections.defaultdict(dict)

                data["MNR"][stat.name][level]["Total"] = found
                data["MNR"][stat.name][level]["Out of"] = maximum
                data["MNR"][stat.name][level]["Proportion"] = proportion
                continue
            else:
                continue

    for tablefmt in args.tablefmt:
        if tablefmt in ("tsv", "html", "csv", "rst"):
            ext = tablefmt
        elif "latex" in tablefmt:
            ext = f"{tablefmt}.tex"
        else:
            ext = f"{tablefmt}.txt"

        for key in ["Sn", "Pr", "F1", "Matches"]:
            if not ("all" in args.levels or key.lower() in args.levels):
                continue
            if key == "Sn":
                name = "Sensitivity"
            elif key == "Pr":
                name = "Precision"
            elif key in ("F1", "Matches"):
                name = key[:]
            else:
                continue
            lname = name.lower()
            with open(f"{args.out}.{lname}.{ext}", "wt") as out:
                accuracy = pd.DataFrame.from_dict(dict(**data[key]), orient="index")
                accuracy.index.name = f"{name}"
                print(tabulate.tabulate(accuracy, showindex="always", headers="keys", tablefmt=tablefmt),
                      file=out)

        # Now the last table
        if not ("all" in args.levels or "missed_novel" in args.levels):
            continue
        dfs = []
        for fname, fdata in data["MNR"].items():
            df = pd.DataFrame.from_dict(fdata, orient="index")
            df.columns.names = [fname]
            dfs.append(df)

        df = pd.concat(dfs, keys=[df.columns.names[0] for df in dfs], axis=1)
        df.columns.names = ["File", "Statistic"]
        h = [df.columns.names[0] + "\n" + df.columns.names[1]] + [
            "\n".join([fname, stat]) if num % 3 == 1 else "\n" + stat for num, (fname, stat) in
            enumerate(df.columns.tolist())]
        with open(f"{args.out}.missed_novel.{ext}", "wt") as out:
            print(tabulate.tabulate(df, showindex="always", headers=h, tablefmt=tablefmt), file=out)


def parser():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", default="mikado_compare_selected_stats",
                        help="Prefix for the output files")
    parser.add_argument("-fmt", "--format", nargs="+", choices=tabulate.tabulate_formats,
                        dest="tablefmt", default=["grid"], metavar="",
                        help="List of formats to print the tables into, separated by a space. Available formats can be "
                             "seen using --available-formats. "
                             "Terminate the list with --. Default: %(default)s")
    parser.add_argument("-l", "--levels", nargs="+", choices=["all", "f1", "sn", "pr", "matches", "missed_novel"],
                        default=["f1", "sn", "pr"],
                        help="Levels to print, separated by a space. Terminate the list with --. Default: %(default)s.")
    parser.add_argument("-avf", "--available-formats", action="store_true", default=False, dest="avf",
                        help="Print out a list of available formats and exit.")
    parser.add_argument("stat", type=argparse.FileType("rt"),
                        nargs="*")
    parser.set_defaults(func=launch)
    return parser


if __name__ == "__main__":
    args = parser()
    launch(args)
