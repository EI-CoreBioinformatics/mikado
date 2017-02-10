import sys
import argparse
import tabulate
from ...scales import class_codes
import textwrap
from itertools import zip_longest


def launch(args):

    rows = []

    if len(args.code) > 0:
        codes = [class_codes.codes[_] for _ in args.code]
    elif len(args.category) > 0:
        codes = [class_codes.codes[_] for _ in class_codes.codes
                 if class_codes.codes[_].category in args.category]
    else:
        codes = [class_codes.codes[_] for _ in class_codes.codes]

    for code in codes:
        definition = textwrap.wrap(code.definition, 30)
        code_rows = zip_longest([code.code],
                                definition,
                                [code.ref_multi],
                                [code.pred_multi],
                                [code.nucl],
                                [code.junc],
                                [code.reverse],
                                code.category.split())
        rows.extend(code_rows)

    __table_format = tabulate._table_formats[args.format]

    if args.format not in ("grid", "fancy_grid"):

        print(tabulate.tabulate(rows,
                                headers=["Class code",
                                         "Definition",
                                         "Reference multiexonic?",
                                         "Prediction multiexonic?",
                                         "Nucleotide: RC, PC, F1",
                                         "Junction: RC, PC, F1",
                                         "Reverse",
                                         "Category"],
                                tablefmt=args.format))

    else:
        out_of_header = False
        separator = None

        for row in tabulate.tabulate(rows,
                                     headers=["Class code",
                                              "Definition",
                                              "Reference multiexonic?",
                                              "Prediction multiexonic?",
                                              "Nucleotide: RC, PC, F1",
                                              "Junction: RC, PC, F1",
                                              "Reverse",
                                              "Category"],
                                     tablefmt=args.format).split("\n"):
            if row[:2] == __table_format.lineabove.begin + __table_format.lineabove.hline:
                separator = row
                if not out_of_header:
                    print(row)
                continue
            if row[:2] == __table_format.linebelowheader.begin + __table_format.linebelowheader.hline:
                out_of_header = True
                print(row, file=args.out)
                continue
            elif out_of_header is False:
                print(row, file=args.out)
            elif row[:2] == __table_format.linebetweenrows[0] + __table_format.linebetweenrows[1]:
                continue
            elif row[0] == __table_format.datarow.begin:
                if row.strip().split(__table_format.datarow.sep)[1].strip() != "":
                    print(separator, file=args.out)
                print(row, file=args.out)
        print(separator, file=args.out)
        print(file=args.out)
    return


def code_parser():

    """
    Command line parser for the utility.
    """

    parser = argparse.ArgumentParser("Script to generate the available class codes.")
    parser.add_argument("-f", "--format", choices=tabulate.tabulate_formats, default="rst")
    parser.add_argument("-c", "--category", nargs="+", default=[],
                        choices=list(set(_.category for _ in class_codes.codes.values())))
    parser.add_argument("-o", "--out", type=argparse.FileType("w"), default=sys.stdout)
    parser.add_argument("code", nargs="*", help="Codes to query.",
                        default=[],
                        choices=[[]] + list(class_codes.codes.keys()))
    parser.set_defaults(func=launch)
    return parser





