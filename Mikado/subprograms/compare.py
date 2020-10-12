#!/usr/bin/env python3

"""
This subprogram is my version of ParsEval/CuffCompare. It produces
class codes for each prediction transcript, like cuffcompare,
and its output statistics are blatantly copied by that program;
however, it is more affine to ParsEval in the sense that its only
purpose is to compare two sets of annotations, with no interest in
assembling or checking the expression.
"""

import argparse
from ..scales.compare import compare
from ..parsers import to_gff
from multiprocessing import cpu_count


__author__ = "Luca Venturini"


def compare_parser():
    """
    The parser for the comparison function

    :return: the argument parser
    """

    def get_procs(arg):

        return max(min(int(arg), cpu_count()), 1)

    parser = argparse.ArgumentParser(
        'Tool to define the spec/sens of predictions vs. references.')
    input_files = parser.add_argument_group(
        'Prediction and annotation files.')
    input_files.add_argument('-r', '--reference',
                             type=to_gff,
                             required=True,
                             help="""Reference annotation file.
                             By default, an index will be crated and saved with the suffix
                             ".midx".""")
    targets = input_files.add_mutually_exclusive_group(required=True)
    targets.add_argument('-p', '--prediction',
                         type=to_gff,
                         help='Prediction annotation file.')
    targets.add_argument("--self", default=False,
                         action="store_true",
                         help="""Flag. If set, the reference will be compared with itself. \
Useful for understanding how the reference transcripts interact with each other. If this option is selected,
the stats file will not be produced.""")
    targets.add_argument("--internal", default=False,
                         action="store_true",
                         help="""Flag. If set, for each gene with more than one transcript isoform
                         each will be compared to the others. Useful for understanding the structural
                         relationships between the transcripts in each gene.""")
    targets.add_argument("--index", default=False,
                         action="store_true",
                         help="""Flag. If set, compare will stop after
                         having generated the GFF index for the reference.""")
    parser.add_argument('--distance', type=int, default=2000,
                        help='''Maximum distance for a transcript to be considered
                        a polymerase run-on. Default: %(default)s''')
    parser.add_argument('-pc', '--protein-coding',
                        dest="protein_coding", action="store_true",
                        default=False,
                        help="""Flag. If set, only transcripts with a CDS
                        (both in reference and prediction) will be considered.""")
    parser.add_argument("-o", "--out", default="mikado_compare", type=str,
                        help="Prefix for the output files. Default: %(default)s")
    parser.add_argument("-fm", "--fuzzy-intron-match", dest="fuzzymatch", default=0,
                        type=int,
                        help="""This parameter controls whether an intron should be considered as "matched" even if
                        its splices are within N bases from the annotated splice junctions. By default this is set to 0
                        (ie only proper matches count).""")
    parser.add_argument("--lenient", action="store_true", default=False,
                        help="""If set, exonic statistics will be calculated leniently
in the TMAP as well - ie they will consider an exon as match even if only
the internal junction has been recovered.""")
    parser.add_argument("-nF", "--do-not-report-fusions", dest="report_fusions", action="store_false",
                        default=True,
                        help="Flag. If invoked, Mikado compare will not report fusions in the input. Useful \
when the reference transcripts have not been clustered properly into genes (e.g. a Mikado prepare run).")
    parser.add_argument("-eu", "--exclude-utr", dest="exclude_utr",
                        default=False, action="store_true",
                        help="""Flag. If set, reference and prediction transcripts
                        will be stripped of their UTRs (if they are coding).""")
    parser.add_argument("-n", "--no-index", "--no-save-index", dest="no_save_index",
                        action="store_true", default=False,
                        help="""Unless this flag is set, compare will save an index of the
                        reference to quicken multiple calls.""")
    parser.add_argument("-erm", "--extended-refmap", action="store_true", default=False,
                        dest="extended_refmap",
                        help="""Flag. If set, the RefMap will also contain recall
                        and precision statistics - not just the F1.""")
    parser.add_argument("-upa", "--use-prediction-alias", action="store_true", default=False,
                        dest="use_prediction_alias",
                        help="""Flag. If set, Mikado Compare will use the alias - rather than the transcript ID -
                        to report the results for prediction transcripts in the TMAP and REFMAP files.""")
    parser.add_argument("-l", "--log", default=None, type=str)
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument("-z", "--gzip",
                        action="store_true",
                        default=False,
                        help="Flag. If set, TMAP and REFMAP files will be GZipped.")
    parser.add_argument("-x", "--processes", default=1,
                        type=get_procs)
    parser.set_defaults(func=compare)

    return parser


if __name__ == '__main__':
    __args__ = compare_parser().parse_args()
    __args__.func(__args__)
