import argparse
from mikado_lib.scales.compare import compare
from mikado_lib.parsers.GTF import GTF
from mikado_lib.parsers.GFF import GFF3


def compare_parser():
    """
    The parser for the comparison function

    :return: the argument parser
    """

    def to_gtf(string):
        """Function to recognize the input file type and create the parser.

        :param string: input file name.
        :type string: str
        """

        if string.endswith(".gtf"):
            return GTF(string)
        elif string.endswith('.gff') or string.endswith('.gff3'):
            return GFF3(string)
        else:
            raise ValueError('Unrecognized file format.')

    parser = argparse.ArgumentParser('Tool to define the spec/sens of predictions vs. references.')
    input_files = parser.add_argument_group('Prediction and annotation files.')
    input_files.add_argument('-r', '--reference', type=to_gtf, help='Reference annotation file.', required=True)
    input_files.add_argument('-p', '--prediction', type=to_gtf, help='Prediction annotation file.', required=True)
    parser.add_argument('--distance', type=int, default=2000,
                        help='''Maximum distance for a transcript to be considered a polymerase run-on.
                        Default: %(default)s''')
    parser.add_argument('-pc', '--protein-coding', dest="protein_coding", action="store_true", default=False,
                        help="""Flag. If set, only transcripts with a CDS (both in reference and prediction)
                        will be considered.""")
    parser.add_argument("-o", "--out", default="mikado_compare", type=str,
                        help="Prefix for the output files. Default: %(default)s")
    parser.add_argument("-eu", "--exclude-utr", dest="exclude_utr", default=False, action="store_true",
                        help="""Flag. If set, reference and prediction transcripts will be stripped
                        of their UTRs (if they are coding)."""
                        )
    parser.add_argument("-l", "--log", default=None, type=str)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.set_defaults(func=compare)

    return parser

if __name__ == '__main__':
    args = compare_parser().parse_args()
    args.func(args)
