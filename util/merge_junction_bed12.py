#!/usr/bin/env python3
# coding: utf-8

"""Script to merge BED12 files *based on the thickStart/End features*.
Necessary for merging junction files such as those produced by TopHat."""

import Mikado.parsers.bed12
import sys
import collections
import operator
import argparse
import io
import multiprocessing


def serialise(filename, tophat=False):
    """Function to serialise the BED12. Useful for multiprocessing.

    :param filename: name of the file to analyse
    :type filename: io.TextIOWrapper | str

    :param tophat: flag, derived from Namespace
    :type tophat: bool

    :rtype : (dict, str)
    """

    juncs = collections.defaultdict(list)

    if type(filename) is str:
        bed = Mikado.parsers.bed12.Bed12Parser(open(filename))
    else:
        if type(filename) is not io.TextIOWrapper:
            raise TypeError("Invalid BED file: {0}".format(type(filename)))
        bed = Mikado.parsers.bed12.Bed12Parser(filename)

    header = ''
    for record in bed:
        if record.header is True:
            if header == '' and record._line.rstrip() != '':
                header = record._line.rstrip()
            continue
        if tophat is True:
            record.thickStart = record.start+record.blockSizes[0]-1
            record.thickEnd = record.start+record.blockStarts[1]-1

        key = (record.chrom, record.strand, record.thickStart, record.thickEnd)
        juncs[key].append(record)

    return juncs, header


def main():
    """
    Main function of the script.
    :return: None
    """

    parser = argparse.ArgumentParser("""Script to merge BED12 files *based on the thickStart/End features*.
    Necessary for merging junction files such as those produced by TopHat""")
    parser.add_argument("--delim", default=";", type=str,
                        help="Delimiter for merged names. Default: %(default)s")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads to use for multiprocessing. Default: %(default)s")
    parser.add_argument("--tophat", default=False, action="store_true",
                        help="""Flag. If set, tophat-like junction style is assumed. This means that junctions
                        are defined using the blockSizes rather than thickStart/End. The script will convert
                        the lines to this latter format. By default, the script assumes that the intron start/end
                        are defined using thickStart/End like in portcullis.
                        Mixed-type input files are not supported.""")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"),
                        help="Output file. Default: stdout")
    parser.add_argument("bed", type=str, nargs="+",
                        help="Input BED files. Use \"-\" for stdin.")
    args = parser.parse_args()

    if args.bed[0] == "-":
        args.bed = [Mikado.parsers.bed12.Bed12Parser(sys.stdin)]
        results = [serialise(sys.stdin, args.tophat)]
    else:
        # The process pool has size which is the minimum of
        # number of available cores
        # number of BED files
        # number of requested cores
        pool = multiprocessing.Pool(processes=min(len(args.bed),
                                                  args.threads,
                                                  multiprocessing.cpu_count()
                                                  ))
        results = pool.starmap(serialise, zip(args.bed, [args.tophat]*len(args.bed)))
        pool.close()
        pool.join()

    header = ''
    juncs = collections.defaultdict(list)

    for res in results:
        if header == '' and res[1] != '':
            header = res[1]
        else:
            pass
        for key in res[0]:
            juncs[key].extend(res[0][key])

    print(header, file=args.output)
    # Sort by chrom, then start, then end, finally strand
    for key in sorted(juncs.keys(), key=operator.itemgetter(0, 2, 3, 1)):
        newrecord = Mikado.parsers.bed12.BED12(None)
        newrecord.header = False
        newrecord.chrom = key[0]
        newrecord.start = min(x.start for x in juncs[key])
        newrecord.end = max(x.end for x in juncs[key])
        newrecord.strand = key[1]
        newrecord.name = args.delim.join(x.name for x in juncs[key])
        newrecord.score = sum(x.score if x.score is not None else 0 for x in juncs[key])
        newrecord.thick_start = key[2]+1
        newrecord.thick_end = key[3]
        newrecord.rgb = "255,0,0"
        newrecord.block_count = 2
        newrecord.block_sizes = [newrecord.thick_start-newrecord.start+1,
                                newrecord.end-newrecord.thick_end+1]
        newrecord.block_starts = [0, newrecord.thick_end - newrecord.start+1]
        print(newrecord, file=args.output)

if __name__ == "__main__":
    main()


__author__ = 'venturil'
__version__ = "0.8"
