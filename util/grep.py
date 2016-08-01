#!/usr/bin/env python3
# coding: utf-8

""" A clone of grep -f."""

import sys
import time
import argparse
import os


def main():
    """Main script class."""

    parser = argparse.ArgumentParser(
        description='''This script is basically an efficient version of the GNU "grep -f" utility for table-like files, and functions with a similar sintax.''')
    parser.add_argument('-v', '--reverse', action='store_true', default=False,
                        help='Equivalent to the "-v" grep option')
    parser.add_argument('-s', '--separator', type=str, default=None,
                        help='The field separator. Default: consecutive whitespace(s)')
    parser.add_argument('-f', '--field', type=int, default=1, help='The field to look in the target file.')
    parser.add_argument('-q', '--quiet', default=True, action='store_false', help='No logging.')
    parser.add_argument('ids', help='The file of patterns to extract', type=str)
    parser.add_argument('target', help='The file to filter', type=argparse.FileType())
    parser.add_argument('out', nargs='?', help='The output file', type=argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()

    start = time.time()
    if not args.quiet:
        print(time.ctime(), 'Begin loading', file=sys.stderr)

    if os.path.exists(args.ids):
        ids_d = set([line.rstrip() for line in open(args.ids).readlines()])
    else:
        ids_d = set(args.ids.split(","))
    if not args.quiet:
        print(len(ids_d), file=sys.stderr)

    if not args.quiet:
        print(time.ctime(), 'Finished loading', file=sys.stderr)
    if args.field > 0:
        args.field -= 1

    try:
        for line in args.target:
            fields = line.rstrip().split(args.separator)
            if len(fields) <= args.field:
                continue
                # print("Too few fields ({0}), you asked for field n. {1}".format(len(fields), args.field), file=sys.stderr)
                # print(fields, file=sys.stderr)
                # raise ValueError("Wrong number of fields")
            else:
                try:
                    idl = fields[args.field]
                except IndexError:
                    raise IndexError(fields)

            if args.reverse:
                if idl not in ids_d:
                    print(line, file=args.out, end='')
            else:
                if idl in ids_d:
                    print(line, file=args.out, end='')
            if not line.strip():
                print(line, file=args.out, end='')
    except BrokenPipeError:
        pass
    except KeyboardInterrupt:
        pass

    if not args.quiet:
        print(time.ctime(), 'Finished', file=sys.stderr)
    stop = time.time()
    tim = stop - start
    if not args.quiet:
        print(time.ctime(), 'The run lasted %f seconds' % tim, file=sys.stderr)

    return 0


if __name__ == '__main__':
    main()
