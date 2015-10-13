#!/usr/bin/env python3
# coding: utf-8

""" Script to check the correctness of the strand for aligned/assembled transcripts. """

import sys
import os
from copy import deepcopy
import mikado_lib.parsers
import mikado_lib.loci_objects
from Bio import SeqIO
import argparse


def main():
    """
    Main function
    """

    def to_gff(string):
        """
        Function to check the input file type (GTF/GFF)
        :param string: filename
        :type string: str
        """
        if not os.path.exists(string) or not os.path.isfile(string) or not os.stat(string).st_size > 0:
            raise ValueError("Invalid input file.")
        if string[-4:] == ".gtf":
            gff_function = mikado_lib.parsers.GTF.GTF
        else:
            gff_function = mikado_lib.parsers.GFF.GFF3

        gff = gff_function(string)
        for ann_record in gff:
            if ann_record.header is False:
                gff.close()
                return gff_function(string)
        raise ValueError("Empty GFF file provided.")

    def to_seqio(string):
        """ Function to load a FASTA file into a SeqIO dictionary.
        :param string: filename
        :type string: str
        """
        if not os.path.exists(string) or not os.path.isfile(string) or not os.stat(string).st_size > 0:
            raise ValueError("Invalid input file.")
        print("Loading reference")
        seqdict = SeqIO.to_dict(SeqIO.parse(open(string), 'fasta'))
        print("Reference loaded")
        return seqdict

    parser = argparse.ArgumentParser("Script to check the correctness of the strand for aligned/assembled transcripts.")
    parser.add_argument("--fasta", type=to_seqio, required=True,
                        help="Genome FASTA file. Required.")
    parser.add_argument("-s", "--strand-specific", dest="strand_specific", action="store_true", default=False,
                        help="""Flag. If set, monoexonic transcripts will be left on their strand
                        rather than being moved to the unknown strand.""")
    parser.add_argument("-l", "--lenient", action="store_true", default=False,
                        help="""Flag. If set, transcripts with mixed +/- splices will not
                        cause exceptions but rather be annotated as problematic.
                        If set, the output will be GFF3, regardless of the input format.""")
    parser.add_argument("gff", type=to_gff, help="Input GFF3/GTF file.")
    parser.add_argument("out", default=sys.stdout, nargs='?', type=argparse.FileType('w'),
                        help="Output file. Default: STDOUT.")
    args = parser.parse_args()

    current_transcripts = dict()
    current_parent = None

    if args.gff.name[-3:] == "gtf":
        is_gff = False
    else:
        is_gff = True

    #   print(is_gff, file=sys.stderr)

    current_seq = None

    reversed_transcripts = 0

    for record in args.gff:
        if record.header is True:
            continue
        if current_seq is None or record.chrom not in current_seq:
            current_seq = dict()
            current_seq[record.chrom] = args.fasta[record.chrom]
        if record.is_parent is True:
            if is_gff is True and record.is_transcript is False:
                to_delete = []
                for tid, tr in current_transcripts.items():
                    try:
                        tr.check_strand()
                        if tr.reversed is True:
                            reversed_transcripts += 1
                    except mikado_lib.exceptions.IncorrectStrandError:
                        to_delete.append(tid)
                    except mikado_lib.exceptions.InvalidTranscript:
                        to_delete.append(tid)
                for tid in to_delete:
                    del current_transcripts[tid]
                if len(current_transcripts) > 0:
                    strands = dict()
                    for tid, tr in current_transcripts.items():
                        if tr.reversed is True:
                            reversed_transcripts += 1
                        if tr.strand not in strands:
                            strands[tr.strand] = []
                        strands[tr.strand].append(tr)
                    for strand in strands:
                        if strand == current_parent.strand:
                            print(current_parent, file=args.out)
                            parent_id = current_parent.id
                        else:
                            new_parent = deepcopy(current_parent)
                            new_parent.strand = strand
                            new_parent.id = "{0}.strand{1}".format(current_parent.id, strand)
                            print(new_parent, file=args.out)
                            parent_id = new_parent.id
                        for tr in strands[strand]:
                            tr.parent = parent_id
                            print(tr, file=args.out)
                current_parent = record
                current_transcripts = dict()

            else:
                tr = current_transcripts[list(current_transcripts.keys())[0]]
                try:
                    tr.check_strand()
                    print(tr, file=args.out)
                except mikado_lib.exceptions.IncorrectStrandError:
                    pass
                except mikado_lib.exceptions.InvalidTranscript:
                    pass
                current_transcripts = dict()
                transcript = mikado_lib.loci_objects.transcriptchecker.TranscriptChecker(
                    record,
                    current_seq,
                    strand_specific=args.strand_specific,
                    lenient=args.lenient
                )

                current_transcripts[transcript.id] = transcript

        elif record.is_transcript is True and is_gff is True:
            try:
                assert current_parent.id in record.parent, record
            except TypeError:
                raise TypeError(str(current_parent), str(record))
            transcript = mikado_lib.loci_objects.transcriptchecker.TranscriptChecker(
                record,
                current_seq,
                strand_specific=args.strand_specific,
                lenient=args.lenient
            )
            current_transcripts[transcript.id] = transcript

        elif record.is_exon is True:
            assert len(current_transcripts) > 0, (str(current_parent), current_transcripts)
            for parent in record.parent:
                try:
                    current_transcripts[parent].add_exon(record)
                except KeyError:
                    print(current_transcripts)
                    raise

    to_delete = []
    for tid, tr in current_transcripts.items():
        try:
            tr.check_strand()
            if tr.reversed is True:
                reversed_transcripts += 1
        except mikado_lib.exceptions.IncorrectStrandError:
            to_delete.append(tid)
        except mikado_lib.exceptions.InvalidTranscript:
            to_delete.append(tid)
    for tid in to_delete:
        del current_transcripts[tid]
    if len(current_transcripts) > 0:
        strands = dict()
        for tid, tr in current_transcripts.items():
            if tr.strand not in strands:
                strands[tr.strand] = []
            strands[tr.strand].append(tr)
        for strand in strands:
            parent_id = None
            if is_gff is True:
                if strand == current_parent.strand:
                    print(current_parent, file=args.out)
                    parent_id = current_parent.id
                else:
                    new_parent = deepcopy(current_parent)
                    new_parent.strand = strand
                    new_parent.id = "{0}.strand{1}".format(current_parent.id, strand)
                    print(new_parent, file=args.out)
                    parent_id = new_parent.id
            for tr in strands[strand]:
                if parent_id is not None:
                    tr.parent = parent_id
                print(tr, file=args.out)


if __name__ == '__main__':
    main()
