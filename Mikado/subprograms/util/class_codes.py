#!/usr/bin/env python3
# coding: utf-8

"""Script to print out the class code definitions."""

import sys
from ...loci import Transcript
import re
import argparse
from collections import defaultdict

__author__ = 'Luca Venturini'


def result():

    table = defaultdict(dict)

    # Perfect matches
    table["="]["definition"] = "Complete intron chain match."
    table["="]["Pred_multiexonic"] = True
    table["="]["Ref_multiexonic"] = True
    table["="]["Defining constraints"] = "Junction F1 = 100.0%"

    table["_"]["definition"] = "Complete match between two monoexonic transcripts"
    table["_"]["Pred_multiexonic"] = False
    table["="]["Ref_multiexonic"] = False
    table["_"]["Defining constraints"] = "Nucleotide F1 >= 80%"

    table["m"]["definition"] = "Match between two monoexonic transcripts"
    table["m"]["Pred_multiexonic"] = False
    table["m"]["Ref_multiexonic"] = False
    table["m"]["Defining constraints"] = "Nucleotide F1 < 80%"

    # Multiexonic vs multiexonic
    table["n"]["definition"] = "Intron chain extension."
    table["n"]["Pred_multiexonic"] = False
    table["n"]["Ref_multiexonic"] = False
    table["n"]["Defining constraints"] = "Junction precision < 100% and junction recall == 100%, all novel splice sites outside of the reference exons."

    table["J"]["definition"] = "Intron chain extension, with novel splice sites within the terminal exons of the reference."
    table["J"]["Pred_multiexonic"] = False
    table["J"]["Ref_multiexonic"] = False
    table["J"]["Defining constraints"] = "Junction precision < 100% and junction recall == 100%, start or end of the first novel splice sites within the reference transcript terminal exons."

    table["j"]["definition"] = "Alternative splicing event."
    table["j"]["Pred_multiexonic"] = False
    table["j"]["Ref_multiexonic"] = False
    table["j"]["Defining constraints"] = "Junction F1 between 0 and 100%, junction precision < 100%."

    table["c"]["definition"] = "Transcript completely contained"
    table["c"]["Pred_multiexonic"] = "NA"
    table["c"]["Ref_multiexonic"] = True
    table["c"]["Defining constraints"] = "Nucleotide precision == 100%, junction precision == 100%, and junction recall < 100%."

    table["C"]["definition"] = "Transcript completely contained"
    table["C"]["Pred_multiexonic"] = True
    table["C"]["Ref_multiexonic"] = True
    table["C"]["Defining constraints"] = "Junction precision == 100% but nucleotide precision lower than 100%."

    table["h"]["definition"] = "Structural match between two multiexonic transcripts, where no splice site is conserved."
    table["h"]["Pred_multiexonic"] = True
    table["h"]["Ref_multiexonic"] = True
    table["h"]["Defining constraints"] = "Junction precision = 0%, nucleotide F1 > 0%, and at least one intron of the reference must overlap partially one intron of the prediction."

    table["o"]["definition"] = "Generic overlap between two multiexonic transcripts, without sharing any intron nucleotide."
    table["o"]["Pred_multiexonic"] = True
    table["o"]["Ref_multiexonic"] = True
    table["o"]["Defining constraints"] = "Junction precision = 0%, nucleotide F1 > 0%, no overlap between any combination of introns."

    table["I"]["definition"] = "Prediction completely contained within the introns of the reference.."
    table["I"]["Pred_multiexonic"] = True
    table["I"]["Ref_multiexonic"] = True
    table["I"]["Defining constraints"] = "Nucleotide F1 = 0%, the prediction is within the boundaries of the reference."

    table["rI"]["definition"] = "Reference completely contained within the introns of the prediction."
    table["rI"]["Pred_multiexonic"] = True
    table["rI"]["Ref_multiexonic"] = True
    table["rI"]["Defining constraints"] = "Nucleotide F1 = 0%, the reference is within the boundaries of the prediction."

    # At least one monoexonic
    table["mo"]["definition"] = "Generic monoexonic match vs a multiexonic reference, excluding cases where the transcript is matching against only one reference exon and has its borders inside the reference introns."
    table["mo"]["Pred_multiexonic"] = False
    table["mo"]["Ref_multiexonic"] = True
    table["mo"]["Defining constraints"] = "Nucleotide F1 between 0 and 1, the transcript must have its ends either within an exon or outside of the transcript - not inside an intron."

    table["e"]["definition"] = "Single exon transcript overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment. The transcript should not be bridging two or more reference exons."
    table["e"]["Pred_multiexonic"] = False
    table["e"]["Ref_multiexonic"] = True
    table["e"]["Defining constraints"] = "Nucleotide F1 between 0 and 1, the transcript must have at least one of its ends within an intron and over 10 bps distant from the its matching exon boundary."

    table["i"]["definition"] = "Monoexonic prediction intronic transcript."
    table["i"]["Pred_multiexonic"] = False
    table["i"]["Ref_multiexonic"] = True
    table["i"]["Defining constraints"] = "Nucleotide F1 = 0, prediction boundaries within one of the introns of the reference transcript."

    table["ri"]["definition"] = "Reverse intronic transcript - reference is contained within the intron of a prediction transcript."
    table["ri"]["Pred_multiexonic"] = True
    table["ri"]["Ref_multiexonic"] = False
    table["ri"]["Defining constraints"] = "Nucleotide F1 = 0, reference boundaries within one of the introns of the prediction transcript."

    table["O"]["definition"] = "Generic match of a multiexonic prediction vs a monoexonic reference."
    table["O"]["Pred_multiexonic"] = False
    table["O"]["Ref_multiexonic"] = True
    table["O"]["Defining constraints"] = "Nucleotide F1 between 0 and 1."

    # Opposite strand
    table["x"]["definition"] = "Monoexonic match on the opposite strand."
    table["x"]["Pred_multiexonic"] = False
    table["x"]["Ref_multiexonic"] = "NA"
    table["x"]["Defining constraints"] = "Nucleotide F1 > 0, opposite strands."

    table["X"]["definition"] = "Monoexonic match on the opposite strand."
    table["X"]["Pred_multiexonic"] = True
    table["X"]["Ref_multiexonic"] = "NA"
    table["X"]["Defining constraints"] = "Nucleotide F1 > 0, opposite strands."

    # Polymerase run-ons
    table["p"]["definition"] = "The prediction is on the same strand of a neighboiring transcript."
    table["p"]["Pred_multiexonic"] = False
    table["p"]["Ref_multiexonic"] = "NA"
    table["p"]["Defining constraints"] = "Nucleotide F1 = 0%, the prediction and the reference transcripts coordinates do not overlap at all, they share the same strand. Probably polymerase run-on."

    # Polymerase run-ons
    table["p"]["definition"] = "The prediction is on the opposite strand of a neighboiring transcript."
    table["p"]["Pred_multiexonic"] = False
    table["p"]["Ref_multiexonic"] = "NA"
    table["p"]["Defining constraints"] = "Nucleotide F1 = 0%, the prediction and the reference transcripts coordinates do not overlap at all, they are located on opposite strands. Possibly a polymerase run-on, but it might be just an unannotated transcript."

    
def launch(args):

    """
    Main caller.

    :param args: the argparse Namespace.
    """

    metric_names = Transcript.get_available_metrics()
    print(file=args.out)
    metrics = ["tid", "parent", "score"]
    metrics.extend(
        sorted(metric for metric in metric_names
               if metric not in ("tid", "parent", "score")))

    for metric in metrics:
        docstring = getattr(Transcript, metric).__doc__
        if docstring is None:
            docstring = ''
        print("- *{0}*:".format(metric), re.sub(" +", " ", re.sub("\n", " ", docstring)),
              sep="\t", file=args.out)

    print(file=args.out)


def metric_parser():

    """
    Command line parser for the utility.
    """

    parser = argparse.ArgumentParser("Script to generate the available metrics")
    parser.add_argument("-o", "--out", type=argparse.FileType("w"), default=sys.stdout)
    parser.set_defaults(func=launch)
    return parser
