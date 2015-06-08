#!/usr/bin/env python3

import sys, os, argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from shanghai_lib.serializers.blast_utils import xmlSerializer
from Bio import SeqIO

def to_seqio(string):
	assert os.path.exists(string) and os.path.isfile(string) and os.stat(string).st_size>0
	return SeqIO.index(string, "fasta")


def main():

	parser=argparse.ArgumentParser("Script to conver the BlastXML output into a tabular format.")
	parser.add_argument("--max_target_seqs", type=int, default=float("Inf"), help="Maximum number of target sequences.")
	parser.add_argument("--maxobjects", type=int, default=10**5, help="Maximum number of objects to cache in memory.")
	parser.add_argument("--definition", action="store_true", default=False, help="Use query def instead of ID for the output.")
	parser.add_argument("--query_seqs", default=None, type=to_seqio, help="Query sequences")
	parser.add_argument("--target_seqs", default=None, type=to_seqio, help="Target sequences")
	parser.add_argument("xml", type=str, help="XML file to parse.")
	parser.add_argument("dbout", type=str, default=":memory:",
			    nargs='?', help="Optional output file. Default: :memory:")

	args=parser.parse_args()
	
	xmlSerializer(args.dbout, args.xml,
				keep_definition=args.definition,
				max_target_seqs=args.max_target_seqs,
				maxobjects=args.maxobjects,
				target_seqs=args.target_seqs,
				query_seqs=args.query_seqs
				)()


if __name__=='__main__': main()
