#!/usr/bin/env python3

import sys, os, subprocess, gzip,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.blast_utils import Base,Query,Hit,Hsp,Target,xmlSerializer
from Bio.Blast.NCBIXML import parse as xparser
import io
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker

def main():

	parser=argparse.ArgumentParser("Script to conver the BlastXML output into a tabular format.")
	parser.add_argument("--max_target_seqs", type=int, default=float("Inf"), help="Maximum number of target sequences.")
	parser.add_argument("--definition", action="store_true", default=False, help="Use query def instead of ID for the output.")
	parser.add_argument("xml", type=str, help="XML file to parse.")
	parser.add_argument("dbout", type=str, default=":memory:",
			    nargs='?', help="Optional output file. Default: :memory:")

	args=parser.parse_args()
	
	xmlSerializer(args.dbout, args.xml, keep_definition=args.definition, max_target_seqs=args.max_target_seqs).serialize()


if __name__=='__main__': main()
