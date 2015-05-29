#!/usr/bin/env python3

import sys, os, subprocess, gzip,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.blast_utils import Base,Query,Hit,Hsp,Target
from Bio.Blast.NCBIXML import parse as xparser
import io
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker

def main():

	parser=argparse.ArgumentParser("Script to conver the BlastXML output into a tabular format.")
	parser.add_argument("--max_target_seqs", type=int, default=float("Inf"), help="Maximum number of target sequences.")
	parser.add_argument("--definition", action="store_true", default=False, help="Use query def instead of ID for the output.")
	parser.add_argument("db", type=str, help="DB file to parse.")
	parser.add_argument("out", type=argparse.FileType("w"), default=sys.stdout,
			    nargs='?', help="Optional output file. Default: %(default)s")

	args=parser.parse_args()

	engine=create_engine("sqlite:///{0}".format(args.db))
	session=sessionmaker()
	session.configure(bind=engine)
	Base.metadata.create_all(engine)
	current_session=session()

	for hit in current_session.query(Hit):
		print(hit)
		for hsp in hit.hsps:
			print("\t", hsp)


if __name__=='__main__': main()
