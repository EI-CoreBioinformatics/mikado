#!/usr/bin/env python3

import sys, os, subprocess, gzip,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.blast_utils import Base,Query,Hit,Hsp,Target
from Bio.Blast.NCBIXML import parse as xparser
import io
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker

def main():

	def to_reader(string):
		if 'gz' in string:
			if ".asn.gz" in string:
				zcat=subprocess.Popen(["zcat", string], shell=False, stdout=subprocess.PIPE) #I cannot seem to make it work with gzip.open
				
				blast_formatter=subprocess.Popen( ['blast_formatter', '-outfmt', '5',
								   '-archive', '-'], shell=False,
								  stdin=zcat.stdout,
								  stdout=subprocess.PIPE)
				reader=io.TextIOWrapper(blast_formatter.stdout, encoding="UTF-8")
			else:
				reader=gzip.open(string, mode="rt")

		elif 'bz2' in string: reader=subprocess.Popen(['bzcat',string],stdout=subprocess.PIPE).stdout
		elif ".asn" in string:
			blast_formatter=subprocess.Popen( ['blast_formatter', '-outfmt', '5',
													'-archive', string], shell=False,
												stdout=subprocess.PIPE)
			reader=io.TextIOWrapper(blast_formatter.stdout, encoding="UTF-8")
		else:
			reader=open(string)
		return reader

	parser=argparse.ArgumentParser("Script to conver the BlastXML output into a tabular format.")
	parser.add_argument("--max_target_seqs", type=int, default=float("Inf"), help="Maximum number of target sequences.")
	parser.add_argument("--definition", action="store_true", default=False, help="Use query def instead of ID for the output.")
	parser.add_argument("xml", type=to_reader, help="XML file to parse.")
	parser.add_argument("dbout", type=str, default=":memory:",
			    nargs='?', help="Optional output file. Default: :memory:")

	args=parser.parse_args()

	engine=create_engine("sqlite:///{0}".format(args.dbout))
	session=sessionmaker()
	session.configure(bind=engine)
	Base.metadata.create_all(engine)
	current_session=session()

	q_mult=1 #Query multiplier: 3 for BLASTX (nucleotide vs. protein), 1 otherwise
	h_mult=1 #Target multiplier: 3 for TBLASTN (protein vs. nucleotide), 1 otherwise

	for record in xparser(args.xml):
		if record.application=="BLASTN":
			q_mult=1
			h_mult=1
		elif record.application=="BLASTX":
			q_mult=3
			h_mult=1
		elif record.application=="TBLASTN":
			q_mult=1
			h_mult=3
		if len(record.descriptions)==0:
			continue
		
		if args.definition==True:
			name=record.query.split()[0]
		else:
			name=record.query_id
		current_query=Query(name, record.query_length)
		current_session.add(current_query)
		current_session.commit()
		
		for ccc,alignment in filter(lambda x: x[0]<=args.max_target_seqs, enumerate(record.alignments)):

			evalue=record.descriptions[ccc].e
			bits = record.descriptions[ccc].bits
			alignment = record.alignments[ccc]
			hit_num=ccc+1
			current_target=current_session.query(Target).filter(Target.name==record.alignments[ccc].accession).all()
			if len(current_target)==0:
				current_target=Target(record.alignments[ccc].accession, record.alignments[ccc].length)
				current_session.add(current_target)
				current_session.commit()
			else:
				current_target=current_target[0]	

			current_hit = Hit(current_query.id, current_target.id, alignment, evalue, bits,
							hit_number=hit_num,
							query_multiplier=q_mult,
							target_multiplier=h_mult)
			current_session.add(current_hit)
			current_session.commit()

			for counter,hsp in enumerate(alignment.hsps):
				current_hsp = Hsp(hsp, counter, current_hit.id)
				current_session.add(current_hsp)
	current_session.commit() 

if __name__=='__main__': main()
