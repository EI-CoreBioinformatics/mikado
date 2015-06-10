import sys, os
import sqlalchemy
import gzip
import subprocess
from Bio import SeqIO
import Bio.File
# import collections
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from scipy import mean
import operator
from sqlalchemy import Column,String,Integer,Float,ForeignKey
from sqlalchemy.sql.schema import PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy.orm import relationship, backref
from Bio.Blast.NCBIXML import parse as xparser
import io
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker
from shanghai_lib.serializers.dbutils import dbBase
import logging

'''This module is used to serialise BLAST objects into a database.'''


#For profiling
if '__builtin__' not in dir() or not hasattr(__builtin__, "profile") or "profile" not in dir(): #@UndefinedVariable
	def profile(function):
		def inner(*args, **kwargs):
			return function(*args, **kwargs)
		return inner


def merge(intervals):
	
	'''This function is used to merge together intervals, which have to be supplied as a list
	of duplexes - (start,stop). The function will then merge together overlapping tuples and
	return a list of non-overlapping tuples.
	If the list is composed by only one element, the function returns immediately. 
	'''
	
	#Assume tuple of the form (start,end)
	#And return 0- and 1-length intervals
	new_intervals=[]
	for interval in intervals:
		new_intervals.append(tuple(sorted(interval)))

	intervals=new_intervals[:]
	if len(intervals)<2: return intervals

	#Sort according to start, end
	intervals=sorted(intervals, key=operator.itemgetter(0,1))
#	print(intervals, file=sys.stderr)
	final_list=[intervals[0]]
	
	for start,end in intervals[1:]:
		if start>final_list[-1][1]:
			final_list.append(tuple([start,end]))
		elif end>final_list[-1][1]:
			final_list[-1]=tuple([final_list[-1][0], end])
#	print(final_list, file=sys.stderr)
	return final_list


class Query(dbBase):
	
	__tablename__="query"
	id=Column(Integer, primary_key=True)
	name=Column(String(200), unique=True, index=True)
	length=Column(Integer, nullable=True) #This so we can load data also from the orf class
	
	@profile
	def __init__(self, name, length):
		self.name=name
		self.length=length

class Target(dbBase):
	__tablename__="target"

	id=Column(Integer, primary_key=True)
	name=Column(String(200), unique=True, index=True)
	length=Column(Integer)
	
	@profile
	def __init__(self, name, length):
		self.name=name
		self.length=length

	
class Hit(dbBase):
	
	'''This class is used to serialise and store in a DB a BLAST hit.
	Stored attributes:
	
	- id				Indexing key
	- query_id			Foreign ID key for the query table
	- target_id			Foreign ID key for the target table
	- qt_constrating	
	
	'''
	
	__tablename__ ="hit"
	query_id=Column(Integer, ForeignKey(Query.id), unique=False)
	target_id=Column(Integer, ForeignKey(Target.id), unique=False)
	qt_constraint = PrimaryKeyConstraint("query_id", "target_id", name="hit_id")
	evalue = Column(Float)
	bits = Column(Float)
	global_identity = Column(Float)
	query_start = Column(Integer)
	query_end = Column(Integer)
	target_start = Column(Integer)
	target_end = Column(Integer)
	hit_number = Column(Integer)
	query_multiplier = Column(Float) # Probably I should move this to a separate table!
	target_multiplier = Column(Float)
	query_aligned_length = Column(Integer)
	target_aligned_length = Column(Integer)
	
	query_object = relationship(Query, uselist=False, lazy="joined", backref=backref("hits" ))
	target_object = relationship(Target, uselist=False, lazy="joined", backref=backref("hits" ) )
	
	__table_args__ = (qt_constraint,)
	
	def __init__(self, query_id, target_id, alignment, evalue, bits, hit_number=1, query_multiplier=1, target_multiplier=1):
		'''This function takes as input the id of a target, the id of the query, and a hit-object from the XML.
		The multiplier keyword is used to calculate the ratio between the query and the target.'''
		
		self.query_id = query_id
		self.target_id = target_id
		self.query_multiplier = query_multiplier
		self.target_multiplier = target_multiplier
		self.hit_number = hit_number
		self.evalue=evalue
		self.bits = bits
		
		self.global_identity = mean([ hsp.identities/hsp.align_length*100 for hsp in alignment.hsps  ]) 
		q_intervals = [tuple([hsp.query_start,hsp.query_end]) for hsp in alignment.hsps]
		q_merged_intervals = sorted(merge(q_intervals), key=operator.itemgetter(0,1))
		q_aligned = sum([tup[1]-tup[0]+1 for tup in q_merged_intervals])
		self.query_aligned_length = q_aligned
		self.query_start = q_merged_intervals[0][0]
		self.query_end = q_merged_intervals[-1][1]
		
		t_intervals = [tuple([hsp.sbjct_start,hsp.sbjct_end]) for hsp in alignment.hsps]
		t_merged_intervals = sorted(merge(t_intervals), key=operator.itemgetter(0,1))
		t_aligned = sum([tup[1]-tup[0]+1 for tup in t_merged_intervals])
		self.target_aligned_length = t_aligned
		self.target_start = t_merged_intervals[0][0]
		self.target_end = t_merged_intervals[-1][1]
	
	def __str__(self):
		line=[]
		line.append(self.query)
		line.append(self.target)
		line.append(self.evalue)
		line.append(self.bits)
		line.append(self.query_start)
		line.append(self.query_end)
		line.append(self.target_start)
		line.append(self.target_end)
		line.append(self.query_len)
		line.append(self.target_len)
		
		return "\t".join(str(x) for x in line)


	def as_dict(self):
		'''Method to return a dict representation of the object.
		Necessary for storing.'''
		
		keys = [
			"evalue",
			"bits",
			"global_identity",
			"query_start",
			"query_end",
			"target_start",
			"target_end",
			"hit_number",
			"query_multiplier",
			"target_multiplier",
			"query_aligned_length",
			"target_aligned_length",
			]

		special_keys = ["query", "target", "hsps","query_len",
			"target_len", "query_hit_ratio", "hit_query_ratio" ] #Keys we want to set directly

		state = dict().fromkeys(keys+special_keys)

		for key in keys:
			state[key]=getattr(self, key)
		
		state["query"]=self.query
		state["target"]=self.target
		state["query_len"]=self.query_len
		state["target_len"]=self.target_len
		state["query_hit_ratio"]=self.query_hit_ratio
		state["hit_query_ratio"]=self.hit_query_ratio
		
		state["hsps"]=[]
		for hsp in self.hsps:
			state["hsps"].append(hsp.as_dict())
			
		return state

	@property
	def query_len(self):
		return self.query_object.length
	@property
	def query(self):
		return self.query_object.name
	
	@property
	def target(self):
		return self.target_object.name
	@property
	def target_len(self):
		return self.target_object.length

	@property
	def query_hit_ratio(self):
		return self.query_len*self.query_multiplier/(self.target_len*self.target_multiplier)

	@property
	def hit_query_ratio(self):
		return self.target_len*self.target_multiplier/(self.query_len*self.query_multiplier)
	

class Hsp(dbBase):

	'''This class serializes and stores into the DB the various HSPs.
	It is directly connected to the Hit table, through the "hit_id" 
	reference key.
	The Hit reference can be accessed through the hit_object attribute;
	back-reference (Hit to Hsps) is given by the "hsps" attribute.
	
	Keys:
	- hit_id 				Reference for the Hit table
	- counter				Indicates the number of the HSP for the hit
	- query_hsp_start		Start position on the query
	- query_hsp_end			End position on the query
	- target_hsp_start		Start position on the target
	- target_hsp_end		End position on the target
	- hsp_evalue			Evalue of the HSP
	- hsp_bits				Bit-score of the HSP
	- hsp_identity			Identity (in %) of the alignment
	- hsp_length			Length of the HSP
	
	An HSP row has the following constraints:
	- Counter,hit_id must be unique (and are primary keys)
	- The combination ("Hit_id","query_hsp_start","query_hsp_end", "target_hsp_start", "target_hsp_end") must be unique
	'''

	__tablename__ = "hsp"
	counter = Column(Integer) #Indicates the number of the HSP inside the hit
	query_id=Column(Integer, ForeignKey(Query.id), unique=False)
	target_id=Column(Integer, ForeignKey(Target.id), unique=False)
	pk_constraint = PrimaryKeyConstraint("counter", "query_id", "target_id", name="hsp_constraint")
	query_hsp_start = Column(Integer)
	query_hsp_end = Column(Integer)
	target_hsp_start = Column(Integer)
	target_hsp_end = Column(Integer)
	uni_constraint = UniqueConstraint("query_id", "target_id", "query_hsp_start", "query_hsp_end", "target_hsp_start", "target_hsp_end")
	hsp_evalue = Column(Float)
	hsp_bits = Column(Float)
	hsp_identity = Column(Float)
	hsp_length = Column(Integer)

	query_object=relationship(Query, uselist=False, lazy="joined")
	target_object=relationship(Target, uselist=False, lazy="joined")

	hit_object=relationship(Hit, uselist=False, lazy="joined", backref=backref("hsps"),
						foreign_keys=[query_id, target_id],
						primaryjoin="and_(Hit.query_id==Hsp.query_id, Hit.target_id==Hsp.target_id)")
	
	__table_args__ = (pk_constraint,)
	
	def __init__(self, hsp, counter, query_id, target_id):
		self.counter = counter
		self.query_hsp_start = hsp.query_start
		self.query_hsp_end = hsp.query_end
		self.target_hsp_start = hsp.sbjct_start
		self.target_hsp_end = hsp.sbjct_end
		self.hsp_identity = float(hsp.identities)/hsp.align_length*100
		self.hsp_length = hsp.align_length
		self.hsp_bits = hsp.bits
		self.hsp_evalue = hsp.expect
		self.query_id = query_id
		self.target_id = target_id
	
	def __str__(self):
		'''Simple printing function.'''
		line=[]
		line.append(self.query)
		line.append(self.target)
		line.append(self.query_hsp_start)
		line.append(self.query_hsp_end)
		line.append(self.target_hsp_start)
		line.append(self.target_hsp_end)
		line.append(self.hsp_evalue)
		return "\t".join([str(x) for x in line])
		
		
	def as_dict(self):
		
		'''Method to return a dict representation of the object. Necessary for storing.'''
		
		keys=[
			"query_hsp_start",
			"query_hsp_end",
			"target_hsp_start",
			"target_hsp_end",
			"hsp_evalue",
			]
		
		special_keys=["query","target","query_hsp_cov","target_hsp_cov"] #Keys which we want to assign directly rather than inside a loop
		
		state=dict().fromkeys(keys+special_keys)
		state["query"]=self.query
		state["target"]=self.target
		state["query_hsp_cov"]=self.query_hsp_cov
		state["target_hsp_cov"]=self.target_hsp_cov
		
		for key in keys:
			state[key]=getattr(self, key)
		
		return state

	@property
	def query(self):
		'''Returns the name of the query sequence, through a nested SQL query.'''
		return self.query_object.name

	@property
	def target(self):
		'''Returns the name of the target sequence, through a nested SQL query.'''
		return self.target_object.name

	
	@property
	def query_hsp_cov(self):
		'''This property returns the percentage of the query which is covered by the HSP.'''
		return (self.query_hsp_end-self.query_hsp_start+1)/(self.query_object.length)
	
	@property
	def target_hsp_cov(self):
		'''This property returns the percentage of the target which is covered by the HSP.'''
		return (self.target_hsp_end-self.target_hsp_start+1)/(self.target_object.length)

class xmlSerializer:
	
	'''This class has the role of taking in input a blast XML file and (partially) serialise it into
	a database. We are using SQLalchemy, so the database type could be any of SQLite, MySQL, PSQL, etc.'''
	
	def __init__(self, db, xml, max_target_seqs=float("Inf"),
				target_seqs=None,
				query_seqs=None,
				keep_definition=False, maxobjects=10000  ):
		'''Initializing method. Arguments:
		-db
		-xml (it can be a handle or a valid file address)
		
		Optional arguments:
		- max_target_seqs 	Maximum number of targets to serialize
		- keep_definition	Keep the definition instead of the ID for queries
		
		'''
		
		self.logger = logging.getLogger("main")
		self.logger.setLevel(logging.INFO)
		self.handler = logging.StreamHandler()
		self.formatter = logging.Formatter("{asctime} - {name} - {message}", style='{')
		self.handler.setFormatter(self.formatter)
		self.logger.addHandler(self.handler)
		
		engine=create_engine("sqlite:///{0}".format(db))
		session=sessionmaker()
		session.configure(bind=engine)
		dbBase.metadata.create_all(engine) #@UndefinedVariable
		self.session=session()
		self.logger.info("Created the session")
		if type(xml) is not xparser:
			if type(xml) is str:
				if not os.path.exists(xml):
					raise OSError("Invalid file address.")
				if xml.endswith("xml"):
					xml=open(xml)
				elif xml.endswith(".xml.gz"):
					xml=gzip.open(xml, mode="rt") #The mode must be textual, otherwise xparser dies 
				elif xml.endswith(".asn.gz"):
					zcat=subprocess.Popen(["zcat", xml], shell=False, stdout=subprocess.PIPE) #I cannot seem to make it work with gzip.open
					blast_formatter=subprocess.Popen( ['blast_formatter', '-outfmt', '5',
								   		'-archive', '-'], shell=False,
								  			stdin=zcat.stdout,
								  			stdout=subprocess.PIPE)
					xml=io.TextIOWrapper(blast_formatter.stdout, encoding="UTF-8")
			assert type(xml) in (gzip.GzipFile,io.TextIOWrapper), type(xml)
			self.xml=xml
			self.xml_parser=xparser(xml)
		else:
			self.xml_parser=xml # This is the BLAST object we will serialize
		#Runtime arguments
		self.keep_definition=keep_definition
		
		if type(query_seqs) is str:
			assert os.path.exists(query_seqs)
			self.query_seqs=SeqIO.index(query_seqs,"fasta")
		elif query_seqs is None:
			self.query_seqs=None
		else:
			assert type(query_seqs) is Bio.File._IndexedSeqFileDict
			self.query_seqs=query_seqs
		
		if type(target_seqs) is str:
			assert os.path.exists(target_seqs)
			self.target_seqs=SeqIO.index(target_seqs,"fasta")
		elif target_seqs is None:
			self.target_seqs = None
		else:
			assert type(target_seqs) is Bio.File._IndexedSeqFileDict
			self.target_seqs=target_seqs
		
		self.max_target_seqs=max_target_seqs
		self.maxobjects = maxobjects 
		
	def serialize(self):
		
		'''Method to serialize the BLAST XML file into a database provided with the __init__ method '''
		
		q_mult=1 #Query multiplier: 3 for BLASTX (nucleotide vs. protein), 1 otherwise
		h_mult=1 #Target multiplier: 3 for TBLASTN (protein vs. nucleotide), 1 otherwise

		objects=[]
		
		targets = dict()
		queries = dict()
		self.logger.info("Loading previous IDs")
		for query in self.session.query(Query):
			queries[query.name] = query.id
		for query in self.session.query(Target):
			targets[query.name] = query.id
		self.logger.info("Loaded previous IDs")
		
		self.logger.info("Started the serialization")
		if self.target_seqs is not None:
			self.logger.info("Started to serialize the targets")
			for record in self.target_seqs:
				if record in targets:
					continue
				objects.append(Target(record, len(self.target_seqs[record])  ))
				if len(objects)>=self.maxobjects:
						self.logger.info("Loading {0} objects into the \"target\" table".format(self.maxobjects))
						self.session.bulk_save_objects(objects, return_defaults=False)
						self.logger.info("Loaded {0} objects into the \"target\" table".format(self.maxobjects))
						objects=[]
			self.logger.info("Loading {0} objects into the \"target\" table".format(self.maxobjects))
			self.session.bulk_save_objects(objects)
			self.logger.info("Loaded {0} objects into the \"target\" table".format(self.maxobjects))
			objects=[]
			self.logger.info("Loaded targets")

		
		if self.query_seqs is not None:
			self.logger.info("Started to serialize the queries")
			for record in self.query_seqs:
				if record in queries:
					continue
				objects.append(Query(record, len(self.query_seqs[record])  ))
				if len(objects)>=self.maxobjects:
					self.logger.info("Loading {0} objects into the \"query\" table".format(self.maxobjects))
					self.session.bulk_save_objects(objects, return_defaults=False)
					self.logger.info("Loaded {0} objects into the \"query\" table".format(self.maxobjects))
					objects=[]
			self.logger.info("Loading {0} objects into the \"query\" table".format(self.maxobjects))
			self.session.bulk_save_objects(objects)
			self.logger.info("Loaded {0} objects into the \"query\" table".format(self.maxobjects))
			objects=[]
			self.logger.info("Queries serialized")

		self.logger.info("Loading all IDs")
		for query in self.session.query(Query):
			queries[query.name] = (query.id, query.length is not None)
		for query in self.session.query(Target):
			targets[query.name] = (query.id, query.length is not None)
		self.logger.info("Loaded all IDs")

		#Memorize the mapping ... it is just faster

		for record in self.xml_parser:
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
			
			if self.keep_definition is True:
				name=record.query.split()[0]
			else:
				name=record.query_id
			self.logger.debug("Started with {0}".format(name))

			if name in queries:
				current_query = queries[name][0]
				if queries[name][1] is False:
					self.session.query(Query).filter(Query.name==name).update({"length": record.query_length})
					self.session.commit()
			else:
				self.logger.info("Adding {0} to the db".format(name))
				current_query=Query(name, record.query_length)
				self.session.add(current_query)
				self.session.commit()
				current_query = current_query.id
				queries[name] = (current_query.id, True)
					
		
			for ccc,alignment in filter(lambda x: x[0]<=self.max_target_seqs, enumerate(record.alignments)):

				self.logger.debug("Started the hit {0}-{1}".format(name,record.alignments[ccc].accession ))
				evalue=record.descriptions[ccc].e
				bits = record.descriptions[ccc].bits
				alignment = record.alignments[ccc]
				hit_num=ccc+1
				if record.alignments[ccc].accession in targets:
					current_target = targets[record.alignments[ccc].accession][0]
					if targets[record.alignments[ccc].accession][1] is False:
						self.session.query(Target).filter(Target.name==record.alignments[ccc].accession).update(
																											{"length": record.query_length}
																											)
					
				else:
					current_target=Target(record.alignments[ccc].accession, record.alignments[ccc].length)
					self.session.add(current_target)
					try:
						self.session.commit()
						assert type(current_target.id) is int 
						targets[record.alignments[ccc].accession] = (current_target.id, True)
						current_target = current_target.id
					except sqlalchemy.exc.IntegrityError:
						self.session.rollback()
						continue

				current_hit = Hit(current_query, current_target, alignment, evalue, bits,
									hit_number=hit_num,
									query_multiplier=q_mult,
									target_multiplier=h_mult)
				objects.append(current_hit)
						
				for counter,hsp in enumerate(alignment.hsps):
					current_hsp = Hsp(hsp, counter, current_query, current_target)
					objects.append(current_hsp)
			
			if len(objects)>=self.maxobjects:
				self.logger.info("Loading {0} objects into the hit, hsp tables".format(len(objects)))
				self.session.bulk_save_objects(objects, return_defaults=False)
				self.logger.info("Loaded {0} objects into the hit, hsp tables".format(len(objects)))
				objects=[]
			
		self.logger.info("Loading {0} objects into the hit, hsp tables".format(len(objects)))
		self.session.bulk_save_objects(objects, return_defaults=False)
		self.logger.info("Loaded {0} objects into the hit, hsp tables".format(len(objects)))
		self.session.commit() 
		self.logger.info("Finished loading blast hits")

	def __call__(self):
		self.serialize()
