import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from scipy import mean
import operator
from sqlalchemy import Column,String,Integer,Float,ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql.schema import PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy.orm import relationship, backref
# from sqlalchemy.engine import create_engine
# from sqlalchemy.orm.session import sessionmaker


#print("#", *sys.argv)
Base=declarative_base()

def merge(intervals):
	
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


class Query(Base):
	
	__tablename__="query"
	id=Column(Integer, primary_key=True)
	name=Column(String(200), unique=True, index=True)
	length=Column(Integer)
	
	def __init__(self, name, length):
		self.name=name
		self.length=length

class Target(Base):
	__tablename__="target"

	id=Column(Integer, primary_key=True)
	name=Column(String(200), unique=True, index=True)
	length=Column(Integer)
	
	def __init__(self, name, length):
		self.name=name
		self.length=length

	
class Hit(Base):
	
	__tablename__ ="hit"
	id=Column(Integer, primary_key=True, autoincrement=True)
	query_id=Column(Integer, ForeignKey(Query.id), unique=False)
	target_id=Column(Integer, ForeignKey(Target.id), unique=False)
	qt_constraint = UniqueConstraint("query_id", "target_id", name="hit")
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
	
	query_object = relationship(Query, uselist=False, lazy="joined" )
	target_object = relationship(Target, uselist=False, lazy="joined" )
	
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
	

class Hsp(Base):

	__tablename__ = "hsp"
	counter = Column(Integer) #Indicates the number of the HSP inside the hit
	hit_id = Column(Integer, ForeignKey(Hit.id)) #Foreign key to recover the hit
	pk_constraint = PrimaryKeyConstraint("counter", "hit_id", name="hsp_constraint")
	query_hsp_start = Column(Integer)
	query_hsp_end = Column(Integer)
	target_hsp_start = Column(Integer)
	target_hsp_end = Column(Integer)
	uni_constraint = UniqueConstraint("hit_id", "query_hsp_start", "query_hsp_end", "target_hsp_start", "target_hsp_end")
	hsp_evalue = Column(Float)
	hsp_bits = Column(Float)
	hsp_identity = Column(Float)
	hsp_length = Column(Integer)

	hit_object=relationship(Hit, uselist=False, lazy="joined", backref=backref("hsps"))
	
	__table_args__ = (pk_constraint,)
	
	def __init__(self, hsp, counter, hit_id):
		self.counter = counter
		self.hit_id = hit_id
		assert self.hit_id is not None
		self.query_hsp_start = hsp.query_start
		self.query_hsp_end = hsp.query_end
		self.target_hsp_start = hsp.sbjct_start
		self.target_hsp_end = hsp.sbjct_end
		self.hsp_identity = float(hsp.identities)/hsp.align_length*100
		self.hsp_length = hsp.align_length
		self.hsp_bits = hsp.bits
		self.hsp_evalue = hsp.expect
	
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
		
			
	def csv_row(self):
		'''This function will return a CSV-dict ready for printing.'''

	@property
	def query(self):
		return self.hit_object.query

	@property
	def target(self):
		return self.hit_object.target

	
	@property
	def query_hsp_cov(self):
		'''This property returns the percentage of the query which is covered by the HSP.'''
		return (self.query_hsp_end-self.query_hsp_start+1)/(self.hit_object.query_len)
	
	@property
	def target_hsp_cov(self):
		'''This property returns the percentage of the target which is covered by the HSP.'''
		return (self.target_hsp_end-self.target_hsp_start+1)/(self.hit_object.target_len)
