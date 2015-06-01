import sys,os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from sqlalchemy import Column,String,Integer,ForeignKey,CHAR,Index,Float,Boolean
from loci_objects import bed12
from sqlalchemy.orm import relationship, backref
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from loci_objects.dbutils import dbBase
from loci_objects.blast_utils import Query

class orf(dbBase):
    __tablename__="orf"
     
    id=Column(Integer, primary_key=True)
    query_id=Column(Integer, ForeignKey(Query.id), unique=False)
    start=Column(Integer, nullable=False)
    end=Column(Integer, nullable=False)
    name=Column(String)
    strand=Column(CHAR)
    thickStart=Column(Integer, nullable=False)
    thickEnd=Column(Integer, nullable=False)
    score=Column(Float)
    has_start_codon=Column(Boolean, nullable=True)
    has_stop_codon=Column(Boolean, nullable=True)
    cds_len = Column(Integer)
    
    __table_args__ = ( Index("index", "query_id", "thickStart", "thickEnd"  ), )
    
    query_object= relationship(Query, uselist=False, backref=backref("orfs"), lazy="joined")
    
    def __init__(self, bed12_object, query_id):
        if type(bed12_object) is not bed12.BED12:
            raise TypeError("Invalid data type!")
        self.query_id = query_id
        self.start=bed12_object.start
        self.end=bed12_object.end
        self.thickStart=bed12_object.thickStart
        self.thickEnd=bed12_object.thickEnd
        self.name=bed12_object.name
        self.strand=bed12_object.strand
        self.score = bed12_object.score
        self.has_start_codon= bed12_object.has_start_codon
        self.has_stop_codon= bed12_object.has_stop_codon
        self.cds_len = bed12_object.cds_len
                
    def __str__(self):
        return "{chrom}\t{start}\t{end}".format(
                                                chrom=self.query,
                                                start=self.start,
                                                end=self.end
                                                )

    @property
    def query(self):
        return self.query_object.name
        
        
class orfSerializer:
        
    def __init__(self, handle, db, fasta_index=None, dbtype="sqlite"):
        
        self.BED12 = bed12.bed12Parser(handle, fasta_index=fasta_index, transcriptomic=True)
        self.engine=create_engine("{dbtype}:///{db}".format(dbtype=dbtype,
                                                       db=db))
        session=sessionmaker()
        session.configure(bind=self.engine)
        dbBase.metadata.create_all(self.engine) #@UndefinedVariable
        self.session=session()
        
    def serialize(self):
        for row in self.BED12:
            if row.header is True:
                continue
            current_query = self.session.query(Query).filter(Query.name==row.id).all()
            if len(current_query) == 0:
                current_query=Query(row.id, None)
                self.session.add(current_query)
                self.session.commit()
            else:
                current_query=current_query[0]
            current_junction = orf( row, current_query.id)
            self.session.add(current_junction)
            self.session.commit()
        self.session.commit()
            
