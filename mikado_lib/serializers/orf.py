import sys,os
from Bio import SeqIO
import Bio.File
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from sqlalchemy import Column,String,Integer,ForeignKey,CHAR,Index,Float,Boolean
from mikado_lib.parsers import bed12
from sqlalchemy.orm import relationship, backref
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from mikado_lib.serializers.dbutils import dbBase,Inspector
from mikado_lib.serializers.blast_utils import Query

class orf(dbBase):
    __tablename__="orf"
     
    id=Column(Integer, primary_key=True)
    query_id=Column(Integer, ForeignKey(Query.id), unique=False)
    start=Column(Integer, nullable=False)
    end=Column(Integer, nullable=False)
    name=Column(String(200))
    strand=Column(CHAR)
    thickStart=Column(Integer, nullable=False)
    thickEnd=Column(Integer, nullable=False)
    score=Column(Float)
    has_start_codon=Column(Boolean, nullable=True)
    has_stop_codon=Column(Boolean, nullable=True)
    cds_len = Column(Integer)
    
    __table_args__ = ( Index("orf_index", "query_id", "thickStart", "thickEnd"  ), Index("query_index", "query_id") )
    
    query_object= relationship(Query, uselist=False, backref=backref("orfs"), lazy="immediate")
    
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


    def as_bed12(self):
        '''Method to transform the mapper into a BED12 object.'''
        
        b=bed12.BED12()
        
        b.header=False
        b.query = b.chrom = self.query
        b.start=self.start
        b.end = self.end
        b.name = self.name
        b.score = self.score
        b.strand = self.strand
        b.thickStart = self.thickStart
        b.thickEnd = self.thickEnd
        b.rgb=0
        b.blockCount=1
        b.blockSizes=[self.end]
        b.blockStarts=[0]
        
        b.has_start_codon = self.has_start_codon
        b.has_stop_codon = self.has_stop_codon

        return b

    @property
    def query(self):
        return self.query_object.name
        
        
class orfSerializer:
        
    def __init__(self, handle, db, fasta_index=None, maxobjects=1000000, json_conf=None):
        
        '''Constructor function. Arguments:
        - handle         the BED12 file
        - db             Output DB
        - fasta_index    A SeqIO-like index of sequence records. Alternatively, the path to the FASTA file. REQUIRED.
        - maxobjects    Integer. Indicates how big should the cache be for objects to be loaded inside the DB
        
        It is HIGHLY RECOMMENDED to provide the fasta index, as it will make the population of the Query
        table much faster.
        '''

        if type(fasta_index) is str:
            assert os.path.exists(fasta_index)
            self.fasta_index=SeqIO.index(fasta_index,"fasta")
        elif fasta_index is None:
            raise ValueError("A fasta index is needed for the serialization!")
        else:
            assert type(fasta_index) is Bio.File._IndexedSeqFileDict
            self.fasta_index=fasta_index

        
        self.BED12 = bed12.bed12Parser(handle, fasta_index=fasta_index, transcriptomic=True)
        if json_conf is not None:
            if json_conf["dbtype"]=="sqlite":
                self.engine=create_engine("sqlite:///{0}".format(json_conf["db"]))
            else:
                self.engine=create_engine("{dbtype}://{dbuser}:{dbpasswd}@{dbhost}/{db}".format(
                                                                                            dbtype=json_conf["dbtype"],
                                                                                            dbuser=json_conf["dbuser"],
                                                                                            dbpasswd=json_conf["dbpasswd"],
                                                                                            dbhost=json_conf["dbhost"],
                                                                                            db=json_conf["db"]))
        else:
            self.engine=create_engine("sqlite:///{0}".format(db))

        session=sessionmaker()
        session.configure(bind=self.engine)
        
        inspector=Inspector.from_engine(self.engine)
        if not orf.__tablename__ in inspector.get_table_names():
            dbBase.metadata.create_all(self.engine) #@UndefinedVariable
        self.session=session()
        self.maxobjects=maxobjects
        
    def serialize(self):
        objects = []
        cache=dict() #Dictionary to hold the data before bulk loading into the database
        
        for record in self.session.query(Query):
            cache[record.name]=record.id
        
        if self.fasta_index is not None:
            for record in self.fasta_index:
                if record in cache: continue
                objects.append(Query(record, len(self.fasta_index[record])))
                if len(objects)>=self.maxobjects:
                    self.session.bulk_save_objects(objects)
                    objects=[]
            
            self.session.bulk_save_objects(objects)
            self.session.commit()
            objects=[]
            
        for record in self.session.query(Query):
            cache[record.name]=record.id
            
        for row in self.BED12:
            if row.header is True or row.invalid is True:
                continue
            if row.id in cache:
                current_query = cache[row.id]
            else:
                current_query=Query(row.id, row.end)
                self.session.add(current_query)
                self.session.commit()
                cache[current_query.name] = current_query.id
                current_query = current_query.id
                
            current_junction = orf( row, current_query)
            objects.append(current_junction)
            if len(objects)>=self.maxobjects:
                self.session.bulk_save_objects(objects)
                objects=[]

        self.session.bulk_save_objects(objects)
        self.session.commit()
    
    def __call__(self):
        self.serialize()
