import sys,os
import io
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from sqlalchemy import Column,String,Integer,ForeignKey,CHAR,Index,Float
from shanghai_lib.parsers import bed12
from sqlalchemy.orm import relationship, backref
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from shanghai_lib.serializers.dbutils import dbBase,Inspector

class Chrom(dbBase):
    
    __tablename__="chrom"
    id=Column(Integer, primary_key=True)
    name=Column(String(200))
    length=Column(Integer, nullable=True)
    
    def __init__(self,name, length=None):
        self.name=name
        if length is not None:
            assert type(length) is int
        self.length=length


class junction(dbBase):
    __tablename__="junctions"
     
    id=Column(Integer, primary_key=True)
    chrom_id=Column(Integer, ForeignKey(Chrom.id), unique=False)
    start=Column(Integer, nullable=False)
    end=Column(Integer, nullable=False)
    name=Column(String)
    strand=Column(CHAR)
    junctionStart=Column(Integer, nullable=False)
    junctionEnd=Column(Integer, nullable=False)
    score=Column(Float)
    __table_args__ = ( Index("junction_index", "chrom_id", "junctionStart", "junctionEnd"  ), )
    
    chrom_object= relationship(Chrom, uselist=False, backref=backref("junctions"), lazy="joined")
    
    def __init__(self, bed12_object, chrom_id):
        if type(bed12_object) is not bed12.BED12:
            raise TypeError("Invalid data type!")
        self.chrom_id = chrom_id
        self.start=bed12_object.start
        self.end=bed12_object.end
        self.junctionStart=bed12_object.thickStart
        self.junctionEnd=bed12_object.thickEnd
        self.name=bed12_object.name
        self.strand=bed12_object.strand
        self.score = bed12_object.score
        
    def __str__(self):
        return "{chrom}\t{start}\t{end}".format(
                                                chrom=self.chrom,
                                                start=self.start,
                                                end=self.end
                                                )

    @property
    def chrom(self):
        return self.chrom_object.name
        
        
class junctionSerializer:
        
    def __init__(self, handle, db, fai=None, dbtype="sqlite", maxobjects=10000):
        
        self.BED12 = bed12.bed12Parser(handle)
        self.engine=create_engine("{dbtype}:///{db}".format(dbtype=dbtype,
                                                       db=db))
        session=sessionmaker()
        session.configure(bind=self.engine)

        inspector=Inspector.from_engine(self.engine)
        if not junction.__tablename__ in inspector.get_table_names():
            dbBase.metadata.create_all(self.engine) #@UndefinedVariable

        self.session=session()
        self.maxobjects=maxobjects
        self.fai = fai
        
 
    def serialize(self):
        
        sequences = dict()
        if self.fai is not None:
            if type(self.fai) is str:
                assert os.path.exists(self.fai)
                fai=open(self.fai)
            else:
                assert type(fai) is io.TextIOWrapper
            for line in fai:
                name, length = line.rstrip().split()[:2]
                current_chrom = Chrom(name, length=int(length))
                self.session.add(current_chrom)
                sequences[current_chrom.name] = current_chrom.id
            self.session.commit()
            for query in self.session.query(Chrom):
                sequences[query.name] = query.id
        
        objects = []
        
        for row in self.BED12:
            if row.header is True:
                continue
            if row.chrom in sequences:
                current_chrom = sequences[row.chrom]
            else:
                current_chrom=Chrom(row.chrom)
                self.session.add(current_chrom)
                self.session.commit()
                sequences[current_chrom.name] = current_chrom.id
                current_chrom = current_chrom.id
                
            current_junction = junction( row, current_chrom)
            objects.append(current_junction)
            if len(objects)>=self.maxobjects:
                self.session.bulk_save_objects(objects)
                objects=[]
        self.session.bulk_save_objects(objects)
        self.session.commit()
            
    def __call__(self):
        self.serialize()

