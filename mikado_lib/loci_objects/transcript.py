import operator
import os.path,sys,re
from collections import OrderedDict
import inspect
import asyncio
from mikado_lib.exceptions import InvalidTranscript
# import logging
# import pickle
# import copy
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

#SQLAlchemy imports
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.sql.expression import desc, asc
from sqlalchemy import and_
from sqlalchemy.ext import baked
from sqlalchemy import bindparam

#mikado imports
from mikado_lib.serializers.junction import junction, Chrom
import mikado_lib.serializers.orf
from mikado_lib.serializers.blast_utils import Query, Hit
from mikado_lib.serializers.orf import orf
from mikado_lib.parsers import bed12
from mikado_lib.loci_objects.abstractlocus import abstractlocus # Needed for the BronKerbosch algorithm ...
from mikado_lib.parsers.GTF import gtfLine
from mikado_lib.parsers.GFF import gffLine
import mikado_lib.exceptions

#from memory_profiler import profile
#import logging

# if "line_profiler" not in dir(): #@UndefinedVariable
#     def profile(function):
#         def inner(*args, **kwargs):
#             return function(*args, **kwargs)
#         return inner

class metric(property):
    '''Simple aliasing of property. All transcript metrics should use this alias, not "property", as a decorator.'''
    pass

class transcript:
    
    __name__ = "transcript"
    
    '''This class defines a transcript, down to its exon/CDS/UTR components. It is instantiated by a transcript
    GT/GFF3 line.
    Key attributes:
    - chrom    The chromosome
    - source
    - feature            mRNA if at least one CDS is defined
    - start
    - end
    - score
    - strand            one of +,-,None
    - phase            set to None
    - id            the ID of the transcripts (or tid)
    - parent        
    - attributes    a dictionary with additional informations from the GFFline
    
    After all exons have been loaded into the instance (see "addExon"), the class must be finalized with the appropriate method.

    CDS locations can be uploaded from the external, using a dictionary of indexed BED12 entries. 
    
    The database queries are baked at the *class* level in order to minimize overhead.
    
    '''
    
    #Query baking to minimize overhead
    bakery = baked.bakery()   
    query_baked=bakery( lambda session: session.query(Query) )
    query_baked += lambda q: q.filter(Query.name == bindparam("query_name"))

    blast_baked = bakery( lambda session: session.query(Hit) )
    blast_baked += lambda q: q.filter(and_(Hit.query_id == bindparam("query_id"),
                                                               Hit.evalue<=bindparam("evalue")),
                                                          )
        
    blast_baked += lambda q: q.order_by(asc(Hit.evalue))
    blast_baked += lambda q: q.limit(bindparam("max_target_seqs"))

    orf_baked = bakery(lambda session: session.query(mikado_lib.serializers.orf.orf))
    orf_baked += lambda q: q.filter(mikado_lib.serializers.orf.orf.query_id==bindparam("query_id"))
    orf_baked += lambda q: q.order_by(desc(mikado_lib.serializers.orf.orf.cds_len))
    
    ######### Class special methods ####################
    
    def __init__(self, *args, source=None, logger=None):
        
        '''Initialise the transcript object, using a mRNA/transcript line.
        Note: I am assuming that the input line is an object from my own "GFF" class.
        The transcript instance must be initialised by a "(m|r|lnc|whatever)RNA" or "transcript" gffLine.'''
        
        self.chrom = None
        self.exons, self.combined_cds, self.combined_utr = [], [], []
        self.set_logger(logger) #Set it to None by default
        self.introns = []
        self.splices = []
        self.finalized = False # Flag. We do not want to repeat the finalising more than once.        
        self.selected_internal_orf_index = None
        self.has_start_codon, self.has_stop_codon = False,False
        self.non_overlapping_cds = None
        self.verified_introns=set()
        
        if len(args)==0:
            return
        else:
            transcript_row=args[0]
            if type(transcript_row) not in (gffLine, gtfLine):
                raise TypeError("Invalid data type: {0}".format(type(transcript_row))) 
                    
        self.chrom = transcript_row.chrom
        assert transcript_row.is_transcript is True
        self.feature=transcript_row.feature
        self.id = transcript_row.id
        self.name = transcript_row.name
        if source is None:
            self.source=transcript_row.source
        else:
            self.source=source
        self.start=transcript_row.start
        self.strand = transcript_row.strand
        self.end=transcript_row.end
        self.score = transcript_row.score

        self.parent = transcript_row.parent
        self.attributes = transcript_row.attributes
        self.blast_hits = []   
        self.json_dict = None
             
        
    def __str__(self, to_gtf=False, print_cds=True):
        '''Each transcript will be printed out in the GFF style.
        This is pretty rudimentary, as the class does not hold any information on the original source, feature, score, etc.'''
        
        self.finalize() #Necessary to sort the exons
        lines = []
        transcript_counter = 0
#         assert self.selected_internal_orf_index > -1

        if self.strand is None:
            strand="."
        else:   
            strand=self.strand
        
        if to_gtf is True:
            parent_line = gtfLine(None)
        else:
            parent_line=gffLine(None)

        if print_cds is True:
            
            for index in range(len(self.internal_orfs)):
                
                if self.number_internal_orfs>1:
                    transcript_counter+=1
                    tid = "{0}.orf{1}".format(self.id, transcript_counter)
                    
                    if index==self.selected_internal_orf_index: self.attributes["maximal"]=True
                    else: self.attributes["maximal"]=False
                else:
                    tid = self.id
                cds_run = self.internal_orfs[index]
                    
                parent_line.chrom=self.chrom
                parent_line.source=self.source
                parent_line.feature=self.feature
                parent_line.start,parent_line.end=self.start,self.end
                parent_line.score=self.score
                parent_line.strand=strand
                parent_line.phase='.'
                parent_line.attributes=self.attributes
                
                parent_line.parent=self.parent
                parent_line.id=tid
                parent_line.name = self.id
            
                exon_lines = []
            
                cds_begin = False
            
                cds_count=0
                exon_count=0
                five_utr_count=0
                three_utr_count=0
    
                for segment in cds_run:
                    if cds_begin is False and segment[0]=="CDS": cds_begin = True
                    if segment[0]=="UTR":
                        if cds_begin is True:
                            if to_gtf is True:
                                if self.strand=="-": feature="5UTR"
                                else: feature="3UTR"
                            else:
                                if self.strand=="-": feature="five_prime_UTR"
                                else: feature="three_prime_UTR"
                        else:
                            if to_gtf is True:
                                if self.strand=="-": feature="3UTR"
                                else: feature="5UTR"
                            else:
                                if self.strand=="-": feature="three_prime_UTR"
                                else: feature="five_prime_UTR"
                        if "five" in feature or "5" in feature:
                            five_utr_count+=1
                            index=five_utr_count
                        else:
                            three_utr_count+=1
                            index=three_utr_count
                    else:
                        if segment[0]=="CDS":
                            cds_count+=1
                            index=cds_count
                        else:
                            exon_count+=1
                            index=exon_count
                        feature=segment[0]
                    if to_gtf is True:
                        exon_line=gtfLine(None)
                    else:
                        exon_line=gffLine(None)
                        
                    exon_line.chrom=self.chrom
                    exon_line.source=self.source
                    exon_line.feature=feature
                    exon_line.start,exon_line.end=segment[1],segment[2]
                    exon_line.strand=strand
                    exon_line.phase=None
                    exon_line.score = None
                    if to_gtf is True:
                        exon_line.gene=self.parent
                        exon_line.transcript=tid
                    else:
                        exon_line.id="{0}.{1}{2}".format(tid, feature,index)
                        exon_line.parent=tid
                        
                    exon_lines.append(str(exon_line))
            
            
                lines.append(str(parent_line))
                lines.extend(exon_lines) 
        else:
            if to_gtf is True:
                parent_line = gtfLine(None)
            else:
                parent_line=gffLine(None)
                    
            parent_line.chrom=self.chrom
            parent_line.source=self.source
            parent_line.feature=self.feature
            parent_line.start,parent_line.end=self.start,self.end
            parent_line.score=self.score
            parent_line.strand=strand
            parent_line.phase='.'
            parent_line.attributes=self.attributes
                
            parent_line.parent=self.parent
            parent_line.id=self.id
            parent_line.name = self.id
            
            lines=[str(parent_line)]
            exon_lines = []
            
            exon_count=0
            for exon in self.exons:
                exon_count+=1
                if to_gtf is True:
                    exon_line = gtfLine(None)
                else:
                    exon_line = gffLine(None)
                exon_line.chrom=self.chrom
                exon_line.source=self.source
                exon_line.feature="exon"
                exon_line.start,exon_line.end=exon[0],exon[1]
                exon_line.score=None
                exon_line.strand=strand
                exon_line.phase=None
                exon_line.attributes=self.attributes
                
                exon_line.id="{0}.{1}{2}".format(self.id, "exon",exon_count)
                exon_line.parent=self.id
                exon_lines.append(str(exon_line))
            
            lines.extend(exon_lines)
        
        return "\n".join(lines)
    
    def __eq__(self, other) -> bool:
        '''Two transcripts are considered identical if they have the same
        start, end, chromosome, strand and internal exons.
        IDs are not important for this comparison; two transcripts coming from different
        methods and having different IDs can still be identical.'''
        
        if not type(self)==type(other): return False
        self.finalize()
        other.finalize()
           
        if self.strand == other.strand and self.chrom == other.chrom and \
            self.start==other.start and self.end == other.end and \
            self.exons == other.exons:
            return True
          
        return False
    
    def __hash__(self):
        '''Returns the hash of the object (call to super().__hash__()).
        Necessary to be able to add these objects to hashes like sets.'''

        return super().__hash__()
    
    def __len__(self) -> int:
        '''Returns the length occupied by the unspliced transcript on the genome.'''
        return self.end-self.start+1

     
    def __lt__(self, other) -> bool:
        '''A transcript is lesser than another if it is on a lexicographic inferior chromosome,
        or if it begins before the other, or (in the case where they begin at the same location)
        it ends earlier than the other.'''
        if self.chrom!=other.chrom:
            return self.chrom<other.chrom
        if self==other:
            return False
        if self.start<other.start:
            return True
        elif self.start==other.start and self.end<other.end:
            return True
        return False
     
    def __gt__(self, other) -> bool:
        return not self<other
     
    def __le__(self, other) -> bool:
        return (self==other) or (self<other)
     
    def __ge__(self, other) -> bool:
        return (self==other) or (self>other) 
        
    def __getstate__(self):
        state=self.__dict__.copy()
        if hasattr(self, "json_dict"):
            if "requirements" in self.json_dict and "compiled" in self.json_dict["requirements"]:
                del state["json_dict"]["requirements"]["compiled"]

        if hasattr(self, "session"):
            state["session"].expunge_all()
            state["session"].close()
            del state["session"]
        if hasattr(self, "sessionmaker"):            
            del state["sessionmaker"]
            del state["engine"]

        if "blast_baked" in state:
            del state["blast_baked"]
            del state["query_baked"]
        if hasattr(self, "logger"):
            del state['logger']
            
        return state

    def __setstate__(self,state):
        self.__dict__.update(state)

    ######### Class instance methods ####################

    def addExon(self, gffLine):
        '''This function will append an exon/CDS feature to the object.'''

        if self.finalized is True:
            raise mikado_lib.exceptions.ModificationError("You cannot add exons to a finalized transcript!")
        
        if self.id not in gffLine.parent:
            raise mikado_lib.exceptions.InvalidTranscript("""Mismatch between transcript and exon:\n
            {0}\n
            {1}\n
            {2}
            """.format(self.id, gffLine.parent,gffLine))
        assert gffLine.is_exon is True, str(gffLine)
            
        if gffLine.feature=="CDS":
            store=self.combined_cds
        elif "combined_utr" in gffLine.feature or "UTR" in gffLine.feature.upper():
            store=self.combined_utr
        elif gffLine.feature=="exon":
            store=self.exons
        elif gffLine.feature=="start_codon":
            self.has_start_codon = True
            return
        elif gffLine.feature=="stop_codon":
            self.has_stop_codon = True
            return
        else:
            raise mikado_lib.exceptions.InvalidTranscript("Unknown feature: {0}".format(gffLine.feature))
            
        start,end=sorted([gffLine.start, gffLine.end])
        store.append((start, end) )

    #@profile
    def split_by_cds(self):
        '''This method is used for transcripts that have multiple ORFs. It will split them according to the CDS information
        into multiple transcripts. UTR information will be retained only if no ORF is down/upstream.
        The minimal overlap is defined inside the JSON at the key ["chimera_split"]["blast_params"]["minimal_hsp_overlap"]
        - basically, we consider a HSP a hit only if the overlap is over a certain threshold and the HSP evalue under a certain threshold. 
        
        The split by CDS can be executed in three different ways - CLEMENT, LENIENT, STRINGENT:
        
        - STRINGENT: split if two CDSs do not have hits in common - even when one or both do not have a hit at all.
        - CLEMENT: split only if two CDSs have hits and none of those is in common between them.
        - LENIENT: split if *both* lack hits, OR *both* have hits and none of those is in common.
        '''
        self.finalize()
        
        if self.json_dict is not None:
            minimal_overlap=self.json_dict["chimera_split"]["blast_params"]["minimal_hsp_overlap"]
        else:
            minimal_overlap=0
        new_transcripts=[]
        if self.number_internal_orfs<2:
            new_transcripts=[self]
        else:
            
            cds_boundaries=OrderedDict()
            for orf in sorted(self.loaded_bed12, key=operator.attrgetter("thickStart", "thickEnd")):
                cds_boundaries[(orf.thickStart,orf.thickEnd)]=[orf]
            
            if self.json_dict is not None and self.json_dict["chimera_split"]["blast_check"] is True:
                                
                cds_hit_dict=OrderedDict().fromkeys(cds_boundaries.keys())
                for key in cds_hit_dict:
                    cds_hit_dict[key]=set()

                for hit in self.blast_hits: #Determine for each CDS which are the hits available 
                    for hsp in filter(lambda hsp: hsp["hsp_evalue"]<=self.json_dict["chimera_split"]["blast_params"]["hsp_evalue"],  hit["hsps"]):
                        for cds_run in cds_boundaries: #If I have a valid hit b/w the CDS region and the hit, add the name to the set 
                            if abstractlocus.overlap(cds_run, (hsp["query_hsp_start"],hsp["query_hsp_end"]))>=minimal_overlap*(cds_run[1]+1-cds_run[0]):
                                cds_hit_dict[cds_run].add(hit["target"])
                
                new_boundaries = []
                for cds_boundary in cds_boundaries:
                    if new_boundaries==[]:
                        new_boundaries.append([cds_boundary])
                    else:
                        old_boundary=new_boundaries[-1][-1]
                        cds_hits = cds_hit_dict[cds_boundary]
                        old_hits = cds_hit_dict[old_boundary]
                        if cds_hits == set() and old_hits==set(): #No hit found for either CDS
                            if self.json_dict["chimera_split"]["blast_params"]["leniency"]=="CLEMENT": #If we are clement, we DO NOT split
                                new_boundaries[-1].append(cds_boundary)
                            else: #Otherwise, we do split
                                new_boundaries.append([cds_boundary])
                        elif cds_hits == set() or old_hits==set(): #We have hits for only one
                            if self.json_dict["chimera_split"]["blast_params"]["leniency"]=="STRINGENT": #If we are stringent, we split
                                new_boundaries.append([cds_boundary])
                            else:
                                new_boundaries[-1].append(cds_boundary)
                        else:
                            if set.intersection(cds_hits,old_hits)==set(): #We do not have any hit in common
                                new_boundaries.append([cds_boundary])
                            else: #We have hits in common
                                new_boundaries[-1].append(cds_boundary)
                
                final_boundaries=OrderedDict()
                for boundary in new_boundaries:
                    if len(boundary)==1:
                        assert len(boundary[0])==2
                        boundary=boundary[0]
                        final_boundaries[boundary]=cds_boundaries[boundary]
                    else:
                        nb = (boundary[0][0],boundary[-1][1])
                        final_boundaries[nb] = []
                        for boun in boundary:
                            final_boundaries[nb].extend(cds_boundaries[boun])
                            
                cds_boundaries=final_boundaries.copy()
                    
            if len(cds_boundaries)==1:
                new_transcripts=[self]
            else:
                spans = [] 
                
                for counter,(boundary, bed12_objects) in enumerate(cds_boundaries.items()):
                    #I *know* that I am moving left to right 
                    new_transcript=self.__class__()
                    new_transcript.feature="mRNA"
                    for attribute in ["chrom", "source", "score","strand", "attributes" ]:
                        setattr(new_transcript,attribute, getattr(self, attribute))
                    #Determine which ORFs I have on my right and left
                    new_transcript.parent = self.parent
                    if counter==0: #leftmost
                        left = None
                        right = cds_boundaries[list(cds_boundaries.keys())[counter+1]]
                    elif 1+counter==len(cds_boundaries): #rightmost
                        left = cds_boundaries[list(cds_boundaries.keys())[counter-1]]
                        right = None
                    else: #somewhere in the middle
                        left = cds_boundaries[list(cds_boundaries.keys())[counter-1]]
                        right = cds_boundaries[list(cds_boundaries.keys())[counter+1]]
                    counter+=1 #Otherwise they start from 0    
                    new_transcript.id = "{0}.split{1}".format(self.id, counter)
                    
                    my_exons = []
                    
                    
                    discarded_exons=[]
                    if self.strand == "-":
                        exons = sorted(self.exons, reverse=True, key=operator.itemgetter(0))
                    else:
                        exons = sorted(self.exons, key=operator.itemgetter(0))

                    tlength = 0
                    tstart = float("Inf")
                    tend = float("-Inf")
                    self.logger.debug( "Transcript {0} counter {1}, left {2} right {3}".format(self.id, counter, left, right) )
                    
                    for exon in exons:
                        #Translate into transcript coordinates
                        elength = exon[1]-exon[0]+1
                        texon = [tlength+1, tlength+elength]
                        tlength += elength
                        #Exon completely contained in the ORF
                        if texon[0]>=boundary[0] and texon[1]<=boundary[1]: 
                            my_exons.append(exon)
                        #Exon on the left of the CDS
                        elif texon[1]<=boundary[0]:
                            if left is None:
                                my_exons.append(exon)
                            else:
                                discarded_exons.append(exon)
                                continue
                        elif texon[0]>=boundary[1]:
                            if right is None:
                                my_exons.append(exon)
                            else:
                                discarded_exons.append(exon)
                                continue
                        #exon with partial UTR
                        else:
                            new_exon = list(exon)
                            if texon[0]<boundary[0]<=texon[1] and texon[1]<=boundary[1]:
                                if left is not None:
                                    if self.strand=="-":
                                        new_exon[1] = exon[0]+(texon[1]-boundary[0])
                                    else:
                                        new_exon[0]=exon[1]-(texon[1]-boundary[0])
                                    self.logger.debug("Tstart shifted for {0}, {1} to {2}".format(self.id, texon[0], boundary[0]))
                                    self.logger.debug("GStart shifted for {0}, {1} to {2}".format(self.id, exon[0], new_exon[1]))
                                    texon[0] = boundary[0]
                            if texon[0]<boundary[1]<texon[1] and texon[0]>=boundary[0]:
                                if right is not None:
                                    if self.strand == "-":
                                        new_exon[0] = exon[1]-(boundary[1]-texon[0])
                                    else:
                                        new_exon[1] = exon[0]+(boundary[1]-texon[0])
                                    self.logger.debug("Tend shifted for {0}, {1} to {2}".format(self.id, texon[1], boundary[1]))
                                    self.logger.debug("Gend shifted for {0}, {1} to {2}".format(self.id, exon[1], new_exon[1]))
                                    texon[1] = boundary[1]
                            elif texon[0]<=boundary[0]<=boundary[1]<=texon[1]: #Monoexonic
                                if self.strand == "-":
                                    if left is not None:
                                        new_exon[1] = exon[0]+(texon[1]-boundary[0])
                                    if right is not None:
                                        new_exon[0] = exon[1]-(boundary[1]-texon[0])
                                else:
                                    if left is not None:
                                        new_exon[0]=exon[1]-(texon[1]-boundary[0])
                                    if right is not None:
                                        new_exon[1] = exon[0]+(boundary[1]-texon[0])
                                
                                texon[0] = boundary[0]
                                texon[1] = boundary[1]
                                self.logger.debug("Tstart shifted for {0}, {1} to {2}".format(self.id, texon[0], boundary[0]))
                                self.logger.debug("GStart shifted for {0}, {1} to {2}".format(self.id, exon[0], new_exon[1]))
                                self.logger.debug("Tend shifted for {0}, {1} to {2}".format(self.id, texon[1], boundary[1]))
                                self.logger.debug("Gend shifted for {0}, {1} to {2}".format(self.id, exon[1], new_exon[1]))
                            
                            my_exons.append(tuple(sorted(new_exon)))
                        tstart=min(tstart, texon[0])
                        tend=max(tend,texon[1])
                           
                    if right is not None:
                        self.logger.debug( "TID {0} TEND {1} Boun[1] {2}".format(self.id, tend, boundary[1]) )
                        #assert tend == boundary[1]
                    elif left is not None:
                        self.logger.debug( "TID {0} TSTART {1} Boun[0] {2}".format(self.id, tstart, boundary[0]) )
                        #assert tstart == boundary[0]
                         
                    assert len(my_exons)>0, (discarded_exons, boundary)        
                    
                    new_transcript.exons=my_exons
                    
                    new_transcript.start=min( exon[0] for exon in new_transcript.exons )
                    new_transcript.end=max( exon[1] for exon in new_transcript.exons )
                    new_transcript.json_dict=self.json_dict
                    #Now we have to modify the BED12s to reflect the fact that we are starting/ending earlier
                    new_bed12s = []
                    for obj in bed12_objects:
                        assert type(obj) is bed12.BED12, (obj, bed12_objects)
                        obj.start=1
                        obj.end=min(obj.end, tend)-tstart
                        obj.thickStart=min(obj.thickStart,tend)-tstart+1 
                        obj.thickEnd=min(obj.thickEnd,tend)-tstart+1
                        obj.blockSizes=[obj.end]
                        new_bed12s.append(obj)
                    
                    new_transcript.finalize()
                    if new_transcript.monoexonic is True:
                        new_transcript.strand=None
                    new_transcript.load_orfs(new_bed12s)
                    
                    if new_transcript.selected_cds_length<=0:
                        err_message="No CDS information retained for {0} split {1}\n".format(self.id, counter)
                        err_message+="BED: {0}".format("\n\t".join([str(x) for x in new_bed12s]))
                        raise InvalidTranscript(err_message)
                        
                    new_transcripts.append(new_transcript)
                    nspan = (new_transcript.start,new_transcript.end)
                    self.logger.debug("Transcript {0} split {1}, discarded exons: {2}".format( self.id, counter, discarded_exons ))
                    for span in spans:
                        overl = abstractlocus.overlap(span, nspan )  
                        
                        self.logger.debug("Comparing start-ends for split of {0}. SpanA: {1} SpanB: {2} Overlap: {3}".format(
                                                                                                                            self.id, span,
                                                                                                                            nspan, overl
                                                                                                                            ))
                        if overl>0:
                            err_message = "Invalid overlap for {0}! T1: {1}. T2: {2}".format( self.id, span, (new_transcript.start,new_transcript.end) )
                            self.logger.error(err_message)
                            raise InvalidTranscript(err_message) 
                    
                    spans.append([new_transcript.start, new_transcript.end])
                    
                    

        assert len(new_transcripts)>0, str(self)
        for nt in new_transcripts:
            yield nt
        
        return

    def remove_utrs(self):
        '''Method to strip a transcript from its UTRs.
        It will not execute anything if the transcript lacks a CDS or
        it has more than one ORF defined.'''
        
        self.finalize()
        if self.selected_cds_length==0:
            return
        elif self.three_utr_length + self.five_utr_length ==0:
            return #No UTR to strip 
        
        elif self.number_internal_orfs>1:
            return
        elif re.search("\.orf[0-9]+$", self.id  ):
            return
        
        self.finalized = False
        exons = []
        cds_start,cds_end = self.combined_cds[0][0], self.combined_cds[-1][1]
        assert type(cds_start) is int
        assert type(cds_end) is int
        if len(self.selected_cds)==1:
            self.exons = self.selected_cds
        else:
            for exon in self.exons:
                if exon in self.combined_utr:
                    continue
                elif exon in self.selected_cds:
                    exons.append(exon)
                elif exon[0]<=cds_start<=exon[1]:
                    exons.append((cds_start, exon[1]))
                elif exon[0]<=cds_end<=exon[1]:
                    exons.append( (exon[0], cds_end ))
                else:
                    raise InvalidTranscript("Exon: {0}; cds_start: {1}; cds_end: {2}; ID: {3}".format( exon, self.selected_cds_start, self.selected_cds_end, self.id))
            assert len(exons)<len(self.exons) or exons[0][0]>self.exons[0][0] or exons[-1][1]<self.exons[-1][1], (exons, self.exons)
            self.exons = exons
        self.start = cds_start
        self.end = cds_end
        self.combined_utr = []
        self.finalize()
        

    def strip_cds(self):
        '''Method to completely remove CDS information from a transcript. Necessary for those cases where
        the input is malformed.'''
        
        self.logger.warning("Stripping CDS from {0}".format(self.id))
        self.finalized = False
        self.combined_cds = []
        self.combined_utr = []
        self.finalize()
        


    def finalize(self):
        '''Function to calculate the internal introns from the exons.
        In the first step, it will sort the exons by their internal coordinates.'''
        
        # We do not want to repeat this step multiple times
        if self.finalized is True:
            return
        self.introns = []
        self.splices=[]
        if len(self.exons)==0:
            raise mikado_lib.exceptions.InvalidTranscript("No exon defined for the transcript {0}. Aborting".format(self.tid))

        if len(self.exons)>1 and self.strand is None:
            raise mikado_lib.exceptions.InvalidTranscript("Multiexonic transcripts must have a defined strand! Error for {0}".format(self.id))

        if self.combined_utr!=[] and self.combined_cds==[]:
            raise mikado_lib.exceptions.InvalidTranscript("Transcript {tid} has defined UTRs but no CDS feature!".format(tid=self.id))

        self.exons = sorted(self.exons, key=operator.itemgetter(0,1) ) # Sort the exons by start then stop
        
        if self.cdna_length > self.combined_utr_length + self.combined_cds_length:
            if self.combined_utr == [] and self.combined_cds!=[]:
                self.combined_cds = sorted(self.combined_cds, key=operator.itemgetter(0,1))
                for exon in self.exons:
                    if exon in self.combined_cds:
                        continue
                    elif exon[1]<self.combined_cds[0][0] or exon[0]>self.combined_cds[-1][1]:
                        self.combined_utr.append(exon)
                    elif exon[0]<self.combined_cds[0][0] and exon[1]==self.combined_cds[0][1]:
                        self.combined_utr.append( (exon[0], self.combined_cds[0][0]-1)  )
                    elif exon[1]>self.combined_cds[-1][1] and exon[0]==self.combined_cds[-1][0]:
                        self.combined_utr.append( (self.combined_cds[-1][1]+1, exon[1]))
                    else:
                        if len(self.combined_cds)==1:
                            self.combined_utr.append( (exon[0], self.combined_cds[0][0]-1)  )
                            self.combined_utr.append( (self.combined_cds[-1][1]+1, exon[1]))
                        else:
                            raise mikado_lib.exceptions.InvalidCDS("Error while inferring the UTR", exon, self.id, self.exons, self.combined_cds
                                                                            ) 
                if not (self.combined_cds_length==self.combined_utr_length==0 or  self.cdna_length == self.combined_utr_length + self.combined_cds_length):
                    raise mikado_lib.exceptions.InvalidCDS("Failed to create the UTR", self.id, self.exons, self.combined_cds, self.combined_utr)
            else:
#                 raise mikado_lib.InvalidCDS("Invalid input", self.id, (self.cdna_length, self.combined_utr_length, self.combined_cds_length),
#                                                           self.exons, self.combined_cds, self.combined_utr)
                pass


#         assert len(self.exons)>0
        try:
            if self.exons[0][0]!=self.start or self.exons[-1][1]!=self.end:
                if self.exons[0][0]>self.start and self.selected_cds[0][0]==self.start:
                    self.exons[0] = (self.start, self.exons[0][0])
                if self.exons[-1][1]<self.end and self.selected_cds[-1][1]==self.end:
                    self.exons[-1] = (self.exons[-1][0], self.end)
                
                if self.exons[0][0]!=self.start or self.exons[-1][1]!=self.end:    
                    raise mikado_lib.exceptions.InvalidTranscript("""The transcript {id} has coordinates {tstart}:{tend},
                    but its first and last exons define it up until {estart}:{eend}!""".format(
                                                                                               id=self.id,
                                                                                               tstart=self.start,
                                                                                               tend=self.end,
                                                                                               estart=self.exons[0][0],
                                                                                               eend=self.exons[-1][1]
                                                                                               ))
        except IndexError as err:
            raise mikado_lib.exceptions.InvalidTranscript(err, self.id, str(self.exons))

        if len(self.exons)>1:
            for index in range(len(self.exons)-1):
                exonA, exonB = self.exons[index:index+2]
                if exonA[1]>=exonB[0]:
                    raise mikado_lib.exceptions.InvalidTranscript("Overlapping exons found!\n{0} {1}/{2}\n{3}".format(self.id, exonA, exonB, self.exons))
                self.introns.append( (exonA[1]+1, exonB[0]-1) ) #Append the splice junction
                self.splices.extend( [exonA[1]+1, exonB[0]-1] ) # Append the splice locations

        self.combined_cds = sorted(self.combined_cds, key=operator.itemgetter(0,1))
        self.combined_utr = sorted(self.combined_utr, key=operator.itemgetter(0,1))
        if len(self.combined_utr)>0 and self.combined_utr[0][0]<self.combined_cds[0][0]:
            if self.strand=="+": self.has_start_codon=True
            elif self.strand=="-": self.has_stop_codon=True
        if len(self.combined_utr)>0 and self.combined_utr[-1][1]>self.combined_cds[-1][1]:
            if self.strand=="+": self.has_stop_codon=True
            elif self.strand=="-": self.has_start_codon=True
        
        self.internal_orfs = []
        self.segments = [ ("exon",e[0],e[1]) for e in self.exons] + \
                    [("CDS", c[0],c[1]) for c in self.combined_cds ] + \
                    [ ("UTR", u[0], u[1]) for u in self.combined_utr ]
        self.segments =  sorted(self.segments, key=operator.itemgetter(1,2,0) )
                
        self.internal_orfs.append(self.segments)
        if self.combined_cds_length>0:
            self.selected_internal_orf_index=0
        
        self.introns = set(self.introns)
        self.splices = set(self.splices)            
        _ = self.selected_internal_orf
#         assert self.selected_internal_orf_index > -1
        if len(self.combined_cds)>0:
            self.feature="mRNA"
        
        self.set_relative_properties()
        
        self.finalized = True
        return

    def set_relative_properties(self):
        '''Function to set to the basic value relative values like e.g. retained_intron'''
        self.retained_introns=[]
        self.retained_fraction=0
        self.exon_fraction=1
        self.intron_fraction=1
        self.cds_intron_fraction=1
        self.selected_cds_intron_fraction=1

    def reverse_strand(self):
        '''Method to reverse the strand'''
        if self.strand=="+":
            self.strand="-"
        elif self.strand=="-":
            self.strand="+"
        elif self.strand is None:
            pass
        return

    #@profile
    def connect_to_db(self):
        
        '''This method will connect to the database using the information contained in the JSON configuration.'''
        
        db = self.json_dict["db"]
        dbtype = self.json_dict["dbtype"]
        
        self.engine = create_engine("{dbtype}:///{db}".format(
                                                              db=db,
                                                              dbtype=dbtype)
                                    )   #create_engine("sqlite:///{0}".format(args.db))
        
        self.sessionmaker = sessionmaker()
        self.sessionmaker.configure(bind=self.engine)
        self.session=self.sessionmaker()
    
    @asyncio.coroutine
    def load_information_from_db(self, json_dict, introns=None, session=None, loop=None):
        '''This method will invoke the check for:
        
        orfs
        blast hits (this necessarily after the ORF check).
        
        Verified introns can be provided from outside using the keyword. Otherwise, they will be extracted from the database directly.
        '''
        
        self.logger.debug("Loading {0}".format(self.id))
        self.load_json(json_dict)
        if session is None:
            self.connect_to_db()
        else:
            self.session=session
        self.load_verified_introns(introns)
        self.query_id = self.query_baked(self.session).params(query_name=self.id).all()
        if len(self.query_id)==0:
            self.logger.warning("Transcript not in database: {0}".format(self.id))
        else:
            self.query_id = self.query_id[0].id
            yield from self.load_orfs_coroutine()
            yield from self.load_blast()
        self.logger.debug("Loaded {0}".format(self.id))
    
    #@profile
    def load_json(self, json_dict):
        self.json_dict = json_dict
    
    #@profile
    def load_verified_introns(self, introns=None):
        
        '''This method will load verified junctions from the external (usually the superlocus class).'''
        
       
        if introns is None:
            chrom_id = self.session.query(Chrom.id).filter(Chrom.name==self.chrom).one().id
            for intron in self.introns:
                if self.session.query(junction).filter(and_(
                                                            junction.chrom_id==chrom_id,
                                                            junction.junctionStart==intron[0],
                                                            junction.junctionEnd==intron[1],
                                                            junction.strand==self.strand
                                                           )
                                                      ).count()==1:
                    self.verified_introns.add(intron)
                    
        else:
            for intron in introns:
                if intron in self.introns:
                    self.verified_introns.add(intron)
        
        return       
        
    #@profile
    @asyncio.coroutine
    def retrieve_orfs(self):
        
        '''This method will look up the ORFs loaded inside the database.
        During the selection, the function will also remove overlapping ORFs.
        '''
        
        if self.query_id is None:
            return []
        
        trust_strand = self.json_dict["orf_loading"]["strand_specific"]
        
        orf_results = self.orf_baked(self.session).params(query_id = self.query_id)
        
        if (self.monoexonic is False) or (self.monoexonic is True and trust_strand is True):
            #Remove negative strand ORFs for multiexonic transcripts, or monoexonic strand-specific transcripts
            candidate_orfs=list(filter(lambda orf: orf.strand!="-", orf_results ))
        else:
            candidate_orfs=orf_results.all()
            
        if len(candidate_orfs)==0: return []
        else:
            return [orf.as_bed12() for orf in candidate_orfs]
    
    #@profile
    @asyncio.coroutine
    def load_orfs_coroutine(self):
        '''Asynchronous coroutine for loading orfs from the database'''
        candidate_orfs = yield from self.retrieve_orfs()
        self.logger.debug("Loading ORF for {0}".format(self.id))
        self.load_orfs(candidate_orfs)
        self.logger.debug("Loaded ORF for {0}".format(self.id))
    
    
    def load_orfs(self, candidate_orfs):

        '''This method replicates what is done internally by the "cdna_alignment_orf_to_genome_orf.pl" utility in the
        TransDecoder suite. It takes as argument "candidate_orfs" i.e. a list of BED12 serialised objects.
        The method expects as argument a dictionary containing BED entries, and indexed by the transcript name.
        The indexed name *must* equal the "id" property, otherwise the method returns immediately. 
        If no entry is found for the transcript, the method exits immediately. Otherwise, any CDS information present in
        the original GFF/GTF file is completely overwritten.
        Briefly, it follows this logic:
        - Finalise the transcript
        - Retrieve from the dictionary (input) the CDS object
        - Sort in decreasing order the CDSs on the basis of:
            - Presence of start/stop codon
            - CDS length (useful for monoexonic transcripts where we might want to set the strand)
        - For each CDS:
            - If the ORF is on the + strand:
                - all good
            - If the ORF is on the - strand:
                - if the transcript is monoexonic: invert its genomic strand
                - if the transcript is multiexonic: skip
            - Start looking at the exons
        
        '''

        #Prepare the transcript
        self.finalize()
        
        #Prepare the ORFs to load
        if self.json_dict is not None:
            minimal_secondary_orf_length=self.json_dict["orf_loading"]["minimal_secondary_orf_length"]
        else:
            minimal_secondary_orf_length=0
            
        candidate_cliques=self.find_overlapping_cds(candidate_orfs)
        new_orfs=[]
        for clique in candidate_cliques:
            new_orfs.append(sorted(clique, reverse=True, key=operator.attrgetter("cds_len"))[0])

        candidate_orfs=[]
        if len(new_orfs)>0:
            candidate_orfs=[new_orfs[0]]
            for orf in filter(lambda x: x.cds_len>minimal_secondary_orf_length, new_orfs[1:]):
                candidate_orfs.append(orf)

        if candidate_orfs is None or len(candidate_orfs)==0:
            self.logger.debug( "No ORF for {0}".format(self.id) )
            return
        
        self.combined_utr = []
        self.combined_cds = []
        self.internal_orfs = []
        self.finalized = False
        primary_orf=True # Token to be set to False after the first CDS is exhausted
        self.loaded_bed12 = [] #This will keep in memory the original BED12 objects
        
        #If we are looking at a multiexonic transcript
        if not (self.monoexonic is True and self.strand is None):
            candidate_orfs = list(filter( lambda co: co.strand == "+", candidate_orfs))
        else:
            #Candidate ORFs are already sorted by CDS length
            candidate_orfs = list(filter( lambda co: co.strand == candidate_orfs[0].strand, candidate_orfs  ))
        
                
        for orf in sorted(candidate_orfs, key=operator.attrgetter( "cds_len"  ), reverse=True):
            #Minimal check
            if primary_orf is True:
                self.has_start_codon, self.has_stop_codon = orf.has_start_codon, orf.has_stop_codon
                primary_orf = False
            
            if not (orf.thickStart>=1 and orf.thickEnd<=self.cdna_length): #Leave leeway for TD
                self.logger.warning("Wrong ORF:\n\t{0}\ncDNA length: {1}\nStart,end: {2}-{3}".format(str(orf), self.cdna_length, self.start, self.end))
                continue

            if self.strand is None:
                self.strand = orf.strand
            
            self.loaded_bed12.append(orf)
            cds_exons = []
            current_start, current_end = 0,0
            if self.strand == "+":
                for exon in sorted(self.exons, key=operator.itemgetter(0,1)):
                    cds_exons.append(("exon", exon[0], exon[1] ) )
                    current_start+=1
                    current_end+=exon[1]-exon[0]+1
                    #Whole UTR
                    if current_end<orf.thickStart or current_start>orf.thickEnd:
                        cds_exons.append( ("UTR", exon[0], exon[1])  )
                    else:
                        c_start = exon[0] + max(0, orf.thickStart-current_start )
                        if c_start > exon[0]:
                            u_end = c_start-1
                            cds_exons.append( ("UTR", exon[0], u_end) )
                        c_end = exon[1] - max(0, current_end - orf.thickEnd )
#                         assert c_end>=exon[0] 
                        if c_start<=c_end:
                            cds_exons.append(("CDS", c_start, c_end))
                        if c_end < exon[1]:
                            cds_exons.append( ("UTR", c_end+1, exon[1]  ) )
                    current_start=current_end
                            
            elif self.strand=="-":
                for exon in sorted(self.exons, key=operator.itemgetter(0,1), reverse=True):
                    cds_exons.append(("exon", exon[0], exon[1] ) )
                    current_start+=1
                    current_end+=exon[1]-exon[0]+1
                    if current_end<orf.thickStart or current_start>orf.thickEnd:
                        cds_exons.append( ("UTR", exon[0], exon[1] ))
                    else:
                        c_end = exon[1] - max(0,orf.thickStart - current_start ) 
#                         assert c_end>=exon[0]
                        if c_end < exon[1]:
                            cds_exons.append(("UTR", c_end+1, exon[1]))
                        c_start = exon[0] + max(0, current_end - orf.thickEnd )
                        cds_exons.append( ("CDS", c_start, c_end) )
                        if c_start>exon[0]:
                            cds_exons.append( ("UTR", exon[0], c_start-1) )
                    current_start=current_end
        
            self.internal_orfs.append( sorted(cds_exons, key=operator.itemgetter(1,2)   ) )

        if len(self.internal_orfs)==1:
            self.combined_cds = sorted(
                              [(a[1],a[2]) for a in filter(lambda x: x[0]=="CDS", self.internal_orfs[0])],
                              key=operator.itemgetter(0,1)
                              
                              )
            self.combined_utr = sorted(
                              [(a[1],a[2]) for a in filter(lambda x: x[0]=="UTR", self.internal_orfs[0])],
                              key=operator.itemgetter(0,1)
                              
                              )
            
            
        elif len(self.internal_orfs)>1:
            
            cds_spans = []
            candidates = []
            for internal_cds in self.internal_orfs:
                candidates.extend([tuple([a[1],a[2]]) for a in filter(lambda tup: tup[0]=="CDS", internal_cds  )])
                              
            candidates=set(candidates)
            #for mc in self.merge_cliques(*self.find_cliques(candidates)):
            for mc in self.find_communities(candidates):
                span=tuple([min(t[0] for t in mc),
                            max(t[1] for t in mc)                        
                            ])
                cds_spans.append(span)
 
            self.combined_cds = sorted(cds_spans, key = operator.itemgetter(0,1))
            
            #This method is probably OBSCENELY inefficient, but I cannot think of a better one for the moment.
            curr_utr_segment = None

            utr_pos = set.difference( 
                                                  set.union(*[ set(range(exon[0],exon[1]+1)) for exon in self.exons]),
                                                  set.union(*[ set(range(cds[0],cds[1]+1)) for cds in self.combined_cds])
                                                  )
            for pos in sorted(list(utr_pos)):
                if curr_utr_segment is None:
                    curr_utr_segment = (pos,pos)
                else:
                    if pos==curr_utr_segment[1]+1:
                        curr_utr_segment = (curr_utr_segment[0],pos)
                    else:
                        self.combined_utr.append(curr_utr_segment)
                        curr_utr_segment = (pos,pos)
                        
            if curr_utr_segment is not None:
                self.combined_utr.append(curr_utr_segment)           
                                   
            assert self.cdna_length == self.combined_cds_length + self.combined_utr_length, (self.cdna_length, self.combined_cds, self.combined_utr)                            
        
        if self.internal_orfs == []:
            self.finalize()
        else:
            self.feature="mRNA"
            self.finalized=True
        return
        

    @asyncio.coroutine
    #@profile
    def load_blast(self):
        
        '''This method looks into the DB for hits corresponding to the desired requirements.
        Hits will be loaded into the "blast_hits" list; we will not store the SQLAlchemy query object,
        but rather its representation as a dictionary (using the Hit.as_dict() method).       
        '''
#         
        if self.query_id is None:
            return
        
        
        if self.json_dict["chimera_split"]["blast_check"] is False:
            return
        
        max_target_seqs = self.json_dict["chimera_split"]["blast_params"]["max_target_seqs"] or float("inf")
        maximum_evalue = self.json_dict["chimera_split"]["blast_params"]["evalue"]
        
        blast_hits_query = self.blast_baked(self.session).params(query_id = self.query_id, evalue = maximum_evalue, max_target_seqs=max_target_seqs )                 
        for hit in blast_hits_query:
            self.blast_hits.append(hit.as_dict())
        
    def load_cds(self, cds_dict, trust_strand=False, minimal_secondary_orf_length=0):
        '''Deprecated'''
        raise NotImplementedError()
                        
    def set_logger(self, logger=None):
        '''Set a logger for the instance.'''
        if logger is None:
            logger = abstractlocus.create_default_logger()
        self.logger = logger
      
    @classmethod
    ####################Class methods#####################################  
    def find_cliques(cls, candidates):
        '''Wrapper for the abstractlocus method. It will pass to the function the class's "is_intersecting" method
        (which would be otherwise be inaccessible from the abstractlocus class method)'''
#         return abstractlocus.find_cliques( candidates, inters=cls.is_intersecting)
        raise NotImplementedError()
    
    @classmethod
    def find_overlapping_cds(cls, candidates: list) -> list:
        '''Wrapper for the abstractlocus method, used for finding overlapping ORFs.
        It will pass to the function the class's "is_overlapping_cds" method
        (which would be otherwise be inaccessible from the abstractlocus class method).
        As we are interested only in the communities, not the cliques, this wrapper discards the cliques
        (first element of the abstractlocus.find_communities results) 
        
        '''
#         graph, cliques = abstractlocus.find_cliques( candidates, inters=cls.is_overlapping_cds)
        d=dict( (x.name, x) for x in candidates )
        communities = abstractlocus.find_communities( abstractlocus.define_graph(d, inters=cls.is_overlapping_cds) )[1]
        final_communities = []
        for comm in communities:
            #Each community is a frozenset of ids
            final_communities.append( set(d[x] for x in comm)  )
        final_communities = [frozenset(x) for x in final_communities]
        return final_communities 
    
    
    @classmethod
    def is_overlapping_cds(cls, first, second):
        if first==second or cls.overlap(
                                        (first.thickStart,first.thickEnd),
                                        (second.thickStart,second.thickEnd))<=0:
            return False
        return True
        
    
    @classmethod
    def is_intersecting(cls, first, second):
        '''Implementation of the is_intersecting method.'''
        if first==second or cls.overlap(first, second)<0:
            return False
        return True

    @classmethod
    def overlap(cls, first,second):
        lend = max(first[0], second[0])
        rend = min(first[1], second[1])
        return rend-lend
        
    @classmethod
    def find_communities(cls, objects: list) -> list:
        '''Wrapper for the abstractlocus method.
        As we are interested only in the communities, not the cliques, this wrapper discards the cliques
        (first element of the abstractlocus.find_communities results) 
        '''
        d=dict((x, x) for x in objects )
        communities = abstractlocus.find_communities( abstractlocus.define_graph(d, inters=cls.is_intersecting) )[1]

        return communities 

    
    @classmethod
    def get_available_metrics(cls) -> list:
        '''This function retrieves all metrics available for the class.'''
        metrics = list(x[0] for x in filter(lambda y: "__" not in y[0] and type(cls.__dict__[y[0]]) is metric, inspect.getmembers(cls)))
        assert "tid" in metrics and "parent" in metrics and "score" in metrics
        final_metrics = ["tid","parent","score"] + sorted(list(filter(lambda x: x not in ["tid","parent","score"], metrics )))  
        return final_metrics 

    ####################Class properties##################################

    @property
    def id(self):
        '''ID of the transcript - cannot be an undefined value.'''
        return self.__id
    
    @id.setter
    def id(self, Id):
        if type(Id) is not str:
            raise ValueError("Invalid value for id: {0}, type {1}".format(
                                                                          Id, type(Id)))
        self.__id = sys.intern(Id)

    @property
    def available_metrics(self) -> list:
        '''Return the list of available metrics, using the "get_metrics" function.'''
        return self.get_metrics()

    @property
    def strand(self):
        return self.__strand
    
    @strand.setter
    def strand(self,strand):
        if strand in ("+", "-"):
            self.__strand=strand
        elif strand in (None,".","?"):
            self.__strand=None
        else:
            raise ValueError("Invalid value for strand: {0}".format(strand))
        
    @property
    def selected_internal_orf(self):
        '''This property will return the tuple of tuples of the ORF selected as "best".
        To avoid memory wasting, the tuple is accessed in real-time using 
        a token (__max_internal_orf_index) which holds the position in the __internal_cds list of the longest CDS.'''
        if len(self.combined_cds)==0: # Non-sense to calculate the maximum CDS for transcripts without it
            self.__max_internal_orf_length=0
            self.selected_internal_orf_index=0
            return tuple([])
        else:
            return self.internal_orfs[self.selected_internal_orf_index]
    
    @property
    def selected_internal_orf_cds(self):
        '''This property will return the tuple of tuples of the CDS segments of the selected ORF
        inside the transcript. To avoid memory wasting, the tuple is accessed in real-time using 
        a token (__max_internal_orf_index) which holds the position in the __internal_cds list of the longest CDS.'''
        if len(self.combined_cds)==0: # Non-sense to calculate the maximum CDS for transcripts without it
            return tuple([])
        else:
            return list(filter(lambda x: x[0]=="CDS", self.internal_orfs[self.selected_internal_orf_index])) 

    @property
    def five_utr(self):
        '''Returns the exons in the 5' UTR of the selected ORF. If the start codon is absent, no UTR is given.'''
        if len(self.combined_cds)==0:
            return []
#         elif self.has_start_codon is False:
#             return []
        if self.strand=="-":
            return list(filter( lambda exon: exon[0]=="UTR" and exon[1]>self.selected_cds_start, self.selected_internal_orf  )  )
        else:
            return list(filter( lambda exon: exon[0]=="UTR" and exon[2]<self.selected_cds_start, self.selected_internal_orf  )  )
            

    @property
    def three_utr(self):
        '''Returns the exons in the 3' UTR of the selected ORF. If the end codon is absent, no UTR is given.'''
        if len(self.combined_cds)==0:
            return []
#         elif self.has_stop_codon is False:
#             return []
        if self.strand=="-":
            return list(filter( lambda exon: exon[0]=="UTR" and exon[2]<self.selected_cds_end, self.selected_internal_orf  )  )
        else:
            return list(filter( lambda exon: exon[0]=="UTR" and exon[1]>self.selected_cds_end, self.selected_internal_orf  )  )

    @property
    def selected_internal_orf_index(self):
        '''Token which memorizes the position in the ORF list of the selected ORF.'''
        return self.__max_internal_orf_index
    
    @selected_internal_orf_index.setter
    def selected_internal_orf_index(self, index):
        if index is None:
            self.__max_internal_orf_index = index
            return
        if type(index) is not int:
            raise TypeError()
        if index<0 or index>=len(self.internal_orfs):
            raise IndexError("No ORF corresponding to this index: {0}".format(index))
        self.__max_internal_orf_index = index
        
    @property
    def internal_orf_lengths(self):
        '''This property returns a list of the lengths of the internal ORFs.'''
        lengths = []
        for internal_cds in self.internal_orfs:
            lengths.append( sum( x[2]-x[1]+1 for x in filter(lambda c: c[0]=="CDS", internal_cds) ) )
        lengths = sorted(lengths, reverse=True)
        return lengths
        
    @property
    def non_overlapping_cds(self):
        '''This property returns a set containing the set union of all CDS segments inside the internal CDSs.
        In the case of a transcript with no CDS, this is empty.
        In the case where there is only one CDS, this returns the combined_cds holder.
        In the case instead where there are multiple CDSs, the property will calculate the set union of all CDS segments.'''
        if self.__non_overlapping_cds is None: 
            self.finalize()
            self.__non_overlapping_cds=set()
            for internal_cds in self.internal_orfs:
                segments = set([(x[1],x[2]) for x in filter(lambda segment: segment[0]=="CDS", internal_cds   )])
                self.__non_overlapping_cds.update(segments)
        return self.__non_overlapping_cds
    
    @non_overlapping_cds.setter
    def non_overlapping_cds(self,arg):
        '''Setter for the non_overlapping_cds property.'''
        self.__non_overlapping_cds = arg
    
    @property
    def exons(self):
        '''This property stores the exons of the transcript as (start,end) tuples.'''
        return self.__exons
    
    @exons.setter
    def exons(self, *args):
        if type(args[0]) not in (set,list):
            raise TypeError(type(args[0]))
        self.__exons=args[0]

    @property
    def combined_cds_introns(self):
        '''This property returns the introns which are located between CDS segments in the combined CDS.'''
        if len(self.combined_cds)<2:
            return set()
        cintrons=[]
        all_cintrons=[]
        for position in range(len(self.combined_cds)-1):
            former=self.combined_cds[position]
            latter=self.combined_cds[position+1]
            junc=(former[1]+1,latter[0]-1)
            if junc in self.introns:
                cintrons.append(junc)
        if len(self.selected_cds_introns)>0:
            assert len(cintrons)>0,(self.id, self.selected_cds_introns,all_cintrons,self.introns) 
        cintrons=set(cintrons)
        return cintrons

    @property
    def selected_cds_introns(self):
        '''This property returns the introns which are located between CDS segments in the selected ORF.'''
        cintrons=[]
        for position in range(len(self.selected_internal_orf_cds)-1):
            cintrons.append(
                            (self.selected_internal_orf_cds[position][1]+1, self.selected_internal_orf_cds[position+1][2]-1)
                            )
        cintrons=set(cintrons)
        return cintrons
    
    @property
    def combined_cds_start(self):
        '''This property returns the location of the start of the combined CDS for the transcript.
        If no CDS is defined, it defaults to the transcript start.'''
        if len(self.combined_cds)==0:
            if self.strand=="+":
                return self.start
            else:
                return self.end
        if self.strand=="+":
            return self.combined_cds[0][0]
        else:
            return self.combined_cds[-1][1]
       
    @property
    def combined_cds(self):
        '''This is a list which contains all the non-overlapping CDS segments inside the cDNA.
        The list comprises the segments as duples (start,end).'''
        return self.__combined_cds


    @combined_cds.setter
    def combined_cds(self, *args):
        if type(args[0]) is not list or (len(args[0])>0 and len(list(filter(
                                                                            lambda x: len(x)!=2 or type(x[0]) is not int or type(x[1]) is not int, args[0])) )>0):
            raise TypeError("Invalid value for combined CDS: {0}".format(args[0]))
        self.__combined_cds = args[0]

    @property
    def combined_utr(self):
        '''This is a list which contains all the non-overlapping UTR segments inside the cDNA.
        The list comprises the segments as duples (start,end).'''
        return self.__combined_utr

    @combined_utr.setter
    def combined_utr(self, *args):
        if type(args[0]) is not list or (len(args[0])>0 and len(list(filter(
                                                                            lambda x: len(x)!=2 or type(x[0]) is not int or type(x[1]) is not int, args[0])) )>0):
            raise TypeError("Invalid value for combined CDS: {0}".format(args[0]))
        self.__combined_utr = args[0]



    @property
    def combined_cds_end(self):
        '''This property returns the location of the end of the combined CDS for the transcript.
        If no CDS is defined, it defaults to the transcript end.'''
        if len(self.combined_cds)==0:
            if self.strand=="+":
                return self.end
            else:
                return self.start
        if self.strand=="-":
            return self.combined_cds[0][0]
        else:
            return self.combined_cds[-1][1]

    @property
    def selected_cds(self):
        '''This property return the CDS exons of the ORF selected as best inside the cDNA, in the form of duplices (start, end)'''
        if len(self.combined_cds)==0:
            self.__selected_cds=[]
        else:
            self.__selected_cds=list((x[1],x[2]) for x in filter(lambda x: x[0]=="CDS", self.selected_internal_orf))
        return self.__selected_cds

    @property
    def selected_cds_start(self):
        '''This property returns the location of the start of the best CDS for the transcript.
        If no CDS is defined, it defaults to the transcript start.'''

        if len(self.combined_cds)==0:
            return None
#             if self.strand=="+":
#                 return self.start
#             else:
#                 return self.end
            
        if self.strand=="-":
            return self.selected_cds[-1][1]
        else:
            return self.selected_cds[0][0]
        
    @property
    def selected_cds_end(self):
        '''This property returns the location of the end of the best CDS for the transcript.
        If no CDS is defined, it defaults to the transcript start.'''

        if len(self.combined_cds)==0:
            return None
#             if self.strand=="+":
#                 return self.end
#             else:
#                 return self.start
        if self.strand=="-":
            return self.selected_cds[0][0]
        else:
            return self.selected_cds[-1][1]


    @property
    def monoexonic(self):
        if len(self.exons)==1:
            return True
        return False


    #################### Class metrics ##################################

    @metric
    def tid(self):
        '''ID of the transcript - cannot be an undefined value. Alias of id.'''
        return self.id
    
    @tid.setter
    def tid(self,tid):
        self.id=tid
        
    @metric
    def parent(self):
        '''Name of the parent feature of the transcript.'''
        return self.__parent
    
    @parent.setter
    def parent(self,parent):
        if type(parent) is list or parent is None:
            self.__parent=parent
        elif type(parent) is str:
            if "," in parent:
                self.__parent=parent.split(",")
            else:
                self.__parent=[parent]
        else:
            raise ValueError("Invalid value for parent: {0}, type {1}".format(
                                                                          parent, type(parent)))
            
    @metric
    def score(self):
        '''Numerical value which summarizes the reliability of the transcript.'''
        return self.__score
        
    @score.setter
    def score(self,score):
        
        if score is not None:
            if type(score) not in (float,int):
                try:
                    score=float(score)
                except:
                    raise ValueError("Invalid value for score: {0}, type {1}".format(
                                                                          score, type(score)))
        self.__score=score
                
        

    @metric
    def combined_cds_length(self):
        '''This property return the length of the CDS part of the transcript.'''
        return sum([ c[1]-c[0]+1 for c in self.combined_cds ])
    
    @metric
    def combined_cds_num(self):
        '''This property returns the number of non-overlapping CDS segments in the transcript.'''
        return len( self.combined_cds )

    @metric
    def combined_cds_num_fraction(self):
        '''This property returns the fraction of non-overlapping CDS segments in the transcript
        vs. the total number of exons'''
        return len( self.combined_cds )/len(self.exons)

    @metric
    def combined_cds_fraction(self):
        '''This property return the percentage of the CDS part of the transcript vs. the cDNA length'''
        return self.combined_cds_length/self.cdna_length
    
    @metric
    def combined_utr_length(self):
        '''This property return the length of the UTR part of the transcript.'''
        return sum([ e[1]-e[0]+1 for e in self.combined_utr ])
    
    @metric
    def combined_utr_fraction(self):
        '''This property returns the fraction of the cDNA which is not coding according
        to any ORF. Complement of combined_cds_fraction'''
        return 1-self.combined_cds_fraction
        
    @metric
    def cdna_length(self):
        '''This property returns the length of the transcript.'''
        return sum([ e[1]-e[0]+1 for e in self.exons ])
    
    @metric
    def number_internal_orfs(self):
        '''This property returns the number of ORFs inside a transcript.'''
        return len(self.internal_orfs)

    @metric
    def selected_cds_length(self):
        '''This property calculates the length of the CDS selected as best inside the cDNA.'''
        if len(self.combined_cds)==0:
            self.__max_internal_orf_length=0
        else:
            self.__max_internal_orf_length=sum(x[2]-x[1]+1 for x in filter(lambda x: x[0]=="CDS", self.selected_internal_orf))
        return self.__max_internal_orf_length
    
    @metric
    def selected_cds_num(self):
        '''This property calculates the number of CDS exons for the selected ORF'''
        return len(list( filter(lambda exon: exon[0]=="CDS", self.selected_internal_orf) ))
    
    @metric
    def selected_cds_fraction(self):
        '''This property calculates the fraction of the selected CDS vs. the cDNA length.'''
        return self.selected_cds_length/self.cdna_length
    
    @metric
    def highest_cds_exons_num(self):
        '''Returns the number of CDS segments in the selected ORF (irrespective of the number of exons involved)'''
        return len(list(filter(lambda x: x[0]=="CDS", self.selected_internal_orf)))
    
    @metric
    def selected_cds_exons_fraction(self):
        '''Returns the fraction of CDS segments in the selected ORF (irrespective of the number of exons involved)'''
        return len(list(filter(lambda x: x[0]=="CDS", self.selected_internal_orf)))/len(self.exons)
    

    @metric
    def highest_cds_exon_number(self):
        '''This property returns the maximum number of CDS segments among the ORFs; this number
        can refer to an ORF *DIFFERENT* from the maximal ORF.'''
        cds_numbers = []
        for cds in self.internal_orfs:
            cds_numbers.append(len(list(filter(lambda x: x[0]=="CDS", cds))))
        return max(cds_numbers)
    
    @metric
    def selected_cds_number_fraction(self):
        '''This property returns the proportion of best possible CDS segments vs. the number of exons.
        See selected_cds_number.'''
        return self.selected_cds_number/self.exon_num

    @metric
    def cds_not_maximal(self):
        '''This property returns the length of the CDS excluded from the selected ORF.'''
        if len(self.internal_orfs)<2:
            return 0
        return self.combined_cds_length-self.selected_cds_length
    
    @metric
    def cds_not_maximal_fraction(self):
        '''This property returns the fraction of bases not in the selected ORF compared to
        the total number of CDS bases in the cDNA.'''
        if self.combined_cds_length==0:
            return 0
        else:
            return self.cds_not_maximal/self.combined_cds_length
    
    @metric
    def five_utr_length(self):
        '''Returns the length of the 5' UTR of the selected ORF.'''
        if len(self.combined_cds)==0:
            return 0
        return sum(x[2]-x[1]+1 for x in self.five_utr)
                            
    @metric
    def five_utr_num(self):
        '''This property returns the number of 5' UTR segments for the selected ORF.'''
        return len(self.five_utr)

    @metric
    def five_utr_num_complete(self):
        '''This property returns the number of 5' UTR segments for the selected ORF, considering only those which are complete exons.'''
        return len(list(filter(lambda utr: (utr[1],utr[2]) in self.exons, self.five_utr)  ))


    @metric
    def three_utr_length(self):
        '''Returns the length of the 5' UTR of the selected ORF.'''
        if len(self.combined_cds)==0:
            return 0
        return sum(x[2]-x[1]+1 for x in self.three_utr)
                            
    @metric
    def three_utr_num(self):
        '''This property returns the number of 3' UTR segments (referred to the selected ORF).'''
        return len(self.three_utr)

    @metric
    def three_utr_num_complete(self):
        '''This property returns the number of 3' UTR segments for the selected ORF, considering only those which are complete exons.'''
        return len(list(filter(lambda utr: (utr[1],utr[2]) in self.exons, self.three_utr)   ))


    @metric
    def utr_num(self):
        '''Returns the number of UTR segments (referred to the selected ORF).'''
        return len(self.three_utr+self.five_utr)

    @metric
    def utr_num_complete(self):
        '''Returns the number of UTR segments which are complete exons (referred to the selected ORF).'''
        return self.three_utr_num_complete+self.five_utr_num_complete


    @metric
    def utr_fraction(self):
        '''This property calculates the length of the UTR of the selected ORF vs. the cDNA length.'''
        return 1-self.selected_cds_fraction

    @metric
    def utr_length(self):
        '''Returns the sum of the 5'+3' UTR lengths'''
        return self.three_utr_length+self.five_utr_length
    
    @metric
    def has_start_codon(self):
        '''Boolean. True if the selected ORF has a start codon.'''
        return self.__has_start
    
    @has_start_codon.setter
    def has_start_codon(self, *args):
        if args[0] not in (None, False,True):
            raise TypeError("Invalid value for has_start_codon: {0}".format(type(args[0])))
        self.__has_start=args[0]
        
    @metric
    def has_stop_codon(self):
        '''Boolean. True if the selected ORF has a stop codon.'''
        return self.__has_stop
    
    @has_stop_codon.setter
    def has_stop_codon(self, *args):
        if args[0] not in (None, False,True):
            raise TypeError("Invalid value for has_stop_codon: {0}".format(type(args[0])))
        self.__has_stop=args[0]

    @metric
    def is_complete(self):
        '''Boolean. True if the selected ORF has both start and end.'''
        return self.has_start_codon and self.has_stop_codon

    @metric
    def exon_num(self):
        '''This property returns the number of exons of the transcript.'''
        return len(self.exons)
    
    @metric
    def exon_fraction(self):
        '''This property returns the fraction of exons of the transcript which are contained in the sublocus.
        If the transcript is by itself, it returns 1. Set from outside.'''
        
        return self.__exon_fraction
    
    @exon_fraction.setter
    def exon_fraction(self, *args):
        if type(args[0]) not in (float,int) or (args[0]<=0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__exon_fraction=args[0]
    
    @metric
    def intron_fraction(self):
        '''This property returns the fraction of introns of the transcript vs. the total number of introns in the locus.
        If the transcript is by itself, it returns 1. Set from outside.'''
        return self.__intron_fraction
    
    @intron_fraction.setter
    def intron_fraction(self, *args):
        if type(args[0]) not in (float,int) or (args[0]<0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        if not self.monoexonic and args[0]==0:
            raise ValueError("It is impossible that the intron fraction is null when the transcript has at least one intron!")
        self.__intron_fraction=args[0]

    @metric
    def max_intron_length(self):
        '''This property returns the greatest intron length for the transcript.'''
        if len(self.introns)==0:
            return 0
        return max(intron[1]+1-intron[0] for intron in self.introns)

    @metric
    def start_distance_from_tss(self):
        '''This property returns the distance of the start of the combined CDS from the transcript start site.
        If no CDS is defined, it defaults to 0.'''
        if len(self.internal_orfs)<2: return self.selected_start_distance_from_tss
        distance=0
        if self.strand=="+" or self.strand is None:
            for exon in self.exons:
                distance+=min(exon[1],self.combined_cds_start-1)-exon[0]+1
                if self.combined_cds_start<=exon[1]:break
        elif self.strand=="-":
            exons=self.exons[:]
            exons.reverse()
            for exon in exons:
                distance+=exon[1]+1-max(self.combined_cds_start+1,exon[0])
                if self.combined_cds_start>=exon[0]:break
        return distance
                
    @metric
    def selected_start_distance_from_tss(self):
        '''This property returns the distance of the start of the best CDS from the transcript start site.
        If no CDS is defined, it defaults to 0.'''
        if len(self.combined_cds)==0: return 0
        distance=0
        if self.strand=="+" or self.strand is None:
            for exon in self.exons:
                distance+=min(exon[1],self.selected_cds_start-1)-exon[0]+1
                if self.selected_cds_start<=exon[1]:break
        elif self.strand=="-":
            exons=self.exons[:]
            exons.reverse()
            for exon in exons:
                distance+=exon[1]+1-max(self.selected_cds_start+1,exon[0])
                if self.selected_cds_start>=exon[0]:break
        return distance

    @metric
    def selected_end_distance_from_tes(self):
        '''This property returns the distance of the end of the best CDS from the transcript end site.
        If no CDS is defined, it defaults to 0.'''
        if len(self.combined_cds)==0: return 0
        distance=0
        if self.strand=="-":
            for exon in self.exons:
                distance+=min(exon[1],self.selected_cds_end-1)-exon[0]+1
                if self.selected_cds_end<=exon[1]:break
        elif self.strand=="+" or self.strand is None:
            exons=self.exons[:]
            exons.reverse()
            for exon in exons:
                distance+=exon[1]+1-max(self.selected_cds_end+1,exon[0])
                if self.selected_cds_end>=exon[0]:break
        return distance

    @metric
    def end_distance_from_tes(self):
        '''This property returns the distance of the end of the combined CDS from the transcript end site.
        If no CDS is defined, it defaults to 0.'''
        if len(self.internal_orfs)<2: return self.selected_end_distance_from_tes
        distance=0
        if self.strand=="-":
            for exon in self.exons:
                distance+=min(exon[1],self.combined_cds_end-1)-exon[0]+1
                if self.combined_cds_end<=exon[1]:break
        elif self.strand=="+" or self.strand is None:
            exons=self.exons[:]
            exons.reverse()
            for exon in exons:
                distance+=exon[1]+1-max(self.combined_cds_end+1,exon[0])
                if self.combined_cds_end>=exon[0]:break
        return distance
    

    @metric
    def combined_cds_intron_fraction(self):
        '''This property returns the fraction of CDS introns of the transcript vs. the total number of CDS introns in the locus.
        If the transcript is by itself, it returns 1.'''
        return self.__combined_cds_intron_fraction
    
    @combined_cds_intron_fraction.setter
    def combined_cds_intron_fraction(self, *args):
        if type(args[0]) not in (float,int) or (args[0]<0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__combined_cds_intron_fraction=args[0]

    @metric
    def selected_cds_intron_fraction(self):
        '''This property returns the fraction of CDS introns of the selected ORF of the transcript vs. the total number of CDS introns in the locus
        (considering only the selected ORF).
        If the transcript is by itself, it should return 1.'''
        return self.__selected_cds_intron_fraction
    
    @selected_cds_intron_fraction.setter
    def selected_cds_intron_fraction(self, *args):
        if type(args[0]) not in (float,int) or (args[0]<0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__selected_cds_intron_fraction=args[0]


    @metric
    def retained_intron_num(self):
        '''This property records the number of introns in the transcripts which are marked as being retained.
        See the corresponding method in the sublocus class.'''
        return len(self.retained_introns)
    
    @metric
    def retained_fraction(self):
        '''This property returns the fraction of the cDNA which is contained in retained introns.'''
        return self.__retained_fraction        

    @retained_fraction.setter
    def retained_fraction(self, *args):
        if type(args[0]) not in (float,int) or (args[0]<0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__retained_fraction=args[0]
    
    @metric
    def proportion_verified_introns(self):
        '''This metric returns, as a fraction, how many of the transcript introns
        are validated by external data.'''
        if self.monoexonic is True:
            return 0
        else:
            return len(self.verified_introns)/len(self.introns)
        
    @metric
    def proportion_verified_introns_inlocus(self):
        '''This metric returns, as a fraction, how many of the verified introns inside the locus
        are contained inside the transcript.'''
        return self.__proportion_verified_introns_inlocus
    
    @proportion_verified_introns_inlocus.setter
    def proportion_verified_introns_inlocus(self, value):
        if value==0:
            assert len(self.verified_introns)==0
        assert 0<=value<=1
        self.__proportion_verified_introns_inlocus=value
        
        
        
         
    
    
