#!/usr/bin/env python3

#Core imports
import sys,os.path
import sqlalchemy
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import asyncio

#SQLAlchemy imports
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.sql.expression import and_
from sqlalchemy import bindparam
from sqlalchemy.ext import baked

#Shanghai imports
from shanghai_lib.serializers.junction import junction, Chrom
from shanghai_lib.loci_objects.abstractlocus import abstractlocus
from shanghai_lib.loci_objects.monosublocus import monosublocus
from shanghai_lib.loci_objects.excluded_locus import excluded_locus
from shanghai_lib.loci_objects.transcript import transcript
from shanghai_lib.loci_objects.sublocus import sublocus
from shanghai_lib.loci_objects.monosublocus_holder import monosublocus_holder
from shanghai_lib.parsers.GFF import gffLine
from shanghai_lib.exceptions import InvalidLocusError,NoJsonConfigError

class superlocus(abstractlocus):
    
    '''The superlocus class is used to define overlapping regions on the genome, and it receives as input
    transcript class instances.'''
    
    __name__ = "superlocus"
    available_metrics = transcript.get_available_metrics()
    
    bakery = baked.bakery()
    db_baked = bakery(lambda session: session.query(Chrom))
    db_baked += lambda q: q.filter(Chrom.name==bindparam("chrom_name"))
    
    junction_baked = bakery(lambda session: session.query(junction))
    junction_baked += lambda q: q.filter(and_(
                                                            junction.chrom_id==bindparam("chrom_id"),
                                                            junction.junctionStart==bindparam("junctionStart"),
                                                            junction.junctionEnd==bindparam("junctionEnd"),
                                                            junction.strand==bindparam("strand"))
                                         )    
    
    ####### Special methods ############
    
    def __init__(self, transcript_instance, stranded=True, json_dict = None, logger=None):
        
        '''The superlocus class is instantiated from a transcript_instance class, which it copies in its entirety.
        
        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - introns - a *set* which contains the positions of each *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus
        
        The constructor method takes the following keyword arguments:
        - stranded    True if all transcripts inside the superlocus are required to be on the same strand
        - json_dict    Required. A dictionary with the coniguration necessary for scoring transcripts.
        - purge        Flag. If True, all loci holding only transcripts with a 0 score will be deleted from further consideration. 
        
        '''
        
        super().__init__()
        self.stranded=stranded
        self.feature=self.__name__
        if json_dict is None or type(json_dict) is not dict:
            raise NoJsonConfigError("I am missing the configuration for prioritizing transcripts!")
        self.json_dict = json_dict
        self.purge = self.json_dict["run_options"]["purge"]
        
        #Dynamically load required modules
        if "modules" in self.json_dict:
            import importlib
            for mod in self.json_dict["modules"]:
                globals()[mod]=importlib.import_module(mod)
                
        #self.__dict__.update(transcript_instance.__dict__)
        self.splices = set(self.splices)
        self.introns = set(self.introns)
        self.transcripts = dict()
        self.loci=[]
        super().add_transcript_to_locus(transcript_instance)
        if self.stranded is True:
            self.strand = transcript_instance.strand
        self.available_monolocus_metrics = []
        self.available_sublocus_metrics = []
        self.set_flags()
        self.set_logger(logger)

    def __str__(self, level=None, print_cds=True):
        
        '''This function will return the desired level of children loci.
        The keyword level accepts the following four values:
        - "None" - print whatever is available.
        - "loci" - print the final loci
        - "monosubloci" - print the monosubloci
        - "subloci" - print the subloci.
        
        The function will then return the desired location in GFF3-compliant format. 
        '''

        if abs(self.start)==float("inf"): return ''

        superlocus_line=gffLine('')
        superlocus_line.chrom=self.chrom
        superlocus_line.feature=self.__name__
        superlocus_line.start,superlocus_line.end,superlocus_line.score=self.start, self.end, "."
        superlocus_line.strand=self.strand
        superlocus_line.phase, superlocus_line.score=None,None
        new_id = "{0}_{1}".format(self.source,self.id)
        superlocus_line.id,superlocus_line.name=new_id, self.name

        lines=[]
        if level not in (None, "loci", "subloci", "monosubloci"):
            raise ValueError("Unrecognized level: {0}".format(level))
        
        if level=="loci" or (level is None and self.loci_defined is True):
            self.define_loci()
            if len(self.loci)>0:
                source="{0}_loci".format(self.source)
                superlocus_line.source=source
                lines.append(str(superlocus_line))
                found=dict()
                
                for locus_instance in self.loci:
                    locus_instance.source=source
                    locus_instance.parent = new_id
                    if locus_instance.id in found:
                        found[locus_instance.id]+=1
                        locus_instance.counter=found[locus_instance.id]
                    else:
                        found[locus_instance.id]=0
                    lines.append(locus_instance.__str__(print_cds=print_cds).rstrip())
        elif level=="monosubloci" or (level is None and self.monosubloci_defined is True):
            self.define_monosubloci()
            if len(self.monosubloci)>0:
                source="{0}_monosubloci".format(self.source)
                superlocus_line.source=source
                lines.append(str(superlocus_line))
                found=dict()
                for monosublocus_instance in self.monosubloci:
                    monosublocus_instance.source=source
                    monosublocus_instance.parent = new_id
                    if monosublocus_instance.id in found:
                        found[monosublocus_instance.id]+=1
                        monosublocus_instance.counter=found[monosublocus_instance.id]
                    else:
                        found[monosublocus_instance.id]=0
                        
                    lines.append(monosublocus_instance.__str__(print_cds=print_cds).rstrip())
        elif level=="subloci" or (level is None and self.monosubloci_defined is False):
            source="{0}_subloci".format(self.source)
            superlocus_line.source=source
            lines=[str(superlocus_line)]
            self.define_subloci()
            found=dict()
            for sublocus_instance in self.subloci:
                sublocus_instance.source=source
                sublocus_instance.parent = new_id
                if sublocus_instance.id in found:
                    found[sublocus_instance.id]+=1
                    sublocus_instance.counter=found[sublocus_instance.id]
                else:
                    found[sublocus_instance.id]=0
                lines.append(sublocus_instance.__str__(print_cds=print_cds).rstrip())
        
        if len(lines)>0:
            lines.append("###")
        return "\n".join(lines)
            

    ############ Class instance methods ############

    def split_strands(self):
        '''This method will divide the superlocus on the basis of the strand.
        The rationale is to parse a GFF file without regard for the strand, in order to find all intersecting loci;
        and subsequently break the superlocus into the different components.
        Notice that each strand might generate more than one superlocus, if genes on a different strand link what are
        two different superloci.
        '''
        
        self.logger.debug("Splitting by strand for {0}".format(self.id))
        if self.stranded is True:
            self.logger.warn("Trying to split by strand a stranded locus, {0}!".format(self.id))
            yield self
        
        else:
            plus, minus, nones = [], [], []
            for cdna_id in self.transcripts:
                cdna=self.transcripts[cdna_id]
                if cdna.strand == "+":
                    plus.append(cdna)
                elif cdna.strand == "-":
                    minus.append(cdna)
                elif cdna.strand is None:
                    nones.append(cdna)

            new_loci = []
            for strand in plus, minus, nones:
                if len(strand)>0:
                    strand = sorted(strand)
                    new_locus = superlocus(strand[0], stranded=True, json_dict=self.json_dict, logger=self.logger)
                    for cdna in strand[1:]:
                        if new_locus.in_locus(new_locus, cdna):
                            new_locus.add_transcript_to_locus(cdna)
                        else:
                            new_loci.append(new_locus)
                            new_locus = superlocus(cdna, stranded=True, json_dict=self.json_dict, logger=self.logger)
                            
                    new_loci.append(new_locus)
                    
            self.logger.debug("Defined {0} loci by splitting by strand at {1}.".format(len(new_loci), self.id))
            for new_locus in iter(sorted(new_loci)):
                yield new_locus

    def set_flags(self):
        '''Method called by __init__ to set basic flags. These are used throughout the program to avoid unnecessary calculations.'''
        self.subloci_defined = False
        self.monosubloci_defined = False
        self.loci_defined = False
        self.monosubloci_metrics_calculated = False

    #@profile
    def connect_to_db(self):
                
        '''This method will connect to the database using the information contained in the JSON configuration.'''
        
        db = self.json_dict["db"]
        dbtype = self.json_dict["dbtype"]
        
        self.engine = create_engine("{dbtype}:///{db}".format(
                                                              db=db,
                                                              dbtype=dbtype),
                                    connect_args={"check_same_thread": False},
                                    poolclass = sqlalchemy.pool.StaticPool,
                                    )   #create_engine("sqlite:///{0}".format(args.db))
        
        self.sessionmaker = sessionmaker()
        self.sessionmaker.configure(bind=self.engine)
        self.session=self.sessionmaker()

    #@profile   
    @asyncio.coroutine 
    def load_transcript_data(self, tid ):
        '''This routine is used to load data for a single transcript.'''
        
        self.transcripts[tid].set_logger(self.logger)
        yield from self.transcripts[tid].load_information_from_db( self.json_dict, introns = self.locus_verified_introns, session = self.session)
        if self.json_dict["chimera_split"]["execute"] is True and self.transcripts[tid].number_internal_orfs>1:
            new_tr = list(self.transcripts[tid].split_by_cds())
            if len(new_tr)>1:
                for tr in new_tr:
                    self.add_transcript_to_locus(tr, check_in_locus=False)
                self.remove_transcript_from_locus(tid)
        return 
    #@profile
    def load_all_transcript_data(self):
        
        '''This method will load data into the transcripts instances, and perform the split_by_cds if required
        by the configuration.'''
        
        if "db" not in self.json_dict or self.json_dict["db"] is None:
            return #No data to load
        self.connect_to_db()
        self.locus_verified_introns=[]
        dbquery = self.db_baked(self.session).params(chrom_name = self.chrom).all()
        if len(dbquery)>0:
            chrom_id = dbquery[0].id
            for intron in self.introns:
                if len(self.junction_baked(self.session).params(
                                                            chrom_id = chrom_id,
                                                            junctionStart = intron[0],
                                                            junctionEnd = intron[1],
                                                            strand = self.strand
                                                            ).all()) == 1:
                    self.locus_verified_introns.append(intron)

        tids = list(self.transcripts.keys())
        loop = asyncio.get_event_loop()
        tasks = []

        for tid in tids:
            tasks.append(asyncio.async(self.load_transcript_data(tid))) 
        loop.run_until_complete(asyncio.wait(tasks))
        self.session.close()

    ###### Sublocus-related steps ######
               
    def define_subloci(self):
        '''This method will define all subloci inside the superlocus.
        Steps:
            - Call the BronKerbosch algorithm to define cliques
            - Call the "merge_cliques" algorithm the merge the cliques.
            - Create "sublocus" objects from the merged cliques and store them inside the instance store "subloci"       
        '''
        
        self.compile_requirements()
        if self.subloci_defined is True:
            return
        self.subloci = []
        not_passing=set()
        self.excluded_transcripts=None

        if "requirements" in self.json_dict:
            for tid in self.transcripts:
                evaluated=dict()
                for key in self.json_dict["requirements"]["parameters"]:
                    name=self.json_dict["requirements"]["parameters"][key]["name"]
                    value=getattr(self.transcripts[tid],name)
                    evaluated[key]=self.evaluate(value, self.json_dict["requirements"]["parameters"][key])
                
                if eval(self.json_dict["requirements"]["compiled"]) is False:
                    not_passing.add(tid)
                    self.transcripts[tid].score=0
                    
        if len(not_passing)>0 and self.purge is True:
            tid=not_passing.pop()
            self.transcripts[tid].score=0
            monosub=monosublocus(self.transcripts[tid], logger=self.logger)
            self.excluded_transcripts=excluded_locus(monosub, json_dict=self.json_dict, logger=self.logger)
            self.excluded_transcripts.__name__ = "excluded_locus"
            self.remove_transcript_from_locus(tid)
            for tid in not_passing:
                self.transcripts[tid].score=0
                self.excluded_transcripts.add_transcript_to_locus(self.transcripts[tid])
                self.remove_transcript_from_locus(tid)

        if len(self.transcripts)==0:
            #we have removed all transcripts from the locus. Set the flag to True and exit.
            self.subloci_defined=True
            return

        candidates = set(self.transcripts.values()) 
        if len(candidates)==0:
            raise InvalidLocusError("This superlocus has no transcripts in it!")
        subloci = self.find_communities(candidates, inters=self.is_intersecting,
                                        cds_only=self.json_dict["run_options"]["subloci_from_cds_only"])

        #Now we should define each sublocus and store it in a permanent structure of the class
                
        for subl in subloci:
            if len(subl)==0:
                continue
            subl=sorted(subl)
            new_sublocus = sublocus(subl[0], json_dict=self.json_dict, logger=self.logger)
            for ttt in subl[1:]:
                new_sublocus.add_transcript_to_locus(ttt)
            new_sublocus.parent = self.id
            self.subloci.append(new_sublocus)
        self.subloci=sorted(self.subloci)
        
        self.subloci_defined = True

    def get_sublocus_metrics(self):
        '''Wrapper function to calculate the metrics inside each sublocus.'''
        
        self.define_subloci()
        self.sublocus_metrics = []
        for sublocus_instance in self.subloci:
            sublocus_instance.get_metrics()

    def define_monosubloci(self):

        '''This is a wrapper method that defines the monosubloci for each sublocus.
        '''
        if self.monosubloci_defined is True:
            return
        
        self.define_subloci()
        self.monosubloci=[]
        #Extract the relevant transcripts
        for sublocus_instance in sorted(self.subloci):
            self.excluded_transcripts=sublocus_instance.define_monosubloci(purge=self.purge, excluded=self.excluded_transcripts)
            for ml in sublocus_instance.monosubloci:
                ml.parent = self.id
                self.monosubloci.append(ml)
        self.monosubloci = sorted(self.monosubloci)
        self.monosubloci_defined = True

    def print_monolocus_metrics(self, rower):
        '''Wrapper function to pass to a csv.DictWriter object the metrics of the transcripts in the monosubloci.'''
        
        raise NotImplementedError()

    def print_monolocus_scores(self, rower):
        '''Wrapper function to pass to a csv.DictWriter object the metrics of the transcripts in the monosubloci.'''
        
        raise NotImplementedError()


    def print_subloci_metrics(self ):
        
        '''Wrapper method to create a csv.DictWriter instance and call the sublocus.print_metrics method
        on it for each sublocus.'''
        
        self.get_sublocus_metrics()
        
        for slocus in self.subloci:
            for row in slocus.print_metrics():
                yield row
        if self.excluded_transcripts is not None:
            for row in self.excluded_transcripts.print_metrics():
                yield row

    def print_subloci_scores(self ):
        
        '''Wrapper method to create a csv.DictWriter instance and call the sublocus.print_metrics method
        on it for each sublocus.'''
        
        self.get_sublocus_metrics()
        
        for slocus in self.subloci:
            for row in slocus.print_scores():
                yield row
        if self.excluded_transcripts is not None:
            for row in self.excluded_transcripts.print_scores():
                yield row


    def print_monoholder_metrics(self ):

        '''Wrapper method to create a csv.DictWriter instance and call the monosublocus_holder.print_metrics method
        on it.'''
        
        
        self.define_loci()

        #self.available_monolocus_metrics = set(self.monoholder.available_metrics)
        if len(self.monoholders)==0:
            yield ''
        for monoholder in self.monoholders:
            for row in monoholder.print_metrics():
                yield row

    def print_monoholder_scores(self ):

        '''Wrapper method to create a csv.DictWriter instance and call the monosublocus_holder.print_scores method
        on it.'''
        
        
        self.define_loci()

        #self.available_monolocus_metrics = set(self.monoholder.available_metrics)
        if len(self.monoholders)==0:
            yield ''
        for monoholder in self.monoholders:
            for row in monoholder.print_scores():
                yield row


            
    def define_loci(self):
        '''This is the final method in the pipeline. It creates a container for all the monosubloci
        (an instance of the class monosublocus_holder) and retrieves the loci it calculates internally.'''
        
        if self.loci_defined is True:
            return
        
        self.define_monosubloci()
        self.calculate_mono_metrics()

        self.loci = []        
        if len(self.monoholders)==0:
            self.loci_defined = True
            return
        for monoholder in self.monoholders:
            monoholder.define_loci(purge=self.purge)
            for locus_instance in monoholder.loci:
                locus_instance.parent = self.id
                self.loci.append(locus_instance)
            
        self.loci=sorted(self.loci)
        self.loci_defined = True
        
        return
    
    def calculate_mono_metrics(self):
        '''Wrapper to calculate the metrics for the monosubloci.'''
        self.monoholders = []
        
        for monosublocus_instance in sorted(self.monosubloci):
            if self.monoholders == []:
                holder = monosublocus_holder(monosublocus_instance, json_dict=self.json_dict, logger=self.logger)
                self.monoholders.append(holder)
            else:
                found_holder=False
                for holder in self.monoholders:
                    if monosublocus_holder.in_locus(holder, monosublocus_instance):
                        holder.add_monosublocus(monosublocus_instance)
                        found_holder=True
                        break
                if found_holder is False:
                    holder = monosublocus_holder(monosublocus_instance, json_dict=self.json_dict, logger=self.logger)
                    self.monoholders.append(holder)
                
    def compile_requirements(self):
        '''Quick function to evaluate the filtering expression, if it is present.'''
        
        if "requirements" in self.json_dict:
            if "compiled" in self.json_dict["requirements"]:
                return
            else:
                self.json_dict["requirements"]["compiled"]=compile(self.json_dict["requirements"]["expression"], "<json>", "eval")
                return
        else:
            return

    ############# Class methods ###########
    
    @classmethod
    def is_intersecting(cls,transcript_instance, other, cds_only=False):
        '''When comparing two transcripts, for the definition of subloci inside superloci we follow these rules:
        If both are multiexonic, the function verifies whether there is at least one intron in common.
        If both are monoexonic, the function verifies whether there is some overlap between them.
        If one is monoexonic and the other is not,  the function will return False by definition.        
        '''
        
        transcript_instance.finalize()
        other.finalize()
        if transcript_instance.id==other.id:
            return False # We do not want intersection with oneself
#         monoexonic_check = len( list(filter(lambda x: x.monoexonic is True, [transcript_instance, other]   )  )   )
        
        if transcript_instance.monoexonic is False and other.monoexonic is False:
            if cds_only is False:
                intersection = set.intersection(transcript_instance.introns, other.introns)
            else:
                intersection = set.intersection(transcript_instance.combined_cds_introns,
                                                other.combined_cds_introns)
            if len(intersection)>0:
                return True
            else:
                return False

        elif transcript_instance.monoexonic is True and other.monoexonic is True:
            if transcript_instance.start==other.start:
                return True
            if transcript_instance.end==other.end:
                return True
            test_result =cls.overlap(
                           (transcript_instance.start, transcript_instance.end),
                           (other.start, other.end)
                           ) 
            if test_result>0: #A simple overlap analysis will suffice
                return True
            else:
                return False
        else:
            return False
    
