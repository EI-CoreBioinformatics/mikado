#!/usr/bin/env python3

from abstractlocus import abstractlocus
import operator
from copy import copy
from sublocus import sublocus

class superlocus(abstractlocus):
    
    '''The superlocus class is used to define overlapping regions on the genome, and it receives as input
    transcript class instances.'''
    
    __name__ = "superlocus"
    
    ####### Special methods ############
    
    def __init__(self, transcript_instance, stranded=True):
        
        '''The superlocus class is instantiated from a transcript_instance class, which it copies in its entirety.
        
        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - junctions - a *set* which contains the positions of each *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus'''
        
        transcript_instance.finalize()
        self.stranded=stranded
        self.feature="superlocus"
        self.__dict__.update(transcript_instance.__dict__)
        self.splices = set(self.splices)
        self.junctions = set(self.junctions)
        self.transcripts = dict()
        super().add_transcript_to_locus(transcript_instance)
        self.available_metrics = []
        self.monosubloci_defined = False
        self.loci_defined = False
        return

    def __str__(self):
        
        '''Before printing, the class calls the define_subloci method. It will then print:
        # a "superlocus" line
        # for each "sublocus":
        ## a "sublocus" line
        ## all the transcripts inside the sublocus (see the transcript class)'''

        if self.strand is not None:
            strand=self.strand
        else:
            strand="."
        
        self.define_subloci()
        superlocus_line = [self.chrom, "locus_pipeline", "superlocus", self.start, self.end, ".", strand, ".", "ID={0}".format(self.id) ]
        superlocus_line = "\t".join(str(s) for s in superlocus_line)
        counter=0
        
        lines=[superlocus_line]        
        for subl in iter(sorted(self.subloci, key=operator.attrgetter("start","end") )):
            counter+=1
            subl.parent = self.id
            lines.append(str(subl).rstrip())

        try:
            return "\n".join(lines)
        except TypeError:
            raise TypeError(lines)    

    ############ Class instance methods ############

    def split_strands(self):
        '''This method will divide the superlocus on the basis of the strand.
        The rationale is to parse a GFF file without regard for the strand, in order to find all intersecting loci;
        and subsequently break the superlocus into the different components.
        '''
        
        if self.stranded is True:
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

            for strand in plus, minus, nones:
                if len(strand)>0:
                    new_locus = superlocus(strand[0], stranded=True)
                    for cdna in strand[1:]:
                        new_locus.add_transcript_to_locus(cdna)
                    yield new_locus
                    
    def define_subloci(self):
        '''This method will define all subloci inside the superlocus.
        Steps:
            - Call the BronKerbosch algorithm to define cliques
            - Call the "merge_cliques" algorithm the merge the cliques.
            - Create "sublocus" objects from the merged cliques and store them inside the instance store "subloci"       
        '''
        
        candidates = set(self.transcripts.values()) # This will order the transcripts based on their position
        if len(candidates)==0:
            raise ValueError("This superlocus has no transcripts in it!")
        
        
        original=copy(candidates)
        
        cliques = set( tuple(clique) for clique in self.BronKerbosch(set(), candidates, set(), original))
        
        subloci = self.merge_cliques(cliques)
        
        #Now we should define each sublocus and store it in a permanent structure of the class
        self.subloci = []
    
        for subl in subloci:
            if len(subl)==0:
                continue
            
            new_sublocus = sublocus(subl[0])
            if len(subl)>1:
                for ttt in subl[1:]:
                    new_sublocus.add_transcript_to_locus(ttt)
                    
            self.subloci.append(new_sublocus)

    def define_monosubloci(self):

        '''This function calculates the monosubloci from the subloci present in the superlocus.'''
        if self.monosubloci_defined is True:
            return
        
        self.define_subloci()
        self.monosubloci = []
        self.monosubloci_transcripts = dict()
        self.monosubloci_exons=set()
        self.monosubloci_splices=set()
        self.monosubloci_junctions=set()
        #Extract the relevant transcripts
        for sublocus_instance in sorted(self.subloci):
            sublocus_instance.define_monosubloci()
            for monosublocus_instance in sublocus_instance.monosubloci:
                tid = monosublocus_instance.tid
                self.monosubloci_transcripts[tid]=self.transcripts[tid]
                self.monosubloci_exons.update(self.transcripts[tid].exons)
                self.monosubloci_splices.update(self.transcripts[tid].splices)
                self.monosubloci_junctions.update(self.transcripts[tid].junctions)
                
            self.monosubloci.extend(sublocus_instance.monosubloci)
        self.monosubloci_defined = True


    def get_monolocus_metrics(self):
        self.define_monosubloci()
        for tid in self.monosubloci_transcripts:
            self.calculate_monolocus_metrics(tid)


            
    def define_loci(self):
        self.get_monolocus_metrics()
        raise NotImplementedError()
            
        #Calculate the best transcripts
            
        
    def calculate_monolocus_metrics(self, tid):
        '''This function will calculate the metrics which will be used to derive a score for a transcript.
        The different attributes of the transcript will be stored inside the transcript class itself,
        to be more precise, into an internal dictionary ("metrics").
         
        This function re-generates only the following metrics, as they are relative:
        The scoring function must consider the following factors:
        - "exon_frac":          % of exons on the total of the exons present in the monosubloci
        - "intron_frac":        % of introns on the total of the introns present in the monosubloci
        - "retained_introns":   no. of retained introns
        - "retained_frac":      % of cdna_length that is in retained introns 
        '''

        transcript_instance = self.monosubloci_transcripts[tid]
        transcript_instance.metrics["exon_frac"] = len(set.intersection( self.monosubloci_exons,transcript_instance.exons   ))/len(self.monosubloci_exons)
        if len(self.monosubloci_junctions)>0:
            transcript_instance.metrics["intron_frac"] = len(set.intersection( self.monosubloci_junctions,transcript_instance.junctions   ))/len(self.monosubloci_junctions)
        else:
            transcript_instance.metrics["intron_frac"]=None
        self.find_retained_introns_in_monosubloci(transcript_instance)
        transcript_instance.metrics["retained_introns"] = len(transcript_instance.retained_introns)
        transcript_instance.metrics["retained_frac"]=sum(e[1]-e[0]+1 for e in transcript_instance.retained_introns)/transcript_instance.cdna_length
        for metric in transcript_instance.metrics:
            if metric not in self.available_metrics:
                self.available_metrics.append(metric)

        self.transcripts[tid]=transcript_instance

    def find_retained_introns_in_monosubloci(self, transcript_instance):
         
        '''This method checks the number of exons that are possibly retained introns for a given transcript.
        To perform this operation, it checks for each non-CDS exon whether it exists a sublocus intron that
        is *completely* contained within a transcript exon.
        CDS exons are ignored because their retention might be perfectly valid.
        The results are stored inside the transcript instance, in the "retained_introns" tuple.'''
         
        transcript_instance.retained_introns=[]
        for exon in filter(lambda e: e not in transcript_instance.cds, transcript_instance.exons):
            #Check that the overlap is at least as long as the minimum between the exon and the intron.
            if any(filter(
                          lambda junction: self.overlap(exon,junction)>=junction[1]-junction[0],
                          self.monosubloci_junctions
                          )) is True:
                    transcript_instance.retained_introns.append(exon)
#                     print("Retained:", transcript_instance.id, exon, file=sys.stderr)
        transcript_instance.retained_introns=tuple(transcript_instance.retained_introns)

    ############# Class methods ###########
    
    @classmethod
    def is_intersecting(cls,transcript, other):
        '''When comparing two transcripts, for the definition of subloci inside superloci we follow these rules:
        If both are multiexonic, the function verifies whether there is at least one intron in common.
        If both are monoexonic, the function verifies whether there is some overlap between them.
        If one is monoexonic and the other is not,  the function will return False by definition.        
        '''
        
        if transcript.id==other.id: return False # We do not want intersection with oneself
        monoexonic_check = len( list(filter(lambda x: x.monoexonic is True, [transcript, other]   )  )   )
        
        flag=False
        if monoexonic_check==0: #Both multiexonic
            for junc in transcript.junctions:
                if junc in other.junctions:
                    flag=True
                    break
                    
            else:
                flag=False
#             return any( filter( lambda j: j in other.junctions, transcript.junctions ) )
        
        elif monoexonic_check==1: #One monoexonic, the other multiexonic: different subloci by definition
            flag=False
        
        elif monoexonic_check==2:
            if cls.overlap(
                           (transcript.start, transcript.end),
                           (other.start, other.end)
                           )>=0: #A simple overlap analysis will suffice
                flag=True
        return flag