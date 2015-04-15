#!/usr/bin/env python3

from abstractlocus import abstractlocus

class superlocus(abstractlocus):
    
    '''The superlocus class is used to define overlapping regions on the genome, and it receives as input
    transcript class instances.'''
    
    def __init__(self, transcript):
        
        '''The superlocus class is instantiated from a transcript class, which it copies in its entirety.
        
        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - junctions - a *set* which contains the positions of each *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus'''
        
        transcript.finalize()
        self.__dict__.update(transcript.__dict__)
        self.splices = set(self.splices)
        self.junctions = set(self.junctions)
        self.transcripts = set()
        self.transcripts.add(transcript)
        return
    
    @classmethod
    def is_intersecting(cls,transcript, other):
        '''I use this class method to verify that two transcripts are intersecting.
        If both are multiexonic, this amounts to checking whether there is at least one intron in common.
        If both are monoexonic, this amounts to checking whether there is some overlap between them.
        If one is monoexonic and the other is not, this function will return False by definition.        
        '''
        
        if transcript.id==other.id: return False # We do not want intersection with oneself
        monoexonic_check = len( list(filter(lambda x: x.monoexonic is True, [transcript, other]   )  )   )
        
        if monoexonic_check==0: #Both multiexonic
            inter=set.intersection(transcript.junctions, other.junctions )
            #print(transcript.id, other.id, inter)
            if len(inter)>0: return True
        
        elif monoexonic_check==1: #One monoexonic, the other multiexonic: different subloci by definition
            return False
        
        elif monoexonic_check==2:
            if cls.overlap(transcript, other)>0: #A simple overlap analysis will suffice
                return True
        
        return False
    
    def define_subloci(self):
        '''This method will define all subloci inside the superlocus.
        The method performs multiple calls to the BronKerbosch class method to define the possible groups of transcripts.        
        '''
        candidates = set(self.transcripts) # This will order the transcripts based on their position
        if len(candidates)==0:
            raise ValueError("This superlocus has no transcripts in it!")
        
        subloci = dict()
        index=0
        
        while len(candidates)>0:
            result=self.BronKerbosch(set(), candidates, set())
            index+=1
            subloci[index]=result
            candidates = set.difference( candidates, result )

        #Now we should define each sublocus and store it in a permanent structure of the class
        self.subloci = []
    
        for sublocus in subloci:
            transcripts = list(subloci[sublocus])
            if len(transcripts)==0:
                continue
            new_superlocus = superlocus(transcripts[0])
            if len(transcripts)>1:
                for ttt in transcripts[1:]:
                    new_superlocus.add_transcript_to_locus(ttt)
                    
            self.subloci.append((new_superlocus, transcripts[0].monoexonic))
    
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
        superlocus_id = "superlocus:{0}{3}:{1}-{2}".format(self.chrom, self.start, self.end,strand)

        superlocus_line = [self.chrom, "locus_pipeline", "superlocus", self.start, self.end, ".", strand, ".", "ID={0}".format(superlocus_id) ]
        superlocus_line = "\t".join(str(s) for s in superlocus_line)
        
        sublocus_lines = []
        counter=0
        for sublocus in self.subloci:
            counter+=1
            sublocus, monoexonic = sublocus
            sublocus_id = "{0}.{1}".format(superlocus_id, counter)
            attr_field = "ID={0};Parent={1};".format(sublocus_id, superlocus_id)
            if monoexonic is True:
                tag="multiexonic=false;"
            else:
                tag="multiexonic=true;"
            attr_field="{0}{1}".format(attr_field, tag)
            if sublocus.strand is None:
                substrand="."
            else:
                substrand=sublocus.strand
            sublocus_line = [ self.chrom, "locus_pipeline", "sublocus", sublocus.start, sublocus.end,
                             ".", substrand, ".", attr_field]
            sublocus_line = "\t".join([str(s) for s in sublocus_line])
            
            sublocus_lines.append(sublocus_line)
            for transcript in sublocus.transcripts:
                transcript.parent=sublocus_id
                sublocus_lines.append(str(transcript).rstrip())
                
        sublocus_lines="\n".join(sublocus_lines)
        lines=[superlocus_line]
        lines.append(sublocus_lines)
        try:
            return "\n".join(lines)
        except TypeError:
            raise TypeError(lines)
            
