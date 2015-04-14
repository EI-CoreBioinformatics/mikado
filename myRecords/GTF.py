#coding: utf_8

from myRecords import Parser,HeaderError
import re
import io

class gtfLine(object):
    '''This class defines a typical GTF line, with some added functionality to make it useful in e.g. parsing cufflinks GTF files or creating GTF lines from scratch.
    Fields:
    - chrom
    - source
    - feature
    - start,stop
    - score
    - phase
    - strand
    - info (a dictionary containing all the annotations contained in the last field)

    - gene: the gene_id
    - transcript: the transcript_id

    For cufflinks files, also:
    - nearest_ref
    - tss_id
    - ccode'''

    _slots=['chrom','source','feature','start',\
    'stop','score','strand','phase','info']

    def __init__(self,line):
        if line==None or line[0]=="#":
            for i in self._slots:
                self.__dict__[i]=None
            self.fields=[]
            self.info={'gene_id': None, 'transcript_id': None}
            self.transcript=None
            self.gene=None


        else:
            assert isinstance(line,str)
            self.fields=line.rstrip().split('\t')
            self.chrom,self.source,self.feature=self.fields[0:3]
            self.start,self.stop=tuple(int(i) for i in self.fields[3:5])
            try: self.score=float(self.fields[5])
            except ValueError:
                if self.fields[5]=='.': self.score=None
            self.strand,self.phase=self.fields[6:8]
            if self.strand in ('.','\x00'):
                self.strand=None
            assert self.strand in ('+',None,'-'), (self.strand,line)
            if self.phase=='.': self.phase=None
            else:
                try:
                    self.phase=int(self.phase)
                    assert self.phase in (0,1,2)
                except: raise

            self._info=self.fields[8].split(';')
            try: self._info.remove('')
            except ValueError: pass
            self.info=dict()
            for info in self._info:
                info=info.lstrip().split(' ') #Divido l'annotazione, e rimuovo le virgolette
                try: self.info[info[0]]=re.sub(r'"','',info[1])
                except IndexError: raise IndexError(self._info, info)

            if 'exon_number' in self.info: self.info['exon_number']=int(self.info['exon_number'])
            assert 'gene_id','transcript_id' in self.info
            self.gene=None #Questa parte serve per le annotazioni di cufflinks.
            self.transcript=None
            self.nearest_ref=None
            self.tss_id=None
            self.ccode=None
            if 'gene_id' in self.info: self.gene=self.info['gene_id']
            if 'transcript_id' in self.info: self.transcript= self.info['transcript_id']
            if 'nearest_ref' in self.info: self.nearest_ref= self.info['nearest_ref']
            if 'tss_id' in self.info: self.tss_id= self.info['tss_id']
            if 'class_code' in self.info: self.ccode=self.info['class_code']
            for tag in [x for x in list(self.info.keys()) if x not in ('gene_id','transcript_id','nearest_ref','tss_id','class_code')]:
                                self.__dict__[tag.lower()]=self.info[tag]
                                
          
    def __str__(self):
        '''Returns the GTF string.'''
        self.fields=[]
        if self.fields==[]: #Devo ricostruire i campi se ho inizializzato da una riga vuota
            self.fields=[self.chrom,self.source,self.feature,
                         self.start,self.stop,self.score,
                         self.strand, self.phase]
            for i in range(len(self.fields)):
                if self.fields[i]==None: self.fields[i]='.' #Correggere per campi vuoti
                self.fields[i]=str(self.fields[i])
            self._info=[]
            assert 'gene_id','transcript_id' in self.info
            if self.info['gene_id']==None: self.info['gene_id']=self.gene
            if self.info['transcript_id']==None: self.info['transcript_id']=self.transcript

            order=['gene_id','transcript_id','exon_number','gene_name','transcript_name'] #Questo è l'ordine originale dei campi nel gtf di umano

            for tag in order:
                try: self._info.append(tag+' "'+str(self.info[tag])+'"')
                except KeyError: pass

            for info in [x for x in list(self.info.keys()) if x not in order]: #Aggiungo eventuali altri campi
                self._info.append(info+' "'+self.info[info]+'"')
            self.fields.append('; '.join(self._info))
            self.fields[-1]+=';' #Fields finito, si può stampare.

        assert self.fields[0]!="", self.fields
        return '\t'.join(self.fields)

    def toGFF3(self, feature_type="gene", source=None):
        '''Converts the GTF line into a GFF3 one.'''
        attributes=[]

        #Redefine source variable if asked
        if source is not None:
            self.source=source
        
        if self.feature=='gene':
            if feature_type=="match": return
            if 'gene_name' not in self.info: self.info['gene_name']=self.info['gene_id']
            if 'description' not in self.info: self.info['description']='NA'
            attributes=['ID='+self.info['gene_id'],'Name='+self.info['gene_name']]

        elif self.feature=='mRNA' or self.feature=="transcript":
            if feature_type=="gene":
                if 'transcript_name' not in self.info: self.info['transcript_name']=self.info['transcript_id']
                attributes=['ID='+self.info['transcript_id'],
                            'Parent='+self.info['gene_id'],
                            'Name='+self.info['transcript_name']]
            elif feature_type=="match":
                self.feature="match"
                attributes=["ID={0}".format(self.info["transcript_id"]),
                            "Name={0}".format(self.info["transcript_id"])]
                
        elif self.feature in ('exon','CDS'):

#            assert 'exon_number' in self.info
            if feature_type=="match":
                self.feature="match_part"

            if "exon_number" in self.info:
                attributes=['ID={0}:exon-{1}'.format(self.transcript, self.info["exon_number"]),
                            'Parent={0}'.format(self.info['transcript_id'])]
            else:
                attributes=['Parent={0}'.format(self.info['transcript_id'])]
            

        elif self.feature in ('UTR', 'five prime UTR', 'three prime UTR'):
            if self.feature=='UTR': raise ValueError('I cannot work with "UTR" only currently! Error in: {0}'.format(self.info['transcript_id'])) #I have to think about a smart way of doing this..
            if self.feature=='five prime UTR': ut='5'
            else: ut='3'
            attributes=['ID=utr.'+ut,
                        'Parent='+self.info['transcript_id']]

        elif self.feature in ('start_codon','stop_codon'):
            attributes=['ID='+self.feature,
                        'Parent='+self.info['transcript_id']]

        if self.score==None: score='.'
        else: score=str(self.score)

        if self.phase==None: phase='.'
        else: phase=str(self.phase)

        if self.strand==None: strand='.'
        else: strand=self.strand

        start=min(self.start,self.stop)
        stop=max(self.start,self.stop)

        line=[self.chrom,self.source,self.feature,str(start),str(stop),score,strand,phase]
        line+=[';'.join(attributes)]
        return '\t'.join(line)

class GTF(Parser):
    '''The parsing class.'''
    def __init__(self,handle): super(GTF, self).__init__(handle)
        
    def __next__(self):
        line=self._handle.readline()
        if line=='': raise StopIteration
        return gtfLine(line)
