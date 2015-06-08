import sys,argparse,os
sys.path.append(
                os.path.dirname(
                                os.path.dirname(__file__)
                                ))
import csv
from shanghai_lib.loci_objects.transcript import transcript
from loci_objects.exceptions import *
from shanghai_lib.loci_objects.abstractlocus import abstractlocus
from shanghai_lib.parsers.GTF import GTF
from shanghai_lib.parsers.GFF import GFF3


def define_class_code(transcript_instance, reference_instance):
    tran_introns=set(transcript_instance.introns)
    other_introns=set(reference_instance.introns)
    if abstractlocus.overlap((transcript_instance.start,transcript_instance.end),(reference_instance.start,reference_instance.end))==0:
        return "p" #Polymerase run-on
    
    if transcript_instance.strand==reference_instance.strand:
        if transcript_instance.monoexonic is True and reference_instance.monoexonic is True:
            if (transcript_instance.start,transcript_instance.end)==(reference_instance.start,reference_instance.end):
                return "="
            elif (transcript_instance.start >= reference_instance.start) and (transcript_instance.end<=reference_instance.end):
                return "c"
            else:
                return "e" 
        elif transcript_instance.monoexonic is True and reference_instance.monoexonic is False:
            boundaries=(transcript_instance.start,transcript_instance.end)
            for exon in reference_instance.exons:
                if abstractlocus.overlap(boundaries, exon)>0:
                    return "e"
            for intron in reference_instance.introns:
                if abstractlocus.overlap(boundaries, intron)>0:
                    return "i"
            raise ValueError("I have not been able to define a ccode for {0}".format(transcript_instance.id))
        elif reference_instance.monoexonic is True and transcript_instance.monoexonic is False:
            boundaries=(transcript_instance.start,transcript_instance.end)
            
        
        
        elif set.intersection(tran_introns, other_introns)!=set():
            if transcript_instance.introns==reference_instance.introns:
                return "="
            elif len(tran_introns)<len(other_introns) and True:
                pass 
             

def ccode_order(ccode):
    ccodes=["=","c","j","e","i","o","p","r","u","x","s"]
    if ccode not in ccodes:
        raise ValueError("Unknown CCODE: {0}".format(ccode))
    return ccodes.index(ccode)
    


def checkTranscript(transcript_instance, reference, args):
    
    transcript_instance.finalize()
    matches = list(filter( lambda key: abstractlocus.overlap((key[0],key[1]), (transcript_instance.start,transcript_instance.end), flank=args.max_distance)  ))
    
    rows=[]
    row=dict().fromkeys(args.tmap.fieldnames)
    if len(matches)==0:
        row["ref_gene_id"]=row["ref_id"]="-"
        row["ccode"]="u"
        row["gene_id"]=transcript_instance.parent
        row["transcript_id"]=transcript_instance.id
        rows.append(row)
        
    
    for row in rows:
        args.tmap.writerow(row)


def main():
    
    def toRower(string):
        assert os.path.exists(string) and os.path.isfile(string) and os.stat(string).st_size>0
        if ".gtf" in string:
            rower=GTF(string)
        else:
            rower=GFF3(string)
        correct=False
        while True:
            try:
                line=next(rower)
                if line.header is False:
                    correct=True
                    break
            except StopIteration:
                break
            except:
                raise TypeError("Invalid input file!")
        if correct is False:
            raise TypeError("Invalid input file!")
        if ".gtf" in string:
            return GTF(string)
        else:
            return GFF3(string)
    
    parser=argparse.ArgumentParser("Clone of cuffcompare.")
    parser.add_argument("-r","--reference", required=True,
                        type=toRower )
    parser.add_argument("--max_distance", default=2000,
                        help="Maximum distance from a transcript to consider it as a polymerase run-on rather than a new locus. Default: %(default)s")
    parser.add_argument("assemblies", type=toRower  )
    parser.add_argument("out", default="compare", nargs="?", type=str)
    args=parser.parse_args()

    #Destroy previous output files.
    args.tmap=csv.DictWriter(open("{0}.tmap".format(args.out),"w"), ["ref_gene_id", "ref_id", "ccode", "gene_id", "transcript_id" ])    
    reference=dict()
    
    for line in args.reference:
        if line.header: continue
        elif line.is_transcript:
            if line.chrom not in reference:
                reference[line.chrom]=dict()
            reference[line.chrom][(line.start,line.end,line.id)]=transcript(line)
            key=(line.start,line.end,line.id)
        elif line.is_exon:
            if key[-1]!=line.parent:
                try:
                    key=next(filter(lambda x: x[-1]==line.parent, reference[line.chrom]  ))
                except:
                    raise KeyError("Unable to find a transcript for this exon: Parent {0}\n{1}".format(line.parent, line))
            reference[line.chrom][(line.start,line.end,line.id)].addExon(line)
        else:
            continue
        
    for chrom in reference:
        for key in reference[chrom]:
            reference[chrom][key].finalize()
            
    currentTranscript=None
            
    for line in args.assemblies:
        if line.header: continue
        elif line.is_transcript:
            if currentTranscript is not None:
                checkTranscript(currentTranscript, reference, args)
            currentTranscript=transcript(line)
        elif line.is_exon:
            currentTranscript.addExon(line)
        else:
            continue
    if currentTranscript is not None:
        checkTranscript(currentTranscript, reference, args)
    
if __name__=="__main__": main()
            
            
            
            
            
            
    
    
    
    