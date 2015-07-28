#!/usr/bin/env python3

import sys,argparse
import mikado_lib.parsers
import mikado_lib.loci_objects

def main():
    
    parser=argparse.ArgumentParser("GFF=>BED12 converter")
    parser.add_argument("gff", type=argparse.FileType("r"))
    parser.add_argument("out", nargs="?", type=argparse.FileType("w"), default=sys.stdout)
    args=parser.parse_args()
    
    currentTranscript=None
    
    for row in mikado_lib.parsers.GFF.GFF3(args.gff):
        if row.is_parent and not row.is_transcript:
            continue
        if row.is_exon:
            currentTranscript.add_exon(row)
        elif row.is_transcript:
            if currentTranscript is not None:
                currentTranscript.finalize()
                assert len(currentTranscript.introns)==1
                introns=sorted(currentTranscript.introns)
                bed12=mikado_lib.parsers.bed12.BED12()
                bed12.chrom=currentTranscript.chrom
                bed12.strand=currentTranscript.strand
                bed12.score=currentTranscript.score
                bed12.start=currentTranscript.start
                bed12.end=currentTranscript.end
                bed12.thickStart=introns[0][0]
                bed12.thickEnd=introns[0][1]
                bed12.name=currentTranscript.id
                bed12.rgb="255,0,0"
                bed12.blockCount=2
                bed12.blockSizes=[(e[1]-e[0]) for e in currentTranscript.exons  ]
                bed12.blockStarts=[]
                bed12.blockStarts=[0, bed12.blockSizes[0]+introns[0][1]-introns[0][0] ] 
                print(bed12, file=args.out)
                                           
            currentTranscript=mikado_lib.loci_objects.transcript.Transcript(row)
            
    if currentTranscript is not None:
        currentTranscript.finalize()
        assert len(currentTranscript.introns)==1
        introns=sorted(currentTranscript.introns)
        bed12=mikado_lib.parsers.bed12.BED12()
        bed12.chrom=currentTranscript.chrom
        bed12.strand=currentTranscript.strand
        bed12.start=currentTranscript.start
        bed12.score=currentTranscript.score
        bed12.end=currentTranscript.end
        bed12.thickStart=introns[0][0]
        bed12.thickEnd=introns[0][1]
        bed12.name=currentTranscript.id
        bed12.rgb="255,0,0"
        bed12.blockCount=2
        bed12.blockSizes=[(e[1]-e[0]) for e in currentTranscript.exons  ]
        bed12.blockStarts=[]
        bed12.blockStarts=[0, bed12.blockSizes[0]+introns[0][1]-introns[0][0] ] 
        print(bed12, file=args.out)

if __name__=="__main__": main()
  
            