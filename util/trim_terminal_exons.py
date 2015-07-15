import sys,argparse
import shanghai_lib.parsers
import shanghai_lib.loci_objects.transcript

def main():
    
    def to_ann(string):
        if string.endswith("gtf"):
            return shanghai_lib.parsers.GTF.GTF(string)
        elif string.endswith("gff") or string.endswith("gff3"):
            return shanghai_lib.parsers.GFF.GFF3(string)
        else:
            raise ValueError("Unrecognized format")
    
    parser=argparse.ArgumentParser("Script to trim down the terminal exons of multiexonic transcripts")
    parser.add_argument("-ml", "--max_length", type=int, default=50, help="Maximmal length of terminal exons")
    parser.add_argument("ann", type=to_ann)
    parser.add_argument("out", nargs="?", default=sys.stdout)
    args=parser.parse_args()
    

    def strip_terminal(currentTranscript, args):
        currentTranscript.strip_cds()
        if currentTranscript.monoexonic is False:
            currentTranscript.finalized = False
            first = list(currentTranscript.exons[0])
            if (first[1]-first[0]+1)>args.max_length:
                first[0]=first[1]-args.max_length
                currentTranscript.start = first[0]
                currentTranscript.exons[0]=tuple(first)
            last =  list(currentTranscript.exons[-1])
            if (last[1]-last[0]+1)>args.max_length:
                last[1]=last[0]+args.max_length
                currentTranscript.end = last[1]
                currentTranscript.exons[-1]=tuple(last)
            currentTranscript.finalize()
        return currentTranscript
    
    currentTranscript=None
    
    for record in args.ann:
        if record.is_transcript is True:
            if currentTranscript is not None:
                print(strip_terminal(currentTranscript, args), file=args.out)
            currentTranscript=shanghai_lib.loci_objects.transcript.transcript(record)
        elif record.is_exon is True:
            currentTranscript.addExon(record)
            
    print(strip_terminal(currentTranscript, args), file=args.out)
    
if __name__ == "__main__": main()