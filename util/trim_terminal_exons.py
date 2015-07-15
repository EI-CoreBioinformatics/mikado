import sys,argparse
import shanghai_lib.parsers
import shanghai_lib.loci_objects.transcript

def strip_terminal(currentTranscript: shanghai_lib.loci_objects.transcript.transcript, args: argparse.Namespace) -> shanghai_lib.loci_objects.transcript.transcript:
    '''This function will take as input a transcript and then:
    - return immediately if the transcript is monoexonic
    - trim the terminal exons to a length of max(args.max_length (if longer), CDS boundary)
    '''
    
    currentTranscript.finalize()
    #Return immediately if the transcript is monoexonic
    if currentTranscript.monoexonic is True:
        return currentTranscript

    currentTranscript.finalized = False
    
    if currentTranscript.selected_cds_length == 0:
        #Non-coding transcript: trim the terminal exons and the return
        first = list(currentTranscript.exons[0])
        if (first[1]-first[0]+1)>args.max_length:
            newfirst=list(first)
            newfirst[0]=first[1]-args.max_length
            currentTranscript.start = newfirst[0]
            currentTranscript.exons[0]=tuple(newfirst)
        last =  currentTranscript.exons[-1]
        if (last[1]-last[0]+1)>args.max_length:
            newlast=list(last[:])
            newlast[1]=last[0]+args.max_length
            currentTranscript.end = newlast[1]
            currentTranscript.exons[-1]=tuple(newlast)
            if currentTranscript.selected_cds_length>0:
                last_utr = list(filter(lambda u: u[1]==last[1], currentTranscript.combined_utr))[0]
                currentTranscript.combined_utr.remove(last_utr)
                currentTranscript.combined_utr.append(tuple(newlast))

    else:
        #Coding transcript
        #Order cds_start and end irrespectively of strand
        cds_start,cds_end = sorted([currentTranscript.selected_cds_end, currentTranscript.selected_cds_start])
        exons = currentTranscript.exons
        first=exons[0]
        last=exons[-1]
        if first[1]<cds_start: #first exon is only UTR
            if first[1]-first[0]+1>args.max_length:
                # First exon longer than max_length; trim it, remove from UTR and exons, substitute with trimmed exon
                newfirst = (first[1]-args.max_length, first[1])
                currentTranscript.combined_utr.remove(first)
                currentTranscript.exons.remove(first)
                currentTranscript.combined_utr.append(newfirst)
                currentTranscript.exons.append(newfirst)
            else:
                #Leave as it is
                newfirst = first
        elif first[0]<cds_start<=first[1]:
            #Partly UTR partly CDS
            if first[1]-first[0]+1>args.max_length:
                newfirst = ( min( cds_start, first[1]-args.max_length   ) , first[1]  )
                u = list(filter(lambda u: u[0]==first[0], currentTranscript.combined_utr))[0] #Retrieve and remove UTR segment
                currentTranscript.combined_utr.remove(u)
                currentTranscript.exons.remove(first)
                currentTranscript.exons.append(newfirst)
                if newfirst[0]<cds_start:
                    #Create new UTR segment
                    newu = (newfirst[0], cds_start-1)
                    currentTranscript.combined_utr.append(newu)
            else:
                newfirst=first
        else: #Beginning of first exon == beginning of CDS
            newfirst=first
        currentTranscript.start = newfirst[0]
        
        if last[0]>cds_end:
            if last[1]-last[0]+1>args.max_length:
                newlast = (last[0], last[0]+args.max_length)
                currentTranscript.combined_utr.remove(last)
                currentTranscript.exons.remove(last)
                currentTranscript.combined_utr.append(newlast)
                currentTranscript.exons.append(newlast)
            else:
                newlast=last
        elif last[0]<=cds_end<last[1]:
            if last[1]-last[0]+1>args.max_length:
                newlast = (last[0],max( cds_end, last[0]+args.max_length  ))
                u = list(filter(lambda u: u[1]==last[1], currentTranscript.combined_utr))[0]
                currentTranscript.combined_utr.remove(u)
                currentTranscript.exons.remove(last)
                currentTranscript.exons.append(newlast)
                if newlast[1]>cds_end:
                    newu = (cds_end+1,newlast[1])
                    currentTranscript.combined_utr.append(newu)
            else:
                newlast=last
        else:
            newlast=last
        
        currentTranscript.end = newlast[1]
            
    currentTranscript.finalize()        
    return currentTranscript


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
    parser.add_argument("--keep_cds", action="store_true", default=False, help="Keep the CDS information.")
    parser.add_argument("ann", type=to_ann)
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType('w'))
    args=parser.parse_args()
    

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