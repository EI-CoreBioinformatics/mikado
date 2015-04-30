#!/usr/bin/env python3

import sys
import argparse,os,re
import json
import multiprocessing

try:
    from Bio import SeqIO
except ImportError as err:
    pass


from loci_objects.superlocus import superlocus
from loci_objects.transcript import transcript
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.bed12 import BED12
from loci_objects.abstractlocus import abstractlocus

def locus_printer( slocus, args, cds_dict=None, lock=None ):
#     if slocus is None:
#         return
    slocus.load_cds(cds_dict, trust_strand = args.strand_specific )
    stranded_loci = sorted(list(slocus.split_strands()))
    assert lock is not None
    
    for stranded_locus in stranded_loci:
        stranded_locus.define_subloci()
        lock.acquire()
        if not os.path.exists(args.sub_metrics ):
            sub_metrics = open(args.sub_metrics,"w")
        else:
            sub_metrics = open(args.sub_metrics,"a")
        stranded_locus.print_subloci_metrics(sub_metrics)
        sub_out=open(args.sub_out, "a")
        print(stranded_locus, file=sub_out)
        sub_metrics.close()
        lock.release()
        stranded_locus.define_monosubloci()
        lock.acquire()
        mono_out=open(args.mono_out, "a")
        print(stranded_locus, file=mono_out)
        mono_out.close()
        lock.release()
        stranded_locus.calculate_mono_metrics()
        lock.acquire()
        if not os.path.exists(args.locus_metrics ):
            locus_metrics = open(args.locus_metrics,"w")
        else:
            locus_metrics = open(args.locus_metrics,"a")
        stranded_locus.print_monoholder_metrics(locus_metrics)
        locus_metrics.close()
        lock.release()
        stranded_locus.define_loci()
        lock.acquire()
        locus_out=open(args.locus_out,"a")
        print(stranded_locus, file=locus_out)
        locus_out.close()
        lock.release()
    return

def main():
    
    def to_json(string):
        with open(string) as json_file:
            json_dict = json.load(json_file)
        return json_dict
    
    def to_index(string):
        try:
            from Bio import SeqIO
        except ImportError as err:
            print("Error importing the Bio module, no indexing performed:\n{0}",format(err) )
            return None
        index_db = "{0}.db".format(string)
        SeqIO.index_db(index_db, format="fasta", filenames=[string])
        return index_db
    
    parser=argparse.ArgumentParser("Quick test utility for the superlocus classes.")
    parser.add_argument("-p", "--procs", type=int, default=1, help="Number of processors to use.")
    parser.add_argument("--json_conf", type=argparse.FileType("r"), required=True, help="JSON configuration file for scoring transcripts.")
    parser.add_argument("--strand_specific", action="store_true", default=False)
    parser.add_argument("--sub_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--mono_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--locus_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--cds", type=argparse.FileType("r"), default=None)
    parser.add_argument("--transcript_fasta", type=to_index, default=None)
    parser.add_argument("gff", type=argparse.FileType("r"))
    
    args=parser.parse_args()

    args.json_conf = to_json(args.json_conf.name)

    currentLocus=None
    currentTranscript=None
    currentChrom=None

    args.sub_out.close()
    args.sub_out=args.sub_out.name
    
    args.mono_out.close()
    args.mono_out=args.mono_out.name

    args.locus_out.close()
    args.locus_out=args.locus_out.name

    args.sub_metrics=re.sub("$",".metrics",  re.sub(".gff3$", "", args.sub_out  ))
    args.locus_metrics = re.sub("$",".metrics",  re.sub(".gff3$", "", args.locus_out  ))
    if os.path.exists(args.sub_metrics):
        os.remove(args.sub_metrics) 
    if os.path.exists(args.locus_metrics):
        os.remove(args.locus_metrics)
    
    cds_dict=None
    
    if args.cds is not None:
        print("Starting to extract CDS data", file=sys.stderr)
        cds_dict = dict()
        if args.transcript_fasta is not None:
            index=SeqIO.index_db(args.transcript_fasta)
        else:
            index=None
        
        for line in BED12(args.cds, fasta_index=index):
            if line.header is True: continue
            if line.chrom not in cds_dict:
                cds_dict[line.chrom]=[]
            to_append = True
            indices_to_remove = []
            for index in range(len(cds_dict[line.chrom])):
                entry=cds_dict[line.chrom][index]
                overl=abstractlocus.overlap( (entry.cdsStart,entry.cdsEnd), (line.cdsStart,line.cdsEnd) )
                if overl==entry.cds_len:
                    indices_to_remove.append(index)
                elif overl==line.cds_len:
                    to_append=False
                    break
            if to_append is True:
                for index in indices_to_remove:
                    del cds_dict[line.chrom][index]
                cds_dict[line.chrom].append(line)
        print("Finished extracting CDS data", file=sys.stderr)

    args.cds.close()
    args.cds=args.cds.name
    args.gff.close()
    args.gff=args.gff.name

    if args.gff[-3:]=="gtf":
        rower=GTF(args.gff)
    else: rower=GFF3(args.gff)

    manager=multiprocessing.Manager()
    lock=manager.RLock()
    pool=multiprocessing.Pool(processes=args.procs)
    first = True    
    for row in rower:
        
        if row.header is True: continue
        if row.chrom!=currentChrom:
            if currentChrom is not None:
                if currentTranscript is None:
                    pass
                elif currentLocus is not None and superlocus.in_locus(currentLocus, currentTranscript):
                    currentLocus.add_transcript_to_locus(currentTranscript)
                else:
                    if currentLocus is not None:
                        if first is True:
                            locus_printer(currentLocus, args, cds_dict=cds_dict)
                            first=False
                        else:
                            pool.apply_async(locus_printer, args=(currentLocus, args), kwds={"cds_dict": cds_dict,
                                                                                             "lock": lock})
                    currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
                    
            currentChrom=row.chrom
            currentTranscript=None
            currentLocus=None
            
        if row.feature=="transcript" or "RNA" in row.feature.upper():
            if currentLocus is not None:
                if currentTranscript is None:
                    pass
                elif superlocus.in_locus(currentLocus, currentTranscript):
                    currentLocus.add_transcript_to_locus(currentTranscript)
                else:
                    if first is True:
                        locus_printer(currentLocus, args, cds_dict=cds_dict, lock=lock)
                        first=False
                    else:
                        pool.apply_async(locus_printer, args=(currentLocus, args), kwds={"cds_dict": cds_dict,
                                                                                         "lock": lock})
                        
                    currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
            elif currentLocus is None:
                if currentTranscript is not None:
                    currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
            currentTranscript=transcript(row)
        elif row.feature in ("exon", "CDS") or "UTR" in row.feature.upper():
            currentTranscript.addExon(row)
        else:
            continue
   
    if currentLocus is not None:
        if currentTranscript is None:
            pass
        elif superlocus.in_locus(currentLocus, currentTranscript):
            currentLocus.add_transcript_to_locus(currentTranscript)
        else:
            pool.apply_async(locus_printer, args=(currentLocus, args),
                             kwds={"cds_dict": cds_dict, "lock": lock})
            currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
        pool.apply_async(locus_printer, args=(currentLocus, args),
                             kwds={"cds_dict": cds_dict, "lock": lock})
        #locus_printer(currentLocus, args, cds_dict=cds_dict, lock=lock)

    pool.close()
    pool.join()
       
if __name__=="__main__":
    main()