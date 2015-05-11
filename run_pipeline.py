#!/usr/bin/env python3

import sys
import pickle
import argparse,re
import multiprocessing
import csv

try:
    from Bio import SeqIO
except ImportError as err:
    pass

from loci_objects.json_utils import check_json
from loci_objects.json_utils import to_json
from loci_objects.superlocus import superlocus
from loci_objects.sublocus import sublocus
from loci_objects.transcript import transcript
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.bed12 import BED12
from loci_objects.abstractlocus import abstractlocus

def locus_printer( slocus, args, cds_dict=None, lock=None ):

    slocus.load_cds(cds_dict, trust_strand = args.strand_specific )
    stranded_loci = sorted(list(slocus.split_strands()))
    
    for stranded_locus in stranded_loci:
        stranded_locus.define_loci()
    
    for stranded_locus in stranded_loci:
        if args.remove_overlapping_fragments is True and len(stranded_loci)>1:
            for final_locus in stranded_locus.loci:
                for other_superlocus in filter(lambda x: x!=stranded_locus, stranded_loci):
                    for other_final_locus in other_superlocus.loci:
                        if other_final_locus.other_is_fragment( final_locus, percentage=0.5 ) is True:
                            try:
                                stranded_locus.loci.remove(final_locus)
                            except ValueError as err:
                                if "not in list" in err: pass
                                else:
                                    raise ValueError(err)
                            finally:
                                break
        sub_lines = stranded_locus.__str__(level="subloci", print_cds=not args.no_cds )
        sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()]
        mono_lines = stranded_locus.__str__(level="monosubloci", print_cds=not args.no_cds)
        locus_metrics_rows=[x for x in stranded_locus.print_monoholder_metrics()]
        locus_lines = stranded_locus.__str__(print_cds=not args.no_cds)

        #Print out
        if lock is not None:
            lock.acquire()
        with open(args.sub_metrics,'a') as out_file:
            csv_out=csv.DictWriter(out_file, superlocus.available_metrics, delimiter="\t")
            for row in sub_metrics_rows: csv_out.writerow(row)
        with open(args.locus_metrics,'a') as out_file:
            csv_out=csv.DictWriter(out_file, superlocus.available_metrics, delimiter="\t")
            for row in locus_metrics_rows: csv_out.writerow(row)
 
        with open(args.sub_out,'a') as sub_out:
            print(sub_lines, file=sub_out)
        if mono_lines!='':
            with open(args.mono_out,'a') as mono_out:
                print(mono_lines, file=mono_out)
        if locus_lines!='':
            with open(args.locus_out,'a') as locus_out:
                print(locus_lines, file=locus_out)
        if lock is not None:
            lock.release()
    return

def main():
    
   
    def to_index(string):
        if "SeqIO" not in globals():
            print("Error importing the Bio module, no indexing performed:\n{0}",format(err) )
            return None

        return SeqIO.index(string, format="fasta") 
    
    parser=argparse.ArgumentParser("Quick test utility for the superlocus classes.")
    parser.add_argument("-p", "--procs", type=int, default=1, help="Number of processors to use.")
    parser.add_argument('-x',  "--remove_overlapping_fragments", action="store_true",
                        default=False, help="""Flag. If set, the program will remove monoexonic loci
                        overlapping non-monoexonic loci on the opposite strand.""")
    parser.add_argument("--json_conf", type=argparse.FileType("r"), required=True, help="JSON configuration file for scoring transcripts.")
    parser.add_argument("--strand_specific", action="store_true", default=False)
    parser.add_argument("--sub_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--mono_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--locus_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--no_cds", action="store_true", default=False,
                        help="Flag. If set, not CDS information will be printed out in the GFF output files.")
    parser.add_argument('--source', type=str, default=None,
                        help='Source field to use for the output files.')
    parser.add_argument('--purge', action='store_true', default=False,
                        help='Flag. If set, the pipeline will suppress any loci whose transcripts do not pass the requirements set in the JSON file.'
                        )
    parser.add_argument("--cds", type=argparse.FileType("r"), default=None)
    parser.add_argument("--transcript_fasta", type=to_index, default=None)
    parser.add_argument("gff", type=argparse.FileType("r"))
    
    args=parser.parse_args()

    args.json_conf = to_json(args.json_conf.name)
    check_json(args.json_conf)
    if ("requirements" in args.json_conf and "compiled" in args.json_conf["requirements"]) or ("compiled" in args.json_conf):
        raise KeyError("Why is compiled here again?")
    
            
    currentLocus=None
    currentTranscript=None
    currentChrom=None

    print('##gff-version 3', file=args.sub_out)
    print('##gff-version 3', file=args.mono_out)
    print('##gff-version 3', file=args.locus_out)

    args.sub_out.close()
    args.sub_out=args.sub_out.name
    
    args.mono_out.close()
    args.mono_out=args.mono_out.name

    args.locus_out.close()
    args.locus_out=args.locus_out.name

    args.sub_metrics=re.sub("$",".metrics.tsv",  re.sub(".gff3$", "", args.sub_out  ))
    args.locus_metrics = re.sub("$",".metrics.tsv",  re.sub(".gff3$", "", args.locus_out  ))
    with open(args.sub_metrics, "w") as out_file:
        csv_out=csv.DictWriter( out_file, superlocus.available_metrics, delimiter="\t" )
        csv_out.writeheader()
    with open(args.locus_metrics, "w") as out_file:
        csv_out=csv.DictWriter( out_file, superlocus.available_metrics, delimiter="\t" )
        csv_out.writeheader()
    
    cds_dict=None
    
    if args.cds is not None:
        print("Starting to extract CDS data", file=sys.stderr)
        cds_dict = dict()
        
        for line in BED12(args.cds, fasta_index=args.transcript_fasta):
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

    if args.cds is not None:
        args.cds.close()    
        args.cds=args.cds.name
    args.gff.close()
    args.gff=args.gff.name
    args.transcript_fasta = None

    if args.gff[-3:]=="gtf":
        rower=GTF(args.gff)
    else: rower=GFF3(args.gff)

    manager=multiprocessing.Manager() # @UndefinedVariable
    lock=manager.RLock()
    pool=multiprocessing.Pool(processes=args.procs) # @UndefinedVariable
    first = True    
    jobs=dict()
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
#                             first=False
                        else:
                            if ("requirements" in args.json_conf and "compiled" in args.json_conf["requirements"]) or ("compiled" in args.json_conf):
                                raise KeyError("Why is compiled here again?")

                            jobs[currentLocus]=pool.apply_async(locus_printer, args=(currentLocus, args), kwds={"cds_dict": cds_dict,
                                                                                             "lock": lock})
                    currentLocus=superlocus(currentTranscript, stranded=False,
                                            json_dict = args.json_conf,
                                            purge=args.purge)
                    
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
#                         first=False
                    else:
                        if ("requirements" in args.json_conf and "compiled" in args.json_conf["requirements"]) or ("compiled" in args.json_conf):
                            raise KeyError("Why is compiled here again?")
                        
                        jobs[currentLocus]=pool.apply_async(locus_printer, args=(currentLocus, args), kwds={"cds_dict": cds_dict,
                                                                                         "lock": lock})
                        
                    currentLocus=superlocus(currentTranscript,
                                            stranded=False, json_dict = args.json_conf,
                                            purge=args.purge)
            elif currentLocus is None:
                if currentTranscript is not None:
                    currentLocus=superlocus(currentTranscript,
                                            stranded=False, json_dict = args.json_conf,
                                            purge=args.purge)
            currentTranscript=transcript(row, source=args.source)
        elif row.feature in ("exon", "CDS") or "UTR" in row.feature.upper():
            currentTranscript.addExon(row)
        else:
            continue
   
    if currentLocus is not None:
        if currentTranscript is not None:
            if superlocus.in_locus(currentLocus, currentTranscript):
                currentLocus.add_transcript_to_locus(currentTranscript)
            else:
                if ("requirements" in args.json_conf and "compiled" in args.json_conf["requirements"]) or ("compiled" in args.json_conf):
                    raise KeyError("Why is compiled here again?")
            
                jobs[currentLocus]=pool.apply_async(locus_printer, args=(currentLocus, args),
                                     kwds={"cds_dict": cds_dict, "lock": lock})

                currentLocus=superlocus(currentTranscript,
                                        stranded=False, json_dict = args.json_conf,
                                        purge=args.purge)
                
        jobs[currentLocus]=pool.apply_async(locus_printer, args=(currentLocus, args), kwds={"cds_dict": cds_dict,
                                                                                         "lock": lock})
        
#         process = multiprocessing.Process(target=locus_printer, # @UndefinedVariable
#                                           args=(currentLocus, args),
#                                           kwargs={"cds_dict": cds_dict, "lock": lock})
#         process.start()
#         process.join()
#         jobs.append(pool.apply_async(locus_printer, args=(currentLocus, args),
#                              kwds={"cds_dict": cds_dict, "lock": lock}))
        
        for job in jobs:
            try:
                jobs[job].get()
            except:
                print("Something has gone wrong with the final process.", file=sys.stderr)
                print(job.chrom, job.start, job.end, ",".join(job.transcripts))
                locus_printer(job, args, cds_dict=cds_dict, lock=lock)
                pickle.dumps(job)
#                 locus_printer(currentLocus, args, cds_dict=cds_dict, lock=lock)

    pool.close()
    pool.join()
       
if __name__=="__main__":
    #Check that the metrics are alright
    metrics_not_found = []
    for metric in filter(lambda m: m not in ("score","parent","tid"), sublocus.available_metrics ):
        if not hasattr(transcript, metric):
            metrics_not_found.append(metric)

    if len(metrics_not_found)>0:
        raise RuntimeError("The following metrics are requested but are not defined in the transcript class:\n\t{0}".format(
                                                                                                                        "\n\t".join(metrics_not_found)
                                                                                                                        ))

    main()