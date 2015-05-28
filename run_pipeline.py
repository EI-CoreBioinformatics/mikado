#!/usr/bin/env python3

#import sys
#from logging import Logger
import argparse,re
import multiprocessing
import csv

try:
    from Bio import SeqIO
except ImportError as err:
    pass

from loci_objects.json_utils import to_json
from loci_objects.superlocus import superlocus
from loci_objects.sublocus import sublocus
from loci_objects.transcript import transcript
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.bed12 import BED12
from loci_objects.abstractlocus import abstractlocus

def analyse_locus( slocus, args, queue, cds_dict=None, lock=None ):

    '''This function takes as input a "superlocus" instance and the pipeline configuration.
    It also accepts as optional keywords a dictionary with the CDS information (derived from a BED12)
    and a "lock" used for avoiding writing collisions during multithreading.
    The function splits the superlocus into its strand components and calls the relevant methods
    to define the loci. It also prints out the results to the requested output files.
    '''

    #Load the CDS information
    slocus.load_cds(cds_dict, trust_strand = args.strand_specific,
                    minimal_secondary_orf_length=args.minimal_secondary_orf_length,
                    split_chimeras=args.split_chimeras, fasta_index = args.transcript_fasta )
    #Split the superlocus in the stranded components
    stranded_loci = sorted(list(slocus.split_strands()))
    #Define the loci
    for stranded_locus in stranded_loci:
        stranded_locus.define_loci()

    #Remove overlapping fragments.
    #This part should be rewritten in order to make it more flexible and powerful.
    for stranded_locus in stranded_loci:
        if args.remove_overlapping_fragments is True and len(stranded_loci)>1:
            for final_locus in stranded_locus.loci:
                for other_superlocus in filter(lambda x: x!=stranded_locus, stranded_loci):
                    for other_final_locus in other_superlocus.loci:
                        if other_final_locus.other_is_fragment( final_locus ) is True:
                            stranded_locus.loci.remove(final_locus)
#                             except ValueError as err:
#                                 if "not in list" in err: pass
#                                 else:
#                                     raise ValueError(err)
#                             finally:
#                                 break
        queue.put(stranded_locus)
    return

def printer(args,queue):

    '''Listener process that will print out the loci recovered by the analyse_locus function.'''
    
    sub_metrics=csv.DictWriter(open(args.sub_metrics,'a'), superlocus.available_metrics, delimiter="\t")
    sub_scores=csv.DictWriter(open(args.sub_scores,'a'), args.score_keys, delimiter="\t")
    locus_metrics=csv.DictWriter(open(args.locus_metrics,'a'), superlocus.available_metrics, delimiter="\t")
    locus_scores=csv.DictWriter(open(args.locus_scores,'a'), args.score_keys, delimiter="\t")
    sub_out=open(args.sub_out,'a')
    mono_out=open(args.mono_out,'a')
    locus_out=open(args.locus_out,'a') 
    
    while True:
        stranded_locus=queue.get()
        if stranded_locus is None:
            break
        sub_lines = stranded_locus.__str__(level="subloci", print_cds=not args.no_cds )
        sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()]
        sub_scores_rows = [x for x in stranded_locus.print_subloci_scores()]
        mono_lines = stranded_locus.__str__(level="monosubloci", print_cds=not args.no_cds)
        locus_metrics_rows=[x for x in stranded_locus.print_monoholder_metrics()]
        locus_scores_rows = [x for x in stranded_locus.print_monoholder_scores()]
        locus_lines = stranded_locus.__str__(print_cds=not args.no_cds)
        for row in sub_metrics_rows: sub_metrics.writerow(row)
        for row in sub_scores_rows: sub_scores.writerow(row)
        for row in locus_metrics_rows: locus_metrics.writerow(row)
        for row in locus_scores_rows: locus_scores.writerow(row)
        print(sub_lines, file=sub_out)
        if mono_lines!='':
            print(mono_lines, file=mono_out)
        if locus_lines!='':
            print(locus_lines, file=locus_out)    
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
    parser.add_argument("-sc", "--split-chimeras", dest="split_chimeras", default=False,
                        action="store_true", help="""Flag. If set, transcripts with multiple ORFs will be split into separate transcripts,
                        trying to retain as much UTR as possible.""" )
    parser.add_argument("--json_conf", type=argparse.FileType("r"), required=True, help="JSON configuration file for scoring transcripts.")
    parser.add_argument("--minimal_secondary_orf_length", type=int, default=200,
                        help="Any secondary ORF shorter than this value will be ignored. Useful to avoid legitimate cases of ORFs in the UTR. Default: %(default)s.")
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
    if args.json_conf["chimera_split"]["blast_check"] is True and args.transcript_fasta is None:
        print("No FASTA, blast check disabled")
        args.json_conf["chimera_split"]["blast_check"]=False
    
#     check_json(args.json_conf)
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
    args.sub_scores=re.sub("$",".scores.tsv",  re.sub(".gff3$", "", args.sub_out  ))
    args.locus_metrics = re.sub("$",".metrics.tsv",  re.sub(".gff3$", "", args.locus_out  ))
    args.locus_scores = re.sub("$",".scores.tsv",  re.sub(".gff3$", "", args.locus_out  ))
    
    args.score_keys = ["tid","parent","score"] + sorted(list(args.json_conf["parameters"].keys()))
    
    #Prepare output files
    with open(args.sub_metrics, "w") as out_file:
        csv_out=csv.DictWriter( out_file, superlocus.available_metrics, delimiter="\t" )
        csv_out.writeheader()
    with open(args.sub_scores,'w') as out_file:
        csv_out=csv.DictWriter( out_file, args.score_keys, delimiter="\t" )
        csv_out.writeheader()
    with open(args.locus_metrics, "w") as out_file:
        csv_out=csv.DictWriter( out_file, superlocus.available_metrics, delimiter="\t" )
        csv_out.writeheader()
    with open(args.locus_scores,'w') as out_file:
        csv_out=csv.DictWriter( out_file, args.score_keys, delimiter="\t" )
        csv_out.writeheader()
        
    cds_dict=None
    #Load the CDS information from the BED12, if one is available
    if args.cds is not None:
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

    if args.cds is not None:
        args.cds.close()    
        args.cds=args.cds.name
    args.gff.close()
    args.gff=args.gff.name
#     args.transcript_fasta = None
    #Determine the type of input file (GTF/GFF3)
    if args.gff[-3:]=="gtf":
        rower=GTF(args.gff)
    else: rower=GFF3(args.gff)

    ctx=multiprocessing.get_context("fork") #@UndefinedVariable
    
    manager=ctx.Manager() # @UndefinedVariable
    queue=manager.Queue()
    lock=manager.RLock()
    pool=ctx.Pool(processes=args.procs) # @UndefinedVariable
    printer_process=multiprocessing.Process(target=printer, args=(args,queue)) # @UndefinedVariable
    printer_process.start()
    first = True    
    for row in rower:
        if row.header is True:
            continue
        if row.chrom!=currentChrom:
            if currentChrom is not None:
                if currentTranscript is None:
                    pass
                elif currentLocus is not None and superlocus.in_locus(currentLocus, currentTranscript):
                    currentLocus.add_transcript_to_locus(currentTranscript)
                else:
                    if currentLocus is not None:
                        if first is True:
                            analyse_locus(currentLocus, args, queue, cds_dict=cds_dict)
                        else:
                            pool.apply_async(analyse_locus, args=(currentLocus, args, queue), kwds={"cds_dict": cds_dict,
                                                                                             "lock": lock})
                    currentLocus=superlocus(currentTranscript, stranded=False,
                                                json_dict = args.json_conf,
                                                purge=args.purge)
                    
            currentChrom=row.chrom
            currentTranscript=None
            currentLocus=None
            
        if row.is_transcript is True:
            if currentLocus is not None:
                if currentTranscript is None:
                    pass
                elif superlocus.in_locus(currentLocus, currentTranscript):
                    currentLocus.add_transcript_to_locus(currentTranscript)
                else:
                    if first is True:
                        analyse_locus(currentLocus, args, queue, cds_dict=cds_dict, lock=lock)
                    else:
                        pool.apply_async(analyse_locus,
                                         args=(currentLocus, args, queue),
                                         kwds={"cds_dict": cds_dict,"lock": lock})
                    currentLocus=superlocus(currentTranscript, stranded=False,
                                            json_dict = args.json_conf,
                                            purge=args.purge)

            elif currentLocus is None:
                if currentTranscript is not None:
                    currentLocus=superlocus(currentTranscript,
                                            stranded=False, json_dict = args.json_conf,
                                            purge=args.purge)
            currentTranscript=transcript(row, source=args.source)
        elif row.feature in ("exon", "CDS") or "UTR" in row.feature.upper() or "codon" in row.feature:
            currentTranscript.addExon(row)
        else:
            continue

    if currentLocus is not None:
        if currentTranscript is not None:
            if superlocus.in_locus(currentLocus, currentTranscript):
                currentLocus.add_transcript_to_locus(currentTranscript)
            else:
                pool.apply_async(analyse_locus, args=(currentLocus, args, queue),
                                kwds={"cds_dict": cds_dict, "lock": lock})
                currentLocus=superlocus(currentTranscript,
                                        stranded=False, json_dict = args.json_conf,
                                        purge=args.purge)
        pool.apply_async(analyse_locus, 
                         args=(currentLocus, args, queue),
                         kwds={"cds_dict": cds_dict, "lock": lock})
        
    pool.close()
    pool.join()
    queue.put(None)
    printer_process.join()
    printer_process.terminate()
       
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
