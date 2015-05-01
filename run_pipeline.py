#!/usr/bin/env python3

import sys
import argparse,re
import json
import multiprocessing
import csv

try:
    from Bio import SeqIO
except ImportError as err:
    pass

from loci_objects.superlocus import superlocus
from loci_objects.sublocus import sublocus
from loci_objects.transcript import transcript
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.bed12 import BED12
from loci_objects.abstractlocus import abstractlocus

def check_json(json_conf):
    '''Quick function to check that the JSON dictionary is well formed.'''
    
    parameters_not_found=[]
    parameters_found=set()
    double_parameters=[]
    for parameter in json_conf["parameters"]:
        if parameter not in superlocus.available_metrics:
            parameters_not_found.append(parameter)
        if parameter in parameters_found:
            double_parameters.add(parameter)
        parameters_found.add(parameter)
    
    import importlib    
    mods_not_found = [] 
    for mod in json_conf["modules"]:
        try:
            importlib.import_module(mod)
        except ImportError:
            mods_not_found.append(mod)

    if len(parameters_not_found)>0 or len(double_parameters)>0 or len(mods_not_found)>0:
        err_message=''
        if len(parameters_not_found)>0:
            err_message="The following parameters, present in the JSON file, are not available!\n\t{0}\n".format("\n\t".join(parameters_not_found))
        if len(double_parameters)>0:
            err_message+="The following parameters have been specified more than once, please correct:\n\t{0}".format("\n\t".join(list(double_parameters)))
        if len(mods_not_found)>0:
            err_message+="The following requested modules are unavailable:\n\t{0}\n".format("\n\t".join(mods_not_found))
        print(err_message, file=sys.stderr)
        sys.exit(1)
           

def locus_printer( slocus, args, cds_dict=None, lock=None ):
#     if slocus is None:
#         return
    slocus.load_cds(cds_dict, trust_strand = args.strand_specific )
    stranded_loci = sorted(list(slocus.split_strands()))
#    assert lock is not None
    
    for stranded_locus in stranded_loci:
        stranded_locus.define_subloci()
        sub_lines = str(stranded_locus)
        sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()] 
        stranded_locus.define_monosubloci()
        mono_lines = str(stranded_locus)
        stranded_locus.calculate_mono_metrics()
        locus_metrics_rows=[x for x in stranded_locus.print_monoholder_metrics()]
        stranded_locus.define_loci()
        locus_lines = str(stranded_locus)

        #Print out
#         lock.acquire()
        with open(args.sub_metrics,'a') as out_file:
            csv_out=csv.DictWriter(out_file, superlocus.available_metrics, delimiter="\t")
            for row in sub_metrics_rows: csv_out.writerow(row)
        with open(args.locus_metrics,'a') as out_file:
            csv_out=csv.DictWriter(out_file, superlocus.available_metrics, delimiter="\t")
            for row in locus_metrics_rows: csv_out.writerow(row)

        with open(args.sub_out,'a') as sub_out:
            print(sub_lines, file=sub_out)
        with open(args.mono_out,'a') as mono_out:
            print(mono_lines, file=mono_out)
        with open(args.locus_out,'a') as locus_out:
            print(locus_lines, file=locus_out)
#         lock.release()
    return

def main():
    
    def to_json(string):
        with open(string) as json_file:
            json_dict = json.load(json_file)
        return json_dict
    
    def to_index(string):
        if "SeqIO" not in globals():
            print("Error importing the Bio module, no indexing performed:\n{0}",format(err) )
            return None

        return SeqIO.index(string, format="fasta") 
    
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
    check_json(args.json_conf)
            

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
#                         first=False
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