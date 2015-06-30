#!/usr/bin/env python3


import sys,argparse,os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from shanghai_lib.parsers import GFF

def main():
    parser=argparse.ArgumentParser('Script to parse and retrieve given features from a GFF file.')
    parser.add_argument('-v', action='store_true', dest='reverse', help="Exclude from the gff all the records in the id file.")
    parser.add_argument("--genes", action="store_true", help="Flag. If set, the program expects as ids only a list of genes, and will exclude/include all the transcripts children of the selected genes.")
    parser.add_argument('ids', type=argparse.FileType('r'), help="ID file (format: mrna_id, gene_id - tab separated)")
    parser.add_argument('gff', type=argparse.FileType('r'), help="The GFF file to parse.")
    parser.add_argument('out', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Optional output file")
    args=parser.parse_args()

    gene_ids, mrna_ids = set(), set()
    for line in args.ids:
        if args.genes is False:
            mrna_id, gene_id = line.rstrip().split()[:2]
            mrna_ids.add(mrna_id)
            gene_ids.add(gene_id)
        else:
            gene_id = line.rstrip()
            gene_ids.add(gene_id)


    curr_gene = None
    curr_transcripts = dict()

    for record in GFF.GFF3(args.gff):
        if record.is_transcript is True: #Potential gene line
            if args.reverse is False and (record.id in mrna_ids or ( args.genes is True and any(lambda p: p in gene_ids, record.parent))): 
                curr_transcripts[record.id]=[record]
            elif args.reverse is True and ((args.genes is False and record.id not in mrna_ids) or (args.genes is True and not any(lambda p: p in gene_ids, record.parent))):
                curr_transcripts[record.id]=[record]
        elif record.is_exon is True:
            for parent in record.parent:
                if parent in curr_transcripts: curr_transcripts[parent].append(record)
        else:
            if curr_gene is not None and len(curr_transcripts)>0:
                print(curr_gene, file=args.out)
                for tid in curr_transcripts:
                    print(curr_transcripts[tid][0], file=args.out)
                    for rec in curr_transcripts[tid][1:]:
                        rec.parent = [tid]
                        print(rec, file=args.out)
                
            curr_gene = None
            curr_transcripts = dict()
            
            if record.id is None:
                continue
            elif args.reverse is True and record.id not in gene_ids:
                curr_gene = record
            elif args.reverse is False and record.id in gene_ids:
                curr_gene = record
        
    if curr_gene is not None and len(curr_transcripts)>0:
        print(curr_gene, file=args.out)
        for tid in curr_transcripts:
            print(curr_transcripts[tid][0], file=args.out)
            for rec in curr_transcripts[tid][1:]:
                rec.parent = [tid]
                print(rec, file=args.out)



                

if __name__=="__main__": main()
