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


    first=True
    for record in GFF.GFF3(args.gff):
        #Add header if absent
        if first:
            first=False
            if "##gff-version" in record._line:
                print(record._line, file=args.out)
                continue
            else:
                print("##gff-version 3", file=args.out)
        elif not first and "##gff-version" in record._line:
            continue #Skip header lines after the first one

#         if record._line[0]=="#":
#             
# 
#             print(record._line.rstrip(), file=args.out)
# #            if args.reverse: print(record._line.rstrip())
#             continue
#         if 'Parent' in record.attributes: parent=record.attributes['Parent']
#         else: parent=None
        if 'Target' in record.attributes: target=record.attributes['Target'].split()[0]
        else: target=None
        if 'Name' in record.attributes: name=record.attributes['Name']
        else: name=None
        if record.header is True:
            continue

        if record.feature=="gene":
            if record.attributes['ID'] in gene_ids and not args.reverse:
                print(record, file=args.out)
            elif record.attributes['ID'] not in gene_ids and args.reverse:
                print(record, file=args.out)
        elif record.feature in ("mRNA", "transcript") or "transcript" in record.feature:
            if args.genes is True and record.parent in gene_ids:
                mrna_ids.add(record.id)

            if (args.reverse):
                if record.id in mrna_ids: continue
                else:
                    print(record, file=args.out)
            elif (not args.reverse) and (record.id in mrna_ids) or (target in mrna_ids) or (name in mrna_ids): print(record, file=args.out)
        elif record.is_exon is True:
#             try: parents = set(parent.split(","))
#             except: continue
            if record.parent is None: continue
            parents = set(record.parent)
            parent_intersection=set.intersection(parents, mrna_ids)
            # if record.id is None:
            #     if record.feature not in current_counter: current_counter[record.feature]=0
                
            #     current_counter[record.feature]+=1
            #     if "," in record.parent:
            #         par = re.sub(",", "_", record.parent)
            #     else:
            #         par = record.parent

                
            #     record.id = "{0}.{1}{2}".format(par, record.feature.lower(), current_counter[record.feature])
            
            if args.reverse:
                if record.id not in mrna_ids and target not in mrna_ids and name not in mrna_ids:
                    if parent_intersection==set(): #No intersection
                        print(record, file=args.out)
                    elif parent_intersection!=parents: #Intersection non-empty
                        good_parents=",".join(set.difference(parents, mrna_ids))
                        record.parent=good_parents
                        print(record, file=args.out)
            else:
                if record.id in mrna_ids or parent_intersection!=set(): #or record.target in mrna_ids or record.name in mrna_ids:
                        if parent_intersection==parents:
                            print(record, file=args.out)
                        else:
                            #print(good_parents, file=sys.stderr)
                            record.parent=parent_intersection
                            print(record, file=args.out)
                

if __name__=="__main__": main()
