import Mikado
import sys
import os
import argparse
from copy import deepcopy


def main():

    parser = argparse.ArgumentParser("Script to tease apart exons in the MAKER GFF files.")
    parser.add_argument("maker", type=Mikado.parsers.to_gff, help="MAKER input GFF.")
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("w"))
    args = parser.parse_args()

    genes = dict()
    tid2gid = dict()
    
    for line in args.maker:
        if line.feature == "gene":
            genes[line.id] = dict()
            genes[line.id]["transcripts"] = dict()
            genes[line.id]["gene_line"] = line
        elif line.feature == "mRNA":
            tid2gid[line.id] = line.parent[0]
            assert tid2gid[line.id] in genes
            genes[tid2gid[line.id]]["transcripts"][line.id] = dict()
            genes[tid2gid[line.id]]["transcripts"][line.id]["exons"] = []
            genes[tid2gid[line.id]]["transcripts"][line.id]["mrna_line"] = line
        elif line.is_exon:
            parents = line.parent
            my_genes = set()
            for parent in parents:
                assert parent in tid2gid, parent
                my_genes.add(tid2gid[parent])
            assert len(my_genes) == 1
            my_gene = my_genes.pop()
            for parent in parents:
                new_line = deepcopy(line)
                new_line.id = None
                new_line.parent = parent
                genes[my_gene]["transcripts"][parent]["exons"].append(new_line)
                continue
        else:
            pass
        continue

    print("###gff-version\t3", file=args.out)
    for gid, gene in genes.items():
        assert len(gene["transcripts"]) > 0, gid
        print(gene["gene_line"])
        for tid, transcript in gene["transcripts"].items():
            print(transcript["mrna_line"])
            print(*transcript["exons"], sep="\n")
            continue
        continue
    return

main()
                
            
