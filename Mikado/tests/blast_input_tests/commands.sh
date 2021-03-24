mkdir -p ./{uniprot,sanitised}
cp uniprot.fasta uniprot/
sanitize_blast_db.py -o sanitised/uniprot.fasta uniprot.fasta
export blast_keys=$(python -c "from Mikado.utilities import blast_keys; print(*blast_keys)")

for folder in uniprot sanitised; do
    mkdir -p ${folder}/{blast,blast_parse_seqids,diamond}/{tsv,xml};
    # mkdir -p ${folder}/{blast,blast_parse_seqids}/asn;
    mkdir -p ${folder}/diamond/daa/;
    cd ${folder};
    # DIAMOND
    diamond makedb --in uniprot.fasta --db diamond/sequences.dnd;
    diamond blastx --query ../queries.fasta --db diamond/sequences.dnd.dmnd \
	    --outfmt 6 ${blast_keys} --compress 1 --out diamond/tsv/blast.tsv.gz;
    diamond blastx --query ../queries.fasta --db diamond/sequences.dnd.dmnd \
	    --outfmt 5 --compress 1 --out diamond/xml/blast.xml.gz
    diamond blastx --query ../queries.fasta --db diamond/sequences.dnd.dmnd \
	    --outfmt 100 --out diamond/daa/blast.daa
    # NCBI BLASTX+
    makeblastdb -dbtype prot -in uniprot.fasta -out blast/sequences
    blastx -query ../queries.fasta -db blast/sequences -outfmt "6 ${blast_keys}" | gzip -c > blast/tsv/blast.tsv.gz
    blastx -query ../queries.fasta -db blast/sequences -outfmt 5 | gzip -c > blast/xml/blast.xml.gz
    # ASN is a *pain* to analyse as it relies on the path of the database. Excluding it for now.
    # blastx -query ../queries.fasta -db blast/sequences -outfmt 11 | gzip -c > blast/asn/blast.asn.gz

    # NCBI BLASTX+ with -parse_seqids
    makeblastdb -parse_seqids -dbtype prot -in uniprot.fasta -out blast_parse_seqids/sequences
    blastx -query ../queries.fasta -db blast_parse_seqids/sequences -outfmt "6 ${blast_keys}" | gzip -c > blast_parse_seqids/tsv/blast.tsv.gz
    blastx -query ../queries.fasta -db blast_parse_seqids/sequences -outfmt 5 | gzip -c > blast_parse_seqids/xml/blast.xml.gz
    # ASN is a *pain* to analyse as it relies on the path of the database. Excluding it for now.
    # blastx -query ../queries.fasta -db blast_parse_seqids/sequences -outfmt 11 | gzip -c > blast_parse_seqids/asn/blast.asn.gz
    # Go back
    cd ../;
done
