all: sample.loci.gff3 trinity.loci.gff3 trinity.loci.from_gff.gff3

clean:
	rm -f ./*loci* ./trinity.diff

sample.loci.gff3: sample.gff3
	/usr/bin/time -v python3 ../run_pipeline.py -p 1 --json_conf scoring.json --cds sample.cds.bed --sub_out sample.subloci.gff3 --mono_out sample.monoloci.gff3 --locus_out sample.loci.gff3 sample.gff3

trinity.loci.gff3: trinity.gtf
	/usr/bin/time -v python3 ../run_pipeline.py -p 1 --json_conf scoring.json --transcript_fasta trinity.transcripts.fasta --cds trinity.cds.bed --sub_out trinity.subloci.gff3 --mono_out trinity.monoloci.gff3 --locus_out trinity.loci.gff3 trinity.gtf

trinity.loci.from_gff.gff3: trinity.gff3
	/usr/bin/time -v python3 ../run_pipeline.py -p 1 --json_conf scoring.json --transcript_fasta trinity.transcripts.fasta --cds trinity.cds.bed --sub_out trinity.subloci.from_gff.gff3 --mono_out trinity.monoloci.from_gff.gff3 --locus_out trinity.loci.from_gff.gff3 trinity.gff3 
	
