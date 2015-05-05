all: sample.loci.gff3 trinity.loci.gff3 cufflinks.loci.gff3 class.loci.gff3 stringtie.loci.gff3 all.loci.gff3 

clean:
	rm -f ./*loci* 

sample.loci.gff3: sample.gff3
	python3 ../run_pipeline.py --source mock --json_conf scoring.json --cds sample.bed --sub_out sample.subloci.gff3 --mono_out sample.monoloci.gff3 --locus_out sample.loci.gff3 sample.gff3

trinity.loci.gff3: trinity.gtf
	python3 ../run_pipeline.py --source trinity --json_conf scoring.json --transcript_fasta trinity.fasta --cds trinity.bed --sub_out trinity.subloci.gff3 --mono_out trinity.monoloci.gff3 --locus_out trinity.loci.gff3 trinity.gtf

cufflinks.loci.gff3: cufflinks.gtf
	python3 ../run_pipeline.py --source cufflinks --json_conf scoring.json --transcript_fasta cufflinks.fasta --cds cufflinks.bed --sub_out cufflinks.subloci.gff3 --mono_out cufflinks.monoloci.gff3 --locus_out cufflinks.loci.gff3 cufflinks.gtf
	
class.loci.gff3: class.gtf
	python3 ../run_pipeline.py --source class --json_conf scoring.json --transcript_fasta class.fasta --cds class.bed --sub_out class.subloci.gff3 --mono_out class.monoloci.gff3 --locus_out class.loci.gff3 class.gtf
	
stringtie.loci.gff3: stringtie.gtf
	python3 ../run_pipeline.py --source stringtie --json_conf scoring.json --transcript_fasta stringtie.fasta --cds stringtie.bed --sub_out stringtie.subloci.gff3 --mono_out stringtie.monoloci.gff3 --locus_out stringtie.loci.gff3 stringtie.gtf

all.loci.gff3: all.gtf
	python3 ../run_pipeline.py --purge --source all --json_conf scoring.json --transcript_fasta all.fasta --cds all.bed --sub_out all.subloci.gff3 --mono_out all.monoloci.gff3 --locus_out all.loci.gff3 all.gtf
