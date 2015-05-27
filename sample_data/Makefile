all: sample.loci.scores.tsv trinity.loci.scores.tsv cufflinks.loci.scores.tsv class.loci.scores.tsv stringtie.loci.scores.tsv all.loci.scores.tsv

clean:
	rm -f ./*loci* 

sample.loci.scores.tsv: sample.loci.gff3
#	python3 ../util/metrics_to_score.py --json_conf scoring.json sample.subloci.metrics.tsv sample.subloci.scores.tsv
#	python3 ../util/metrics_to_score.py --json_conf scoring.json sample.loci.metrics.tsv sample.loci.scores.tsv 

sample.loci.gff3: sample.gff3
	python3 ../run_pipeline.py -sc --source mock --json_conf scoring.json --cds sample.bed --sub_out sample.subloci.gff3 --mono_out sample.monoloci.gff3 --locus_out sample.loci.gff3 sample.gff3

trinity.loci.scores.tsv: trinity.loci.gff3
#	python3 ../util/metrics_to_score.py --json_conf scoring.json trinity.subloci.metrics.tsv trinity.subloci.scores.tsv
#	python3 ../util/metrics_to_score.py --json_conf scoring.json trinity.loci.metrics.tsv trinity.loci.scores.tsv 

trinity.loci.gff3: trinity.gtf
	python3 ../run_pipeline.py -x -sc --source trinity --json_conf scoring.json --transcript_fasta trinity.fasta --cds trinity.bed --sub_out trinity.subloci.gff3 --mono_out trinity.monoloci.gff3 --locus_out trinity.loci.gff3 trinity.gtf

cufflinks.loci.scores.tsv: cufflinks.loci.gff3
#	python3 ../util/metrics_to_score.py --json_conf scoring.json cufflinks.subloci.metrics.tsv cufflinks.subloci.scores.tsv
#	python3 ../util/metrics_to_score.py --json_conf scoring.json cufflinks.loci.metrics.tsv cufflinks.loci.scores.tsv 

cufflinks.loci.gff3: cufflinks.gtf
	python3 ../run_pipeline.py -x -sc --source cufflinks --json_conf scoring.json --transcript_fasta cufflinks.fasta --cds cufflinks.bed --sub_out cufflinks.subloci.gff3 --mono_out cufflinks.monoloci.gff3 --locus_out cufflinks.loci.gff3 cufflinks.gtf

class.loci.scores.tsv: class.loci.gff3
#	python3 ../util/metrics_to_score.py --json_conf scoring.json class.subloci.metrics.tsv class.subloci.scores.tsv
#	python3 ../util/metrics_to_score.py --json_conf scoring.json class.loci.metrics.tsv class.loci.scores.tsv 
class.loci.gff3: class.gtf
	python3 ../run_pipeline.py -x -sc --source class --json_conf scoring.json --transcript_fasta class.fasta --cds class.bed --sub_out class.subloci.gff3 --mono_out class.monoloci.gff3 --locus_out class.loci.gff3 class.gtf

stringtie.loci.scores.tsv: stringtie.loci.gff3
#	python3 ../util/metrics_to_score.py --json_conf scoring.json stringtie.subloci.metrics.tsv stringtie.subloci.scores.tsv
#	python3 ../util/metrics_to_score.py --json_conf scoring.json stringtie.loci.metrics.tsv stringtie.loci.scores.tsv 
stringtie.loci.gff3: stringtie.gtf
	python3 ../run_pipeline.py -x -sc --source stringtie --json_conf scoring.json --transcript_fasta stringtie.fasta --cds stringtie.bed --sub_out stringtie.subloci.gff3 --mono_out stringtie.monoloci.gff3 --locus_out stringtie.loci.gff3 stringtie.gtf

all.loci.scores.tsv: all.loci.gff3
#	python3 ../util/metrics_to_score.py --json_conf scoring.json all.subloci.metrics.tsv all.subloci.scores.tsv
#	python3 ../util/metrics_to_score.py --json_conf scoring.json all.loci.metrics.tsv all.loci.scores.tsv 

all.loci.gff3: all.gtf
	python3 ../run_pipeline.py -x -sc --source all --json_conf scoring.json --transcript_fasta all.fasta --cds all.bed --sub_out all.subloci.gff3 --mono_out all.monoloci.gff3 --locus_out all.loci.gff3 all.gtf
