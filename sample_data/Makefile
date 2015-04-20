all: sample.loci.gff3 trinity.diff

clean:
	rm -f ./*loci* ./trinity.diff

trinity.diff: trinity.loci.gff3 trinity.loci.from_gff.gff3
	diff --suppress-common-lines trinity.loci.gff3 trinity.loci.from_gff.gff3 > trinity.diff

sample.loci.gff3: sample.subloci.gff3 sample.scores
	python3 ../sublocus_test.py -s sample.scores sample.subloci.gff3 sample.loci.gff3
sample.subloci.gff3: sample.gff3
	python3 ../superlocus_test.py sample.gff3 sample.subloci.gff3

trinity.loci.gff3: trinity.subloci.gff3 trinity.scores
	python3 ../sublocus_test.py -s trinity.scores trinity.subloci.gff3 trinity.loci.gff3

trinity.subloci.gff3: trinity.gtf
	python3 ../superlocus_test.py trinity.gtf trinity.subloci.gff3

trinity.loci.from_gff.gff3: trinity.subloci.from_gff.gff3 trinity.scores
	python3 ../sublocus_test.py -s trinity.scores trinity.subloci.from_gff.gff3 trinity.loci.from_gff.gff3

trinity.subloci.from_gff.gff3: trinity.gff3
	python3 ../superlocus_test.py trinity.gff3 trinity.subloci.from_gff.gff3
