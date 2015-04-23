all: sample.loci.gff3 trinity.loci.gff3 trinity.loci.from_gff.gff3

clean:
	rm -f ./*loci* ./trinity.diff

sample.loci.gff3: sample.gff3
	python3 ../superlocus_test_whole.py --sub_out sample.subloci.gff3 --mono_out sample.monoloci.gff3 --locus_out sample.loci.gff3 sample.gff3

trinity.loci.gff3: trinity.gtf
	python3 ../superlocus_test_whole.py --sub_out trinity.subloci.gff3 --mono_out trinity.monoloci.gff3 --locus_out trinity.loci.gff3 trinity.gtf

trinity.loci.from_gff.gff3: trinity.gff3 
	python3 ../superlocus_test_whole.py --sub_out trinity.subloci.from_gff.gff3 --mono_out trinity.monoloci.from_gff.gff3 --locus_out trinity.loci.from_gff.gff3 trinity.gff3
