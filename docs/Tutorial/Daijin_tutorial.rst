.. _ERR588044: http://www.ebi.ac.uk/ena/data/view/ERR588044
.. _PRJEB7093: http://www.ebi.ac.uk/ena/data/view/PRJEB7093
.. _ERR588044_1.fastq.gz: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR588/ERR588044/ERR588044_1.fastq.gz
.. _ERR588044_2.fastq.gz: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR588/ERR588044/ERR588044_2.fastq.gz
.. _Daijin-Tutorial:

Tutorial for Daijin
===================

This tutorial will guide you through the task of configuring and running the whole Daijin pipeline on a *A. thaliana* dataset, using one aligner (HISAT) and two assemblers (Stringtie and CLASS2) as chosen methods. A modern desktop computer with a multicore processor and 4GB of RAM or more should suffice to execute the pipeline.


Overview
~~~~~~~~

The tutorial will guide you through the following tasks:

#. Configuring Daijin to analyse the chosen data
#. Running the alignment and assemblies using ``daijin assemble``
#. Running the Mikado analysis using ``daijin mikado``
#. Comparing the results against the reference transcriptome

Required software
~~~~~~~~~~~~~~~~~

Mikado should be installed and configured properly (see our :ref:`installation instructions <Installation>`). Additionally, you should have the following software tools at disposal (between brackets is indicated the version used at the time of the writing):

* BLAST+ (v2.3.0)
* TransDecoder (v3.0.0)
* Portcullis (v0.17.2)
* HISAT2 (v2.0.4)
* Stringtie (v1.2.4)
* CLASS2 (v2.12)
* SAMtools (v1.2)


Input data
~~~~~~~~~~

Throughout this tutorial, we will use data coming from release 32 of EnsEMBL and from the PRJEB7093_ experiment on ENA. In particular, we will need:

 * the `genome FASTA file <ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz>`_ of *Arabidopsis thaliana*
 * its `relative genome annotation <ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.32.gff3.gz>`_
 * RNA-Seq from one sample, ERR588044_, of the PRJEB7093_ study: left (`ERR588044_1.fastq.gz`_) and right (`ERR588044_2.fastq.gz`_)
 * the SwissProt plants database (downloaded from `here <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz>`_, using Uniprot release **2016-07**)

Uncompress the genome, the annotation GFF3 and the Uniprot DAT files using GUNZip. The protein database must be converted into FASTA format, and we should remove the *A. thaliana* proteins to avoid biasing our results. To do so, use the script ``remove_from_embl.py`` (included in Mikado; please include the quotes in the command, or the script will fail!)::

    remove_from_embl.py -o "Arabidopsis thaliana" uniprot_sprot_plants.dat uniprot_sprot_plants.fasta

.. hint:: it is not necessary nor advisable to remove proteins from the database in a real experiment. We are doing so here only to avoid biasing the Mikado results. You can convert the DAT files to FASTA using the same script but using "NULL" in place of "Arabidopsis thaliana".

It is possible to have a feel for the annnotation of this species - the size of its genes and transcripts, the average number of exons per transcript, etc - by using ``mikado util stats``; just issue the following command::

    mikado util stats Arabidopsis_thaliana.TAIR10.32.gff3 Arabidopsis_thaliana.TAIR10.32.gff3.stats

The file *Arabidopsis_thaliana.TAIR10.32.gff3.stats* contains the following information:

+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Stat                             | Total    | Average   | Mode    | Min     | 1%     | 5%   | 10%   | 25%   | Median   | 75%   | 90%   | 95%   | 99%    | Max     |
+==================================+==========+===========+=========+=========+========+======+=======+=======+==========+=======+=======+=======+========+=========+
| Number of genes                  | 33602    | NA        | NA      | NA      | NA     | NA   | NA    | NA    | NA       | NA    | NA    | NA    | NA     | NA      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Number of genes (coding)         | 27397    | NA        | NA      | NA      | NA     | NA   | NA    | NA    | NA       | NA    | NA    | NA    | NA     | NA      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Number of monoexonic genes       | 11079    | NA        | NA      | NA      | NA     | NA   | NA    | NA    | NA       | NA    | NA    | NA    | NA     | NA      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Transcripts per gene             | 41671    | 1.24      | 1       | 1       | 1      | 1    | 1     | 1     | 1        | 1     | 2     | 2     | 4      | 10      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Coding transcripts per gene      | 35359    | 1.05      | 1       | 0       | 0      | 0    | 0     | 1     | 1        | 1     | 2     | 2     | 4      | 10      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDNA lengths                     | 64867032 | 1,556.65  | 72      | 22      | 73     | 255  | 463   | 835   | 1,361    | 1,963 | 2,835 | 3,614 | 5,517  | 16,347  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDNA lengths (mRNAs)             | 54248365 | 1,534.22  | 357     | 22      | 144    | 371  | 555   | 897   | 1,383    | 1,919 | 2,664 | 3,268 | 4,826  | 16,347  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDS lengths                      | 43498974 | 1,043.87  | 0       | 0       | 0      | 0    | 0     | 386   | 903      | 1,449 | 2,154 | 2,724 | 4,238  | 16,182  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDS lengths (mRNAs)              | NA       | 1,043.87  | 0       | 0       | 0      | 0    | 0     | 386   | 903      | 1,449 | 2,154 | 2,724 | 4,238  | 16,182  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDS/cDNA ratio                   | NA       | 78.44     | 100.0   | 2       | 35     | 48   | 56    | 68    | 80       | 92    | 100   | 100   | 100    | 100     |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Monoexonic transcripts           | 11494    | 1,377.97  | 72      | 22      | 71     | 77   | 136   | 465   | 972      | 1,841 | 3,110 | 4,335 | 5,949  | 15,195  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| MonoCDS transcripts              | 7367     | 862.50    | 357     | 22      | 99     | 138  | 192   | 357   | 675      | 1,206 | 1,791 | 2,136 | 2,822  | 6,885   |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Exons per transcript             | 217183   | 5.21      | 1       | 1       | 1      | 1    | 1     | 1     | 3        | 7     | 12    | 15    | 24     | 79      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Exons per transcript (mRNAs)     | 2699     | 5.86      | 1       | 1       | 1      | 1    | 1     | 2     | 4        | 8     | 13    | 16    | 25     | 79      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Exon lengths                     | NA       | 298.67    | 72      | 1       | 32     | 51   | 63    | 89    | 150      | 314   | 615   | 990   | 2,400  | 15,195  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Exon lengths (mRNAs)             | NA       | 261.83    | 72      | 1       | 32     | 51   | 63    | 89    | 147      | 300   | 559   | 857   | 1,736  | 7,761   |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Intron lengths                   | NA       | 165.84    | 87      | 8       | 67     | 74   | 78    | 86    | 100      | 168   | 334   | 461   | 849    | 57,631  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Intron lengths (mRNAs)           | NA       | 164.69    | 87      | 8       | 67     | 75   | 78    | 86    | 100      | 168   | 333   | 458   | 838    | 11,602  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDS exons per transcript         | 2308     | 4.73      | 1       | 0       | 0      | 0    | 0     | 1     | 3        | 7     | 11    | 14    | 24     | 78      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDS exons per transcript (mRNAs) | 2308     | 5.57      | 1       | 1       | 1      | 1    | 1     | 2     | 4        | 8     | 12    | 15    | 24     | 78      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDS exon lengths                 | 43498974 | 220.92    | 72      | 1       | 19     | 45   | 58    | 81    | 127      | 228   | 458   | 735   | 1,572  | 7,761   |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| CDS Intron lengths               | 25183403 | 155.90    | 84      | 7       | 66     | 73   | 77    | 84    | 98       | 155   | 307   | 427   | 789    | 10,233  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| 5'UTR exon number                | 35359    | 0.98      | 1       | 0       | 0      | 0    | 0     | 1     | 1        | 1     | 2     | 2     | 3      | 12      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| 3'UTR exon number                | 35359    | 0.87      | 1       | 0       | 0      | 0    | 0     | 1     | 1        | 1     | 1     | 2     | 2      | 15      |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| 5'UTR length                     | 4119344  | 116.50    | 0       | 0       | 0      | 0    | 0     | 11    | 81       | 163   | 280   | 378   | 609    | 3,214   |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| 3'UTR length                     | 6630047  | 187.51    | 0       | 0       | 0      | 0    | 0     | 77    | 181      | 257   | 348   | 434   | 755    | 3,164   |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Stop distance from junction      | NA       | 7.62      | 0       | 0       | 0      | 0    | 0     | 0     | 0        | 0     | 0     | 10    | 213    | 2,286   |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Intergenic distances             | NA       | 1,429.33  | 254;260 | -57,916 | -1,251 | -42  | 68    | 279   | 822      | 1,899 | 3,619 | 4,996 | 8,734  | 72,645  |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+
| Intergenic distances (coding)    | NA       | 2,163.73  | 187     | -10,136 | -571   | -34  | 70    | 284   | 866      | 2,143 | 4,500 | 6,878 | 20,134 | 490,690 |
+----------------------------------+----------+-----------+---------+---------+--------+------+-------+-------+----------+-------+-------+-------+--------+---------+

From this summary it is quite apparent that the *A. thaliana* genome preferentially encodes multiexonic transcripts with a median CDS/UTR ratio of approximately 80%. Each transcript has on average 5 exons (median 3), and intron lengths are generally short - 165 bps on average and 100 as median, with a maximum value of ~11,000 bps for mRNAs. It is also a compact genome, with a median intergenic distance under 1 kbps and a modal distance of only 250 bps.

Step 1: configuring Daijin
~~~~~~~~~~~~~~~~~~~~~~~~~~

The first task is to create a configuration file for Daijin using ``daijin configure``. On the command line, we have to configure the following:

* name of the configuration file (daijin.yaml)
* number of threads per process (5)
* reference genome and the name of the species (**without spaces or non-ASCII/special characters**)
* reads to use, with their strandedness and the name of the sample
* aligners (HISAT2) and assemblers (Stringtie, CLASS2) to use
* output directory (athaliana)
* the scoring file to use for Mikado pick; we will ask to copy it in-place to have a look at it (plants.yaml)
* the protein database for homology searches for Mikado (uniprot_sprot_plants.fasta)
* flank: as *A. thaliana* has a quite compact genome, we should decrease the maximum distance for grouping together transcripts. We will decrease from 1kbps (default) to 250.
* (optional) the scheduler to use for the cluster (we will presume SLURM, but we also support LSF and PBS)
* (optional) name of the cluster configuration file, which will have to be edited manually.

This can all be specified with the following command line::

  daijin configure \
    -o daijin.yaml \
    --flank 250 \
    --threads 5 \
    --genome Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --name "Athaliana" \
    -od athaliana \
    -r1 ERR588044_1.fastq.gz -r2 ERR588044_2.fastq.gz --strandedness fr-secondstrand -s ERR588044 \
    -al hisat -as stringtie class \
    --scoring plants.yaml --copy-scoring plants.yaml \
    --prot-db uniprot_sprot_plants.fasta \
    --scheduler SLURM \
    -c slurm.yaml;

This will produce a file, *daijin.yaml*, which should specify the same parameters as :download:`this preformed example <daijin.yaml>`.

If you are executing this program on a cluster, you will have to modify the "load" section and specify for each software which operations should be performed to make the tool available to the submitted script (eg. "source mikado-1.0.0" to make Mikado available). For example, this is how our load section looks like after editing:

.. code-block:: yaml

  load:
      #  Commands to use to load/select the versions of the programs to use. Leave an empty
      #  string if no loading is necessary.
      blast: 'source blast-2.3.0'
      class: 'source class-2.12'
      cufflinks: ''
      gmap: ''
      hisat: 'source HISAT-2.0.4'
      mikado: 'source mikado-devel'
      portcullis: 'source portcullis-0.17.2'
      samtools: 'source samtools-1.2'
      star: ''
      stringtie: 'source stringtie-1.2.4'
      tophat: ''
      transdecoder: 'source transdecoder-3.0.0'
      trinity: ''

Also, if you are operating on a cluster, you might want to edit the cluster configuration file "slurm.conf". This is how it appears:

.. literalinclude:: ../Usage/hpc.yaml

The most important parameter to modify is the queue to use - set it according to the configuration of your own cluster.


Step 2: running the assemble part
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created a proper configuration file, it is time to launch Daijin assemble and inspect the results. Issue the command::

    daijin assemble -c slurm.yaml daijin.yaml

After checking that the configuration file is valid, Daijin will start the alignment and assembly of the dataset. On a normal desktop computer, this should take ~2 hours. Before launching the pipeline, you can obtain a graphical representation of the steps with::

  daijin assemble -c slurm.yaml --dag daijin.yaml | dot -Tsvg > assemble.svg

You can also ask Daijin to display the steps to be executed, inclusive of their command lines, by issuing the following command::

  daijin assemble -c slurm.yaml --dryrun daijin.yaml

When Daijin is finished, have a look inside the folder athaliana/3-assemblies/output/; you will find the following two GTF files:

* class-0-hisat--0.gtf
* stringtie-0-hisat--0.gtf

These are standard GTF files reporting the assembled transcripts for each method. We can have a feel for how they compare with our reference annotation by, again, using ``mikado util stats``. Conveniently, Daijin has already performed this analysis for us, and the files will be present in the same folder:

* class-0-hisat-ERR588044-0.gtf**.stats**
* stringtie-0-hisat-ERR588044-0.gtf**.stats**

From this quick analysis, it seems that Stringtie has assembled - correctly - mostly monoexonic transcripts, but has probably massively overestimated the number of genes (11,062, almost 50% more than those found in the reference annotation). CLASS2 has been more conservative, but contrary to Stringtie it has assembled mostly multiexonic transcripts. Stringtie transcripts are also much shorter than CLASS2 transcripts, suggesting that this tool might have reported many fragmentary loci.

It is possible to compare them to the current annotation by using :ref:`Mikado compare <compare>`. To do so, issue the following commands:

.. code-block:: bash

  # First index the reference GFF3
  mikado compare -r Arabidopsis_thaliana.TAIR10.32.gff3 --index;
  # Compare the CLASS2 assembly against the reference
  mikado compare -r Arabidopsis_thaliana.TAIR10.32.gff3 -p athaliana/3-assemblies/output/class-0-hisat-ERR588044-0.gtf -o athaliana/3-assemblies/output/class-0-hisat-ERR588044-0.compare -l athaliana/3-assemblies/output/class-0-hisat-ERR588044-0.log
  # Compare the Stringtie assembly against the reference
  mikado compare -r Arabidopsis_thaliana.TAIR10.32.gff3 -p athaliana/3-assemblies/output/stringtie-0-hisat-ERR588044-0.gtf -o athaliana/3-assemblies/output/stringtie-0-hisat-ERR588044-0.compare -l athaliana/3-assemblies/output/stringtie-0-hisat-ERR588044-0.compare.log

The analysis will produce *TMAP*, *REFMAP* and *STATS* files for each of the two assemblies. This is the report of the stats file for the CLASS2 assembly, ::

  41671 reference RNAs in 33602 genes
  35276 predicted RNAs in  29029 genes
  --------------------------------- |   Sn |   Pr |   F1 |
                          Base level: 40.73  80.36  54.06
              Exon level (stringent): 43.80  56.22  49.24
                Exon level (lenient): 65.81  74.31  69.80
                        Intron level: 68.00  86.79  76.25
                  Intron chain level: 28.37  26.88  27.61
        Transcript level (stringent): 0.00  0.00  0.00
    Transcript level (>=95% base F1): 2.30  2.70  2.48
    Transcript level (>=80% base F1): 16.33  19.19  17.64
           Gene level (100% base F1): 0.00  0.00  0.00
          Gene level (>=95% base F1): 2.77  3.20  2.97
          Gene level (>=80% base F1): 19.43  22.37  20.80

  #   Matching: in prediction; matched: in reference.

              Matching intron chains: 8549
               Matched intron chains: 8575
     Matching monoexonic transcripts: 870
      Matched monoexonic transcripts: 883
          Total matching transcripts: 9419
           Total matched transcripts: 9458

            Missed exons (stringent): 95133/169282  (56.20%)
             Novel exons (stringent): 57739/131888  (43.78%)
              Missed exons (lenient): 52381/153221  (34.19%)
               Novel exons (lenient): 34867/135707  (25.69%)
                      Missed introns: 40932/127896  (32.00%)
                       Novel introns: 13234/100198  (13.21%)

                  Missed transcripts: 15490/41671  (37.17%)
                   Novel transcripts: 99/35276  (0.28%)
                        Missed genes: 14849/33602  (44.19%)
                         Novel genes: 70/29029  (0.24%)

and this instead for Stringtie::

  /tgac/software/testing/bin/core/../..//mikado/devel/x86_64/bin/mikado compare -r Arabidopsis_thaliana.TAIR10.32.gff3 -p athaliana/3-assemblies/output/stringtie-0-hisat-ERR588044-0.gtf -o athaliana/3-assemblies/output/stringtie-0-hisat-ERR588044-0.compare -l athaliana/3-assemblies/output/stringtie-0-hisat-ERR588044-0.compare.log
  41671 reference RNAs in 33602 genes
  55055 predicted RNAs in  40566 genes
  --------------------------------- |   Sn |   Pr |   F1 |
                          Base level: 55.21  59.42  57.23
              Exon level (stringent): 45.31  38.81  41.81
                Exon level (lenient): 68.12  57.04  62.09
                        Intron level: 72.94  64.23  68.31
                  Intron chain level: 42.54  29.94  35.15
        Transcript level (stringent): 0.00  0.00  0.00
    Transcript level (>=95% base F1): 21.26  18.02  19.50
    Transcript level (>=80% base F1): 34.77  32.56  33.63
           Gene level (100% base F1): 0.00  0.00  0.00
          Gene level (>=95% base F1): 24.76  22.49  23.57
          Gene level (>=80% base F1): 40.20  40.28  40.24

  #   Matching: in prediction; matched: in reference.

              Matching intron chains: 12817
               Matched intron chains: 12854
     Matching monoexonic transcripts: 1827
      Matched monoexonic transcripts: 1836
          Total matching transcripts: 14644
           Total matched transcripts: 14690

            Missed exons (stringent): 92576/169282  (54.69%)
             Novel exons (stringent): 120959/197665  (61.19%)
              Missed exons (lenient): 49435/155043  (31.88%)
               Novel exons (lenient): 79549/185157  (42.96%)
                      Missed introns: 34606/127896  (27.06%)
                       Novel introns: 51958/145248  (35.77%)

                  Missed transcripts: 12822/41671  (30.77%)
                   Novel transcripts: 155/55055  (0.28%)
                        Missed genes: 12288/33602  (36.57%)
                         Novel genes: 138/40566  (0.34%)

These comparisons suggest the following:

#. Most of the transcripts reconstructed by both assemblers are near enough known genes so as not to be categorised as completely novel. Only 70 genes for CLASS2 and 130 genes for Stringtie are in non-annotated genomic loci.
#. CLASS2 has performed worse than Stringtie for this particular species and samples - the number of reconstructed transcript is lower, and so is the F1 at the transcript and gene level.
#. Although the number of novel loci is low for both CLASS and Stringtie, the high number of loci


Step 3: running the Mikado steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created the input assemblies, it is time to run Mikado. First of all, let us have a look at the *spombe.yaml* file, which we have conveniently copied to our current directory:

.. literalinclude:: plants.yaml

With this file, we are telling Mikado that we are looking for transcripts with a low number of exons ("exon_num: {rescaling: min}"), with a ratio of CDS/UTR of ~80% ("selected_cds_fraction: {rescaling: target, value: 0.8}") and good homology support ("blast_score: {rescaling: max}"). We are also rewarding long cDNAs ("cdna_length: {rescaling: max}") with a good homology to known proteins ("blast_score: {rescaling: max}") and only one ORF ("cds_not_maximal: {rescaling: min}"; "number_internal_orfs: {rescaling: target, value: 1}"). This set of rules is tailored for compact fungal genomes, where this kind of models is expected at a much higher rate than in other eukaryotes (eg plants or mammals).

Before launching ``daijin mikado``, we can have a look at the BLAST and TransDecoder sections. It is advisable to modify the number of chunks for BLASTX to a high number, eg. 100, if you are using a cluster. Please note that now we are going to use a different configuration file, **athaliana/mikado.yaml**, which Daijin created as last step in the assembly pipeline.

Issue the command::

  daijin mikado -c slurm.yaml athaliana/mikado.yaml

This part of the pipeline should be quicker than the previous stage. After the pipeline is finished, Daijin will have created the final output files in athaliana/5-mikado/pick/. As we requested only for the *permissive* mode, we only have one output - *athaliana/5-mikado/pick/mikado-permissive.loci.gff3*. These are basic statistics on this annotation:

+------------------------------+---------+--------------------+---------------+------------------------+-----------------------+--------------------------+---------+------------------------+-----------------+
| File                         |   genes |   monoexonic_genes |   transcripts |   transcripts_per_gene |   transcript_len_mean |   monoexonic_transcripts |   exons |   exons_per_transcript |   exon_len_mean |
+==============================+=========+====================+===============+========================+=======================+==========================+=========+========================+=================+
| mikado-permissive.loci.stats |    5633 |               2915 |          6054 |                   1.07 |               1642.08 |                     3095 |   11850 |                   1.96 |          838.92 |
+------------------------------+---------+--------------------+---------------+------------------------+-----------------------+--------------------------+---------+------------------------+-----------------+

Mikado has created an annotation that is in between those produced by Stringtie and CLASS2. Compared to CLASS2, the Mikado assembly is probably more comprehensive, given the higher number of genes. Half the genes are monoexonic, an improvement compared to CLASS2, and at the same time the number of genes is much lower than in Stringtie, with a higher average cDNA length. These statistics suggest that the Mikado filtered annotation has been able to retain most of the real genes, while discarding many of the fragments present in either assembly. We can verify this by comparing the Mikado results against the reference annotation::

    mikado compare -r Schizosaccharomyces_pombe.ASM294v2.32.gff3 -p athaliana/5-mikado/pick/permissive/mikado-permissive.loci.gff3 -o athaliana/5-mikado/mikado_prepared.compare -l athaliana/5-mikado/mikado_prepared.compare.log

A cursory look at the STATS file confirms the impression; the Mikado annotation has removed many of the spurious transcripts and is quite more precise than either of the input assemblies, while retaining their sensitivity::

    Command line:
    /tgac/software/testing/bin/core/../..//mikado/devel/x86_64/bin/mikado compare -r Arabidopsis_thaliana.TAIR10.32.gff3 -p athaliana/5-mikado/pick/permissive/mikado-permissive.loci.gff3 -o athaliana/5-mikado/pick/permissive/mikado-permissive.loci.compare -l athaliana/5-mikado/pick/permissive/mikado-permissive.loci.compare.log
    41671 reference RNAs in 33602 genes
    30806 predicted RNAs in  24730 genes
    --------------------------------- |   Sn |   Pr |   F1 |
                            Base level: 53.43  87.95  66.48
                Exon level (stringent): 45.44  58.77  51.26
                  Exon level (lenient): 68.71  86.45  76.57
                          Intron level: 73.13  93.43  82.05
                    Intron chain level: 43.73  53.48  48.11
          Transcript level (stringent): 0.01  0.01  0.01
      Transcript level (>=95% base F1): 18.87  25.49  21.69
      Transcript level (>=80% base F1): 33.71  45.78  38.83
             Gene level (100% base F1): 0.01  0.02  0.01
            Gene level (>=95% base F1): 22.36  30.35  25.75
            Gene level (>=80% base F1): 39.52  53.95  45.62

    #   Matching: in prediction; matched: in reference.

                Matching intron chains: 13173
                 Matched intron chains: 13213
       Matching monoexonic transcripts: 1702
        Matched monoexonic transcripts: 1710
            Total matching transcripts: 14875
             Total matched transcripts: 14923

              Missed exons (stringent): 92353/169282  (54.56%)
               Novel exons (stringent): 53962/130891  (41.23%)
                Missed exons (lenient): 48471/154914  (31.29%)
                 Novel exons (lenient): 16682/123125  (13.55%)
                        Missed introns: 34362/127896  (26.87%)
                         Novel introns: 6576/100110  (6.57%)

                    Missed transcripts: 13373/41671  (32.09%)
                     Novel transcripts: 132/30806  (0.43%)
                          Missed genes: 12824/33602  (38.16%)
                           Novel genes: 120/24730  (0.49%)



Moreover, Mikado models have an ORF assigned to them. We can ask Mikado compare to consider only the coding component of transcripts with the following command line::

    mikado compare -eu -r Arabidopsis_thaliana.TAIR10.32.gff3 -p athaliana/5-mikado/pick/permissive/mikado-permissive.loci.gff3 -o athaliana/5-mikado/pick/permissive/mikado-permissive.loci.compare.eu -l athaliana/5-mikado/pick/permissive/mikado-permissive.loci.compare.eu.log

The statistics file looks as follows::

    Command line:
    mikado compare -eu -r Arabidopsis_thaliana.TAIR10.32.gff3 -p athaliana/5-mikado/pick/permissive/mikado-permissive.loci.gff3 -o athaliana/5-mikado/pick/permissive/mikado-permissive.loci.compare.eu -l athaliana/5-mikado/pick/permissive/mikado-permissive.loci.compare.eu.log
    41671 reference RNAs in 33602 genes
    30806 predicted RNAs in  24730 genes
    --------------------------------- |   Sn |   Pr |   F1 |
                            Base level: 53.00  93.07  67.54
                Exon level (stringent): 60.90  78.70  68.67
                  Exon level (lenient): 69.57  88.37  77.85
                          Intron level: 73.44  94.86  82.79
                    Intron chain level: 46.39  57.65  51.41
          Transcript level (stringent): 25.84  34.31  29.48
      Transcript level (>=95% base F1): 35.32  47.28  40.43
      Transcript level (>=80% base F1): 39.45  53.03  45.25
             Gene level (100% base F1): 26.78  36.40  30.86
            Gene level (>=95% base F1): 37.63  51.16  43.37
            Gene level (>=80% base F1): 42.29  57.54  48.75

    #   Matching: in prediction; matched: in reference.

                Matching intron chains: 13988
                 Matched intron chains: 14097
       Matching monoexonic transcripts: 2687
        Matched monoexonic transcripts: 2693
            Total matching transcripts: 16675
             Total matched transcripts: 16790

              Missed exons (stringent): 61443/157150  (39.10%)
               Novel exons (stringent): 25902/121609  (21.30%)
                Missed exons (lenient): 44733/146988  (30.43%)
                 Novel exons (lenient): 13456/115711  (11.63%)
                        Missed introns: 32026/120580  (26.56%)
                         Novel introns: 4795/93349  (5.14%)

                    Missed transcripts: 13627/41671  (32.70%)
                     Novel transcripts: 141/30806  (0.46%)
                          Missed genes: 13060/33602  (38.87%)
                           Novel genes: 127/24730  (0.51%)

The similarity is quite higher, suggesting that for many models the differences between the Mikado annotation and the reference lies in the UTR component.

We suggest to visualise assemblies with one of the many tools currently at disposal, such as eg `WebApollo <http://genomearchitect.org/>`_ [Apollo]_. Mikado files are GFF3-compliant and can be loaded directly into Apollo or similar tools. GTF files can be converted into proper GFF3 files using the convert utility::

  mikado util convert <input gtf> <output GFF3>

