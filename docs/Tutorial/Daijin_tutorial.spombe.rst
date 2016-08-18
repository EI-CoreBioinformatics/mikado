.. _Daijin-Tutorial:

Tutorial for Daijin
===================

This tutorial will guide you through the task of configuring and running the whole Daijin pipeline on a *S. pombe* dataset, using one aligner (HISAT) and two assemblers (Stringtie and CLASS2) as chosen methods. A modern desktop computer with a multicore processor and 4GB of RAM or more should suffice to execute the pipeline.


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

Throughout this tutorial, we will use data coming from release 32 of EnsEMBL and from the SRR1617247 experiment on ENA. In particular, we will need:

 * the `genome FASTA file <ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz>`_ of *Schizosaccharomyces pombe*
 * its `relative genome annotation <ftp://ftp.ensemblgenomes.org/pub/fungi/release-32/gff3/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.32.gff3.gz>`_
 * RNA-Seq from the SRR1617247 experiment: `left <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/007/SRR1617247/SRR1617247_1.fastq.gz>`_ and `right <ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR161/007/SRR1617247/SRR1617247_2.fastq.gz>`_.
 * the SwissProt fungi database (downloaded from `here <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_fungi.dat.gz>`_, using Uniprot release **2016-07**)

Uncompress the genome, the annotation GFF3 and the Uniprot DAT files using GUNZip. The protein database must be converted into FASTA format, and we should remove the *Sc. pombe* proteins to avoid biasing our results. To do so, use the script ``remove_from_embl.py`` (included in Mikado)::

    remove_from_embl.py --format fasta -o "Schizosaccharomyces pombe" uniprot_sprot_fungi.dat uniprot_sprot_fungi.fasta

It is possible to have a feel for the annnotation of this species - the size of its genes and transcripts, the average number of exons per transcript, etc - by using ``mikado util stats``; just issue the following command::

    mikado util stats Schizosaccharomyces_pombe.ASM294v2.32.gff3 Schizosaccharomyces_pombe.ASM294v2.32.gff3.stats

The file *Schizosaccharomyces_pombe.ASM294v2.32.gff3.stats* contains the following information:

+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Stat                             | Total    | Average   | Mode         | Min    | 1%     | 5%     | 10%    | 25%   | Median   | 75%   | 90%   | 95%   | 99%   | Max    |
+==================================+==========+===========+==============+========+========+========+========+=======+==========+=======+=======+=======+=======+========+
| Number of genes                  | 7014     | NA        | NA           | NA     | NA     | NA     | NA     | NA    | NA       | NA    | NA    | NA    | NA    | NA     |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Number of genes (coding)         | 5145     | NA        | NA           | NA     | NA     | NA     | NA     | NA    | NA       | NA    | NA    | NA    | NA    | NA     |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Number of monoexonic genes       | 4461     | NA        | NA           | NA     | NA     | NA     | NA     | NA    | NA       | NA    | NA    | NA    | NA    | NA     |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Transcripts per gene             | 7015     | 1.00      | 1            | 1      | 1      | 1      | 1      | 1     | 1        | 1     | 1     | 1     | 1     | 2      |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Coding transcripts per gene      | 5146     | 0.73      | 1            | 0      | 0      | 0      | 0      | 0     | 1        | 1     | 1     | 1     | 1     | 2      |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDNA lengths                     | 12586389 | 1,794.21  | 72           | 47     | 72     | 185    | 448    | 922   | 1,549    | 2,363 | 3,400 | 4,151 | 5,894 | 15,022 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDNA lengths (mRNAs)             | 10592917 | 2,058.48  | 1315         | 75     | 283    | 546    | 753    | 1,201 | 1,792    | 2,613 | 3,675 | 4,425 | 6,404 | 15,022 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDS lengths                      | 7178717  | 1,023.34  | 0            | 0      | 0      | 0      | 0      | 0     | 795      | 1,500 | 2,339 | 3,045 | 4,994 | 14,775 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDS lengths (mRNAs)              | NA       | 1,023.34  | 0            | 0      | 0      | 0      | 0      | 0     | 795      | 1,500 | 2,339 | 3,045 | 4,994 | 14,775 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDS/cDNA ratio                   | NA       | 67.48     | 100.0        | 5      | 15     | 28     | 37     | 53    | 70       | 84    | 93    | 100   | 100   | 100    |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Monoexonic transcripts           | 4461     | 1,745.01  | 72           | 47     | 72     | 100    | 336    | 811   | 1,489    | 2,344 | 3,432 | 4,179 | 6,019 | 14,362 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| MonoCDS transcripts              | 2753     | 1,485.99  | 4002;375;432 | 93     | 218    | 332    | 418    | 732   | 1,191    | 1,806 | 2,849 | 3,784 | 5,875 | 14,154 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Exons per transcript             | 12379    | 1.76      | 1            | 1      | 1      | 1      | 1      | 1     | 1        | 2     | 4     | 4     | 7     | 16     |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Exons per transcript (mRNAs)     | 3081     | 2.03      | 1            | 1      | 1      | 1      | 1      | 1     | 1        | 3     | 4     | 5     | 7     | 16     |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Exon lengths                     | NA       | 1,016.75  | 72           | 2      | 25     | 61     | 83     | 184   | 561      | 1,502 | 2,548 | 3,326 | 5,015 | 14,362 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Exon lengths (mRNAs)             | NA       | 1,013.39  | 106          | 3      | 24     | 59     | 86     | 175   | 499      | 1,516 | 2,605 | 3,414 | 5,117 | 14,362 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Intron lengths                   | NA       | 83.69     | 49           | 1      | 17     | 38     | 41     | 46    | 56       | 85    | 162   | 226   | 411   | 2,526  |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Intron lengths (mRNAs)           | NA       | 84.16     | 46           | 1      | 36     | 39     | 41     | 46    | 56       | 86    | 162   | 227   | 412   | 2,526  |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDS exons per transcript         | 2091     | 1.46      | 1            | 0      | 0      | 0      | 0      | 0     | 1        | 2     | 3     | 4     | 7     | 16     |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDS exons per transcript (mRNAs) | 2091     | 1.99      | 1            | 1      | 1      | 1      | 1      | 1     | 1        | 3     | 4     | 5     | 7     | 16     |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDS exon lengths                 | 7178717  | 702.76    | 69;99        | 1      | 8      | 29     | 49     | 106   | 307      | 990   | 1,797 | 2,474 | 4,370 | 14,154 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| CDS Intron lengths               | 414501   | 81.77     | 44;45        | 0      | 35     | 38     | 40     | 45    | 55       | 84    | 157   | 220   | 395   | 2,525  |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| 5'UTR exon number                | 5146     | 0.95      | 1            | 0      | 0      | 0      | 1      | 1     | 1        | 1     | 1     | 1     | 2     | 3      |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| 3'UTR exon number                | 5146     | 0.93      | 1            | 0      | 0      | 0      | 1      | 1     | 1        | 1     | 1     | 1     | 2     | 3      |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| 5'UTR length                     | 1372304  | 266.67    | 0            | 0      | 0      | 0      | 17     | 71    | 154      | 309   | 586   | 932   | 1,935 | 4,397  |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| 3'UTR length                     | 2041896  | 396.79    | 0            | 0      | 0      | 0      | 46     | 126   | 243      | 441   | 865   | 1,386 | 2,644 | 5,911  |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Stop distance from junction      | NA       | 7.58      | 0            | 0      | 0      | 0      | 0      | 0     | 0        | 0     | 0     | 0     | 23    | 3,385  |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Intergenic distances             | NA       | -59.13    | -66          | -9,461 | -3,779 | -2,220 | -1,425 | -188  | 80       | 383   | 857   | 1,303 | 2,734 | 31,961 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+
| Intergenic distances (coding)    | NA       | 297.72    | -66          | -7,815 | -3,477 | -1,598 | -440   | -32   | 176      | 600   | 1,302 | 1,913 | 3,924 | 78,421 |
+----------------------------------+----------+-----------+--------------+--------+--------+--------+--------+-------+----------+-------+-------+-------+-------+--------+

From this summary it is quite apparent that the *Sc. pombe* genome preferentially encodes single-exon transcripts, and that the multiexonic transcripts have relatively short introns (average of 84 bps, median of 56 bps).

Step 1: configuring Daijin
~~~~~~~~~~~~~~~~~~~~~~~~~~

The first task is to create a configuration file for Daijin using ``daijin configure``. On the command line, we have to configure the following:

* name of the configuration file (daijin.yaml)
* number of threads per process (5)
* reference genome and the name of the species (**without spaces or non-ASCII/special characters**)
* reads to use, with their strandedness and the name of the sample
* aligners (HISAT2) and assemblers (Stringtie, CLASS2) to use
* output directory (daijin_spombe)
* the scoring file to use for Mikado pick; we will ask to copy it in-place to have a look at it (spombe.yaml)
* the protein database for homology searches for Mikado (uniprot_sprot.fungi.fasta)
* (optional) the scheduler to use for the cluster (we will presume SLURM, but we also support LSF and PBS)
* (optional) name of the cluster configuration file, which will have to be edited manually.

This can all be specified with the following command line::

    daijin configure \
        -o daijin.yaml \
        --threads 5 \
        --genome Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa --name "Sc_pombe" \
        -od daijin_spombe \
        -r1 SRR1617247_1.fastq.gz -r2 SRR1617247_2.fastq.gz --strandedness fr-secondstrand -s SRR1617247 \
        -al hisat -as stringtie class \
        --scoring spombe.yaml --copy-scoring spombe.yaml \
        --prot-db uniprot_sprot_fungi.fasta \
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

When Daijin is finished, have a look inside the folder daijin_spombe/3-assemblies/output/; you will find the following two GTF files:

* class-0-hisat-SRR1617247-0.gtf
* stringtie-0-hisat-SRR1617247-0.gtf

These are standard GTF files reporting the assembled transcripts for each method. We can have a feel for how they compare with our reference annotation by, again, using ``mikado util stats``. Conveniently, Daijin has already performed this analysis for us, and the files will be present in the same folder:

* class-0-hisat-SRR1617247-0.gtf**.stats**
* stringtie-0-hisat-SRR1617247-0.gtf**.stats**

+------------------------------------------+---------+--------------------+---------------+------------------------+-----------------------+--------------------------+---------+------------------------+-----------------+
| File                                     |   genes |   monoexonic_genes |   transcripts |   transcripts_per_gene |   transcript_len_mean |   monoexonic_transcripts |   exons |   exons_per_transcript |   exon_len_mean |
+==========================================+=========+====================+===============+========================+=======================+==========================+=========+========================+=================+
| class-0-hisat-SRR1617247-0.gtf.stats     |    4148 |                781 |          4415 |                   1.05 |               1604.63 |                      781 |   10944 |                   2.48 |          647.34 |
+------------------------------------------+---------+--------------------+---------------+------------------------+-----------------------+--------------------------+---------+------------------------+-----------------+
| stringtie-0-hisat-SRR1617247-0.gtf.stats |   11062 |               7938 |         12625 |                   1.11 |               1188.7  |                     8236 |   21427 |                   1.7  |          700.39 |
+------------------------------------------+---------+--------------------+---------------+------------------------+-----------------------+--------------------------+---------+------------------------+-----------------+

From this quick analysis, it seems that Stringtie has assembled - correctly - mostly monoexonic transcripts, but has probably massively overestimated the number of genes (11,062, almost 50% more than those found in the reference annotation). CLASS2 has been more conservative, but contrary to Stringtie it has assembled mostly multiexonic transcripts. Stringtie transcripts are also much shorter than CLASS2 transcripts, suggesting that this tool might have reported many fragmentary loci.

It is possible to compare them to the current annotation by using :ref:`Mikado compare <compare>`. To do so, issue the following commands:

.. code-block:: bash

  mikado compare -r Schizosaccharomyces_pombe.ASM294v2.32.gff3 -p daijin_spombe/3-assemblies/output/stringtie-0-hisat-SRR1617247-0.gtf -o daijin_spombe/3-assemblies/output/stringtie-0-hisat-SRR1617247-0.compare -l daijin_spombe/3-assemblies/output/stringtie-0-hisat-SRR1617247-0.compare.log;
  mikado compare -r Schizosaccharomyces_pombe.ASM294v2.32.gff3 -p daijin_spombe/3-assemblies/output/class-0-hisat-SRR1617247-0.gtf -o daijin_spombe/3-assemblies/output/class-0-hisat-SRR1617247-0.compare -l daijin_spombe/3-assemblies/output/class-0-hisat-SRR1617247-0.compare.log;

The analysis will produce *TMAP*, *REFMAP* and *STATS* files for each of the two assemblies. This is the report of the stats file for the CLASS2 assembly::

  /tgac/software/testing/bin/core/../..//mikado/devel/x86_64/bin/mikado compare -r Schizosaccharomyces_pombe.ASM294v2.32.gff3 -p daijin_spombe/3-assemblies/output/class-0-hisat-SRR1617247-0.gtf -o daijin_spombe/3-assemblies/output/class-0-hisat-SRR1617247-0.compare -l daijin_spombe/3-assemblies/output/class-0-hisat-SRR1617247-0.compare.log
  7015 reference RNAs in 7014 genes
  4415 predicted RNAs in  4148 genes
  --------------------------------- |   Sn |   Pr |   F1 |
                          Base level: 34.43  69.23  45.98
              Exon level (stringent): 19.11  22.61  20.72
                Exon level (lenient): 81.19  65.65  72.60
                        Intron level: 86.12  74.36  79.81
                  Intron chain level: 68.68  48.21  56.65
        Transcript level (stringent): 0.00  0.00  0.00
    Transcript level (>=95% base F1): 3.12  4.92  3.82
    Transcript level (>=80% base F1): 17.01  26.66  20.77
           Gene level (100% base F1): 0.00  0.00  0.00
          Gene level (>=95% base F1): 3.12  5.21  3.90
          Gene level (>=80% base F1): 17.01  28.28  21.24

  #   Matching: in prediction; matched: in reference.

              Matching intron chains: 1752
               Matched intron chains: 1753
     Matching monoexonic transcripts: 310
      Matched monoexonic transcripts: 326
          Total matching transcripts: 2062
           Total matched transcripts: 2079

            Missed exons (stringent): 10012/12378  (80.89%)
             Novel exons (stringent): 8097/10463  (77.39%)
              Missed exons (lenient): 1490/7923  (18.81%)
               Novel exons (lenient): 3366/9799  (34.35%)
                      Missed introns: 744/5361  (13.88%)
                       Novel introns: 1592/6209  (25.64%)

                  Missed transcripts: 2937/7015  (41.87%)
                   Novel transcripts: 3/4415  (0.07%)
                        Missed genes: 2936/7014  (41.86%)
                         Novel genes: 3/4148  (0.07%)

and this instead for Stringtie::

  Command line:
  /tgac/software/testing/bin/core/../..//mikado/devel/x86_64/bin/mikado compare -r Schizosaccharomyces_pombe.ASM294v2.32.gff3 -p daijin_spombe/3-assemblies/output/stringtie-0-hisat-SRR1617247-0.gtf -o daijin_spombe/3-assemblies/output/stringtie-0-hisat-SRR1617247-0.compare -l daijin_spombe/3-assemblies/output/stringtie-0-hisat-SRR1617247-0.compare.log
  7015 reference RNAs in 7014 genes
  12625 predicted RNAs in  11062 genes
  --------------------------------- |   Sn |   Pr |   F1 |
                          Base level: 70.12  67.71  68.89
              Exon level (stringent): 19.88  12.74  15.53
                Exon level (lenient): 68.66  60.10  64.10
                        Intron level: 88.73  70.16  78.36
                  Intron chain level: 73.81  42.99  54.34
        Transcript level (stringent): 0.00  0.00  0.00
    Transcript level (>=95% base F1): 25.69  14.51  18.55
    Transcript level (>=80% base F1): 47.68  27.84  35.16
           Gene level (100% base F1): 0.00  0.00  0.00
          Gene level (>=95% base F1): 25.69  16.14  19.82
          Gene level (>=80% base F1): 47.69  30.69  37.35

  #   Matching: in prediction; matched: in reference.

              Matching intron chains: 1890
               Matched intron chains: 1884
     Matching monoexonic transcripts: 1739
      Matched monoexonic transcripts: 1746
          Total matching transcripts: 3629
           Total matched transcripts: 3630

            Missed exons (stringent): 9917/12378  (80.12%)
             Novel exons (stringent): 16857/19318  (87.26%)
              Missed exons (lenient): 3027/9658  (31.34%)
               Novel exons (lenient): 4402/11033  (39.90%)
                      Missed introns: 604/5361  (11.27%)
                       Novel introns: 2023/6780  (29.84%)

                  Missed transcripts: 523/7015  (7.46%)
                   Novel transcripts: 9/12625  (0.07%)
                        Missed genes: 523/7014  (7.46%)
                         Novel genes: 9/11062  (0.08%)

These comparisons suggest the following:

#. Most of the transcripts reconstructed by both assemblers are near enough known genes so as not to be categorised as completely novel. Only 3 genes for CLASS2 and 9 genes for Stringtie are in non-annotated genomic loci.
#. CLASS2 has performed worse than Stringtie for this particular species - the number of reconstructed transcript is lower, and so is the F1 at the transcript and gene level.
#. Gene precision is pretty low for Stringtie - 30% - suggesting that the assembler has probably detected many spurious loci near annotated ones.


Step 3: running the Mikado steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created the input assemblies, it is time to run Mikado. First of all, let us have a look at the *spombe.yaml* file, which we have conveniently copied to our current directory:

.. literalinclude:: spombe.yaml

With this file, we are telling Mikado that we are looking for transcripts with a low number of exons ("exon_num: {rescaling: min}"), with a ratio of CDS/UTR of ~80% ("selected_cds_fraction: {rescaling: target, value: 0.8}") and good homology support ("blast_score: {rescaling: max}"). We are also rewarding long cDNAs ("cdna_length: {rescaling: max}") with a good homology to known proteins ("blast_score: {rescaling: max}") and only one ORF ("cds_not_maximal: {rescaling: min}"; "number_internal_orfs: {rescaling: target, value: 1}"). This set of rules is tailored for compact fungal genomes, where this kind of models is expected at a much higher rate than in other eukaryotes (eg plants or mammals).

Before launching ``daijin mikado``, we can have a look at the BLAST and TransDecoder sections. It is advisable to modify the number of chunks for BLASTX to a high number, eg. 100, if you are using a cluster. Please note that now we are going to use a different configuration file, **daijin_spombe/mikado.yaml**, which Daijin created as last step in the assembly pipeline.

Issue the command::

  daijin mikado -c slurm.yaml daijin_spombe/mikado.yaml

This part of the pipeline should be quicker than the previous stage. After the pipeline is finished, Daijin will have created the final output files in daijin_spombe/5-mikado/pick/. As we requested only for the *permissive* mode, we only have one output - *daijin_spombe/5-mikado/pick/mikado-permissive.loci.gff3*. These are basic statistics on this annotation:

+------------------------------+---------+--------------------+---------------+------------------------+-----------------------+--------------------------+---------+------------------------+-----------------+
| File                         |   genes |   monoexonic_genes |   transcripts |   transcripts_per_gene |   transcript_len_mean |   monoexonic_transcripts |   exons |   exons_per_transcript |   exon_len_mean |
+==============================+=========+====================+===============+========================+=======================+==========================+=========+========================+=================+
| mikado-permissive.loci.stats |    5633 |               2915 |          6054 |                   1.07 |               1642.08 |                     3095 |   11850 |                   1.96 |          838.92 |
+------------------------------+---------+--------------------+---------------+------------------------+-----------------------+--------------------------+---------+------------------------+-----------------+

Mikado has created an annotation that is in between those produced by Stringtie and CLASS2. Compared to CLASS2, the Mikado assembly is probably more comprehensive, given the higher number of genes. Half the genes are monoexonic, an improvement compared to CLASS2, and at the same time the number of genes is much lower than in Stringtie, with a higher average cDNA length. These statistics suggest that the Mikado filtered annotation has been able to retain most of the real genes, while discarding many of the fragments present in either assembly. We can verify this by comparing the Mikado results against the reference annotation::

    mikado compare -r Schizosaccharomyces_pombe.ASM294v2.32.gff3 -p daijin_spombe/5-mikado/pick/permissive/mikado-permissive.loci.gff3 -o daijin_spombe/5-mikado/mikado_prepared.compare -l daijin_spombe/5-mikado/mikado_prepared.compare.log

A cursory look at the STATS file confirms the impression; the Mikado annotation has removed many of the spurious transcripts and is both more sensitive and more precise than either of the input assemblies::

  /tgac/software/testing/bin/core/../..//mikado/devel/x86_64/bin/mikado compare -r Schizosaccharomyces_pombe.ASM294v2.32.gff3 -p daijin_spombe/5-mikado/pick/permissive/mikado-permissive.loci.gff3 -o daijin_spombe/5-mikado/mikado_prepared.compare -l daijin_spombe/5-mikado/mikado_prepared.compare.log
  7015 reference RNAs in 7014 genes
  6054 predicted RNAs in  5633 genes
  --------------------------------- |   Sn |   Pr |   F1 |
                          Base level: 64.04  85.97  73.40
              Exon level (stringent): 20.37  21.96  21.13
                Exon level (lenient): 70.61  81.36  75.61
                        Intron level: 88.21  87.03  87.61
                  Intron chain level: 79.58  68.60  73.68
        Transcript level (stringent): 0.09  0.10  0.09
    Transcript level (>=95% base F1): 24.55  28.51  26.38
    Transcript level (>=80% base F1): 47.04  54.64  50.56
           Gene level (100% base F1): 0.09  0.11  0.09
          Gene level (>=95% base F1): 24.55  30.64  27.26
          Gene level (>=80% base F1): 47.05  58.73  52.24

  #   Matching: in prediction; matched: in reference.

              Matching intron chains: 2030
               Matched intron chains: 2030
     Matching monoexonic transcripts: 1688
      Matched monoexonic transcripts: 1690
          Total matching transcripts: 3718
           Total matched transcripts: 3720

            Missed exons (stringent): 9857/12378  (79.63%)
             Novel exons (stringent): 8960/11481  (78.04%)
              Missed exons (lenient): 2822/9602  (29.39%)
               Novel exons (lenient): 1553/8333  (18.64%)
                      Missed introns: 632/5361  (11.79%)
                       Novel introns: 705/5434  (12.97%)

                  Missed transcripts: 1745/7015  (24.88%)
                   Novel transcripts: 1/6054  (0.02%)
                        Missed genes: 1745/7014  (24.88%)
                         Novel genes: 1/5633  (0.02%)

We suggest to visualise assemblies with one of the many tools currently at disposal, such as eg `WebApollo <http://genomearchitect.org/>`_ [Apollo]_. Mikado files are GFF3-compliant and can be loaded directly into Apollo or similar tools. GTF files can be converted into proper GFF3 files using the convert utility::

  mikado util convert <input gtf> <output GFF3>

