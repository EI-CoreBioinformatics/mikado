.. _ERR588044: http://www.ebi.ac.uk/ena/data/view/ERR588044
.. _ERX1174638: http://www.ebi.ac.uk/ena/data/view/ERX1174638
.. _ERX1174639: http://www.ebi.ac.uk/ena/data/view/ERX1174639
.. _PRJEB11515: http://www.ebi.ac.uk/ena/data/view/PRJEB11515
.. _ERR588044_1.fastq.gz: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR588/ERR588044/ERR588044_1.fastq.gz
.. _ERR588044_2.fastq.gz: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR588/ERR588044/ERR588044_2.fastq.gz
.. _Daijin-Tutorial:

Tutorial for Daijin
===================

This tutorial will guide you through the task of configuring and running the whole Daijin pipeline on a *S. cerevisiae* dataset comprising two different samples, using one aligner (HISAT) and two assemblers (Stringtie and CLASS2) as chosen methods. A modern desktop computer with a multicore processor and 4GB of RAM or more should suffice to execute the pipeline within two hours.


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

.. * BLAST+ (v2.3.0)
.. * TransDecoder (v3.0.0)

* DIAMOND (v0.8.22 or later)
* Prodigal (v2.6.3 or later)
* Portcullis (v1.0.2 or later)
* HISAT2 (v2.0.4)
* Stringtie (v1.2.4)
* CLASS2 (v2.12)
* SAMtools (v1.1 or later)


Input data
~~~~~~~~~~

Throughout this tutorial, we will use data coming from EnsEMBL v89, and from the PRJEB11515_ experiment on ENA. In particular, we will need:

 * the `genome FASTA file <ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz>`_ of *Saccharomyces cerevisiae*
 * its `relative genome annotation <ftp://ftp.ensembl.org/pub/release-89/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.89.gtf.gz>`_
 * RNA-Seq from two samples of the PRJEB11515_ study:
     * ERX1174638_, `left <ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1633_1_1_ATCACG_L006_R1.fastq.gz>`_ and `right <ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1633_1_1_ATCACG_L006_R2.fastq.gz>`_ reads
     * ERX1174639_, `left <ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1637__2_2_ACAGTG_L006_R1.fastq.gz>`_ and `right <ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1637__2_2_ACAGTG_L006_R2.fastq.gz>`_ reads
 * protein sequences for the related species *Pichia pastoris*, `downloaded from Uniprot <http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=organism:%22Pichia%20pastoris%20[644223]%22&fil=&format=fasta&force=yes>`_

Preparation of the input data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First of all, let us set up a folder with the reference data::

  mkdir -p Reference;
  cd Reference;
  wget ftp://ftp.ensembl.org/pub/release-89/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz;
  wget ftp://ftp.ensembl.org/pub/release-89/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.89.gtf.gz;
  wget -O Pichia_pastoris.fasta.gz "http://www.uniprot.org/uniprot/?sort=score&desc=&compress=yes&query=organism:%22Pichia%20pastoris%20[644223]%22&fil=&format=fasta&force=yes"
  gunzip *gz;
  cd ../;

The snippet of the bash script above will create a "Reference" directory, and download the genome of *S. cerevisiae* in FASTA file, the corresponding GTF, and the protein sequences for *Pichia pastoris*. It will also decompress all files.

It is possible to have a feel for the annnotation of this species - the size of its genes and transcripts, the average number of exons per transcript, etc - by using ``mikado util stats``; just issue the following command::

  mikado util stats Reference/Saccharomyces_cerevisiae.R64-1-1.89.gtf Reference/Saccharomyces_cerevisiae.R64-1-1.89.stats.txt

These are the results:

+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Stat                             | Total   | Average  | Mode  | Min    | 1%     | 5%   | 10% | 25% | Median | 75%   | 90%   | 95%   | 99%   | Max    |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Number of genes                  | 7126    | NA       | NA    | NA     | NA     | NA   | NA  | NA  | NA     | NA    | NA    | NA    | NA    | NA     |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Number of genes (coding)         | 6692    | NA       | NA    | NA     | NA     | NA   | NA  | NA  | NA     | NA    | NA    | NA    | NA    | NA     |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Number of monoexonic genes       | 6730    | NA       | NA    | NA     | NA     | NA   | NA  | NA  | NA     | NA    | NA    | NA    | NA    | NA     |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Transcripts per gene             | 7126    | 1.00     | 1     | 1      | 1      | 1    | 1   | 1   | 1      | 1     | 1     | 1     | 1     | 1      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Coding transcripts per gene      | 6692    | 0.94     | 1     | 0      | 0      | 0    | 1   | 1   | 1      | 1     | 1     | 1     | 1     | 1      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDNA lengths                     | 9153986 | 1,284.59 | 72    | 51     | 72     | 96   | 237 | 453 | 1,012  | 1,722 | 2,655 | 3,441 | 5,313 | 14,733 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDNA lengths (mRNAs)             | 9050724 | 1,352.47 | 1323  | 51     | 117    | 237  | 327 | 534 | 1,078  | 1,770 | 2,697 | 3,498 | 5,313 | 14,733 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDS lengths                      | 9031176 | 1,267.36 | 0     | 0      | 0      | 0    | 216 | 441 | 1,004  | 1,715 | 2,642 | 3,424 | 5,310 | 14,730 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDS lengths (mRNAs)              | NA      | 1,349.55 | 1320  | 48     | 114    | 234  | 324 | 531 | 1,076  | 1,767 | 2,694 | 3,495 | 5,310 | 14,730 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDS/cDNA ratio                   | NA      | 99.57    | 100.0 | 94     | 97     | 99   | 99  | 99  | 100    | 100   | 100   | 100   | 100   | 100    |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Monoexonic transcripts           | 6730    | 1,284.75 | 72    | 51     | 72     | 112  | 255 | 474 | 1,035  | 1,731 | 2,631 | 3,342 | 5,052 | 14,733 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| MonoCDS transcripts              | 6363    | 1,344.21 | 1320  | 48     | 111    | 234  | 324 | 549 | 1,095  | 1,770 | 2,664 | 3,374 | 5,192 | 14,730 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Exons per transcript             | 7553    | 1.06     | 1     | 1      | 1      | 1    | 1   | 1   | 1      | 1     | 1     | 2     | 2     | 8      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Exons per transcript (mRNAs)     | 1520    | 1.05     | 1     | 1      | 1      | 1    | 1   | 1   | 1      | 1     | 1     | 1     | 2     | 8      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Exon lengths                     | NA      | 1,211.97 | 72    | 1      | 14     | 72   | 132 | 390 | 948    | 1,668 | 2,594 | 3,364 | 4,990 | 14,733 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Exon lengths (mRNAs)             | NA      | 1,283.79 | 1323  | 1      | 12     | 132  | 265 | 468 | 1,026  | 1,722 | 2,649 | 3,426 | 5,052 | 14,733 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Intron lengths                   | NA      | 271.86   | 1     | 1      | 1      | 1    | 1   | 38  | 100    | 384   | 525   | 1,010 | 2,448 | 2,483  |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Intron lengths (mRNAs)           | NA      | 312.61   | 1     | 1      | 1      | 1    | 1   | 82  | 123    | 401   | 550   | 1,404 | 2,463 | 2,483  |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDS exons per transcript         | 1520    | 0.99     | 1     | 0      | 0      | 0    | 1   | 1   | 1      | 1     | 1     | 1     | 2     | 8      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDS exons per transcript (mRNAs) | 1520    | 1.05     | 1     | 1      | 1      | 1    | 1   | 1   | 1      | 1     | 1     | 1     | 2     | 8      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDS exon lengths                 | 9031176 | 1,281.02 | 1320  | 1      | 12     | 131  | 264 | 465 | 1,023  | 1,719 | 2,646 | 3,423 | 5,049 | 14,730 |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| CDS Intron lengths               | 111558  | 311.61   | 0     | 0      | 0      | 0    | 0   | 81  | 122    | 400   | 549   | 1,403 | 2,462 | 2,482  |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| 5'UTR exon number                | 6692    | 0.00     | 0     | 0      | 0      | 0    | 0   | 0   | 0      | 0     | 0     | 0     | 0     | 0      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| 3'UTR exon number                | 6692    | 0.97     | 1     | 0      | 0      | 1    | 1   | 1   | 1      | 1     | 1     | 1     | 1     | 1      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| 5'UTR length                     | 0       | 0.00     | 0     | 0      | 0      | 0    | 0   | 0   | 0      | 0     | 0     | 0     | 0     | 0      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| 3'UTR length                     | 19548   | 2.92     | 3     | 0      | 0      | 3    | 3   | 3   | 3      | 3     | 3     | 3     | 3     | 3      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Stop distance from junction      | NA      | 0.00     | 0     | 0      | 0      | 0    | 0   | 0   | 0      | 0     | 0     | 0     | 0     | 0      |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Intergenic distances             | NA      | 403.20   | -1322 | -9,349 | -1,322 | -292 | -7  | 180 | 315    | 548   | 955   | 1,367 | 2,693 | 9,140  |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+
| Intergenic distances (coding)    | NA      | 444.63   | -1322 | -9,349 | -1,322 | -293 | -5  | 186 | 324    | 581   | 1,028 | 1,526 | 3,279 | 9,140  |
+----------------------------------+---------+----------+-------+--------+--------+------+-----+-----+--------+-------+-------+-------+-------+--------+

From this summary it is quite apparent that the *S. cerevisiae* genome preferentially encodes monoexonic transcripts, which are almost always completely coding. Intron lengths are generally short - 272 bps on average and 100 as median, with a maximum value of ~2,500 bps for mRNAs. It is also an extremely compact genome, with a median intergenic distance of 324 bps and a negative modal distance - indicating that most often genes are overlapping.

Next, we download the reads that we will use for this example::

  mkdir -p Reads;
  cd Reads;
  wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1633_1_1_ATCACG_L006_R1.fastq.gz;
  wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1633_1_1_ATCACG_L006_R2.fastq.gz;
  wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1637__2_2_ACAGTG_L006_R1.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA526/ERA526598/fastq/1637__2_2_ACAGTG_L006_R2.fastq.gz;
  cd ../;

These files have a total file size of approximately 4GB, so they might take five to ten minutes to download, depending on your connection speed.

Step 1: configuring Daijin
~~~~~~~~~~~~~~~~~~~~~~~~~~

The first task is to create a configuration file for Daijin using ``daijin configure``. On the command line, we have to configure the following:

* name of the configuration file (daijin.yaml)
* number of threads per process (2)
* reference genome and the name of the species (**without spaces or non-ASCII/special characters**)
* reads to use, with their strandedness and the name of the sample
* aligners (HISAT2) and assemblers (Stringtie, CLASS2)
* output directory (Scerevisiae)
* the scoring file to use for Mikado pick; we will ask to copy it in-place to have a look at it (scerevisiae.yaml)
* the protein database for homology searches for Mikado (Pichia_pastoris.fasta)
* flank: as *A. thaliana* has a quite compact genome, we should decrease the maximum distance for grouping together transcripts. We will decrease from 1kbps (default) to 200.
* (optional) the scheduler to use for the cluster (we will presume that the job is being executed locally)
* (optional) name of the cluster configuration file, which will have to be edited manually.

First, we will create a sample sheet, containing the information of the sample that we are going to use. This is a tab-delimited text file, where each line defines a single sample, and with up to 5 columns per line. The first three columns are mandatory, while the last two are optional. The columns are as follows:

 * **Read1**: required. Location of the left reads for the sample.
 * **Read2**: optional, location of the right reads for the sample if it is paired.
 * **Sample**: name of the sample. Required.
 * **Strandedness**: strandedness of the sample. It can be one of:
    * fr-unstranded (Unstranded data)
    * fr-firststrand (Stranded data, first read forward, second read reverse)
    * fr-secondstrand (Stranded data, second read forward, first read reverse)
    * f (Forward, single read only)
    * r (Reverse, single read only)
 * **Long read sample**: Boolean flag. If set to "True", the sample will be considered as coming from a non-second generation sequencing platform (eg. Sanger ESTs or PacBio IsoSeq) and for its reads we would therefore consider only the alignment, without performing any assembly.

For our example, therefore, the sample sheet will look like this::

  Reads/1633_1_1_ATCACG_L006_R1.fastq.gz	Reads/1633_1_1_ATCACG_L006_R1.fastq.gz	ERX1174638	fr-unstranded
  Reads/1637__2_2_ACAGTG_L006_R1.fastq.gz	Reads/1637__2_2_ACAGTG_L006_R2.fastq.gz	ERX1174639	fr-unstranded

Write this into a text file called "sample_sheet.tsv".
Now we will configure Daijin for the run::

  daijin configure \
       --scheduler "" \
       --scoring scerevisiae.yaml \
       --copy-scoring scerevisiae.yaml \
       -m permissive \
       --sample-sheet reads.tsv \
       --flank 200 -i 10 2500 --threads 2 \
       --genome Reference/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
       -al hisat -as class stringtie -od Scerevisiae --name Scerevisiae \
       -o daijin.yaml  --prot-db Reference/Pichia_pastoris.fasta;

This will create three files in the working directory:

  * **daijin.yaml**, the main configuration file. This file is in `YAML format <http://www.yaml.org/spec/1.2/spec.html>`_.
  * *scerevisiae.yaml*: this file is a copy of the scoring configuration file. Please refer to the :ref:`dedicated section <scoring_files>` for details.
  * *daijin_exe.yaml*: this small configuration file contains the instruction to load the necessary software into the working environment. Ignore it if you are working on a local machine. If you are working within a cluster environment, please modify this file with the normal commands you would use in a cluster script to load necessary software. For example, this is how this file looks like on our own cluster system:

.. code-block:: yaml

    blast: 'source blast-2.3.0'
    class: 'source class-2.12'
    cufflinks: ''
    gmap: ''
    hisat: 'source HISAT-2.0.4'
    mikado: 'source mikado-1.1'
    portcullis: 'source portcullis-0.17.2'
    samtools: 'source samtools-1.2'
    star: ''
    stringtie: 'source stringtie-1.2.4'
    tophat: ''
    transdecoder: 'source transdecoder-3.0.0'
    trinity: ''
..

.. important::

  If you are operating on a cluster, instead of a local machine, you will need to specify the scheduler type. Currently we support **SLURM, PBS and LSF**. Add the following switch to the configure command above::

    --scheduler <One of SLURM, PBS or LSF>

  Adding this switch will also create a default *cluster configuration file*, specifying the number of resources per job and the submission queue. This is an example of how it appears on our system:
.. literalinclude:: ../Usage/hpc.yaml
..

Step 2: running the assemble part
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created a proper configuration file, it is time to launch Daijin assemble and inspect the results. Issue the command::

    daijin assemble --cores <Number of maximum cores> daijin.yaml


After checking that the configuration file is valid, Daijin will start the alignment and assembly of the dataset. On a normal desktop computer, this should take less than 2 hours. Before launching the pipeline, you can obtain a graphical representation of the steps with::

  daijin assemble --dag daijin.yaml | dot -Tsvg > assemble.svg

..
.. image:: assemble_pipeline.png

You can also ask Daijin to display the steps to be executed, inclusive of their command lines, by issuing the following command::

  daijin assemble --dryrun daijin.yaml

When Daijin is finished, have a look inside the folder Scerevisiae/3-assemblies/output/; you will find the following GTF files:

* Scerevisiae/3-assemblies/output/class-0-hisat-ERX1174638-0.gtf
* Scerevisiae/3-assemblies/output/class-0-hisat-ERX1174639-0.gtf
* Scerevisiae/3-assemblies/output/stringtie-0-hisat-ERX1174638-0.gtf
* Scerevisiae/3-assemblies/output/stringtie-0-hisat-ERX1174639-0.gtf

These are standard `GTF files <http://www.ensembl.org/info/website/upload/gff.html>`_ reporting the assembled transcripts for each method. We can have a feel for how they compare with our reference annotation by, again, using ``mikado util stats``. Conveniently, Daijin has already performed this analysis for us, and the files will be present in the same folder:

* Scerevisiae/3-assemblies/output/class-0-hisat-ERX1174638-0.gtf.stats
* Scerevisiae/3-assemblies/output/class-0-hisat-ERX1174639-0.gtf.stats
* Scerevisiae/3-assemblies/output/stringtie-0-hisat-ERX1174638-0.gtf.stats
* Scerevisiae/3-assemblies/output/stringtie-0-hisat-ERX1174639-0.gtf.stats

Daijin has also created a summary of these statistics in Scerevisiae/3-assemblies/assembly.stats:

+------------------------------------------+--------+------------------+-------------+----------------------+---------------------+------------------------+--------+----------------------+---------------+
| File                                     | genes  | monoexonic_genes | transcripts | transcripts_per_gene | transcript_len_mean | monoexonic_transcripts | exons  | exons_per_transcript | exon_len_mean |
+------------------------------------------+--------+------------------+-------------+----------------------+---------------------+------------------------+--------+----------------------+---------------+
| class-0-hisat-ERX1174638-0.gtf.stats     | 123.0  | 122.0            | 123.0       | 1.0                  | 1483.59             | 122.0                  | 124.0  | 1.01                 | 1471.63       |
+------------------------------------------+--------+------------------+-------------+----------------------+---------------------+------------------------+--------+----------------------+---------------+
| class-0-hisat-ERX1174639-0.gtf.stats     | 482.0  | 121.0            | 507.0       | 1.05                 | 1563.39             | 121.0                  | 911.0  | 1.8                  | 870.08        |
+------------------------------------------+--------+------------------+-------------+----------------------+---------------------+------------------------+--------+----------------------+---------------+
| stringtie-0-hisat-ERX1174638-0.gtf.stats | 3385.0 | 3090.0           | 3505.0      | 1.03                 | 3028.65             | 3152.0                 | 3881.0 | 1.11                 | 2735.23       |
+------------------------------------------+--------+------------------+-------------+----------------------+---------------------+------------------------+--------+----------------------+---------------+
| stringtie-0-hisat-ERX1174639-0.gtf.stats | 2627.0 | 2336.0           | 2773.0      | 1.05                 | 3973.53             | 2424.0                 | 3145.0 | 1.13                 | 3503.53       |
+------------------------------------------+--------+------------------+-------------+----------------------+---------------------+------------------------+--------+----------------------+---------------+

From this quick analysis, it looks like CLASS under-assembled the datasets, with a very low number of genes per prediction. Stringtie seems to have performed better, but the length of the predicted transcripts is much higher than the average length of annotated transcripts in this species - suggesting that the tool might have concatenated multiple transcripts in long read-through events. We can obtain a clearer picture of the situation by using :ref:`Mikado compare <compare>`. To do so, issue the following commands:

.. code-block:: bash

  # First index the reference GFF3
  mikado compare -r Reference/Saccharomyces_cerevisiae.R64-1-1.89.gtf --index --log Reference/index.log;
  mkdir -p Comparisons;
  # Compare the CLASS2 assemblies against the reference
  mikado compare -r Reference/Saccharomyces_cerevisiae.R64-1-1.89.gtf -p Scerevisiae/3-assemblies/output/class-0-hisat-ERX1174639-0.gtf -o Comparisons/class_ERX1174639 -l Comparisons/class_ERX1174639.log;
  mikado compare -r Reference/Saccharomyces_cerevisiae.R64-1-1.89.gtf -p Scerevisiae/3-assemblies/output/class-0-hisat-ERX1174638-0.gtf -o Comparisons/class_ERX1174638 -l Comparisons/class_ERX1174638.log;
  # Compare the StringTie assemblies against the reference:
  mikado compare -r Reference/Saccharomyces_cerevisiae.R64-1-1.89.gtf -p Scerevisiae/3-assemblies/output/stringtie-0-hisat-ERX1174638-0.gtf -o Comparisons/stringtie_ERX1174638 -l Comparisons/stringtie_ERX1174638.log;
  mikado compare -r Reference/Saccharomyces_cerevisiae.R64-1-1.89.gtf -p Scerevisiae/3-assemblies/output/stringtie-0-hisat-ERX1174639-0.gtf -o Comparisons/stringtie_ERX1174639 -l Comparisons/stringtie_ERX1174639.log;

The analysis will produce *TMAP*, *REFMAP* and *STATS* files for each of the assemblies. As an example, this is the statistics file for the ERX1174639 StringTie assembly::

    Command line:
    /usr/users/ga002/venturil/miniconda3/envs/py360/bin/mikado compare -r Reference/Saccharomyces_cerevisiae.R64-1-1.89.gtf -p Scerevisiae/3-assemblies/output/stringtie-0-hisat-ERX1174639-0.gtf -o Comparisons/stringtie_ERX1174639 -l Comparisons/stringtie_ERX1174639.log
    7126 reference RNAs in 7126 genes
    2773 predicted RNAs in  2627 genes
    --------------------------------- |   Sn |   Pr |   F1 |
                            Base level: 32.73  27.44  29.85
                Exon level (stringent): 0.12  0.29  0.17
                  Exon level (lenient): 55.89  68.89  61.71
                          Intron level: 58.39  67.99  62.83
                    Intron chain level: 56.57  64.74  60.38
          Transcript level (stringent): 0.00  0.00  0.00
      Transcript level (>=95% base F1): 2.55  6.56  3.68
      Transcript level (>=80% base F1): 11.69  30.04  16.83
             Gene level (100% base F1): 0.00  0.00  0.00
            Gene level (>=95% base F1): 2.55  6.93  3.73
            Gene level (>=80% base F1): 11.69  31.40  17.04

    #   Matching: in prediction; matched: in reference.

                Matching intron chains: 226
                 Matched intron chains: 224
       Matching monoexonic transcripts: 774
        Matched monoexonic transcripts: 774
            Total matching transcripts: 1000
             Total matched transcripts: 998

              Missed exons (stringent): 7522/7531  (99.88%)
               Novel exons (stringent): 3098/3107  (99.71%)
                Missed exons (lenient): 367/832  (44.11%)
                 Novel exons (lenient): 210/675  (31.11%)
                        Missed introns: 171/411  (41.61%)
                         Novel introns: 113/353  (32.01%)

                    Missed transcripts: 987/7126  (13.85%)
                     Novel transcripts: 92/2773  (3.32%)
                          Missed genes: 987/7126  (13.85%)
                           Novel genes: 89/2627  (3.39%)


Step 3: running the Mikado steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have created the input assemblies, it is time to run Mikado. First of all, let us have a look at the *scerevisiae.yaml* file, which we have conveniently copied to our current directory:

.. literalinclude:: scerevisiae.yaml

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

