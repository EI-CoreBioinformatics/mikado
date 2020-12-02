.. _Tutorial:

Tutorial
========

This tutorial will guide you through a simple analysis of Mikado, using a small amount of data coming from an experiment on *Arabidopsis thaliana*. RNA-Seq data was obtained from `study PRJEB7093 on ENA <http://www.ebi.ac.uk/ena/data/view/PRJEB7093>`_, aligned with STAR [STAR]_ against the `TAIR10 <http://www.arabidopsis.org>`_ reference genome, and assembled with four different programs. For this small example, we are going to focus on a small genomic region: Chr5, from 26,575,364 to 26,614,205.

During this analysis, you will require the following files:

* :download:`chr5.fas.gz <chr5.fas.gz>`: a FASTA file containing the Chr5 of *A. thaliana*.
* :download:`reference.gff3`: a GFF3 file with the annotation of the genomic slice we are interested in, for comparison purposes.
* :download:`junctions.bed <junctions.bed>`: a BED12 file of reliable splicing junctions in the region, identified using Portcullis [Portcullis]_
* :download:`class.gtf`: a GTF file of transcripts assembled using CLASS [Class2]_
* :download:`cufflinks.gtf`: a GTF file of transcripts assembled using Cufflinks [Cufflinks]_
* :download:`stringtie.gtf`: a GTF file of transcripts assembled using Stringtie [StringTie]_
* :download:`trinity.gff3`: a GFF3 file of transcripts assembled using Trinity [Trinity]_ and aligned using GMAP [GMAP]_
* :download:`orfs.bed`: a BED12 file containing the ORFs of the above transcripts, derived using TransDecoder [Trinity]_
* :download:`uniprot_sprot_plants.fasta.gz`: a FASTA file containing the plant proteins released with SwissProt [Uniprot]_

All of this data can also be found in the ``sample_data`` directory of the `Mikado source <http://www.github.com/EI-CoreBioinformatics/Mikado>`_.

You will also require the following software:

* a functioning installation of SQLite.
* a functioning version of BLAST+ [Blastplus]_.

Creating the configuration file for Mikado
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the first step, we need to create a configuration file to drive Mikado. To do so, we will first create a tab-delimited file describing our assemblies (class.gtf, cufflinks.gtf, stringtie.gtf, trinity.gff3)::

    class.gtf	cl	True
    cufflinks.gtf	cuff	True
    stringtie.gtf	st	True    1
    trinity.gff3	tr	False   -0.5

In this file, the three fields define the following:

#. The file location and name (if no folder is specified, Mikado will look for each file in the current working directory)
#. An alias associated with the file, which has to be unique
#. A binary flag (``True`` / ``False``) indicating whether the assembly is strand-specific or not
#. An optional fourth field which defines a score associated with that sample. All transcripts associated with the label will have their score corrected by the value on this field. So eg. in this example all Stringtie models will receive an additional point, and all Trinity models will be penalised by half a point. Class and Cufflinks have no special bonus or malus associated with them.

Then, we will decompress the chromosome FASTA file:

.. code-block:: bash

    gunzip chr5.fas.gz  # This will create the file chr5.fas

Finally, we will create the configuration file itself using ``mikado configure``:

.. code-block:: bash

    mikado configure --list list.txt --reference chr5.fas --mode permissive --scoring plants.yaml  --copy-scoring plants.yaml --junctions junctions.bed -bt uniprot_sprot_plants.fasta configuration.yaml

This will create a configuration.yaml file with the parameters that were specified on the command line. This is :ref:`simplified configuration file <conf_anatomy>`, containing all the necessary parameters for the Mikado run. It will also copy the ``plants.yaml`` file from the Mikado installation to your current working directory.
The parameters we used for the command line instruct Mikado in the following ways:

Alternatively, we could use a BGZIP-compressed file as input for ``mikado configure``:

.. code-block:: bash
    bgzip chr5.fas
    samtools faidx chr5.fas.gz
    mikado configure --list list.txt --reference chr5.fas.gz --mode permissive --scoring plants.yaml  --copy-scoring plants.yaml --junctions junctions.bed -bt uniprot_sprot_plants.fasta configuration.yaml

* *--list list.txt*: this part of the command line instructs Mikado to read the file we just created to understand where the input files are and how to treat them.
* *--reference chr5.fas*: this part of the command line instructs Mikado on the location of the genome file.
* *--mode permissive*: the mode in which Mikado will treat cases of chimeras. See the :ref:`documentation <chimera_splitting_algorithm>` for details.
* *--junctions junctions.bed*: this part of the command line instructs Mikado to consider this file as the source of reliable splicing junctions.
* *-bt uniprot_sprot_plants.fasta*: this part of the command line instructs Mikado to consider this file as the BLAST database which will be used for deriving homology information.

.. hint:: The *--copy-scoring* argument is usually not necessary, however, it allows you to easily inspect the :ref:`scoring file <scoring_files>` we are going to use  during this run.

.. hint:: Mikado provides a handful of pre-configured scoring files for different species. However, we do recommend inspecting and tweaking your scoring file to cater to your species. We provide a guide on how to create your own configuration files :ref:`here <configure-scoring-tutorial>`.

Mikado prepare
~~~~~~~~~~~~~~

The subsequent step involves running ``mikado prepare`` to create a :ref:`sorted, non-redundant GTF with all the input assemblies <prepare>`. As we have already created a configuration file with all the details regarding the input files, this will require us only to issue the command:

.. code-block:: bash

    mikado prepare --json-conf configuration.yaml

This command will create three files:

#. *mikado_prepared.gtf*: one of the two main output files. This is a sorted, non-redundant GTF containing the transcripts from the four input GTFs
#. *mikado_prepared.fasta*: a FASTA file of the transcripts present in *mikado_prepared.gtf*.
#. *prepare.log*: the log of this step. This should look like the following, minus the timestamps::

    2016-08-10 13:53:58,443 - prepare - prepare.py:67 - INFO - setup - MainProcess - Command line: /usr/users/ga002/venturil/py351/bin/mikado prepare --json-conf configuration.yaml
    2016-08-10 13:53:58,967 - prepare - prepare.py:206 - INFO - perform_check - MainProcess - Finished to analyse 95 transcripts (93 retained)
    2016-08-10 13:53:58,967 - prepare - prepare.py:405 - INFO - prepare - MainProcess - Finished

At the end of this phase, you should have 93 candidate transcripts, as 2 were redundant.

BLAST of the candidate transcripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although it is not strictly necessary, Mikado benefits from integrating homology data from BLAST. Mikado requires this data to be provided either in XML or ASN format (in the latter case, ``blast_formatter`` will be used to convert it in-memory to XML).

To create this file, we will proceed as follows:

#. Uncompress the SwissProt database:

    .. code-block:: bash

        gzip -dc uniprot_sprot_plants.fasta.gz > uniprot_sprot_plants.fasta

#. Prepare the database for the BLAST:

    .. code-block:: bash

        makeblastdb -in uniprot_sprot_plants.fasta -dbtype prot -parse_seqids > blast_prepare.log

#. Execute the BLAST, asking for XML output, and compress it to limit space usage.

    .. code-block:: bash

        blastx -max_target_seqs 5 -num_threads 10 -query mikado_prepared.fasta -outfmt 5 -db uniprot_sprot_plants.fasta -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz

This will produce the ``mikado.blast.xml.gz`` file, which contains the homology information for the run.

Mikado serialise
~~~~~~~~~~~~~~~~

This step involves running ``mikado serialise`` to create a SQLite database with all the information that mikado needs to perform its analysis. As most of the parameters are already specified inside the configuration file, the command line is quite simple:

.. code-block:: bash

    mikado serialise --json-conf configuration.yaml --xml mikado.blast.xml.gz --orfs mikado.bed --blast_targets

After mikado serialise has run, it will have created two files:

#. ``mikado.db``, the SQLite database that will be used by ``pick`` to perform its analysis.
#. ``serialise.log``, the log of the run.

If you inspect the SQLite database ``mikado.db``, you will see it contains six different tables::

    $ sqlite3 mikado.db
    SQLite version 3.8.2 2013-12-06 14:53:30
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .tables
    chrom      hit        hsp        junctions  orf        query      target

These tables contain the information coming, in order, from the genome FAI, the BLAST XML, the junctions BED file, the ORFs BED file, and finally the input transcripts and the proteins. For more details on the database structure, please refer to the section on :ref:`this step <serialise>` in this documentation.

Mikado pick
~~~~~~~~~~~

Finally, during this step ``mikado pick`` will integrate the data present in the database with the positional and structural data present in the GTF file :ref:`to select the best transcript models <pick>`. The command line to be issued is the following:

.. code-block:: bash

    mikado pick --json-conf configuration.yaml --subloci_out mikado.subloci.gff3

At this step, we have to specify only some parameters for ``pick`` to function:

* *json-conf*: the configuration file. This is the only compulsory option.
* *subloci_out*: the partial results concerning the *subloci* step during the selection process will be written to ``mikado.subloci.gff3``.

``mikado pick`` will produce the following output files:

* ``mikado.loci.gff3``, ``mikado.loci.metrics.tsv``, ``mikado.loci.scores.tsv``: the proper output files. These contain the location of the selected transcripts, their metrics, and their scores. Please see :ref:`this section for details <pick-output>`.
* ``mikado.subloci.gff3``, ``mikado.subloci.metrics.tsv``, ``mikado.subloci.scores.tsv``: these files contain the same type of information as those above, but for the *subloci* stage. As such, all the transcripts in the input files are represented, not just those that are going to be selected as the best.
* *mikado_pick.log*: the log file for this operation.

Comparing files with the reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, we can compare our files to the original reference annotation, and see how our results are compared to those. To do so, we will use :ref:`Mikado compare <compare>`.
The first step is to index the reference annotation to make the comparisons faster:

.. code-block:: bash

    mikado compare -r reference.gff3 --index

This will create a new file, ``reference.gff3.midx``. If you inspect with eg. ``zless``, you will notice it is a SQLite database, describing the locations and components of each gene on the annotation. Now that we have indexed the reference, we can perform the comparisons we are interested in:

#. Reference vs. the input transcripts:

.. code-block:: bash

    mikado compare -r reference.gff3 -p mikado_prepared.gtf -o compare_input -l compare_input.log;

#. Reference vs. the subloci stage:

.. code-block:: bash

    mikado compare -r reference.gff3 -p mikado.subloci.gff3 -o compare_subloci -l compare_subloci.log;

#. Reference vs the final output:

.. code-block:: bash

    mikado compare -r reference.gff3 -p mikado.loci.gff3 -o compare -l compare.log;

Each of these comparisons will produce three files:

* a *tmap* file, detailing the best match in the reference for each of the query transcripts;
* a *refmap* file, detailing the best match among the query transcripts for each of the reference transcripts;
* a *stats* file, summarising the comparisons.

The *stats* file for the input GTF should look like this::

    Command line:
    /usr/local/bin/mikado compare -r reference.gff3 -p mikado_prepared.gtf -o compare_input -l compare_input.log
    7 reference RNAs in 5 genes
    93 predicted RNAs in  64 genes
    --------------------------------- |   Sn |   Pr |   F1 |
                            Base level: 95.97  29.39  45.00
                Exon level (stringent): 68.09  18.60  29.22
                  Exon level (lenient): 90.91  31.25  46.51
                          Intron level: 94.74  45.57  61.54
                    Intron chain level: 16.67  1.59  2.90
          Transcript level (stringent): 0.00  0.00  0.00
      Transcript level (>=95% base F1): 14.29  1.08  2.00
      Transcript level (>=80% base F1): 14.29  1.08  2.00
             Gene level (100% base F1): 0.00  0.00  0.00
            Gene level (>=95% base F1): 20.00  1.56  2.90
            Gene level (>=80% base F1): 20.00  1.56  2.90

    #   Matching: in prediction; matched: in reference.

                Matching intron chains: 1
                 Matched intron chains: 1
       Matching monoexonic transcripts: 0
        Matched monoexonic transcripts: 0
            Total matching transcripts: 1
             Total matched transcripts: 1

              Missed exons (stringent): 15/47  (31.91%)
               Novel exons (stringent): 140/172  (81.40%)
                Missed exons (lenient): 4/44  (9.09%)
                 Novel exons (lenient): 88/128  (68.75%)
                        Missed introns: 2/38  (5.26%)
                         Novel introns: 43/79  (54.43%)

                    Missed transcripts: 0/7  (0.00%)
                     Novel transcripts: 24/93  (25.81%)
                          Missed genes: 0/5  (0.00%)
                           Novel genes: 21/64  (32.81%)

For the *subloci* file, where we still have all the transcripts but we have split obvious chimeras, it should look like this::

    Command line:
    /usr/local/bin/mikado compare -r reference.gff3 -p mikado.subloci.gff3 -o compare_subloci -l compare_subloci.log
    7 reference RNAs in 5 genes
    105 predicted RNAs in  26 genes
    --------------------------------- |   Sn |   Pr |   F1 |
                            Base level: 95.96  29.24  44.83
                Exon level (stringent): 70.21  19.08  30.00
                  Exon level (lenient): 88.89  32.00  47.06
                          Intron level: 94.74  46.75  62.61
                    Intron chain level: 33.33  3.17  5.80
          Transcript level (stringent): 0.00  0.00  0.00
      Transcript level (>=95% base F1): 28.57  9.52  14.29
      Transcript level (>=80% base F1): 42.86  11.43  18.05
             Gene level (100% base F1): 0.00  0.00  0.00
            Gene level (>=95% base F1): 40.00  7.69  12.90
            Gene level (>=80% base F1): 60.00  11.54  19.35

    #   Matching: in prediction; matched: in reference.

                Matching intron chains: 3
                 Matched intron chains: 2
       Matching monoexonic transcripts: 9
        Matched monoexonic transcripts: 1
            Total matching transcripts: 12
             Total matched transcripts: 3

              Missed exons (stringent): 14/47  (29.79%)
               Novel exons (stringent): 140/173  (80.92%)
                Missed exons (lenient): 5/45  (11.11%)
                 Novel exons (lenient): 85/125  (68.00%)
                        Missed introns: 2/38  (5.26%)
                         Novel introns: 41/77  (53.25%)

                    Missed transcripts: 0/7  (0.00%)
                     Novel transcripts: 24/105  (22.86%)
                          Missed genes: 0/5  (0.00%)
                           Novel genes: 13/26  (50.00%)

A marked improvement can already be seen - we have now 105 transcripts instead of 93, and the total number of matching transcripts has increased from 1 to 3. Precision is still poor, however, as we have not discarded any transcript yet. Moreover, we have redundancy - 9 transcripts match the same monoexonic gene, and 3 transcripts match 2 intron chains in the reference.
Finally, the comparison against the proper output (``mikado.loci.gff3``) should look like this::

    Command line:
    /usr/local/bin/mikado compare -r reference.gff3 -p mikado.loci.gff3 -o compare -l compare.log
    7 reference RNAs in 5 genes
    15 predicted RNAs in  8 genes
    --------------------------------- |   Sn |   Pr |   F1 |
                            Base level: 85.74  64.73  73.77
                Exon level (stringent): 63.83  42.86  51.28
                  Exon level (lenient): 80.00  52.94  63.72
                          Intron level: 89.47  59.65  71.58
                    Intron chain level: 33.33  14.29  20.00
          Transcript level (stringent): 0.00  0.00  0.00
      Transcript level (>=95% base F1): 28.57  13.33  18.18
      Transcript level (>=80% base F1): 42.86  20.00  27.27
             Gene level (100% base F1): 0.00  0.00  0.00
            Gene level (>=95% base F1): 40.00  25.00  30.77
            Gene level (>=80% base F1): 60.00  37.50  46.15

    #   Matching: in prediction; matched: in reference.

                Matching intron chains: 2
                 Matched intron chains: 2
       Matching monoexonic transcripts: 1
        Matched monoexonic transcripts: 1
            Total matching transcripts: 3
             Total matched transcripts: 3

              Missed exons (stringent): 17/47  (36.17%)
               Novel exons (stringent): 40/70  (57.14%)
                Missed exons (lenient): 9/45  (20.00%)
                 Novel exons (lenient): 32/68  (47.06%)
                        Missed introns: 4/38  (10.53%)
                         Novel introns: 23/57  (40.35%)

                    Missed transcripts: 0/7  (0.00%)
                     Novel transcripts: 6/15  (40.00%)
                          Missed genes: 0/5  (0.00%)
                           Novel genes: 2/8  (25.00%)


After selecting the best transcripts in each locus, Mikado has discarded most of the incorrect transcripts while retaining most of the correct information; this can be seen in the increase in precision at eg. the nucleotide level (from 30% to 65%). The number of genes has also decreased, as Mikado has discarded many loci whose transcripts are just UTR fragments of neighbouring correct genes.

Analysing the tutorial data with Snakemake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The workflow described in this tutorial can be executed automatically using Snakemake [Snake]_ with :download:`this Snakefile <Snakefile>`. Just execute:

.. code-block:: bash

    snakemake

in the directory where you have downloaded all of the tutorial files. In graph representation, this is how the pipeline looks like:

   .. figure:: snakemake_dag.svg
        :align: center
        :scale: 50%
        :figwidth: 100%


