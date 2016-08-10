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

Additionally, you will require the following software:

* a functioning installation of SQLite.
* a functioning version of BLAST+.

.. important:: all the available files are at disposal in the **sample_data** directory of the `Mikado source <http://www.github.com/lucventurini/Mikado>`. The steps in this tutorial are codified in a Snakefile and can be executed automatically using Snakemake [Snake]_.

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

    mikado configure --list list.txt --reference chr5.fas --mode permissive  --copy-scoring plants.yaml --scoring plants.yaml --junctions junctions.bed -bt uniprot_sprot_plants.fasta configuration.yaml

This will create a :download:`configuration.yaml` file with the parameters that were specified on the command line. This is :ref:`simplified configuration file <conf_anatomy>`, containing all the necessary parameters for the Mikado run. It will also copy the ``plants.yaml`` file from the Mikado installation to your current working directory.

.. hint:: The *--copy-scoring* argument is usually not necessary, however, it allows you to easily inspect the :ref:`scoring file <scoring-files>` we are going to use  during this run.

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

At the end of this phase, you should have 93 candidate transcripts.

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

Mikado pick
~~~~~~~~~~~




.. output of the Snakefile:
    Creating the configuration file
    mikado configure --list list.txt --reference chr5.fas --mode permissive         --scoring plants.yaml --junctions junctions.bed -bt uniprot_sprot_plants.fasta configuration.yaml
    Provided cores: 1
    Rules claiming more threads will be scaled down.
    Job counts:
        count	jobs
        1	blast_complete
        1	compare
        1	compare_input
        1	compare_subloci
        1	complete
        1	pick
        1	prepare
        1	prepare_blast
        1	serialise
        1	test_json
        1	uncompress_blast
        11
    gunzip -c chr5.fas.gz > chr5.fas
    Touching output file configuration.yaml.ok.
    1 of 11 steps (9%) done
    mikado prepare --json-conf configuration.yaml
    2 of 11 steps (18%) done
    gzip -dc uniprot_sprot_plants.fasta.gz > uniprot_sprot_plants.fasta
    3 of 11 steps (27%) done
    makeblastdb -in uniprot_sprot_plants.fasta -dbtype prot -parse_seqids > blast_prepare.log
    4 of 11 steps (36%) done
    blastx -max_target_seqs 5 -num_threads 10 -query mikado_prepared.fasta -outfmt 5 -db uniprot_sprot_plants.fasta -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz
    5 of 11 steps (45%) done
    mikado serialise -p 1 --json-conf configuration.yaml --xml mikado.blast.xml.gz         --orfs mikado.bed --blast_targets uniprot_sprot_plants.fasta --force
    Touching output file serialised.ok.
    6 of 11 steps (55%) done
    mikado pick --json-conf configuration.yaml -lv INFO --subloci_out mikado.subloci.gff3 -p 1
    7 of 11 steps (64%) done
    mikado compare -r reference.gff3 -p mikado.subloci.gff3 -o compare_subloci -l compare_subloci.log
    8 of 11 steps (73%) done
    mikado compare -r reference.gff3 -p mikado_prepared.gtf -o compare_input -l compare_input.log
    9 of 11 steps (82%) done
    mikado compare -r reference.gff3 -p mikado.loci.gff3 -o compare -l compare.log
    10 of 11 steps (91%) done
    localrule complete:
        input: compare.stats, compare_subloci.stats, compare_input.stats
        output: finished.ok
    Touching output file finished.ok.
    11 of 11 steps (100%) done

