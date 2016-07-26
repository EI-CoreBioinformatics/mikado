.. _Snakemake: https://bitbucket.org/snakemake/snakemake/wiki/Home
.. _YAML: http://www.yaml.org/spec/1.2/spec.html
.. _TransDecoder: https://github.com/TransDecoder/TransDecoder

.. _assemble_pipeline:

.. |snake_badge| image:: https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)]
   :target: http://snakemake.bitbucket.org

.. _Daijin:

The Daijin pipeline for driving Mikado
======================================

|snake_badge|

No emperor or empress can lead its nation without a trusty chancellor to help him or her in organising the bureaucracy. Daijin, the Japanese minister, has the same role in Mikado - it smooths the path to go from a collection of read inputs (both RNA-Seq or long reads) to a polished transcriptome assembly. The pipeline is based on Snakemake_ [Snake]_; while Snakemake can support any scheduling system, our pipeline manager currently supports only three (SLURM, PBS and LSF), plus any DRMAA-compliant batch submission system. Other schedulers can be added upon request.

.. hint:: It is possible to launch the two steps of the pipeline directly with Snakemake, using the snakefiles located in Mikado.daijin: :download:`tr.snakefile <tr.snakefile>` for the first step, and :download:`mikado.snakefile` for the second.


.. _daijin-configure:

Configure
~~~~~~~~~

This utility creates the configuration file that will drive Daijin, in YAML_ format. The file will need to be edited to

.. _daijin-assemble:

Assemble
~~~~~~~~

In the first step of the pipeline, Daijin will perform the following operations for each of the read datasets provided:

#. Create the necessary indices for each of the aligner programs requested.
#. Align the read dataset using all the different tools requested, in all the possible combinations of parameters requested.
   * For example, it is possible to ask each dataset to be aligned twice with TopHat2 - once with the "micro-exon" mode activated, the second time without. Both alignments will be run independently.
   * It is possible to specify which datasets are strand-specific and which are not, and moreover, it is possible to specify the kind of strand-specificity (fr-secondstrand, fr-firststrand).
#. Call all the reliable junctions across the alignments using `PortCullis <https://github.com/maplesond/portcullis>`_
#. Create the statistics for the assembly using ``samtools stat``, and merge them together in a single file.
#. Assemble each alignment with all the tools requested, in all the parameter combinations desired.
#. Call the statistics on each assembly using :ref:`mikado util stats <stat-command>`, and merge them together in a single file.
#. Create the configuration file for Mikado.

So during this first step Daijin will go from raw reads files to multiple assemblies, and configure Mikado for the second step.

.. topic:: Assembly pipeline, as driven by Daijin

    .. figure:: daijin_assemble.svg
        :align: center
        :scale: 50%
        :figwidth: 100%


        Example of a pipeline to assemble a single paired-end read dataset using three aligners (STAR [STAR]_, Hisat [Hisat]_, TopHat2 [TopHat2]_ ) and two different RNA-Seq assemblers (StringTie [StringTie]_, CLASS2 [Class2]_ ).

Usage::

    $ daijin assemble --help
    usage: daijin assemble [-h] [-c HPC_CONF] [--jobs N] [--cores [N]]
                           [--threads N] [--no_drmaa] [--rerun-incomplete]
                           [--forcerun TARGET [TARGET ...]] [--detailed-summary]
                           [--list] [--dag]
                           config

    positional arguments:
      config                Configuration file to use for running the transcript
                            assembly pipeline.

    optional arguments:
      -h, --help            show this help message and exit
      -c HPC_CONF, --hpc_conf HPC_CONF
                            Configuration file that allows the user to override
                            resource requests for each rule when running under a
                            scheduler in a HPC environment.
      --jobs N, -J N        Maximum number of cluster jobs to execute
                            concurrently.
      --cores [N], -C [N]   Use at most N cores in parallel (default: 1000).
      --threads N, -t N     Maximum number of threads per job. Default: None (set
                            in the configuration file)
      --no_drmaa, -nd       Use this flag if you wish to run without DRMAA, for
                            example, if running on a HPC and DRMAA is not
                            available, or if running locally on your own machine
                            or server.
      --rerun-incomplete, --ri
                            Re-run all jobs the output of which is recognized as
                            incomplete.
      --forcerun TARGET [TARGET ...], -R TARGET [TARGET ...]
                            Force the re-execution or creation of the given rules
                            or files. Use this option if you changed a rule and
                            want to have all its output in your workflow updated.
      --detailed-summary, -D
                            Print detailed summary of all input and output files
      --list, -l            List resources used in the workflow
      --dag                 Do not execute anything and print the redirected
                            acylic graph of jobs in the dot language.





.. _daijin-mikado:

Mikado
~~~~~~



In this step, the Daijin manager will execute all the steps necessary to perform Mikado on the desired inputs. The manager will execute the following steps:

#. Merge all the input assemblies together using :ref:`Mikado prepare <prepare>`
#. Execute TransDecoder_ [Trinity]_ on the transcript sequences, to retrieve their ORFs.
#. Split the FASTA file in as many chunks as specified during configuration, and analyse them separately
#. Execute `BLASTX+ <https://www.ncbi.nlm.nih.gov/books/NBK279690/>`_ [Blastplus]_ on the splitted FASTAs, creating BLAST XML outputs.
#. Run :ref:`Mikado serialise <serialise>` to load the BLAST results, TransDecoder ORFs, and portcullis junctions into a single database.
#. Run :ref:`Mikado pick <pick>` on the data, in the selected modes.
#. Collate and collapse the statistics for each of the filtered assemblies.


.. topic:: Mikado pipeline, as driven by Daijin

    .. figure:: daijin_mikado.svg
        :align: center
        :scale: 50%


        Example of a typical Mikado pipeline. In this case the number of chunks for BLAST is limited - 10 - but we advise to increase this number for big datasets.


Command line usage::

    $ daijin mikado --help
    usage: daijin mikado [-h] [-c HPC_CONF] [--jobs N] [--cores [N]] [--threads N]
                         [--no_drmaa] [--rerun-incomplete]
                         [--forcerun TARGET [TARGET ...]] [--detailed-summary]
                         [--list] [--dag]
                         config

    positional arguments:
      config                Configuration file to use for running the Mikado step
                            of the pipeline.

    optional arguments:
      -h, --help            show this help message and exit
      -c HPC_CONF, --hpc_conf HPC_CONF
                            Configuration file that allows the user to override
                            resource requests for each rule when running under a
                            scheduler in a HPC environment.
      --jobs N, -J N        Maximum number of cluster jobs to execute
                            concurrently.
      --cores [N], -C [N]   Use at most N cores in parallel (default: 1000).
      --threads N, -t N     Maximum number of threads per job. Default: None (set
                            in the configuration file)
      --no_drmaa, -nd       Use this flag if you wish to run without DRMAA, for
                            example, if running on a HPC and DRMAA is not
                            available, or if running locally on your own machine
                            or server.
      --rerun-incomplete, --ri
                            Re-run all jobs the output of which is recognized as
                            incomplete.
      --forcerun TARGET [TARGET ...], -R TARGET [TARGET ...]
                            Force the re-execution or creation of the given rules
                            or files. Use this option if you changed a rule and
                            want to have all its output in your workflow updated.
      --detailed-summary, -D
                            Print detailed summary of all input and output files
      --list, -l            List resources used in the workflow
      --dag                 Do not execute anything and print the redirected
                            acylic graph of jobs in the dot language.


.. The part regarding being able to use directly Mikado configure is not true yet. Gotta work on it!
.. tip:: If you have already created some assemblies and wish to analyse them with Daijin, it is also possible to :ref:`configure Mikado externally <configure>` and use the resulting configuration file to guide Daijin.
