.. _Snakemake: https://bitbucket.org/snakemake/snakemake/wiki/Home
.. _YAML: http://www.yaml.org/spec/1.2/spec.html
.. _TransDecoder: https://github.com/TransDecoder/TransDecoder
.. _Portcullis: https://github.com/maplesond/portcullis

.. _assemble_pipeline:

.. |snake_badge| image:: https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)]
   :target: http://snakemake.bitbucket.org

.. _Daijin:

The Daijin pipeline for driving Mikado
======================================

|snake_badge|

No emperor or empress can lead its nation without a trusty chancellor to help him or her in organising the bureaucracy. Daijin, the Japanese minister, has the same role in Mikado
- it smooths the path to go from a collection of read inputs (both RNA-Seq or long reads) to a polished transcriptome assembly.
The pipeline is based on Snakemake_ [Snake]_; while Snakemake can support any scheduling system, our pipeline manager currently supports only three
(SLURM, PBS and LSF), plus any DRMAA-compliant batch submission system. Other schedulers can be added upon request.

This page contains a detailed explanation on how to use Daijin. We also provide a :ref:`tutorial <Daijin-Tutorial>`
that will guide you through using the manager for analysing an RNA-Seq sample for *Sc. pombe*.

.. hint:: It is possible to launch the two steps of the pipeline directly with Snakemake, using the snakefiles located in
Mikado.daijin: :download:`assemble.smk <assemble.smk>` for the first step, and :download:`mikado.smk` for the second.

.. warning:: Starting from :TODO:


.. _daijin-configure:

Configuring Daijin
~~~~~~~~~~~~~~~~~~

``daijin configure`` creates the configuration file that will drive Daijin, in YAML_ format. Most options can be specified by command line. Available parameters for the command line are:

* *out*: Output file, ie the configuration file for Daijin.
* *od*: Directory that Daijin will use as master for its jobs.
* *genome*: Genome FASTA file. Required.
* *transcriptome*: Optional reference transcriptome, to be used for alignment and assembly. Please note that Mikado will **not** use the reference annotation, it will only be used during the alignment and assembly steps.
* *name*: Name to be applied to the genome, eg the species.

.. warning:: The name must avoid containing spaces or any non-standard character. This string will be used among other things to determine the name of the index to be used by the aligners; if the string contains non-valid characters it will cause them - and on cascade Daijin - to fail.

* *prot-db*: this parameter specifies a protein database to be used for BLAST. If none is specified, this step will be omitted.
* *aligners*: aligner(s) to be used during the run. Currently, Daijin supports the following aligners:

    * *gsnap*
    * *star*
    * *tophat2*
    * *hisat2*
* *assemblers*: assembler(s) to be used during the run. Currently, Daijin supports the following RNA-Seq assemblers:

    * *cufflinks*
    * *class2*
    * *stringtie*
    * *trinity*
* *threads*: Number of threads to be requested for parallelisable steps.
* *modes*: Daijin can run Mikado in multiple modes regarding the :ref:`handling of putative chimeras <chimera_splitting_algorithm>`. Specify those you desire to execute here.
* *flank*: Amount of flanking that Mikado will allow while looking for fragments around good gene loci. Default 1000 bps. It is advisable to reduce it in species with compact genomes.
* *scheduler*: if Daijin has to execute the pipeline on a cluster (potentially using DRMAA), it is necessary to specify the scheduler here. At the moment we support the following widely used schedulers: PBS, LSF, SLURM.
* *r1*, *r2*, *samples*: RNA-Seq reads 1, RNA-Seq reads 2, and sample name. At least one of each is required.
* *strandedness*: if not specified, all samples will be assumed to be unstranded. Specify it as you would with HISAT or TopHat2.
* *cluster_config*: if specified, a sample cluster configuration file will be copied in the specified location. A cluster configuration file has the following structure:

.. literalinclude:: hpc.yaml

As it can be seen, it is a YAML format with two fields: __default__ and asm_trinitygg. The former specifies common parameters for all the steps: the queue to be used, the number of threads to request, and the memory per job. The "asm_trinitygg" field specifies particular parameters to be applied to the asm_trinitygg rule - in this case, increasing the requested memory (Trinity is a *de novo* assembler and uses much more memory than most reference-guided assemblers). Please note that **the number of threads specified on the command line, or in the configuration file proper, takes precedence over the equivalent parameter in the cluster configuration file**.

.. code-block:: bash

    $ daijin configure --help
    usage: daijin configure [-h] [-c CLUSTER_CONFIG] [--threads N] [-od OUT_DIR]
                            [-o OUT] [--scheduler {,SLURM,LSF,PBS}] [--name NAME]
                            --genome GENOME [--transcriptome TRANSCRIPTOME]
                            [-r1 R1 [R1 ...]] [-r2 R2 [R2 ...]]
                            [-s SAMPLES [SAMPLES ...]]
                            [-st {fr-unstranded,fr-secondstrand,fr-firststrand} [{fr-unstranded,fr-secondstrand,fr-firststrand} ...]]
                            -al
                            [{gsnap,star,hisat,tophat2} [{gsnap,star,hisat,tophat2} ...]]
                            -as
                            [{class,cufflinks,stringtie,trinity} [{class,cufflinks,stringtie,trinity} ...]]
                            [--scoring {insects.yaml,human.yaml,plants.yaml,worm.yaml,spombe.yaml}]
                            [--copy-scoring COPY_SCORING]
                            [-m {nosplit,split,permissive,stringent,lenient} [{nosplit,split,permissive,stringent,lenient} ...]]
                            [--flank FLANK] [--prot-db PROT_DB [PROT_DB ...]]


    optional arguments:
      -h, --help            show this help message and exit
      -al [{gsnap,star,hisat,tophat2} [{gsnap,star,hisat,tophat2} ...]], --aligners [{gsnap,star,hisat,tophat2} [{gsnap,star,hisat,tophat2} ...]]
                            Aligner(s) to use for the analysis. Choices: gsnap,
                            star, hisat, tophat2
      -as [{class,cufflinks,stringtie,trinity} [{class,cufflinks,stringtie,trinity} ...]], --assemblers [{class,cufflinks,stringtie,trinity} [{class,cufflinks,stringtie,trinity} ...]]
                            Assembler(s) to use for the analysis. Choices: class,
                            cufflinks, stringtie, trinity

    Options related to how to run Daijin - threads, cluster configuration, etc.:
      -c CLUSTER_CONFIG, --cluster_config CLUSTER_CONFIG
                            Cluster configuration file to write to.
      --threads N, -t N     Maximum number of threads per job. Default: 4
      -od OUT_DIR, --out-dir OUT_DIR
                            Output directory. Default if unspecified: chosen name.
      -o OUT, --out OUT     Output file. If the file name ends in "json", the file
                            will be in JSON format; otherwise, Daijin will print
                            out a YAML file. Default: STDOUT.
      --scheduler {,SLURM,LSF,PBS}
                            Scheduler to use. Default: None - ie, either execute
                            everything on the local machine or use DRMAA to submit
                            and control jobs (recommended).

    Arguments related to the reference species.:
      --name NAME           Name of the species under analysis.
      --genome GENOME, -g GENOME
                            Reference genome for the analysis, in FASTA format.
                            Required.
      --transcriptome TRANSCRIPTOME
                            Reference annotation, in GFF3 or GTF format.

    Arguments related to the input paired reads.:
      -r1 R1 [R1 ...], --left_reads R1 [R1 ...]
                            Left reads for the analysis. Required.
      -r2 R2 [R2 ...], --right_reads R2 [R2 ...]
                            Right reads for the analysis. Required.
      -s SAMPLES [SAMPLES ...], --samples SAMPLES [SAMPLES ...]
                            Sample names for the analysis. Required.
      -st {fr-unstranded,fr-secondstrand,fr-firststrand} [{fr-unstranded,fr-secondstrand,fr-firststrand} ...], --strandedness {fr-unstranded,fr-secondstrand,fr-firststrand} [{fr-unstranded,fr-secondstrand,fr-firststrand} ...]
                            Strandedness of the reads. Specify it 0, 1, or number
                            of samples times. Choices: fr-unstranded, fr-
                            secondstrand, fr-firststrand.

    Options related to the Mikado phase of the pipeline.:
      --scoring {insects.yaml,human.yaml,plants.yaml,worm.yaml}
                            Available scoring files.
      --copy-scoring COPY_SCORING
                            File into which to copy the selected scoring file, for
                            modification.
      -m {nosplit,split,permissive,stringent,lenient} [{nosplit,split,permissive,stringent,lenient} ...], --modes {nosplit,split,permissive,stringent,lenient} [{nosplit,split,permissive,stringent,lenient} ...]
                            Mikado pick modes to run. Choices: nosplit, split,
                            permissive, stringent, lenient
      --flank FLANK         Amount of flanking for grouping transcripts in
                            superloci during the pick phase of Mikado.
      --prot-db PROT_DB [PROT_DB ...]
                            Protein database to compare against, for Mikado.

.. warning:: if DRMAA is requested and no scheduler is specified, Daijin will fail. For this reason, Daijin *requires* a scheduler name. If one is not provided, Daijin will fall back to local execution.

Tweaking the configuration file
-------------------------------

Daijin can be configured so to run each assembler and/or aligner multiple times, with a different set of parameters each. This is achieved by specifying the additional command line arguments in the array for each program. For example, in this section:

.. code-block:: yaml

    align_methods:
      hisat:
      - ''
      - "-k 10"
    asm_methods:
      class:
      - ''
      - "--set-cover"
      stringtie:
      - ''
      - "-T 0.1"

we are asking Daijin to execute each program twice:

* HISAT2: once with default parameters, once reporting the 10 best alignments for each pair
* CLASS2: once with default parameters, once using the "set-cover" algorithm
* Stringtie: once with default parameters, once lowering the minimum expression to 0.1 TPM.

If you are running Daijin on a cluster and the software tools have to be sourced or loaded, you can specify the versions and command to load in the configuration file itself. Edit the file at this section:

.. code-block:: yaml

    load:
      #  Commands to use to load/select the versions of the programs to use. Leave an empty
      #  string if no loading is necessary.
      blast: 'source blast-2.3.0'
      class: 'source class-2.12'
      cufflinks: ''
      gmap: ''
      hisat: 'source HISAT-2.0.4'
      mikado: 'souce mikado-devel'
      portcullis: 'source portcullis-0.17.2'
      samtools: 'source samtools-1.2'
      star: ''
      stringtie: 'source stringtie-1.2.4'
      tophat: ''
      transdecoder: 'source transdecoder-3.0.0'
      trinity: ''

In this case, the cluster requires to load software on-demand using the source command, and we specify the versions of the programs we wish to use.

Regarding read alignment, it is also possible to specify the minimum and maximum permissible intron lengths:

.. code-block:: yaml

 short_reads:
  #  Parameters related to the reads to use for the assemblies. Voices:
  #  - r1: array of left read files.
  #  - r2: array of right read files. It must be of the same length of r1; if one
  #    one or more of the samples are single-end reads, add an empty string.
  #  - samples: array of the sample names. It must be of the same length of r1.
  #  - strandedness: array of strand-specificity of the samples. It must be of the
  #    same length of r1. Valid values: fr-firststrand, fr-secondstrand, fr-unstranded.
  max_intron: 10000
  min_intron: 20
  r1:
  - SRR1617247_1.fastq.gz
  r2:
  - SRR1617247_2.fastq.gz
  samples:
  - SRR1617247
  strandedness:
  - fr-secondstrand

Regarding the Mikado stage of Daijin, the configuration file contains all the fields that can be found in a :ref:`normal Mikado configuration file <configure>`. All mikado-specific parameters are stored under the "mikado" field. It is also possible to modify the following:

* *blastx*: these are parameters regarding the running of BLASTX. This field contains the following variables:

   * prot-db: this is an **array** of FASTA files to use as protein databases for the BLASTX step.
   * chunks: number of chunks to divide the BLASTX into. When dealing with big input FASTA files and with a cluster at disposal, it is more efficient to chunk the input FASTA in multiple smaller files and execute BLASTX on them independently. The default number of chunks is 10. Increase it as you see fit - it often makes sense, especially on clusters, to launch a multitude of small jobs rather than a low number of big jobs.
   * evalue: maximum e-value for the searches.
   * max_target_seqs: maximum number of hits to report.
* *transdecoder*: parameters related to transdecoder. At the moment only one is available:

    * *min_protein_len*: minimum protein length that TransDecoder should report. The default value set by Mikado, 30, is much lower than the default (100) and this is intentional, as the chimera splitting algorithm relies on the capability of TransDecoder of finding incomplete short ORFs at different ends of a transcript.

Structure of the output directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Daijin will organise the output directory in 5 major sections, plus the configuration file for the Mikado step:

#. *1-reads*: this folder contains links to the original read files.
#. *2-alignments*: this folder stores the indices built for each aligner, and the results. The structure is as follows:

    * *output*: this folder contains the final BAM files, sorted and indexed.
    * *alignments.stats*: this file contains a table reporting the most salient parameters derived from ``samtools stats`` calls onto the multiple BAM files.
    * One folder per aligner, inside which are present:

        * *index* folder, containing the genome indices for the tool
        * One folder per sample, containing the output BAM file.
#. *3-assemblies*: this folder contains the RNA-Seq assemblies which will be used as input by Mikado. The structure is as follows:

    * *output*: this folder contains the final GTF/GFF3 assemblies. For each combination of aligner/assembler that has been requested, Daijin will here provide:

        * the GTF file
        * a statistics file derived using ``mikado util stats`` (see the :ref:`section on this utility <stat-command>` for details)
    * *assembly.stats*: a tabular file collecting the most salient data from the statistics files generated for each assembly
    * One folder per assembly, containing tool-specific files and the final assembly.
#. *4-portcullis*: this folder contains the results of Portcullis_, if its execution has been requested. The folder will contain the following:

    * *output*: folder which contains a merged BED file of reliable junctions, creating by merging all junctions from all alignments.
    * One folder per alignment analysed. We redirect you to the `documentation of the tool <http://portcullis.readthedocs.io/en/latest/>`_ for more details.
#. *mikado.yaml*: final output file of the `assemble` part of the pipeline. This file will act both as the configuration for Daijin and for Mikado; for a description of the Mikado specific fields, we remand to the :ref:`section on the configuration of the tool <configure>`.
#. *5-mikado*: this folder contains the results for mikado. It is organised as follows:

    #. a link to the genome FASTA, and corresponding FAI file (generated with samtools)
    #. Files created by the :ref:`prepare step <prepare>`:

        * mikado_prepared.fasta
        * mikado_prepared.gtf
        * prepare.log
    #. *transdecoder*: this folder contains the results of the TransDecoder_ run against the mikado_prepared.fasta file. Mikado will use the file *transcripts.fasta.transdecoder.bed* as source for the ORF location.
    #. *blast*: this folder contains the BLAST data. In particular:

        * *index*: this folder contains the FASTA file of the protein database, and the BLAST database files.
        * *fastas*: Daijin will split mikado_prepared.fasta into multiple files, for easier runs onto a cluster. This folder contains the splitted FASTAs.
        * *xmls*: this folder contains the XML files corresponding to the BLASTs of the files present in *fastas*
        * *logs*: this folder contains the log files corresponding to the BLASTs of the files present in *fastas*
    #. *pick*: this folder contains the results of :ref:`Mikado pick <pick>`. It is organissed as follows:

        * One folder per requested :ref:`Mikado chimera-splitting mode <chimera_splitting_algorithm>`. Inside each folder, it is possible to find:

            * mikado-{mode}.loci.gff3: Final GFF3 output file.
            * mikado-{mode}.metrics.gff3: Final metrics output file, containing the metrics of the transcripts that have been selected.
            * mikado-{mode}.scores.gff3: Final metrics output file, containing the scores associated to the evaluated metrics, for each of the selected transcripts.
            * mikado-{mode}.loci.stats: statistics file derived using ``mikado util stats`` (see the :ref:`section on this utility <stat-command>` for details)
        * *comparison.stats*: this tabular file collects the most salient data from the statistics files generated for each pick mode.

Running the pipeline
~~~~~~~~~~~~~~~~~~~~

Daijin executes the pipeline in two distinct phases, *assemble* and *mikado*. Both commands have the same command line interface, namely::

     $ daijin assemble --help
    usage: daijin assemble [-h] [-c HPC_CONF] [-d] [--jobs N] [--cores [N]]
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
      -d, --dryrun          Do a dry run for testing.
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


The available command parameters are:

* *config*: the configuration file.
* *hpc_conf*: cluster configuration file.
* *jobs*: Maximum number of jobs that can be executed (if Daijin is in local mode) or be present in the submission queue (if Daijin is in DRMAA/cluster mode) at any one time.
* *dryrun*: do not execute, just list all the commands that will be executed. Useful also for listing the rules that have to be executed.
* *cores*: Maximum number of cores that Daijin can claim at any one time.
* *threads*: Maximum number of cores/threads that can be assigned to any step of the pipeline.
* *rerun-incomplete*: Ask Snakemake to check which steps have produced empty or incomplete output files, and re-execute them and all the downstream commands.
* *forcerun*: force the re-run of a


.. _daijin-assemble:

Assemble
~~~~~~~~

In the first step of the pipeline, Daijin will perform the following operations for each of the read datasets provided:

#. Create the necessary indices for each of the aligner programs requested.
#. Align the read dataset using all the different tools requested, in all the possible combinations of parameters requested.

   * For example, it is possible to ask each dataset to be aligned twice with TopHat2 - once with the "micro-exon" mode activated, the second time without. Both alignments will be run independently.
   * It is possible to specify which datasets are strand-specific and which are not, and moreover, it is possible to specify the kind of strand-specificity (fr-secondstrand, fr-firststrand).
#. Call all the reliable junctions across the alignments using Portcullis_.
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


        Example of a pipeline to assemble a single paired-end read dataset using one aligner (Hisat [Hisat]_) and two different RNA-Seq assemblers (StringTie [StringTie]_ and CLASS2 [Class2]_ ). Reliable junctions from the three alignments are called and merged together using Portcullis_.


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

``daijin mikado`` by default should use as config **the configuration file created by** ``daijin assemble``, which will be located in <output directory>/mikado.yaml.

.. topic:: Mikado pipeline, as driven by Daijin

    .. figure:: daijin_mikado.svg
        :align: center
        :scale: 50%


        Example of a typical Mikado pipeline. In this case the number of chunks for BLAST is limited - 10 - but we advise to increase this number for big datasets.

.. The part regarding being able to use directly Mikado configure is not true yet. Gotta work on it!
.. hint:: If you have already created some assemblies and wish to analyse them with Daijin, it is also possible to :ref:`configure Mikado externally <configure>` and use the resulting configuration file to guide Daijin. At the time of this writing, this is also the recommended protocol for including eg Pacbio or EST alignments.
