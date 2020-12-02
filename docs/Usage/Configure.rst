.. _SQLAlchemy: http://www.sqlalchemy.org/
.. _Portcullis: https://github.com/maplesond/portcullis
.. _BED12: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _GffRead: https://github.com/gpertea/gffread
.. _JSON: https://www.json.org/
.. _TOML: https://toml.io/
.. _YAML: https://yaml.org/
.. _configparser: https://docs.python.org/3/library/configparser.html
.. _configure:

Mikado configure
================

This utility prepares the configuration file that will be used throughout the pipeline stages.
While the most important options can be set at runtime through the command line, many algorithmic details can be accessed and intervened upon only through the file produced through this command.

.. important::

  Please note that any value absent from the configuration at runtime **will be imputed to the default value for Mikado, as specified internally**.


.. _input_file_list:

Input annotation files
~~~~~~~~~~~~~~~~~~~~~~

The preferred method for providing Mikado configure with the input annotation files is a **<TAB>** delimited file, where each line specifies a different file name.
The fields in this file are as follows, for each row:

.. <file> <label> <strandedness(def. False)> <score(optional, def. 0)> <is_reference(optional, def. False)>
                        <exclude_redundant(optional, def. True)> strandedness, is_reference and exclude_redundant must be boolean values (True, False) score must be a valid floating number.

- File name of the input file. **Mandatory**
    - The file name must specify a valid path (either absolute or relative) from the folder where Mikado configure is launched to the file itself.
- Label for the input file. **Mandatory**
    - Each label must be **unique** within the run.
- Strandedness of the annotation file. Boolean (either True or False, capitalization is ignored). **Mandatory**.
    - Reference files will be considered as stranded even if this value is set to False.
- Score associated to the input file. Default 0. *Optional*. Floating numbers only.
    - The score will be used to determine tie winners during Mikado prepare in case of redundancy. It will also be applied to transcripts during the pick stage. See :ref:`<this section later on <source_score>`
- *is_reference*: boolean flag (either True or False, capitalization is ignored), default False.. "Reference" transcripts are treated slightly differently from the others in Mikado. Specifically:
    - During the prepare stage, reference transcripts are never discarded on account of redundancy, short size, etc.; moreover, their strand will never be reversed on account of the internal Mikado checks.
    - During the pick stage, reference transcripts can be given an additional boost in their scoring. Moreover, it is possible to instruct Mikado to only look for alternative splicing events of the reference transcripts.
- *exclude_redundant*: boolean flag (either True or False, capitalization is ignored), default False. If set to True, Mikado prepare will remove any transcript in this annotation that is considered redundant (ie identical or contained within another transcript).
    - By default, Mikado will remove only ***identical*** copies of transcripts across annotations. Datasets with this flag set will be treated as if they were going through GffRead_.
- *skip_split*: boolean flag, default False. If set, Mikado will *not* break suspected chimeras, ie transcripts with more than one internal ORF, during the pick stage.


Usage
~~~~~

This command will generate a configuration file (in either JSON or YAML format), with the correct configuration for the parameters set on the command line. See :ref:`the in-depth section on the structure of the configuration file <conf_anatomy>` for details.

Selected command line parameters:

- *--full*: By default, Mikado configure will output a stripped-down configuration file, with only some of the fields explicitly present. Use this flag to show all the available fields.
- *--seed*: random seed for full reproducibility of Mikado runs.
- *--list*: the argument to this option specifies the input annotation files for Mikado. See above.

Command line help:

.. code-block::

    $ mikado configure --help
    usage: Mikado configure [-h] [--full] [--seed SEED] [--minimum-cdna-length MINIMUM_CDNA_LENGTH] [--max-intron-length MAX_INTRON_LENGTH] [--scoring SCORING] [--copy-scoring COPY_SCORING]
                            [-i INTRON_RANGE INTRON_RANGE] [--subloci-out SUBLOCI_OUT] [--monoloci-out MONOLOCI_OUT] [--no-pad] [--only-reference-update] [-eri] [-kdc] [--check-references]
                            [-mco MIN_CLUSTERING_CDNA_OVERLAP] [-mcso MIN_CLUSTERING_CDS_OVERLAP] [--strand-specific] [--no-files | --gff GFF | --list LIST] [--reference REFERENCE] [--junctions JUNCTIONS]
                            [-bt BLAST_TARGETS] [--strand-specific-assemblies STRAND_SPECIFIC_ASSEMBLIES] [--labels LABELS] [--codon-table CODON_TABLE] [--external EXTERNAL] [--daijin] [-bc BLAST_CHUNKS]
                            [--use-blast] [--use-transdecoder] [--mode {nosplit,stringent,lenient,permissive,split} [{nosplit,stringent,lenient,permissive,split} ...]] [--scheduler {local,SLURM,LSF,PBS}]
                            [--exe EXE] [-c CLUSTER_CONFIG] [-t THREADS] [--skip-split SKIP_SPLIT [SKIP_SPLIT ...]] [-j | -y | --toml] [-od OUT_DIR]
                            [out]

    Configuration utility for Mikado

    positional arguments:
      out

    optional arguments:
      -h, --help            show this help message and exit
      --full
      --seed SEED           Random seed number.
      --strand-specific     Boolean flag indicating whether all the assemblies are strand-specific.
      --no-files            Remove all files-specific options from the printed configuration file. Invoking the "--gff" option will disable this flag.
      --gff GFF             Input GFF/GTF file(s), separated by comma
      --list LIST           Tab-delimited file containing rows with the following format: <file> <label> <strandedness(def. False)> <score(optional, def. 0)> <is_reference(optional, def. False)>
                            <exclude_redundant(optional, def. True)>  <skip_split(optional, def. False)> strandedness, is_reference, exclude_redundant and skip_split must be boolean values (True, False) score must be a valid floating number.
      --reference REFERENCE, --genome REFERENCE
                            Fasta genomic reference.
      --strand-specific-assemblies STRAND_SPECIFIC_ASSEMBLIES
                            List of strand-specific assemblies among the inputs.
      --labels LABELS       Labels to attach to the IDs of the transcripts of the input files, separated by comma.
      --codon-table CODON_TABLE
                            Codon table to use. Default: 0 (ie Standard, NCBI #1, but only ATG is considered a valid start codon.
      --external EXTERNAL   External configuration file to overwrite/add values from. Parameters specified on the command line will take precedence over those present in the configuration file.
      -t THREADS, --threads THREADS
      --skip-split SKIP_SPLIT [SKIP_SPLIT ...]
                            List of labels for which splitting will be disabled (eg long reads such as PacBio)
      -j, --json            Output will be in JSON (default: inferred by filename, with TOML as fallback).
      -y, --yaml            Output will be in YAML (default: inferred by filename, with TOML as fallback).
      --toml                Output will be in TOML (default: inferred by filename, with TOML as fallback).
      -od OUT_DIR, --out-dir OUT_DIR
                            Destination directory for the output.

    Options related to the prepare stage.:
      --minimum-cdna-length MINIMUM_CDNA_LENGTH
                            Minimum cDNA length for transcripts.
      --max-intron-length MAX_INTRON_LENGTH
                            Maximum intron length for transcripts.

    Options related to the scoring system:
      --scoring SCORING     Scoring file to use. Mikado provides the following: mammalian.yaml, plant.yaml, HISTORIC/athaliana_scoring.yaml, HISTORIC/celegans_scoring.yaml, HISTORIC/dmelanogaster_scoring.yaml,
                            HISTORIC/hsapiens_scoring.yaml, HISTORIC/human.yaml, HISTORIC/insects.yaml, HISTORIC/plants.yaml, HISTORIC/scerevisiae.yaml, HISTORIC/worm.yaml
      --copy-scoring COPY_SCORING
                            File into which to copy the selected scoring file, for modification.

    Options related to the picking:
      -i INTRON_RANGE INTRON_RANGE, --intron-range INTRON_RANGE INTRON_RANGE
                            Range into which intron lengths should fall, as a couple of integers. Transcripts with intron lengths outside of this range will be penalised. Default: (60, 900)
      --subloci-out SUBLOCI_OUT
                            Name of the optional subloci output. By default, this will not be produced.
      --monoloci-out MONOLOCI_OUT
                            Name of the optional monoloci output. By default, this will not be produced.
      --no-pad              Disable transcript padding. On by default.
      --only-reference-update
                            Flag. If switched on, Mikado will only keep loci where at least one of the transcripts is marked as "reference". CAUTION: new and experimental. If no transcript has been marked as
                            reference, the output will be completely empty!
      -eri, --exclude-retained-introns
                            Exclude all retained intron alternative splicing events from the final output. Default: False. Retained intron events that do not dirsupt the CDS are kept by Mikado in the final
                            output.
      -kdc, --keep-disrupted-cds
                            Keep in the final output transcripts whose CDS is most probably disrupted by a retained intron event. Default: False. Mikado will try to detect these instances and exclude them from
                            the final output.
      --check-references    Flag. If switched on, Mikado will also check reference models against the general transcript requirements, and will also consider them as potential fragments. This is useful in the
                            context of e.g. updating an *ab-initio* results with data from RNASeq, protein alignments, etc.
      -mco MIN_CLUSTERING_CDNA_OVERLAP, --min-clustering-cdna-overlap MIN_CLUSTERING_CDNA_OVERLAP
                            Minimum cDNA overlap between two transcripts for them to be considered part of the same locus during the late picking stages. NOTE: if --min-cds-overlap is not specified, it will be
                            set to this value! Default: 20%.
      -mcso MIN_CLUSTERING_CDS_OVERLAP, --min-clustering-cds-overlap MIN_CLUSTERING_CDS_OVERLAP
                            Minimum CDS overlap between two transcripts for them to be considered part of the same locus during the late picking stages. NOTE: if not specified, and --min-cdna-overlap is
                            specified on the command line, min-cds-overlap will be set to this value! Default: 20%.

    Options related to the serialisation step:
      --junctions JUNCTIONS
      -bt BLAST_TARGETS, --blast_targets BLAST_TARGETS

    Options related to configuring a Daijin run.:
      --daijin              Flag. If set, the configuration file will be also valid for Daijin.
      -bc BLAST_CHUNKS, --blast-chunks BLAST_CHUNKS
                            Number of parallel DIAMOND/BLAST jobs to run. Default: 10.
      --use-blast           Flag. If switched on, Mikado will use BLAST instead of DIAMOND.
      --use-transdecoder    Flag. If switched on, Mikado will use TransDecoder instead of Prodigal.
      --mode {nosplit,stringent,lenient,permissive,split} [{nosplit,stringent,lenient,permissive,split} ...]
                            Mode(s) in which Mikado will treat transcripts with multiple ORFs. - nosplit: keep the transcripts whole. - stringent: split multi-orf transcripts if two consecutive ORFs have both
                            BLAST hits and none of those hits is against the same target. - lenient: split multi-orf transcripts as in stringent, and additionally, also when either of the ORFs lacks a BLAST hit
                            (but not both). - permissive: like lenient, but also split when both ORFs lack BLAST hits - split: split multi-orf transcripts regardless of what BLAST data is available. If multiple
                            modes are specified, Mikado will create a Daijin-compatible configuration file.
      --scheduler {local,SLURM,LSF,PBS}
                            Scheduler to use. Default: None - ie, either execute everything on the local machine or use DRMAA to submit and control jobs (recommended).
      --exe EXE             Configuration file for the executables.
      -c CLUSTER_CONFIG, --cluster_config CLUSTER_CONFIG
                            Cluster configuration file to write to.

.. _conf_anatomy:

Anatomy of the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format of the configuration file
--------------------------------

The configuration files accepted by Mikado can be in any of three dialects:

- TOML_, the default choice. TOML is an intuitive configuration file format, similar to the INI files :ref:`preferred by Python <configparser>`.
- YAML_, a human-readable configuration file format based on indentation. Less preferred because of the unreadability of deeply-nested values.
- JSON_, a less human-readable file format that is commonly used to pass data across processes / programs.

We leave freedom to the user to select their preferred file format. In this section, we will use TOML_ to explain the different sections of the file.

Global options
--------------

The following options apply to all programs in the Mikado pipeline, and they refer to general parameters such as logging verbosity, number of threads, etc.

Parameters:
- *threads*: this is the number of processes/threads that will be requested by the Mikado programs. This parameter can be overridden on the command line.
- *seed*: random seed specification, to ensure maximum reproducibility of the run.
.. _start-methods:
- *multiprocessing_method*: this specifies the way that :ref:`Python will start children processes <https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods>`. The possible choices are "spawn" (default), "fork" and "fork-server". See the sidebar for a more complete explanation.

.. _scheduler-multiprocessing:
.. sidebar:: "Python, multiprocessing, and cluster schedulers"

    Some schedulers, in particular SLURM, are not capable to understand that the processes *forked* by Python are still sharing the same memory with the main process, and think instead that each process is using that memory in isolation. As a result, they might think that the Mikado process is using its memory multiplied by the number of processes - depending on when the forking happens - and therefore shut down the program as it *appears* to be using much more memory than needed. For this reason, :ref:`Daijin <Daijin>` forces Mikado to run in **spawn** mode. Although spawning is slower than forking, it happens only once per run, and it has therefore a limited cost in terms of runtime - while greatly reducing the chances of the program being shut down because of spurious "Out of memory" reasons.

.. code-block:: toml

    threads = 4
    seed = 0
    multiprocessing_method = "spawn"

Log settings
~~~~~~~~~~~~

It is possible to set high-level settings for the logs in the ``log_settings`` section:

- log_level: level of the logging for Mikado. Options: *DEBUG, INFO, WARNING, ERROR, CRITICAL*. By default, Mikado will be quiet and output log messages of severity *WARNING* or greater.
- sql_level: level of the logging for messages regarding the database connection (through `SQLAlchemy`_). By default, SQLAlchemy will be set in quiet mode and asked to output only messages of severity *WARNING* or greater.

.. warning:: Mikado and SQLAlchemy can be greatly verbose if asked to output *DEBUG* or *INFO* messages, to the point of slowing down the program significantly due to the amount of writing to disk. Please consider setting the level to *DEBUG* only when there is a real problem to debug, not otherwise!

.. code-block:: toml

    [log_settings]
    # Settings related to the logs. Keys:
    # - sql_level: verbosity for SQL calls. Default: WARNING. In decreasing order: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    # - log_level: verbosity. Default: INFO. In decreasing order: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    log_level = "INFO"
    sql_level = "WARNING"
    log = ""

.. _db-settings:

Database settings
-----------------

This section deals with the database settings that will be necessary for the :ref:`serialisation <serialise>` and :ref:`picking <pick>` phases of the pipeline. By default, Mikado will use a `SQLite database <https://www.sqlite.org/>`_, but it currently also supports `MySQL <http://www.mysql.com/>`_ and `PostgreSQL <https://www.postgresql.org/>`_ through SQLAlchemy_. Fields:

- db: name of the database to use. In case the database is SQLite, this will be the database file, otherwise it will be the database *name*.
- dbtype: one of:
  * sqlite
  * mysql
  * postgresql
- dbhost: host where the database is located. **Required with MySQL and PostgreSQL**.
- dbuser: User of the database. **Required with MySQL and PostgreSQL**.
- dbpasswd: Database password. **Required with MySQL and PostgreSQL**.
- dbport: Port to access to the database. It defaults to the normal ports for the selected database.

.. code-block:: toml

    [db_settings]
    # Settings related to DB connection. Parameters:
    # db: the DB to connect to. Required. Default: mikado.db
    # dbtype: Type of DB to use. Choices: sqlite, postgresql, mysql. Default: sqlite.
    db = "/c/Users/lucve/PycharmProjects/EICore/mikado/sample_data/mikado.db"
    dbtype = "sqlite"
    dbhost = "localhost"
    dbuser = ""
    dbpasswd = ""
    dbport = 0

.. _ref-settings:

Reference settings
------------------

This section of the configuration file deals with the reference genome. It specifies two fields:

- genome: the genome FASTA file. **Required**.
- genome_fai: FAI index of the genome. Used by :ref:`Mikado serialise <serialise>`, it can be inferred if left null.
- transcriptome: optional annotation file for the genome. Mikado currently ignores this field, but it is used by :ref:`Daijin <Daijin>` to guide some of the RNA-Seq assemblies.

.. code-block:: yaml

    [reference]
    genome = "chr5.fas.gz"
    genome_fai = ""
    transcriptome = ""

.. _prep-settings:

Settings for the prepare stage
------------------------------

This section of the configuration file deals with the :ref:`prepare stage of Mikado <prepare>`. It specifies the input files, their labels, and which of them are strand specific. The available fields are the following:

.. _canonical-configuration:

- *exclude_redundant*: if set to true, Mikado will only keep one copy of transcripts that are identical or contained into a different transcripts.
  - please note that this *global* values, if set to true, overrides the label-specific

- *canonical*: this voice specifies the splice site donors and acceptors that are considered canonical for the species. By default, Mikado uses the canonical splice site (GT/AG) and the two semi-canonical pairs (GC/AG and AT/AC). Type: Array of two-element arrays, composed by two-letter strings.
- *lenient*: boolean value. If set to *false*, transcripts that either only have non-canonical splice sites or have a mixture of canonical junctions on *both* strands will be **removed** from the output. Otherwise, they will left in, be properly tagged.
- *minimum_cdna_length*: minimum length of the transcripts to be kept.
- *max_intron_length*: Transcripts with introns greater than this will be **discarded**. The default is one million base pairs (effectively disabling the option).
- *strand_specific*: boolean. If set to *true*, **all** input assemblies will be treated as strand-specific, therefore keeping the strand of monoexonic fragments as it was. Multiexonic transcripts will not have their strand reversed even if doing that would mean making some or all non-canonical junctions canonical.
- *strip_cds*: boolean. If set to *true*, the CDS features will be stripped off the input transcripts. This might be necessary for eg transcripts obtained through alignment with `GMAP <http://research-pub.gene.com/gmap/>`_ [GMAP]_.
- *single*: boolean. For debug purposes only. If set to *true*, Mikado will disable multiprocessing.

.. code-block:: toml
    [prepare]
    # Options related to the input data preparation.
    # - procs: Number of processes to use.
    # - strand_specific: if set to True, transcripts will be assumed to be in the correct orientation, no strand flipping or removal
    # - strip_cds: Boolean. It indicates whether to remove the CDS from the predictions during preparation.
    exclude_redundant = false
    minimum_cdna_length = 200
    max_intron_length = 1000000
    strip_cds = false
    single = false
    lenient = false
    strand_specific = false
    canonical = [["GT", "AG"], ["GC", "AG"], ["AT", "AC"]]


Settings for the prepare stage: files settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. important:: As this section contains multiple linked lists, it is recommended to not edit this part of the configuration file directly, but rather, to rely on the mikado configure utility / mikado prepare interface to set it up. Specifically, setting up this section through the use of a :ref:`file of file names <input_file_list>` is highly recommended.

This sub-section is the most important for `prepare`, as it contains among other things the locations and labels for the input files.

  - *output_dir*: destination folder for the output files and the log. It will be created automatically, if it does not already exist on disk.
  - *out*: name of the output GTF file. Default: *mikado_prepared.gtf*.
  - *out_fasta*: name of the output GTF file. Default: *mikado_prepared.fasta*.
  - *log*: name of the log file. Default: *prepare.log*.
  - *gff*: list of filenames of the input files.
  - *labels*: list of labels associated with the input files.
  - *reference*: list of boolean values, indicating whether each input assembly is to be considered of "reference" quality.
  - *strand_specific_assemblies*: list of boolean values, indicating whether each input assembly is to be considered having a trustworthy strand information, or not.
  - *strip_cds*: list of boolean values, indicating whether the CDS of a given assembly should be ignored.
  - *exclude_redudant*: list of boolean values, indicating whether redundant models prsent in this assembly should be discarded or not.
  - *source_score*: dictionary linking the scores of each different assembly to a specific score, **using the label as key**, which will be applied in two different points:
    + during the prepare stage itself, in order to give an order priority for transcripts that come from different assemblies.
    + during the picking stage, this score will be added to each model from this assembly, therefore influencing the picked models.

.. code-block:: toml

    [prepare.files]
    # Options related to the input and output files.
    # - out: output GTF file
    # - out_fasta: output transcript FASTA file
    # - gff: array of input predictions for this step.
    # - labels: labels to be associated with the input GFFs. Default: None.
    # - reference: these files are treated as reference-like, ie, these transcripts will never get discarded
    #   during the preparation step.
    output_dir = "."
    out = "mikado_prepared.gtf"
    out_fasta = "mikado_prepared.fasta"
    log = "prepare.log"
    gff = ["class.gtf", "cufflinks.gtf", "stringtie.gtf", "trinity.gff3", "reference.gff3"]
    labels = ["cl", "cuff", "st", "tr", "at"]
    strand_specific_assemblies = ["class.gtf", "cufflinks.gtf", "stringtie.gtf", "reference.gff3"]
    reference = [false, false, false, false, true]
    exclude_redundant = [false, false, true, false, true]
    strip_cds = [false, false, false, false, false]

    [prepare.files.source_score]
    cl = 0
    cuff = 0
    st = 1.0
    tr = -0.5
    at = 5.0

.. _serialise-settings:

Settings for the serialisation stage
------------------------------------

This section of the configuration file deals with the :ref:`serialisation stage of Mikado <serialise>`. It specifies the location of the ORF BED12 files from TransDecoder, the location of the XML files from BLAST, the location of portcullis junctions, and other details important at run time. It has the following voices:

- *substitution_matrix*: the matrix used by BLAST or DIAMOND. Default is the standard BLOSUM62.
- *force*: whether the database should be truncated and rebuilt, or just updated.
- *max_objects*: this parameter is quite important when running with a SQLite database. SQLite does not support caching on the disk before committing the changes, so that every change has to be kept in memory. This can become a problem for RAM quite quickly. On the other hand, committing is an expensive operation, and it makes sense to minimise calls as much as possible. This parameter specifies the maximum number of objects Mikado will keep in memory before committing them to the database. The default number of 10 million privileges speed over RAM parsimony.

.. _max-regression::

- *max_regression*: this parameter is a float comprised between 0 and 1. Prodigal and TransDecoder will sometimes output open ORFs even in the presence of an in-frame start codon. Mikado can try to "regress" along the ORF until it finds one such start codon. This parameter imposes how much Mikado will regress, in percentage of the cDNA length.

.. note:: Recent versions of TransDecoder perform by default an analogous process. As such, we advise to keep this switch off if TransDecoder is used.

- *codon_table*: this parameter indicates the codon table to use. We use the `NCBI nomenclature <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_, with a variation:

  - the code "0" is added to indicate a variation on the standard code (identifier "1"), which differs only in that
    only "ATG" is considered as a valid start codon. This is because *in silico* ORF predictions tend to over-predict
    the presence of non-standard "ATG" codons, which are rare in nature.
- *max_target_seqs*: equivalent to the BLAST+ parameter of the same name - it indicates the maximum number of discrete hits
  that can be assigned to one sequence in the database.
- *single_thread*: boolean, if set to *true* it will forcibly disable multi-threading. Useful mostly for debugging purposes.

.. code-block:: toml

    [serialise]
    # Options related to serialisation
    # - force: whether to drop and reload everything into the DB
    # - files: options related to input files
    # - max_objects: Maximum number of objects to keep in memory while loading data into the database
    # - max_regression: if the ORF lacks a valid start site, this percentage indicates how far
    #   along the sequence Mikado should look for a good start site. Eg. with a value of 0.1,
    #   on a 300bp sequence with an open ORF Mikado would look for an alternative in-frame start codon
    #   in the first 30 bps (10% of the cDNA).
    # - max_target_seqs: equivalently to BLAST, it indicates the maximum number of targets to keep
    #   per blasted sequence.
    # - discard_definition: Boolean. **Deprecated**, it was used for specifying how to load BLAST files.
    # - single_thread: if true, Mikado prepare will force the usage of a single thread in this step.
    # - codon_table: codon table to use for verifying/modifying the ORFs. Default: 0, ie
    #  the universal codon table but enforcing as only valid start codon ATG.
    substitution_matrix = "blosum62"
    max_objects = 10000000
    max_regression = 0.2
    start_adjustment = true
    max_target_seqs = 100000
    force = false
    single_thread = false
    codon_table = 0


Settings for the serialisation stage: files settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This sub-section of the configuration file codifies the location of the input and output files for `serialise`. It contains
the following voices:

.. _reliable_junctions:
- junctions: array of locations of reliable junction files. These must be in BED12 format. The preferred source for this
  is :ref:`Portcullis` [Portcullis]_.
- log: log file.
- orfs: array of locations of ORFs location on the cDNA, as created by eg TransDecoder [Trinity]_.
- output_dir: output directory where the log file and the SQLite database will be written to (if SQLite has been chosen as the database type)
- transcripts: input transcripts. This should be set to be equal to the output of :ref:`Mikado prepare <prepare>`,
  ie the "out_fasta" field of the :ref:`prepare section of the configuration file <prep-settings>`.
- external_scores: this field indicates the location of a tabular file containing additional numeric values to be added to Mikado.
- xml: this array indicates the location of the BLAST output file(s). Please see the :ref:`section on serialisation <serialise_input_blast>` for details. Elements of the array can be:
.. _input_xml:
  + A **custom-formatted** BLAST tabular output file
  + BLAST+ XML files (optionally compressed with gzip)
  + BLAST+ ASN files (optionally compressed with gzip), which will be converted in-memory using ``blast_formatter``
  + a folder containing files of the above types.

.. code-block:: toml

    [serialise.files]
    junctions = ["junctions.bed"]
    xml = []
    blast_loading_debug = false
    external_scores = ""
    orfs = []
    transcripts = "mikado_prepared.fasta"
    log = "serialise.log"
    blast_targets = ["uniprot_sprot_plants.fasta"]
    output_dir = "."

.. _misc-settings:

Settings for the pick stage
---------------------------

This section of the configuration file deals with the :ref:`picking stage of Mikado <pick>`. It specifies details on how to handle BLAST and ORF data, which alternative splicing events are considered as valid during the final stages of the picking, and other important algorithmic details. The section comprises the following subsections:

- alternative_splicing: Options related to which AS events are considered as valid for the primary transcript in a locus.
- chimera_split: Options related to how to handle transcripts with multiple valid ORFs.
- files: Input and output files.
- orf_loading: Options related to how to decide which ORFs to load onto each transcript.
- output_format: options related to how to format the names of the transcripts, the source field of the GFFs, etc.
- run_options: Generic options related either to the general algorithm or to the number of resources requested.
.. _scoring_file_conf:
- scoring_file: This value specifies the :ref:`scoring file <scoring_files>` to be used for Mikado. These can be found in Mikado.configuration.scoring_files.
.. hint:: It is possible to ask for the configuration file to be copied in-place for customisation when calling ``mikado configure``.

Each subsection of the pick configuration will be explained in its own right.

.. _source_score:
Giving different priorities to transcripts from different assemblies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to specify boni/mali to be assigned to specific labels. Eg, it might be possible to assign a bonus of 1 to any transcript coming from PacBio reads, or a malus to any transcript coming from a given assembler. Example of such a configuration:

.. code-block:: toml

    [prepare.files.source_score]
    cl = 0
    cuff = 0
    st = 1.0
    tr = -0.5
    at = 5.0

In this example, we are prioritising the reference annotation ("at") by five points, the StringTie assembly by 1, and slightly penalising the Trinity assembly with a malus of half a point.

.. _configure-alternative-splicing:

Parameters regarding the alternative splicing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After selecting the best model for each locus, Mikado will backtrack and try to select valid alternative splicing events. This section deals with how Mikado will operate the selection. In order to be considered as valid potential AS events, transcripts have to satisfy the minimum :ref:`requirements specified in the scoring file <requirements-section>`. These are the available parameters:

- *report*: boolean. Whether to calculate and report possible alternative splicing events at all. By default this is set to true; ***setting this parameter to false will inactivate all the options in this section***.
- *keep_retained_introns*: boolean. It specifies whether transcripts with :ref:`retained introns <retained_intron_definition>` will be accepted as potentially valid AS events. By default, they are.
- *keep_retained_introns*: boolean. It specifies whether transcripts with :ref:`a CDS disrupted by their retained intron <retained_intron_disrupted_cds>` will be accepted as potentially valid AS events. By default, Mikado will exclude them.
- *min_cdna_overlap*: minimum cDNA overlap between the primary transcript and the AS candidate. By default, this is set to 0.5 (50%). It must be a number between 0 and 1.
- *min_cds_overlap*: minimum CDS overlap between the primary transcript and the AS candidate. By default this is set to 0.6, ie 60%. It must be a number between 0 and 1.
- *min_score_perc*: Minimum percentage of the score of the primary transcript that any candidate AS must have to be considered. By default, this is set to 0.5 (50%). It must be a number between 0 and 1.
- *only_confirmed_introns*: boolean. If set to true (default), Mikado will consider as potential AS events only transcripts whose introns *not shared with the primary transcript* are confirmed :ref:`in the dataset of reliable junctions <reliable_junctions>`.
- *redundant_ccodes*: any candidate AS will be :ref:`compared <Compare>` against all the transcripts already retained in the locus. If any of these comparisons returns one of the :ref:`class codes <ccodes>` specified in this array, **the transcript will be ignored**. The rationale is to avoid bringing back multiple minor variations of the same transcript. Default class codes: c, m, _, =, n.
- valid_ccodes: any candidate AS will be :ref:`compared <Compare>` against *the primary transcript* to determine the type of AS event. If the :ref:`class code <ccodes>` is one of those specified in this array, the transcript will be considered further. Valid class codes are within the categories "Alternative splicing", "Extension" with junction F1 lower than 100%, and Overlap (with the exclusion of "m"). Default class codes: j, J, G, h.

.. _pad-configuration:
- pad: boolean option. If set to True, Mikado will try to pad transcripts so that they share the same 5'. Please :ref:`see this section for further information <padding>`.
- ts_max_splices: numerical. When padding is activated, at *most* how many splice junctions can be introduced?
- ts_distance: numerical. When padding is activated, at *most* of how many base pairs can a transcript be extended?

.. warning:: the AS transcript event does not need to be a valid AS event for *all* transcripts in the locus, only against the *primary* transcript.
.. note:: when padding transcripts, Mikado will consider also transcripts with the same intron structure but differing end points (so "=" or "_"). These will be used to expand the UTRs of other transcripts; however, only *one* of these transcripts with identical structures will be reported in the end.

.. code-block:: toml

    [pick.alternative_splicing]
    # Parameters related to alternative splicing reporting.
    # - report: whether to report at all or not the AS events.
    # - min_cds_overlap: minimum overlap between the CDS of the primary transcript and any AS event. Default: 60%.
    # - min_cdna_overlap: minimum overlap between the CDNA of the primary transcript and any AS event.
    # Default: 0% i.e. disabled, we check for the CDS overlap.
    # - keep_retained_introns: Whether to consider as valid AS events where one intron
    # - max_isoforms: Maximum number of isoforms per locus. 1 implies no AS reported. Default: 5
    # is retained compared to the primary or any other valid AS. Default: false.
    # - valid_ccodes: Valid class codes for AS events. Valid codes are in categories
    # 'Alternative splicing', 'Extension' (with junction F1 lower than 100%), and Overlap (exluding m). Default: j, J, g, G, C, h
    # - max_utr_length: Maximum length of the UTR for AS events. Default: 10e6 (i.e. no limit)
    # - max_fiveutr_length: Maximum length of the 5'UTR for AS events. Default: 10e6 (i.e. no limit)
    # - max_threeutr_length: Maximum length of the 5'UTR for AS events. Default: 10e6 (i.e. no limit)
    # - min_score_perc: Minimum score threshold for subsequent AS events. Only transcripts with a score at least (best) * value are retained.
    # - only_confirmed_introns: bring back AS events only when their introns are either present in the primary transcript or in the set of confirmed introns.
    # - pad: boolean switch. If true, Mikado will pad all the transcript in a gene so that their ends are the same
    # - ts_distance: if padding, this is the maximum distance in base-pairs between the starts of transcripts to be considered to be padded together.
    # - ts_max_splices: if padding, this is the maximum amount of splicing junctions that the transcript to pad is allowed to cross. If padding would lead to cross more than this number, the transcript will not be padded.
    report = true
    min_cds_overlap = 0.5
    min_cdna_overlap = 0.6
    keep_retained_introns = true
    keep_cds_disrupted_by_ri = false
    max_isoforms = 10
    valid_ccodes = ["j", "J", "G", "h"]
    redundant_ccodes = ["c", "m", "_", "=", "n"]
    min_score_perc = 0.5
    only_confirmed_introns = true
    ts_distance = 2000
    pad = true
    ts_max_splices = 2

.. _clustering_specifics:

Parameters regarding the clustering of transcripts in loci
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section influences how Mikado clusters transcripts in its multi-stage selection. The available parameters are:

- *flank*: numerical. When constructing :ref:`Superloci <superloci>`, Mikado will use this value as the maximum distance
between transcripts for them to be integrated within the same superlocus.
- *cds_only*: boolean. If set to true, during the :ref:`picking stage <pick-algo>` Mikado will consider only the **primary ORF** to evaluate whether two transcripts intersect. Transcripts which eg. share introns in their UTR but have completely unrelated CDSs will be clustered separately. Disabled by default.
- *purge*: boolean. If true, any transcript failing the :ref:`specified requirements <requirements-section>` will be purged out. Otherwise, they will be assigned a score of 0 and might potentially appear in the final output, if no other transcript is present in the locus.
- *simple_overlap_for_monoexonic*: boolean. During the :ref:`second clustering <monosubloci>`, by default monoexonic transcripts are clustered together even if they have a very slight overlap with another transcript. Manually setting this flag to *false* will cause Mikado to cluster monoexonic transcripts only if they have a minimum amount of cDNA and CDS overlap with the other transcripts in the holder.
- *min_cdna_overlap*: numerical, between 0 and 1. Minimum cDNA overlap between two multiexonic transcripts for them to be considered as intersecting, if all other conditions fail.
- *min_cdna_overlap*: numerical, between 0 and 1. Minimum CDS overlap between two multiexonic transcripts for them to be considered as intersecting, if all other conditions fail.

.. code-block:: toml

    [pick.clustering]
    # Parameters related to the clustering of transcripts into loci.
    # - cds_only: boolean, it specifies whether to cluster transcripts only according to their CDS (if present).
    # - min_cds_overlap: minimal CDS overlap for the second clustering.
    # - min_cdna_overlap: minimal cDNA overlap for the second clustering.
    # - flank: maximum distance for transcripts to be clustered within the same superlocus.
    # - remove_overlapping_fragments: boolean, it specifies whether to remove putative fragments.
    # - purge: boolean, it specifies whether to remove transcripts which fail the minimum requirements check - or whether to ignore those requirements altogether.
    # - simple_overlap_for_monoexonic: boolean. If set to true (default), then any overlap mean inclusion in a locus for or against a monoexonic transcript. If set to false, normal controls for the percentage of overlap will apply.
    # - max_distance_for_fragments: maximum distance from a valid locus for another to be considered a fragment.
    cds_only = false
    min_cds_overlap = 0.2
    min_cdna_overlap = 0.2
    purge = true
    flank = 200
    simple_overlap_for_monoexonic = true

.. _fragment_options:

Parameters regarding the detection of putative fragments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section determines how Mikado treats :ref:`potential fragments in the output <fragments>`. Available options:

- *remove*: boolean, default true. If set to true, fragments will be excluded from the final output; otherwise, they will be printed out, but properly tagged.

- *max_distance*: numerical. For non-overlapping fragments, this value determines the maximum distance from the valid gene. Eg. with the default setting of 2000, a putative fragment at the distance of 1000 will be tagged and dealt with as a fragment; an identical model at a distance of 3000 will be considered as a valid gene and left untouched.

- *valid_class_codes*: valid :ref:`class codes <ccodes>` for potential fragments. Only Class Codes in the categories Overlap, Intronic, Fragment, with the addition of "_", are considered as valid choices.

.. code-block:: toml

    [pick.fragments]
    # Parameters related to the handling of fragments.
    # - remove: boolean. Whether to remove fragments or leave them, properly tagged.
    # - max_distance: maximum distance of a putative fragment from a valid gene.
    # - valid_class_codes: which class codes will be considered as fragments. Default: (p, P, x, X, i, m, _). Choices: '_' plus any class code with category 'Intronic', 'Fragment', or 'Overlap'.
    remove = true
    max_distance = 2000
    valid_class_codes = ["p", "P", "x", "X", "i", "m", "_", "e", "o"]


.. _orf_loading:

Parameters regarding assignment of ORFs to transcripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section of the configuration file deals with how to determine valid ORFs for a transcript from those present in the database. The parameters to control the behaviour of Mikado are the following:

- *minimal_orf_length*: minimal length of the *primary* ORF to be loaded onto the transcript. By default, this is set at 50 **bps** (not aminoacids)
- *minimal_secondary_orf_length*: minimal length of any ORF that can be assigned to the transcript after the first. This value should be set at a **higher setting** than minimal_orf_length, in order to avoid loading uORFs [uORFs]_ into the transcript, leading to :ref:`spurious break downs of the UTRs <chimera_splitting_algorithm>`. Default: 200 bps.
- *strand_specific*: boolean. If set to *true*, only ORFs on the plus strand (ie the same of the cDNA) will be considered. If set to *false*, monoexonic transcripts mihgt have their strand flipped.


.. code-block:: toml

    [pick.orf_loading]
    # Parameters related to ORF loading.
    # - minimal_secondary_orf_length: Minimum length of a *secondary* ORF to be loaded after the first, in bp. Default: 200 bps
    # - minimal_orf_length: Minimum length in bps of an ORF to be loaded, as the primary ORF, onto a transcript. Default: 50 bps
    # - strand_specific: Boolean flag. If set to true, monoexonic transcripts will not have their ORF reversed even if they would have an ORF on the opposite strand.
    minimal_secondary_orf_length = 200
    minimal_orf_length = 50
    strand_specific = true

.. _chimera_splitting_configuration:

Parameters regarding splitting of chimeras
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section determines how Mikado will deal with :ref:`chimeras <chimera_splitting_algorithm>`. These are the relevant parameters:

- *execute*: boolean. If set to *false*, Mikado will operate in the *nosplit* mode. If set to *true*, the choice of the mode will be determined by the other parameters.
- *skip*: this is list of input assemblies (identified by the label in prepare, :ref:`above <input_file_list>`) that will **never** have the transcripts split.

.. hint:: cDNAs, reference transcripts, and the like should end up in the "skip" category. These are, after all, transcripts
that are presupposed to be originated from a single RNA molecule and therefore without fusions.

- *blast_check*: boolean. Whether to execute the check on the BLAST hits. If set to *false*, Mikado will operate in the *split* mode, unless *execute* is set to *false* (execute takes precedence over the other parameters).
- *blast_params*: this section contains the settings relative to the *permissive*, *lenient* and *stringent* mode.

   * *evalue*: maximum evalue of a hit to be assigned to the transcript and therefore be considered.
   * *hsp_evalue*: maximum evalue of a hsp inside a hit to be considered for the analysis.
   * *leniency*: one of **LENIENT, PERMISSIVE, STRINGENT**. See above for definitions.
   * *max_target_seqs*: integer. when loading BLAST hits from the database, only the first N will be considered for analysis.
   * *minimal_hsp_overlap*: number between 0 and 1. This indicates the overlap that must exist between the HSP and the ORF for the former to be considered for the split.
   .. code section: splitting.py, lines ~152-170

   * *min_overlap_duplication*: in the case of tandem duplicated genes, a chimera will have two ORFs that share the same hits, but possibly in a peculiar way - the HSPs will insist on the same region of the *target* sequence. This parameter controls how much overlap counts as a duplication. The default value is of 0.9 (90%).

.. code-block:: toml

    [pick.chimera_split]
    # Parameters related to the splitting of transcripts in the presence of
    # two or more ORFs. Parameters:
    # - execute: whether to split multi-ORF transcripts at all. Boolean.
    # - blast_check: whether to use BLAST information to take a decision. See blast_params for details.
    # - blast_params: Parameters related to which BLAST data we want to analyse.
    blast_check = true
    execute = true
    skip = [false, false, false, false, false]

    [pick.chimera_split.blast_params]
    # Parameters for the BLAST check prior to splitting.
    # - evalue: Minimum evalue for the whole hit. Default: 1e-6
    # - hsp_evalue: Minimum evalue for any HSP hit (some might be discarded even if the whole hit is valid). Default: 1e-6
    # - leniency: One of 'STRINGENT', 'LENIENT', 'PERMISSIVE'. Default: STRINGENT
    # - max_target_seqs: maximum number of hits to consider. Default: 3
    # - minimal_hsp_overlap: minimum overlap of the ORF with the HSP (*not* reciprocal). Default: 0.8, i.e. 80%
    # - min_overlap_duplication: minimum overlap (in %) for two ORFs to consider them as target duplications. This means that if two ORFs have no HSPs in common, but the coverage of their disjoint HSPs covers more than this % of the length of the *target*, they represent most probably a duplicated gene.
    evalue = 1e-06
    hsp_evalue = 1e-06
    leniency = "PERMISSIVE"
    max_target_seqs = 3
    minimal_hsp_overlap = 0.5
    min_overlap_duplication = 0.8

Parameters regarding input and output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The "files" and "output_format" sections deal respectively with input files for the pick stage and with some basic settings for the GFF output. Options:

- *input*: input GTF file for the run. It should be the one generated by the prepare stage, ie the :ref:`out file of the prepare stage <prep-settings>`.
- *loci_out*: main output file. It contains the winning transcripts, separated in their own gene loci, in GFF3 format. It will also determine the prefix of the *metrics* and *scores* files for this step. See the :ref:`pick manual page for details on the output <pick-output>`.
- *log*: name of the log file. Default: mikado_pick.log
- *monoloci_out*: this optional output file will contain the transcripts that have been passed to the :ref:`monoloci phase <introduction>`. It will also determine the prefix of the *metrics* and *scores* files for this step. See the :ref:`pick manual page for details on the output <pick-output>`.
- *subloci_out*: this optional output file will contain the transcripts that have been passed to the :ref:`subloci phase <introduction>`. It will also determine the prefix of the *metrics* and *scores* files for this step. See the :ref:`pick manual page for details on the output <pick-output>`.

.. code-block:: toml

   [pick.files]
    # Input and output files for Mikado pick.
    # - gff: input GTF/GFF3 file. Default: mikado_prepared.gtf
    # - loci_out: output GFF3 file from Mikado pick. Default: mikado.loci.gff3
    # - subloci_out: optional GFF file with the intermediate subloci. Default: no output
    # - monoloci_out: optional GFF file with the intermediate monoloci. Default: no output
    # - log: log file for this step.
    output_dir = "."
    input = "mikado_prepared.gtf"
    loci_out = "mikado.loci.gff3"
    subloci_out = ""
    monoloci_out = ""
    log = "pick.log"

Parameters regarding the output format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Available parameters:

- *id_prefix*: prefix for all the final Mikado models. The ID will be <prefix>.<chromosome>G<progressive ID>.
- *report_all_orfs*: some Mikado models will have more than one ORF (unless pick is operating in the *split* mode). If this option is set to ``true``, Mikado will report the transcript multiple times, one for each ORF, using different progressive IDs (<model name>.orf<progressive ID>). By default, this option is set to False, and only the primary ORF is reported.
- *source*: prefix for the source field in the output files. Loci GFF3 will have "<prefix>_loci", subloci GFF3s will have "<prefix>_subloci", and monoloci will have "<prefix>_monoloci".

.. code-block:: toml
    [pick.output_format]
    # Parameters related to the output format.
    #   - source: prefix for the source field in the mikado output.
    #   - id_prefix: prefix for the ID of the genes/transcripts in the output
    source = "Mikado"
    id_prefix = "mikado"
    report_all_orfs = false


Generic parameters on the pick run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section deals with other parameters necessary for the run, such as the number of processors to use, but also more important algorithmic parameters such as how to recognise fragments.

Parameters:

- *exclude_cds*: whether to remove CDS/UTR information from the Mikado output. Default: *false*.
- *intron_range*: tuple that indicates the range of lengths in which most introns should fall. Transcripts with introns either shorter or longer than this interval will be potentially penalised, depending on the scoring scheme. For the paper, this parameter was set to a tuple of integers in which *98%* of the introns of the reference annotation were falling (ie cutting out the 1st and 99th percentiles).
- *shm*: boolean. In certain cases, especially when disk access is a severely limiting factor, it might make sense to copy a SQLite database into RAM before querying. If this parameter is set to *true*, Mikado will copy the SQLite database into a temporary file in RAM, and query it from there.
- *only_reference_update*:
- *check_references*:
- *single_thread*: boolean. If set to true, Mikado will completely disable multiprocessing. Useful mostly for debugging reasons.

.. warning:: the shared-memory options are available only on Linux platforms.

.. code-block:: toml

    [pick.run_options]
    # Generic run options.
    # - shm: boolean flag. If set and the DB is sqlite, it will be copied onto the /dev/shm faux partition
    # - exclude_cds: boolean flag. If set, the CDS information will not be printed in Mikado output. Default: false
    # - single_thread: boolean flag. If set, multithreading will be disabled - useful for profiling and debugging.
    shm = false
    exclude_cds = false
    intron_range = [60, 10000]
    only_reference_update = false
    check_references = false
    single_thread = false

Technical details
~~~~~~~~~~~~~~~~~

The configuration file obeys a specific JSON schema which can be found at :download:`Mikado/configuration/configuration_blueprint.json <configuration_blueprint.json>`. Every time a Mikado utility is launched, it checks the configuration file against the schema to validate it. The schema contains non-standard "Comment" and "SimpleComment" string arrays which are used at runtime to generate the comment strings in the YAML output.
