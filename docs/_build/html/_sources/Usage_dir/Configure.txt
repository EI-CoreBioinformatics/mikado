.. _SQLAlchemy: http://www.sqlalchemy.org/

.. _configure:

Mikado configure
================

This utility prepares the configuration file that will be used


.. _conf_anatomy:

Anatomy of the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The guide here describes all voices of the configuration file. However, the configuration created by default by ``mikado configure`` is much simpler

.. _db-settings:

Database settings
-----------------

This section deals with the database settings that will be necessary for the :ref:`serialisation <serialise>` and :ref:`picking <pick>` phases of the pipeline. By default, Mikado will use a `SQLite database <https://www.sqlite.org/>`_, but it currently also supports `MySQL <http://www.mysql.com/>`_ and `PostgreSQL <https://www.postgresql.org/>`_ through SQLAlchemy_. Fields:

* db: name of the database to use. In case the database is SQLite, this will be the database file, otherwise it will be the database *name*.
* dbtype: one of:
  * sqlite
  * mysql
  * postgresql
* dbhost: host where the database is located. **Required with MySQL and PostgreSQL**.
* dbuser: User of the database. **Required with MySQL and PostgreSQL**.
* dbpasswd: Database password. **Required with MySQL and PostgreSQL**.
* dbport: Port to access to the database. It defaults to the normal ports for the selected database.

.. code-block:: yaml

    db_settings:
      #  Settings related to DB connection. Parameters:
      #  db: the DB to connect to. Required. Default: mikado.db
      #  dbtype: Type of DB to use. Choices: sqlite, postgresql, mysql. Default: sqlite.
      #  dbhost: Host of the database. Unused if dbtype is sqlite. Default: localhost
      #  dbuser: DB user. Default:
      #  dbpasswd: DB password for the user. Default:
      #  dbport: Integer. It indicates the default port for the DB.
      db: /usr/users/ga002/venturil/workspace/mikado/docs/mikado.db
      dbhost: localhost
      dbpasswd: ''
      dbport: 0
      dbtype: sqlite
      dbuser: ''

.. _ref-settings:

Reference settings
------------------

This section of the configuration file deals with the reference genome. It specifies two fields:

* genome: the genome FASTA file. **Required**.
* genome_fai: FAI index of the genome. Used by :ref:`Mikado serialise <serialise>`, it can be inferred if left null.
* transcriptome: optional annotation file for the genome. Mikado currently ignores this field, but it is used by :ref:`Daijin <Daijin>` to guide some of the RNA-Seq assemblies.

.. code-block:: yaml

    reference:
      #  Options related to the reference genome.
      genome: ''
      genome_fai: ''
      transcriptome: ''

.. _prep-settings:

Settings for the prepare stage
------------------------------

This section of the configuration file deals with the :ref:`prepare stage of Mikado <prepare>`. It specifies the input files, their labels, and which of them are strand specific. The available fields are the following:

* canonical: this voice specifies the splice site donors and acceptors that are considered canonical for the species. By default, Mikado uses the canonical splice site (GT/AG) and the two semi-canonical pairs (GC/AG and AT/AC). Type: Array of two-element arrays, composed by two-letter strings.
* lenient: boolean value. If set to *false*, transcripts that only have non-canonical splice sites will be **removed** from the output.
* minimum_length: minimum length of the transcripts to be kept.
* procs: number of processors to be used.
* strand_specific: boolean. If set to *true*, **all** input assemblies will be treated as strand-specific, therefore keeping the strand of monoexonic fragments as it was.
* strip_cds: boolean. If set to *true*, the CDS features will be stripped off the input transcripts. This might be necessary for eg transcripts obtained through alignment with `GMAP <http://research-pub.gene.com/gmap/>`_ [GMAP]_.
* files: this sub-section is the most important, as it contains among other things the locations and labels for the input files. Voices:
    * gff: array of the input files, in GFF or GTF format. Please note that only CDS/exon/UTR features will be considered from these files.
    * labels: optional array of the labels to be assigned to the input files. If non-empty, *it must be of the same order and length of the gff array*, and be composed of unique elements.
    * log: name of the log file.
    * out: name of the GTF output file.
    * out_fasta: name of the corresponding output FASTA file.
    * output_dir: output directory. It will be created if it does not exist already.
    * strand_specific_assemblies: array of the names of the GFF/GTF files that are strand specific. **All the file names in this array must also appear in the gff array as well.**.

.. code-block:: yaml

    prepare:
      #  Options related to the input data preparation.
      #  - files: options relative to the input/output files.
      #  - procs: Number of processes to use.
      #  - strip_cds: whether to remove the CDS from the predictions during preparation.
      #  - lenient: if set to True, invalid transcripts will be only flagged and not removed.
      #  EXPERIMENTAL.
      #  - strand_specific: if set to True, transcripts will be assumed to be in the correct
      #  orientation, no strand flipping or removal
      #  - strand_specific_assemblies: array of input predictions which are to be considered
      #  as strand-specific.
      #    Predictions not in this list will be considered as non-strand-specific.
      #  - canonical: canonical splice sites, to infer the correct orientation.
      canonical:
      - - GT
        - AG
      - - GC
        - AG
      - - AT
        - AC
      files:
        #  Options related to the input and output files.
        #  - out: output GTF file
        #  - out_fasta: output transcript FASTA file
        #  - gff: array of input predictions for this step.
        #  - log: output log. Default: prepare.log
        #  - labels: labels to be associated with the input GFFs. Default: None.
        gff: []
        labels: []
        log: prepare.log
        out: mikado_prepared.gtf
        out_fasta: mikado_prepared.fasta
        output_dir: .
        strand_specific_assemblies: []
      lenient: false
      minimum_length: 200
      procs: 1
      single: false
      strand_specific: false
      strip_cds: false

.. _serialise-settings:

Settings for the serialisation stage
------------------------------------

This section of the configuration file deals with the :ref:`serialisation stage of Mikado <serialise>`. It specifies the location of the ORF BED12 files from TransDecoder, the location of the XML files from BLAST, the location of portcullis junctions, and other details important at run time. It has the following voices:

* discard_definition: boolean. This is used to specify whether we will use the ID or the definition of the sequences when parsing BLAST results. This is important when BLAST data might have a mock, local identifier for the sequence ("lcl|1") rather than its original ID.
* force: whether the database should be truncated and rebuilt, or just updated.
* max_objects: this parameter is quite important when running with a SQLite database. SQLite does not support caching on the disk before committing the changes, so that every change has to be kept in memory. This can become a problem for RAM quite quickly. On the other hand, committing is an expensive operation, and it makes sense to minimise calls as much as possible. This parameter specifies the maximum number of objects Mikado will keep in memory before committing them to the database. The default number, 100,000, should ensure that Mikado runs with less than 1GB memory. Increase it to potentially increase speed at the price of greater memory usage.
* max_regression: this parameter is a float comprised between 0 and 1. TransDecoder will sometimes output open ORFs even in the presence of an in-frame start codon. Mikado can try to "regress" along the ORF until it finds one such start codon. This parameter imposes how much Mikado will regress, in percentage of the cDNA length.
* max_target_seqs: equivalent to the BLAST+ parameter of the same name - it indicates the maximum number of discrete hits that can be assigned to one sequence in the database.
* procs: number of processors to use. Most important for serialising BLAST+ files.
* single_thread: boolean, if set to *true* it will forcibly disable multi-threading. Useful mostly for debugging purposes.
* files: this sub-section codifies the location of the input files for serialise. It contains the following voices:
    * junctions: array of locations of reliable junction files. These must be in BED12 format.
    * log: log file.
    * orfs: array of locations of ORFs location on the cDNA, as created by eg TransDecoder [Trinity]_.
    * output_dir: output directory where the log file and the SQLite database will be written to (if SQLite has been chosen as the database type)
    * transcripts: input transcripts. This should be set to be equal to the output of :ref:`Mikado prepare <prepare>`, ie the "out_fasta" field of the :ref:`prepare section of the configuration file <prep-settings>`.
    * xml: this array indicates the location of the BLAST output file. Elements of the array can be:
       * BLAST+ XML files (optionally compressed with gzip)
       * BLAST+ ASN files (optionally compressed with gzip), which will be converted in-memory using ``blast_formatter``
       * a folder containing files of the above types.

.. code-block:: yaml

    serialise:
      #  Options related to serialisation
      #  - force: whether to drop and reload everything into the DB
      #  - files: options related to input files
      #  - max_objects: Maximum number of objects to keep in memory while loading data
      #  into the database
      #  - max_regression: if the ORF lacks a valid start site, this percentage indicates
      #  how far
      #    along the sequence Mikado should look for a good start site. Eg. with a value
      #  of 0.1,
      #    on a 300bp sequence with an open ORF Mikado would look for an alternative in-frame
      #  start codon
      #    in the first 30 bps (10% of the cDNA).
      #  - max_target_seqs: equivalently to BLAST, it indicates the maximum number of
      #  targets to keep
      #    per blasted sequence.
      #  - discard_definition: Boolean. Used to indicate whether Mikado should use the
      #  definition
      #    rather than the ID for BLAST sequences. Necessary as in some instances BLAST
      #  XMLs will have
      #    a mock identifier rather than the original sequence ID (eg lcl|1). Default:
      #  false.
      #  - procs: Number of processors to use. Default: 1.
      #  - single_thread: if true, Mikado prepare will force the usage of a single thread
      #  in this step.
      discard_definition: false
      files:
        blast_targets:
        - ''
        junctions: []
        log: serialise.log
        orfs:
        - ''
        output_dir: .
        transcripts: mikado_prepared.fasta
        xml:
        - ''
      force: false
      max_objects: 100000
      max_regression: 0
      max_target_seqs: 100000
      procs: 1
      single_thread: false

.. hint:: The most expensive operation in a "Mikado serialise" run is by far the serialisation of the BLAST files. Splitting the input files in multiple chunks, and analysing them separately, allows Mikado to parallelise the analysis of the BLAST results. If a single monolythic XML/ASN file is produced, by contrast, Mikado will be quite slow as it will have to parse it all.

.. _misc-settings:

Miscellanea
-----------

It is possible to set high-level settings for the logs in the ``log_settings`` section:

* log_level: level of the logging for Mikado. Options: *DEBUG, INFO, WARNING, ERROR, CRITICAL*. By default, Mikado will be quiet and output log messages of severity *WARNING* or greater.
* sql_level: level of the logging for messages regarding the database connection (through `SQLAlchemy`_). By default, SQLAlchemy will be set in quiet mode and asked to output only messages of severity *WARNING* or greater.

.. warning:: Mikado and SQLAlchemy can be greatly verbose if asked to output *DEBUG* or *INFO* messages, to the point of slowing down the program significantly due to the amount of writing to disk. Please consider setting the level to *DEBUG* only when there is a real problem to debug, not otherwise!

.. code-block:: yaml

    log_settings:
      #  Settings related to the logs. Keys:
      #  - sql_level: verbosity for SQL calls. Default: WARNING.
      #    In decreasing order: DEBUG, INFO, WARNING, ERROR, CRITICAL
      #  - log_level: verbosity. Default: WARNING.
      #    In decreasing order: DEBUG, INFO, WARNING, ERROR, CRITICAL
      log_level: WARNING
      sql_level: WARNING

It is also possible to set the type of multiprocessing method that should be used by Python3. The possible choices are "fork", "spawn", and "fork-server".

.. code-block:: yaml

    multiprocessing_method: spawn

.. warning:: Some schedulers, in particular SLURM, are not capable to understand that the processes *forked* by Python are still sharing the same memory with the main process, and think instead that each process is using that memory in isolation. As a result, they might think that the Mikado process is using its memory multiplied by the number of processes - depending on when the forking happens - and therefore shut down the program as it *appears* to be using much more memory than needed. For this reason, :ref:`Daijin <Daijin>` forces Mikado to run in **spawn** mode. Although spawning is slower than forking, it happens only once per run, and it has therefore a limited cost in terms of runtime - while greatly reducing the chances of the program being shut down because of "Out of memory" reasons.

Technical details
~~~~~~~~~~~~~~~~~

The configuration file obeys a specific JSON schema which can be found at :download:`Mikado/configuration/configuration_blueprint.json <configuration_blueprint.json>`. Every time a Mikado utility is launched, it checks the configuration file against the schema to validate it. The schema contains non-standard "Comment" and "SimpleComment" string arrays which are used at runtime to generate the comment strings in the YAML output.
