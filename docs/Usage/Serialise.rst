.. _Portcullis: https://github.com/maplesond/portcullis
.. _TopHat2: http://ccb.jhu.edu/software/tophat/index.shtml
.. _TransDecoder: http://transdecoder.github.io/
.. _BED12: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _STAR: https://github.com/alexdobin/STAR
.. _SQLalchemy: http://www.sqlalchemy.org/
.. _Prodigal: https://github.com/hyattpd/Prodigal

.. _serialise:

Mikado serialise
================

Mikado integrates data from multiple sources to select the best transcripts. During this step, these sources are brought together inside a common database, simplifying the process of retrieving them at runtime. Currently, Mikado integrates three or more different types of data:

.. _BED12-sidebar:
.. sidebar:: BED12_ files

    During serialisation, Mikado interprets BED12 files as having the relevant information in the **thick-start**, **thick-end** fields, ie columns 7 and 8. For ORF files, this means that the CDS should start at the 7th column and end at the 8th (stop-codon exclusive, like in TransDecoder_) - or viceversa for ORFs on the negative strand. For junction files, it means that the reliable junction, ie the reliable *intron*, starts at the thick-start position and ends at the thick-end position. This format requirement is the opposite of what happens for example in the junctions produced by `TopHat <http://ccb.jhu.edu/software/tophat/index.shtml>`_, where the BED12 file lists the two *exonic* fragments, rather than the intron. STAR_ instead provides a tabular non-standard file indicating the most reliable junctions in the alignment. Portcullis_ provides utilities to convert such files into the format used by Mikado, and to merge together multiple junction files. However, we recommend running Portcullis directly on the alignments, rather than using the junctions indicated by the programs themselves.

#. reliable junctions, as detected by Portcullis_, in BED12_ format
#. ORF data, currently identified using TransDecoder_, but any program capable of generating it in BED12_ format is suitable.
#. BLASTX [Blastplus]_ data in XML format
#. Miscellaneous numeric scores can also be given as input to Mikado as tabular file, with one row per transcript.

After serialisation in the database, these sources will be available to use for any subsequent Mikado pick run. Having the data present in the database allows to run Mikado with multiple configurations and little overhead in terms of pre-loading data; this feature is taken directly advantage of in :ref:`Daijin`, where it is possible to run Mikado using multiple modes.

Mikado serialise can use three different SQL databases as backends - SQLite, MySQL and PostgreSQL - thanks to SQLAlchemy_.
This step, together with the creation of the TransDecoder and BLAST data, is the most time consuming of the pipeline. In particular, although Mikado serialise will try to analyse the XML data in a parallelised fashion if so instructed, the insertion of the data in the database will still happen in a single thread and will therefore be of limited speed. If using SQLite as database (the default option), it is possible to decrease the runtime by modifying the "max_objects" parameters, at the cost however of increased RAM usage.

.. important:: The schema of Mikado databases changed with version 1.3. Any database created prior to this version **should be regenerated**, otherwise Mikado pick will fail.

Transdecoder ORFs
~~~~~~~~~~~~~~~~~

When Mikado analyses ORFs produced by TransDecoder_, Prodigal_ or equivalent program, it performs additionally the following checks:

#. Check the congruence between the length of the transcript in the BED12 file and that found in the FASTA file
#. Check that the ORF does not contain internal stop codons
#. Check that the CDS length is valid, ie a multiple of 3, if the ORF is complete
#. Optionally, if the ORF is open on the 5' side, Mikado can try to find an internal start codon. See :ref:`this section <max-regression>` for details.


Usage
~~~~~

``mikado serialise`` allows to override some of the parameters present in the configuration file through command line options, eg the input files. Notwithstanding, in the interest of reproducibility we advise to configure everything through the configuration file and supply it to Mikado prepare without further modifications.

Available parameters:

* Parameters related to performance:

    - *start-method*: one of fork, spawn, forkserver. It determines the multiprocessing start method. By default, Mikado will use the default for the system (fork on UNIX, spawn on Windows).
    - *procs*: Number of processors to use.
    - *single-thread*: flag. If set, Mikado will disable all multithreading.
    - *max_objects*: Maximum number of objects to keep in memory before committing to the database. See :ref:`this section of the configuration <max-objects>` for details.
* Basic input data and settings:

    - *output-dir*: directory where the SQLite database and the log will be written to.
    - *transcripts*: these are the input transcripts that are present on the GTF file considered by Mikado. Normally this should be the output of Mikado prepare.
    - *genome_fai*: FAIDX file of the genome FASTA. If not given, serialise will derive it from the "reference: genome" field of the configuration.
    - *force*: flag. If set, and the database is already present, it will be truncated rather than updated.
    - **json-conf**: this is the configuration file created with :ref:`Mikado configure <configure>`.
    - *db*: if the database is specified on the command line, ``mikado serialise`` will interpret it as a **SQLite** database. This will overwrite any setting present in the configuration file.
* Parameters related to logging:

    - *log*: log file. It defaults to ``serialise.log``.
    - *log_level*: verbosity of the logging. Please be advised that excessive verbosity can negatively impact the performance of the program - the debug mode is extremely verbose.
* Parameters related to reliable junctions:

    - *junctions*: a BED12_ file of reliable junctions. This can be obtained using Portcullis_. Please see the relative :ref:`sidebar <BED12-sidebar>`.
* Parameters related to the treatment of ORF data:

    - *orfs*: ORF BED12 files, separated by comma.
    - *max-regression*: A percentage, expressed as a number between 0 and 1, which indicates how far can Mikado regress along the ORF to find a valid start codon. See the :ref:`relative section in the configuration <max-regression>` for details.

    - *codon-table*: this parameter specifies the codon table to use for the project. Mikado by default uses the NCBI codon table 1 (standard with eukaryotes) with the modification that only ATG is considered as a valid start codon, as ORF predictions usually inflate the number of non-standard starts.
* Parameters related to BLAST data:

    - *blast_targets*: BLAST FASTA database.
    - *discard-definition*: Flag. Depending on how the database has been created, sometimes BLAST will substitute the ID of the sequence with "lcl|" ids. Mikado circumvents this by looking for the definition field in the XML file. Using this flag will disable this behaviour and force Mikado to use the ID - with the potential of having a mismatch between the sequences in the BLAST DB and the sequences in the BLAST files.
    - *xml*: BLAST files to parse. This can be one of the following:

        + A list of XML BLAST files, optionally compressed with GZip or BZip2, comma separated (suffix .xml)
        + A list of ASN BLAST files, optionally compressed with GZip or BZip2, comma separated (suffix .asn)
        + A list of folders, comma separated, where it is possible to find files of the former 2 types
        + A mixture of the three above types.
    - *max-target-seqs*: maximum number of BLAST targets that can be loaded per sequence, for each BLAST alignment. Please note that if you align against multiple databases, this threshold will be applied once per file.

.. hint:: Mikado will parallelise only the reading of multiple XML files. As such, this part of the pipeline is less performing than the other steps.

.. warning:: It is advised to set this parameter to *spawn* even on UNIX. See :ref:`the dedicated sidebar for details <scheduler-multiprocessing>`.

Usage::

    $ mikado serialise --help
    usage: Mikado serialise [-h] [--start-method {fork,spawn,forkserver}]
                            [--orfs ORFS] [--transcripts TRANSCRIPTS]
                            [-mr MAX_REGRESSION] [--codon-table CODON_TABLE]
                            [--max_target_seqs MAX_TARGET_SEQS]
                            [--blast_targets BLAST_TARGETS] [--xml XML] [-p PROCS]
                            [--single-thread] [--genome_fai GENOME_FAI]
                            [--junctions JUNCTIONS]
                            [--external-scores EXTERNAL_SCORES] [-mo MAX_OBJECTS]
                            [-f] --json-conf JSON_CONF [-l [LOG]] [-od OUTPUT_DIR]
                            [-lv {DEBUG,INFO,WARN,ERROR}]
                            [db]

    optional arguments:
      -h, --help            show this help message and exit
      --start-method {fork,spawn,forkserver}
                            Multiprocessing start method.
      -od OUTPUT_DIR, --output-dir OUTPUT_DIR
                            Output directory. Default: current working directory

      --orfs ORFS           ORF BED file(s), separated by commas
      --transcripts TRANSCRIPTS
                            Transcript FASTA file(s) used for ORF calling and
                            BLAST queries, separated by commas. If multiple files
                            are given, they must be in the same order of the ORF
                            files. E.g. valid command lines are:
                            --transcript_fasta all_seqs1.fasta --orfs all_orfs.bed
                            --transcript_fasta seq1.fasta,seq2.fasta --orfs
                            orfs1.bed,orf2.bed --transcript_fasta all_seqs.fasta
                            --orfs orfs1.bed,orf2.bed These are invalid instead: #
                            Inverted order --transcript_fasta
                            seq1.fasta,seq2.fasta --orfs orfs2.bed,orf1.bed #Two
                            transcript files, one ORF file --transcript_fasta
                            seq1.fasta,seq2.fasta --orfs all_orfs.bed
      -mr MAX_REGRESSION, --max-regression MAX_REGRESSION
                            "Amount of sequence in the ORF (in %) to backtrack in
                            order to find a valid START codon, if one is absent.
                            Default: None
      --codon-table CODON_TABLE
                            Codon table to use. Default: 0 (ie Standard, NCBI #1,
                            but only ATG is considered a valid stop codon.

      --max_target_seqs MAX_TARGET_SEQS
                            Maximum number of target sequences.
      --blast_targets BLAST_TARGETS
                            Target sequences
      --xml XML             XML file(s) to parse. They can be provided in three
                            ways: - a comma-separated list - as a base folder -
                            using bash-like name expansion (*,?, etc.). In this
                            case, you have to enclose the filename pattern in
                            double quotes. Multiple folders/file patterns can be
                            given, separated by a comma.
      -p PROCS, --procs PROCS
                            Number of threads to use for analysing the BLAST
                            files. This number should not be higher than the total
                            number of XML files.
      --single-thread       Force serialise to run with a single thread,
                            irrespective of other configuration options.

      --genome_fai GENOME_FAI
      --junctions JUNCTIONS

      --external-scores EXTERNAL_SCORES
                            Tabular file containing external scores for the
                            transcripts. Each column should have a distinct name,
                            and transcripts have to be listed on the first column.

      -mo MAX_OBJECTS, --max-objects MAX_OBJECTS
                            Maximum number of objects to cache in memory before
                            committing to the database. Default: 100,000 i.e.
                            approximately 450MB RAM usage for Drosophila.
      -f, --force           Flag. If set, an existing databse will be deleted
                            (sqlite) or dropped (MySQL/PostGreSQL) before
                            beginning the serialisation.
      --json-conf JSON_CONF
      -l [LOG], --log [LOG]
                            Optional log file. Default: stderr
      -lv {DEBUG,INFO,WARN,ERROR}, --log_level {DEBUG,INFO,WARN,ERROR}
                            Log level. Default: derived from the configuration; if
                            absent, INFO
      db                    Optional output database. Default: derived from
                            json_conf


Technical details
~~~~~~~~~~~~~~~~~

The schema of the database is quite simple, as it is composed only of 9 discrete tables in two groups. The first group, *chrom* and *junctions*, serialises the information pertaining to the reliable junctions - ie information which is not relative to the transcripts but rather to their genomic locations.
The second group serialises the data regarding ORFs, BLAST files and external arbitrary data. The need of using a database is mainly driven by the latter, as querying a relational database is faster than retrieving the information from the XML files themselves at runtime.

.. database figure generated with `SchemaCrawler <https://github.com/schemacrawler/SchemaCrawler>`_, using the following command line:
    schemacrawler --command=schema --url=jdbc:sqlite:sample_data/daijin/5-mikado/mikado.db -o docs/Usage/database_schema.png --outputformat=png --info-level=maximum

.. topic:: Database schema used by Mikado.

    .. figure:: database_schema.png
        :align: center
        :scale: 50%
