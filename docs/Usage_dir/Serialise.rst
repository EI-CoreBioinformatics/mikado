.. _Portcullis: https://github.com/maplesond/portcullis
.. _TopHat2: http://ccb.jhu.edu/software/tophat/index.shtml
.. _TransDecoder: http://transdecoder.github.io/
.. _BED12: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _STAR: https://github.com/alexdobin/STAR

.. _serialise:

Mikado serialise
================

Mikado integrates data from multiple sources to select the best transcripts. During this step, these sources are brought together inside a common database, simplifying the process of retrieving them at runtime. Currently, Mikado integrates three different types of data:

.. sidebar:: BED12_ files

    During serialisation, Mikado interprets BED12 files as having the relevant information in the **thick-start**, **thick-end** fields, ie columns 7 and 8. For ORF files, this means that the CDS should start at the 7th column and end at the 8th (stop-codon exclusive, like in TransDecoder_) - or viceversa for ORFs on the negative strand. For junction files, it means that the reliable junction, ie the reliable *intron*, starts at the thick-start position and ends at the thick-end position. This format requirement is the opposite of what happens for example in the junctions produced by `TopHat <http://ccb.jhu.edu/software/tophat/index.shtml>`_, where the BED12 file lists the two *exonic* fragments, rather than the intron. STAR_ instead provides a tabular non-standard file indicating the most reliable junctions in the alignment. Portcullis_ provides utilities to convert such files into the format used by Mikado, and to merge together multiple junction files. However, we recommend running Portcullis directly on the alignments, rather than using the junctions indicated by the programs themselves.

#. reliable junctions, as detected by Portcullis_, in BED12_ format (
#. ORF data, currently identified using TransDecoder_, but any program capable of generating it in BED12_ format is suitable.
#. BLASTX [Blastplus]_ data in XML format

After serialisation in the database, these sources will be available to use for any subsequent Mikado pick run. Having the data present in the database allows to run Mikado with multiple configurations and little overhead in terms of pre-loading data; this feature is taken directly advtange of in :ref:`Daijin`, where it is possible to run Mikado using multiple modes.

Mikado serialise can use three different SQL databases as backends - SQLite, MySQL and PostgreSQL - thanks to SQLAlchemy_.
This step, together with the creation of the TransDecoder and BLAST data, is the most time consuming of the pipeline. In particular, although Mikado serialise will try to analyse the XML data in a parallelised fashion if so instructed, the insertion of the data in the database will still happen in a single thread and will therefore be of limited speed. If using SQLite as database (the default option), it is possible to decrease the runtime by modifying the "max_objects" parameters, at the cost however of increased RAM usage.


Usage
~~~~~

``mikado serialise`` allows to override some of the parameters present in the configuration file through command line options, eg the input files. Notwithstanding, in the interest of reproducibility we advise to configure everything through the configuration file and supply it to Mikado prepare without further modifications.

Available parameters:

* *start-method*: one of fork, spawn, forkserver. It determines the multiprocessing start method. By default, Mikado will use the default for the system (fork on UNIX, spawn on Windows).
* *output-dir*: directory where the SQLite database and the log will be written to.
* *orfs*: ORF BED12 files, separated by comma.
* *transcripts*: these are the input transcripts that are present on the GTF file considered by Mikado. Normally this should be the output of Mikado prepare.
* *max-regression*: TransDecoder will greedily try to find the longest ORF in any transcript, even when this means creating an incomplete ORF at the 5' site ... inclusive of a valid internal start codon. Mikado can try to identify those cases by backtracking along the ORF and truncating it on the 5' side when it finds a valid START position. This percentage value, expressed as a number between 0 and 1, indicates how much bps Mikado is allowed to backtrack before considering the ORF as truly open on the 5'. By default, this behaviour is disabled.
* *max-target-seqs*: maximum number of targets 


.. warning:: It is advised to set this parameter to *spawn* even on UNIX. See :ref:`the dedicated sidebar for details <scheduler-multiprocessing>`.

Usage::

    $ mikado serialise --help
    usage: Mikado serialise [-h] [--start-method {fork,spawn,forkserver}]
                            [--orfs ORFS] [--transcripts TRANSCRIPTS]
                            [-mr MAX_REGRESSION]
                            [--max_target_seqs MAX_TARGET_SEQS]
                            [--blast_targets BLAST_TARGETS] [--discard-definition]
                            [--xml XML] [-p PROCS] [--single-thread]
                            [--genome_fai GENOME_FAI] [--junctions JUNCTIONS]
                            [-mo MAX_OBJECTS] [-f] --json-conf JSON_CONF
                            [-l [LOG]] [-od OUTPUT_DIR]
                            [-lv {DEBUG,INFO,WARN,ERROR}]
                            [db]

    Mikado serialise creates the database used by the pick program. It handles
    Junction and ORF BED12 files as well as BLAST XML results.

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

      --max_target_seqs MAX_TARGET_SEQS
                            Maximum number of target sequences.
      --blast_targets BLAST_TARGETS
                            Target sequences
      --discard-definition  Flag. If set, the sequences IDs instead of their
                            definition will be used for serialisation.
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
                            Log level. Default: INFO
      db                    Optional output database. Default: derived from
                            json_conf



Technical details
~~~~~~~~~~~~~~~~~

The schema of the database is quite simple, as it is composed only of 7 discrete tables in two groups. The first group, *chrom* and *junctions*, serialises the information pertaining to the reliable junctions - ie information which is not relative to the transcripts but rather to their genomic locations.
The second group serialises the data regarding ORFs and BLAST files. The need of using a database is mainly driven by the latter, as querying a relational database is faster than retrieving the information from the XML files themselves at runtime.

.. database figure generated with `SchemaCrawler <http://sualeh.github.io/SchemaCrawler/>`_, using the following command line:
    schemacrawler -c graph -url=jdbc:sqlite:sample_data/mikado.db -o docs/Usage_dir/database_schema.png --outputformat=png -infolevel=maximum

.. topic:: Database schema used by Mikado.

    .. figure:: database_schema.png
        :align: center
        :scale: 50%
