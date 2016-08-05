.. _pick:

Mikado pick
===========

This is the final stage of the pipeline, and the one which actually applies the :ref:`core algorithms <Introduction>` to identify the loci and select the best transcripts.
As such, it is the




Usage
~~~~~

``mikado pick`` allows to modify some of the parameters regarding the run at runtime. However, some sections - such as most of the settings regarding alternative splicing - are left untouched by the utility, and are best modified by editing the :ref:`configuration file itself <configure>`. The available parameters are as follows:

* *json-conf*: required. This is the configuration file created in the :ref:`first step <configure>` of the pipeline.
* *gff*; optionally, it is possible to point Mikado prepare to the GTF it should use here on the command line. This file should be the output of :ref:```mikado prepare`` <prepare>`. Please note that this file should be in GTF format, sorted by chromosome and position; if that is not the case, Mikado will fail.
* Options regarding the resources to be used:
    * *procs*: number of processors to use.
    * *start-method*: multiprocessing start method. See the :ref:`explanation on Python multiprocessing <scheduler-multiprocessing>`
    * *single*: flag. If present, multiprocessing will be disabled.
    * *shared-memory*: flag, available on Unix systems only. If set, Mikado will try to copy the SQLite database in RAM. It might provide a small boost if disk access is a limitation. Ineffective with databases other than SQLite.
    * *shared-memory-db*: if the database has already been copied in memory, its new location can be given with this argument. Useful to have multiple Mikado picks share the same DB.
    * *preload*: flag. If present, Mikado will load the database in memory before execution. **Discouraged unless the database size is quite small.**
* Options regarding logging:
    * *log*: name of the log file. By default, "pick.log"
    * *verbose*: sets the log level to DEBUG. Please be advised that the debug mode is **extremely** verbose and is bestly invoked only for real, targeted debugging sessions.
    * *noverbose*: sets the log level to ERROR. In most cases, the log file will be empty.
    * *log-level*: directly sets the log level. One of the logging module available mode:



Usage::

    $ mikado pick --help
    usage: Mikado pick [-h] [--start-method {fork,spawn,forkserver}] [-p PROCS]
                       --json-conf JSON_CONF [-i INTRON_RANGE INTRON_RANGE]
                       [--subloci_out SUBLOCI_OUT] [--monoloci_out MONOLOCI_OUT]
                       [--loci_out LOCI_OUT] [--prefix PREFIX] [--no_cds]
                       [--source SOURCE] [--flank FLANK] [--purge] [-shm]
                       [-shmdb SHM_DB] [--preload] [-db SQLITE_DB]
                       [-od OUTPUT_DIR] [--single] [-l LOG] [-v | -nv]
                       [-lv {DEBUG,INFO,WARN,ERROR,CRITICAL}]
                       [--mode {nosplit,stringent,lenient,permissive,split}]
                       [gff]

    Mikado pick analyses a sorted GTF/GFF files in order to identify its loci and
    choose the best transcripts according to user-specified criteria. It is
    dependent on files produced by the "prepare" and "serialise" components.

    positional arguments:
      gff

    optional arguments:
      -h, --help            show this help message and exit
      --start-method {fork,spawn,forkserver}
                            Multiprocessing start method. (default: None)
      -p PROCS, --procs PROCS
                            Number of processors to use. Default: look in the
                            configuration file (1 if undefined) (default: None)
      --json-conf JSON_CONF
                            JSON/YAML configuration file for scoring transcripts.
                            (default: None)
      -i INTRON_RANGE INTRON_RANGE, --intron-range INTRON_RANGE INTRON_RANGE
                            Range into which intron lengths should fall, as a
                            couple of integers. Transcripts with intron lengths
                            outside of this range will be penalised. Default: (60,
                            900) (default: None)
      --subloci_out SUBLOCI_OUT
      --monoloci_out MONOLOCI_OUT
      --loci_out LOCI_OUT   This output file is mandatory. If it is not specified
                            in the configuration file, it must be provided here.
                            (default: None)
      --prefix PREFIX       Prefix for the genes. Default: Mikado (default: None)
      --no_cds              Flag. If set, not CDS information will be printed out
                            in the GFF output files. (default: None)
      --source SOURCE       Source field to use for the output files. (default:
                            None)
      --flank FLANK         Flanking distance (in bps) to group non-overlapping
                            transcripts into a single superlocus. Default:
                            determined by the configuration file. (default: None)
      --purge               Flag. If set, the pipeline will suppress any loci
                            whose transcripts do not pass the requirements set in
                            the JSON file. (default: False)
      -shm, --shared-memory
                            Flag. If set, the DB will be copied into memory.
                            (default: False)
      -shmdb SHM_DB, --shared-memory-db SHM_DB
                            Name of the shared memory DB. WARNING: if set, the DB
                            copy will be persistently copied into memory, so that
                            multiple pickers can share. (default: None)
      --preload             Flag. If set, the Mikado DB will be pre-loaded into
                            memory for faster access. WARNING: this option will
                            increase memory usage and the preloading might be
                            quite slow. (default: False)
      -db SQLITE_DB, --sqlite-db SQLITE_DB
                            Location of an SQLite database to overwrite what is
                            specified in the configuration file. (default: None)
      -od OUTPUT_DIR, --output-dir OUTPUT_DIR
                            Output directory. Default: current working directory
                            (default: None)
      --single              Flag. If set, Creator will be launched with a single
                            process. Useful for debugging purposes only. (default:
                            False)
      --mode {nosplit,stringent,lenient,permissive,split}
                            Mode in which Mikado will treat transcripts with
                            multiple ORFs. - nosplit: keep the transcripts whole.
                            - stringent: split multi-orf transcripts if two
                            consecutive ORFs have both BLAST hits and none of
                            those hits is against the same target. - lenient:
                            split multi-orf transcripts as in stringent, and
                            additionally, also when either of the ORFs lacks a
                            BLAST hit (but not both). - permissive: like lenient,
                            but also split when both ORFs lack BLAST hits - split:
                            split multi-orf transcripts regardless of what BLAST
                            data is available. (default: None)

    Log options:
      -l LOG, --log LOG     File to write the log to. Default: decided by the
                            configuration file. (default: None)
      -v, --verbose         Flag. If set, the debug mode will be activated.
                            (default: False)
      -nv, --noverbose      Flag. If set, the log will report only errors and
                            critical events. (default: False)
      -lv {DEBUG,INFO,WARN,ERROR,CRITICAL}, --log-level {DEBUG,INFO,WARN,ERROR,CRITICAL}
                            Logging level. Default: retrieved by the configuration
                            file. (default: None)

.. block end


