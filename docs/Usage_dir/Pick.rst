.. _pick:

Mikado pick
===========

This is the final stage of the pipeline, and the one which actually applies the :ref:`Mikado algorithm <Introduction>` to identify the loci and select the best transcripts.


Usage
~~~~~

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


