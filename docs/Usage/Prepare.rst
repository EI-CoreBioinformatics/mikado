.. _prepare:

Mikado prepare
==============

This is the first executive step of the Mikado pipeline. It will accomplish the following goals:

#. Collect annotations from disparate annotation files.
#. Remove redundant assemblies, ie, assemblies that are *identical* across the various input files.
#. Determine the strand of the transcript junctions.
#. Ensure uniqueness of the transcript names.
#. Order the transcript by locus.
#. Extract the transcript sequences.

Usage
~~~~~

``Mikado prepare`` allows to override some of the parameters present in the configuration file through command line options, eg. the input files. Notwithstanding, in the interest of reproducibility we advise to configure everything through the configuration file and supply it to Mikado prepare without further modifications.

Available parameters:

* *json-conf*: the most important parameter. This is the configuration file created through :ref:`Mikado configure <configure>`.
* *fasta*: reference genome. Required, either through the command line or through the configuration file.
* *out*: Output GTF file, with the collapsed transcripts.
* *out_fasta*: Output FASTA file of the collapsed transcripts.
* *start-method*: :ref:`multiprocessing start method <start-methods>`.
* *verbose*, *quiet*: flags to set the verbosity of Mikado prepare. It is generally not advised to turn the verbose mode on, unless there is a problem to debug, given the verbosity of the output.
* *strand-specific*: If set, all assemblies will be treated as strand-specific.
* *strand-specific-assemblies*: comma-separated list of strand specific assemblies.
* *list*: in alternative to specifying all the information on the command line, it is possible to give to Mikado a *tab-separated* file with the following contents:

   #. Location of the file
   #. label for the file
   #. whether that assembly is strand-specific or not (write True/False)
   #. Optionally, a bonus/malus to be associated to transcripts coming from that assembly.
* *log*: log file. Optional, by default Mikado will print to standard error.
* *lenient*: flag. If set, transcripts without any canonical splice site will be output as well. By default, they would be discarded.
* *single*: flag that disables multiprocessing. Mostly useful for debug purposes.
* *strip-cds*: some aligners (eg GMAP) will try calculate a CDS on the fly for alignments. Use this flag to discard such CDS sections and retain only the cDNA information.
* *minimum_length*: minimum length of the transcripts to be kept.

Command line usage:

.. code-block:: bash

    $ mikado prepare --help
    usage: Mikado prepare [-h] [--fasta FASTA] [-v | -q]
                          [--start-method {fork,spawn,forkserver}]
                          [-s | -sa STRAND_SPECIFIC_ASSEMBLIES] [--list LIST]
                          [-l LOG] [--lenient] [-m MINIMUM_LENGTH] [-p PROCS]
                          [-scds] [--labels LABELS] [--single] [-od OUTPUT_DIR]
                          [-o OUT] [-of OUT_FASTA] [--json-conf JSON_CONF] [-k]
                          [gff [gff ...]]
    
    positional arguments:
      gff                   Input GFF/GTF file(s).
    
    optional arguments:
      -h, --help            show this help message and exit
      --fasta FASTA         Genome FASTA file. Required.
      -v, --verbose
      -q, --quiet
      --start-method {fork,spawn,forkserver}
                            Multiprocessing start method.
      -s, --strand-specific
                            Flag. If set, monoexonic transcripts will be left on
                            their strand rather than being moved to the unknown
                            strand.
      -sa STRAND_SPECIFIC_ASSEMBLIES, --strand-specific-assemblies STRAND_SPECIFIC_ASSEMBLIES
                            Comma-delimited list of strand specific assemblies.
      --list LIST           Tab-delimited file containing rows with the following format:
<file>  <label> <strandedness, def. False> <score(optional, def. 0)> <is_reference(optional, def. False)> <keep_redundant(optional, def. True)>
strandedness, is_reference and keep_redundant must be boolean values (True, False)
score must be a valid floating number.
      -l LOG, --log LOG     Log file. Optional.
      --lenient             Flag. If set, transcripts with only non-canonical
                            splices will be output as well.
      -m MINIMUM_LENGTH, --minimum_length MINIMUM_LENGTH
                            Minimum length for transcripts. Default: 200 bps.
      -p PROCS, --procs PROCS
                            Number of processors to use (default None)
      -scds, --strip_cds    Boolean flag. If set, ignores any CDS/UTR segment.
      --labels LABELS       Labels to attach to the IDs of the transcripts of the
                            input files, separated by comma.
      --single              Disable multi-threading. Useful for debugging.
      -od OUTPUT_DIR, --output-dir OUTPUT_DIR
                            Output directory. Default: current working directory
      -o OUT, --out OUT     Output file. Default: mikado_prepared.gtf.
      -of OUT_FASTA, --out_fasta OUT_FASTA
                            Output file. Default: mikado_prepared.fasta.
      --json-conf JSON_CONF
                            Configuration file.
      -k, --keep-redundant  Boolean flag. If invoked, Mikado prepare will retain
                            redundant models.


Collection of transcripts from the annotation files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Different assemblers will produce data in different formats, typically in GFF or GTF format, and not necessarily in the same order (if any is present). Mikado will serialise the transcripts from these files and port them all into a standard GTF format. Moreover, it will ensure that each transcript ID appears only once across the input files. The optional labels provided for each file will be attached to the transcript names as prefixes, and used as the source field in the output GTF, to ensure the uniqueness of each transcript name.
If two or more transcripts are found to be identical, only one will be retained, chosen at random among all the possibilities.
In addition to this, Mikado prepare will also sort the transcripts by coordinate, irrespective of strand, so that they are suitably displayed for the divide-et-impera algorithm of :ref:`Mikado pick <pick>`.

When two or more identical transcripts are present in a locus, Mikado will use the (optionally provided) source score to select the *a priori* best assembly amongst the choices.
For example, if a mikado prepare run comprises both PacBio reads and Illumina assemblies and the experimenter has given a score of 1 or more to the former dataset but not the latter, if a PacBio read is present together with a stringtie assembly, the PacBio read will always be selected over the StringTie.
Please note that this "score-based" selection ***only operates for transcripts that are identical**. No other selection is performed at this stage.

.. warning:: To be considered *identical*, two transcripts must match down to the last base pair. A simple match or containment of the intron chain will not suffice. This is because using the cDNA data alone it is difficult to understand whether the longer form(s) is the correct assembly rather than a chimera or a trans-splice event.
.. note:: From version 1.3 onwards, Mikado considers the CDS as well when performing the redundancy check. So, two transcripts having the same coordinates but different CDS (because of non-overlapping ORFs or disagrement on the frame and/or start codon position) will be kept as non-redundant.
.. note:: Transcripts that are considered to come from a "reference" assembly are never going to be excluded, and will always be prioritised over other assemblies.


Check on strand correctness
---------------------------

During its run, Mikado prepare will also check the correctness of the transcripts. In particular:

* Unless the assembly is marked as strand-specific, any monoexonic transcript will have its strand *removed*.
* If a transcript contains canonical splice junctions on **both** strands, it will be completely removed
* If a transcript contains only non-canonical splice junctions, it will be removed *unless* the "lenient" option is specified either at the command line or in the configuration file.

The couples of splice acceptors and donors which are considered as canonical :ref:`can be specified in the configuration file <canonical-configuration>`. By default, Mikado will consider as canonical both properly canonical splicing event (GT-AG) as well as the semi-canonical events (GC-AG, AT-AC). Any other couple will be considered as non-canonical.

.. warning:: Mikado will check the strand of each junction inside a transcript *independently*. Therefore, if a transcript with 9 junctions on the plus strand is found to have a non-canonical splicing junction **which happens to be the reverse of a canonical one** (eg. CT-AC), it will deem this junction as misassigned to the wrong strand and flip it to the minus strand. In this example, the transcript will therefore be **considered as an error** as it contains both + and - junctions, and discarded.

.. note:: Starting from Mikado version **1.3**, transcripts can be tagged as being from an assembly of "reference" quality. This implies that:

* A transcript which is marked as “reference” will never have its CDS stripped
* A transcript which is marked as “reference” will never be marked for removal due to redundancy, even if there are multiple copies of it, or if other assemblies with a higher score have identical transcripts (normally only one transcript would be retained, and that would be chosen amongst the highest scoring assemblies)
* A transcript which is marked as reference will never have its strand removed or flipped.

Please see the :ref:`configuration help page <configure>` for details.


Output files
------------

Mikado prepare will produce two files:

* a *sorted* GTF file, containing all the transcripts surviving the checks
* a FASTA file of the transcripts, in the proper cDNA orientation.

.. warning:: contrary to other tools such as eg gffread from Cufflinks [Cufflinks]_, Mikado prepare will **not** try to calculate the loci for the transcripts. This task will be performed later in the pipeline. As such, the GTF file is formally incorrect, as multiple transcripts in the same locus but coming from different assemblies will *not* have the same gene_id but rather will have kept their original one. Moreover, if two gene_ids were identical but discrete in the input files (ie located on different sections of the genome), this error will not be corrected. If you desire to use this GTF file for any purpose, please use a tool like gffread to calculate the loci appropriately.
