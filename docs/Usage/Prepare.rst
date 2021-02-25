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

* *configuration*: the most important parameter. This is the configuration file created through :ref:`Mikado configure <configure>`.
* *fasta*: reference genome. Required, either through the command line or through the configuration file.
* *out*: Output GTF file, with the collapsed transcripts.
* *out_fasta*: Output FASTA file of the collapsed transcripts.
* *start-method*: :ref:`multiprocessing start method <start-methods>`.
* *verbose*, *quiet*: flags to set the verbosity of Mikado prepare. It is generally not advised to turn the verbose mode on, unless there is a problem to debug, given the verbosity of the output.
* *strand-specific*: If set, all assemblies will be treated as strand-specific.
* *strand-specific-assemblies*: comma-separated list of strand specific assemblies.
* *strip-cds*: some aligners (eg GMAP) will try calculate a CDS on the fly for alignments. Use this flag to discard all CDS information from all input transcripts.
* *exclude-redundant*: if set, this flag instructs Mikado to look for and simplify redundant intron chains. By default, this option is disabled, or enabled on a per-sample basis. See :ref:`this section for an explanation of redundancy removal in Mikado <redundant-transcripts-in-prepare>`, and the section on the list input files below for an explanation on how to set this value on a per-sample basis (recommended).
* *codon-table*: Mikado prepare will check the validity of the ORFs of input models. This value indicates which codon table Mikado should use for this purpose. See the section on :ref:`the checks on CDSs <orf-check-prepare>`.
* *strip-faulty-cds*: when encountering a transcript with an invalid ORF due to e.g. in-frame stop codons, Mikado will usually discard the whole transcript. If this flag is set, Mikado will instead remove the CDS information and leave the transcript in place.
* *list*: in alternative to specifying all the information on the command line, it is possible to give to Mikado a *tab-separated* file with details of the files to use. See :ref:`this section for details <input_file_list>`.
* *log*: log file. Optional, by default Mikado will print to standard error.
* *lenient*: flag. If set, multiexonic transcripts without any canonical splice site will be output as well. By default, they would be discarded.
* *minimum-cdna-length*: minimum length of the transcripts to be kept, default 200 bps.
* *max-intron-size*:
* *seed*: integer seed to use for reproducibility. By default, Mikado will use the seed set in the configuration file.
* *single*: flag that disables multiprocessing. Mostly useful for debug purposes.


Command line usage:

.. code-block:: bash

    $ mikado prepare --help
    usage: Mikado prepare [-h] [--fasta REFERENCE] [-v | -q] [--start-method {fork,spawn,forkserver}] [-s | -sa STRAND_SPECIFIC_ASSEMBLIES] [--list LIST] [-l LOG] [--lenient] [-m MINIMUM_CDNA_LENGTH]
                          [-MI MAX_INTRON_LENGTH] [-p PROCS] [-scds] [--labels LABELS] [--codon-table CODON_TABLE] [--single] [-od OUTPUT_DIR] [-o OUT] [-of OUT_FASTA] [--configuration CONFIGURATION] [-er]
                          [--strip-faulty-cds] [--seed SEED]
                          [gff [gff ...]]

    positional arguments:
      gff                   Input GFF/GTF file(s).

    optional arguments:
      -h, --help            show this help message and exit
      --fasta REFERENCE, --reference REFERENCE
                            Genome FASTA file. Required.
      -v, --verbose
      -q, --quiet
      --start-method {fork,spawn,forkserver}
                            Multiprocessing start method.
      -s, --strand-specific
                            Flag. If set, monoexonic transcripts will be left on their strand rather than being moved to the unknown strand.
      -sa STRAND_SPECIFIC_ASSEMBLIES, --strand-specific-assemblies STRAND_SPECIFIC_ASSEMBLIES
                            Comma-delimited list of strand specific assemblies.
      --list LIST           Tab-delimited file containing rows with the following format: <file> <label> <strandedness(def. False)> <score(optional, def. 0)> <is_reference(optional, def. False)>
                            <exclude_redundant(optional, def. True)> <strip_cds(optional, def. False)> <skip_split(optional, def. False)> "strandedness", "is_reference", "exclude_redundant", "strip_cds" and
                            "skip_split" must be boolean values (True, False) "score" must be a valid floating number.
      -l LOG, --log LOG     Log file. Optional.
      --lenient             Flag. If set, transcripts with only non-canonical splices will be output as well.
      -m MINIMUM_CDNA_LENGTH, --minimum-cdna-length MINIMUM_CDNA_LENGTH
                            Minimum length for transcripts. Default: 200 bps.
      -MI MAX_INTRON_LENGTH, --max-intron-size MAX_INTRON_LENGTH
                            Maximum intron length for transcripts. Default: 1,000,000 bps.
      -p PROCS, --procs PROCS
                            Number of processors to use (default None)
      -scds, --strip_cds    Boolean flag. If set, ignores any CDS/UTR segment.
      --labels LABELS       Labels to attach to the IDs of the transcripts of the input files, separated by comma.
      --codon-table CODON_TABLE
                            Codon table to use. Default: 0 (ie Standard, NCBI #1, but only ATG is considered a valid start codon.
      --single, --single-thread
                            Disable multi-threading. Useful for debugging.
      -od OUTPUT_DIR, --output-dir OUTPUT_DIR
                            Output directory. Default: current working directory
      -o OUT, --out OUT     Output file. Default: mikado_prepared.gtf.
      -of OUT_FASTA, --out_fasta OUT_FASTA
                            Output file. Default: mikado_prepared.fasta.
      --configuration CONFIGURATION, --json-conf CONFIGURATION
                            Configuration file.
      -er, --exclude-redundant
                            Boolean flag. If invoked, Mikado prepare will exclude redundant models,ignoring the per-sample instructions.
      --strip-faulty-cds    Flag. If set, transcripts with an incorrect CDS will be retained but with their CDS stripped. Default behaviour: the whole transcript will be considered invalid and discarded.
      --seed SEED           Random seed number.


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


.. _redundant-transcripts-in-prepare:

Removal of redundant transcripts
--------------------------------

Many third-party tools, e.g. gffread [GffRead]_, try to simplify transcript assemblies by grouping transcripts according to their intron chains and then keeping only one transcript per group, usually the longest. This removes transcripts with identical intron chains as well as transcripts whose intron chain is completely contained within another one in the group.
In most cases, Mikado explicitly does **not** take this approach because, especially with RNASeq assemblies, longer transcripts might not necessarily be the most correct; rather, in a non-negligible portion of cases, longer transcripts might have originated by an artefactual fusion of two different, neighbouring transcripts. The implicit assumption made by e.g. gffread (that shorter transcripts are the result of fragmentation of the longer transcripts) would therefore lead to incorrect assemblies.
The default approach taken by Mikado, therefore, is to identify cases where transcripts are **completely** identical (both in terms of cDNA and CDS, if kept), and only remove redundancies in those rare, specific cases.

In certain situations, however, a strategy based on intron chain redundancy like in gffread might be warranted. Specifically:

- in long read datasets (e.g. PacBio or ONT alignments) the implicit assumption made by gffread is valid: in these cases it is safe to assume that fragmentation during RNA extraction and library preparation would constitute the main origin of redundancy.
- when dealing with massive transcript datasets (>=5-10 million transcripts), removing excess transcripts might be necessary to keep the analysis manageable, at the cost of slightly reduced accuracy.

Mikado allows to perform this more extensive redundancy removal either on a *per-analysis* or *per-sample* (recommended) basis.
When scanning the transcript assemblies, Mikado will look for intron chains that are completely contained within another. When such an occurence arises, *if and only if Mikado has been instructed to remove redundant cases*, Mikado will do the following:

- if one of the two transcripts comes from a sample for which the redundancy removal is disabled (including, automatically, all "reference" samples), it will always be kept.
- if the transcript marked for redundancy check and removal has a lower baseline score *or* is contained within the other transcript, it will be marked for removal.


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

.. _orf-check-prepare:

Check on ORF correctness
------------------------

Mikado will check that the input transcripts have a formally valid ORF both in terms of the structure as well as of its sequence.

By "formally correct ORF structure", Mikado means that:

- CDS segments must be contained within declared exons
- CDS segments should not overlap each other or any declared UTR segment
- there should be no gap between CDS segments on a transcript's cDNA

By "formally correct ORF sequence", Mikado means that:

- there are no internal in-frame stop codons
- including the initial phase, the length of the CDS should be a multiple of three.

To perform the latter check, Mikado will use the specified codon table (by default "0", ie the NCBI Standard codon table but only considering "ATG" as a valid start codon).

Finally, Mikado will infer whether the transcript has a start and/or stop codon and tag it appropriately.

Output files
------------

Mikado prepare will produce two files:

* a *sorted* GTF file, containing all the transcripts surviving the checks
* a FASTA file of the transcripts, in the proper cDNA orientation.

.. warning:: contrary to other tools such as eg gffread [GffRead]_, Mikado prepare will **not** try to calculate the loci for the transcripts. This task will be performed later in the pipeline. As such, the GTF file is formally incorrect, as multiple transcripts in the same locus but coming from different assemblies will *not* have the same gene_id but rather will have kept their original one. Moreover, if two gene_ids were identical but discrete in the input files (ie located on different sections of the genome), this error will not be corrected. If you desire to use this GTF file for any purpose, please use a tool like gffread to calculate the loci appropriately.
