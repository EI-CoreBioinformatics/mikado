.. _utils:

Mikado miscellaneous scripts
============================

All these utilities can be accessed with the ``mikado util`` CLI. They perform relatively minor tasks.

awk_gtf
~~~~~~~

This utility is used to retrieve specific regions from a GTF file, without breaking any transcript feature. Any transcript falling even partly within the specified coordinates will be retained in its entirety.

Usage::

    $ mikado util awk_gtf --help
    usage: mikado.py util awk_gtf [-h] (-r REGION | --chrom CHROM) [-as]
                                  [--start START] [--end END]
                                  gtf [out]


    positional arguments:
      gtf
      out

    optional arguments:
      -h, --help            show this help message and exit
      -r REGION, --region REGION
                            Region defined as a string like <chrom>:<start>..<end>
      --chrom CHROM
      -as, --assume-sorted
      --start START
      --end END``

.. _class-codes-command:

class_codes
~~~~~~~~~~~

This utility is used to obtain information about any class code or category thereof.

Usage::

    $ mikado util class_codes --help
    usage: mikado util class_codes [-h]
                                   [-f {fancy_grid,grid,html,jira,latex,latex_booktabs,mediawiki,moinmoin,orgtbl,pipe,plain,psql,rst,simple,textile,tsv}]
                                   [-c {Intronic,Match,Alternative splicing,Unknown,Fragment,Overlap,Extension,Fusion} [{Intronic,Match,Alternative splicing,Unknown,Fragment,Overlap,Extension,Fusion} ...]]
                                   [-o OUT]
                                   [{,=,_,n,J,c,C,j,h,g,G,o,e,m,i,I,ri,rI,f,x,X,p,P,u} [{,=,_,n,J,c,C,j,h,g,G,o,e,m,i,I,ri,rI,f,x,X,p,P,u} ...]]

    Script to print out the class codes.

    positional arguments:
      {[],=,_,n,J,c,C,j,h,g,G,o,e,m,i,I,ri,rI,f,x,X,p,P,u}
                            Codes to query.

    optional arguments:
      -h, --help            show this help message and exit
      -f {fancy_grid,grid,html,jira,latex,latex_booktabs,mediawiki,moinmoin,orgtbl,pipe,plain,psql,rst,simple,textile,tsv}, --format {fancy_grid,grid,html,jira,latex,latex_booktabs,mediawiki,moinmoin,orgtbl,pipe,plain,psql,rst,simple,textile,tsv}
      -c {Intronic,Match,Alternative splicing,Unknown,Fragment,Overlap,Extension,Fusion} [{Intronic,Match,Alternative splicing,Unknown,Fragment,Overlap,Extension,Fusion} ...], --category {Intronic,Match,Alternative splicing,Unknown,Fragment,Overlap,Extension,Fusion} [{Intronic,Match,Alternative splicing,Unknown,Fragment,Overlap,Extension,Fusion} ...]
      -o OUT, --out OUT

convert
~~~~~~~

This utility is used to convert between GTF and GFF3 files, with the possibility of giving as output BED12 files as well. It is limited to converting transcript features, and will therefore ignore any other feature present (transposons, loci, etc.). The output of the conversion to GFF3 is completely GFF3 compliant.

Usage::

    $ mikado util convert --help
    usage: mikado.py util convert [-h] [-of {bed12,gtf,gff3}] gf [out]

    positional arguments:
      gf
      out

    optional arguments:
      -h, --help            show this help message and exit
      -of {bed12,gtf,gff3}, --out-format {bed12,gtf,gff3}


.. _grep-command:

grep
~~~~

This utility extracts specific transcripts and genes from an input GTF/GFF3 file. As input, it requires a text file of either the format "<transcript id><tab><gene id>", or simply gene per line (in which case the "--genes" switch has to be invoked). If only some of the transcripts of a gene are included in the text file, the gene feature will be shrunk accordingly. The name is an obvious homage to the invaluable UNIX command that we all love.

Usage::

    $ mikado util grep --help
    usage: mikado.py util grep [-h] [-v] [--genes] ids gff [out]

    Script to extract specific models from GFF/GTF files.

    positional arguments:
      ids         ID file (format: mrna_id, gene_id - tab separated)
      gff         The GFF file to parse.
      out         Optional output file

    optional arguments:
      -h, --help  show this help message and exit
      -v          Exclude from the gff all the records in the id file.
      --genes     Flag. If set, the program expects as ids only a list of genes,
                  and will exclude/include all the transcripts children of the
                  selected genes.

.. _merge-blast-command:

merge_blast
~~~~~~~~~~~

This script merges together various XML BLAST+ files into a single entity. It might be of use when the input data has been chunked into different FASTA files for submission to a cluster queue. It is also capable of converting from ASN files and of dealing with GZipped files.

Usage::

    $ mikado util merge_blast --help
    usage: mikado.py util merge_blast [-h] [-v] [-l LOG] [--out [OUT]]
                                      xml [xml ...]

    positional arguments:
      xml

    optional arguments:
      -h, --help         show this help message and exit
      -v, --verbose
      -l LOG, --log LOG
      --out [OUT]

.. _metrics-command:

metrics
~~~~~~~

This command generates the documentation regarding the available transcript metrics. It is generated dynamycally by inspecting the code. The documentation in the :ref:`introduction <Metrics>` is generated using this utility.

Usage::

    $ mikado util metrics --help
    usage: mikado util metrics [-h]
                               [-f {fancy_grid,grid,html,jira,latex,latex_booktabs,mediawiki,moinmoin,orgtbl,pipe,plain,psql,rst,simple,textile,tsv}]
                               [-o OUT]
                               [-c {CDS,Descriptive,External,Intron,Locus,UTR,cDNA} [{CDS,Descriptive,External,Intron,Locus,UTR,cDNA} ...]]
                               [metric [metric ...]]

    Simple script to obtain the documentation on the transcript metrics.

    positional arguments:
      metric

    optional arguments:
      -h, --help            show this help message and exit
      -f {fancy_grid,grid,html,jira,latex,latex_booktabs,mediawiki,moinmoin,orgtbl,pipe,plain,psql,rst,simple,textile,tsv}, --format {fancy_grid,grid,html,jira,latex,latex_booktabs,mediawiki,moinmoin,orgtbl,pipe,plain,psql,rst,simple,textile,tsv}
                            Format of the table to be printed out.
      -o OUT, --out OUT     Optional output file
      -c {CDS,Descriptive,External,Intron,Locus,UTR,cDNA} [{CDS,Descriptive,External,Intron,Locus,UTR,cDNA} ...], --category {CDS,Descriptive,External,Intron,Locus,UTR,cDNA} [{CDS,Descriptive,External,Intron,Locus,UTR,cDNA} ...]
                            Available categories to select from.


.. _stat-command:

stats
~~~~~

This command generates a statistics file for GFF3/GTF files. The output is a table including Average, Mode, and various quantiles for different features present in a typical GFF file (genes, introns, exons, cDNAs, etc.). The operation can be quite time consuming for large files, in which case it is advisable to ask for multiple processors.

Usage::

    $ mikado util stats --help
    usage: mikado.py util stats [-h] [--only-coding] [-p PROCS] gff [out]

    GFF/GTF statistics script. It will compute median/average length of RNAs,
    exons, CDS features, etc.

    positional arguments:
      gff                   GFF file to parse.
      out

    optional arguments:
      -h, --help            show this help message and exit
      --only-coding
      -p PROCS, --processors PROCS

A typical example statistics file can be found :download:`here, for the TAIR10 annotation <./TAIR10.stats>`.

.. _trim-command:

trim
~~~~

This utility trims down the terminal exons of multiexonic transcripts, until either shrinking them to the desired maximum length or meeting the beginning/end of the CDS. It has been used for generating the "trimmed" annotations for the analysis of the original Mikado paper.

Usage::

    $ mikado util trim --help
    usage: mikado.py util trim [-h] [-ml MAX_LENGTH] [--as-gtf] ann [out]

    positional arguments:
      ann                   Reference GTF/GFF output file.
      out

    optional arguments:
      -h, --help            show this help message and exit
      -ml MAX_LENGTH, --max_length MAX_LENGTH
                            Maximal length of trimmed terminal exons
      --as-gtf              Flag. If set, the output will be in GTF rather than
                            GFF3 format.


.. _included_scripts:

Included scripts
================

All the following scripts are included in the "util" folder in the source code, and will be included on the PATH after installation. Some of this scripts are used by the :ref:`Daijin` pipeline to produce statistics or perform other intermediate steps.

add_transcript_feature_to_gtf.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script is needed to add a top-level transcript feature to GTFs that lack it, eg. those produced by CuffMerge [CuffMerge]_.

Usage::

    $ add_transcript_feature_to_gtf.py --help
    usage: Script to add a transcript feature to e.g. Cufflinks GTFs
           [-h] gtf [out]

    positional arguments:
      gtf         Input GTF
      out         Output file. Default: stdout.

    optional arguments:
      -h, --help  show this help message and exit

align_collect.py
~~~~~~~~~~~~~~~~

This script is used to collect statistics from `samtools stat <www.htslib.org/doc/samtools.html>`_.
Usage::

    $ align_collect.py  --help
    usage: Script to collect info from multiple samtools stats files
           [-h] input [input ...]

    positional arguments:
      input       The list of samtools stats file to process

    optional arguments:
      -h, --help  show this help message and exit

asm_collect.py
~~~~~~~~~~~~~~

This script is used to collect statistics obtained with from the :ref:`mikado util stats <stat-command>` utility. Output is printed directly to the screen. Usage::

    $ asm_collect.py -h
    usage: Script to collect info from multiple mikado util stats files
           [-h] input [input ...]

    positional arguments:
      input       The list of mikado util stats file to process

    optional arguments:
      -h, --help  show this help message and exit

bam2gtf.py
~~~~~~~~~~

This script will use PySam to convert read alignments into a GTF file. Mostly useful to convert from BAM alignment of long reads (eg. PacBio) into a format which Mikado can interpret and use.

Usage::

    $ bam2gtf.py --help
    usage: Script to convert from BAM to GTF, for PB alignments [-h] bam [out]

    positional arguments:
      bam         Input BAM file
      out         Optional output file

    optional arguments:
      -h, --help  show this help message and exit


class_run.py
~~~~~~~~~~~~

Python3 wrapper for the CLASS [Class2]_ assembler. It will perform the necessary operations for the assembler (depth and call of the splicing junctions), and launch the program itself. Usage::

    $ class_run.py --help
    usage: Quick utility to rewrite the wrapper for CLASS. [-h] [--clean]
                                                           [--force]
                                                           [-c CLASS_OPTIONS]
                                                           [-p PROCESSORS]
                                                           [--class_help] [-v]
                                                           [bam] [out]

    positional arguments:
      bam                   Input BAM file.
      out                   Optional output file.

    optional arguments:
      -h, --help            show this help message and exit
      --clean               Flag. If set, remove tepmorary files.
      --force               Flag. If set, it forces recalculation of all
                            intermediate files.
      -c CLASS_OPTIONS, --class_options CLASS_OPTIONS
                            Additional options to be passed to CLASS. Default: no
                            additional options.
      -p PROCESSORS, --processors PROCESSORS
                            Number of processors to use with class.
      --class_help          If called, the wrapper will ask class to display its
                            help and exit.
      -v, --verbose

getFastaFromIds.py
~~~~~~~~~~~~~~~~~~

Script to extract a list of sequences from a FASTA file, using the `pyfaidx <https://pypi.python.org/pypi/pyfaidx>`_ [PyFaidx]_ module. Usage::

    $ getFastaFromIds.py -h
    usage: getFastaFromIds.py [-h] [-v] list fasta [out]

    A simple script that retrieves the FASTA sequences from a file given a list of
    ids.

    positional arguments:
      list           File with the list of the ids to recover, one by line.
                     Alternatively, names separated by commas.
      fasta          FASTA file.
      out            Optional output file.

    optional arguments:
      -h, --help     show this help message and exit
      -v, --reverse  Retrieve entries which are not in the list, as in grep -v (a
                     homage).

gffjunc_to_bed12.py
~~~~~~~~~~~~~~~~~~~

Script to convert a GFF junction file to a BED12 file. Useful to format the input for Mikado serialise.

Usage::

    $ gffjunc_to_bed12.py --help
    usage: GFF=>BED12 converter [-h] gff [out]

    positional arguments:
      gff
      out

    optional arguments:
      -h, --help  show this help message and exit

grep.py
~~~~~~~

A script to extract data from *column* files, using a list of targets. More efficient than a standard "grep -f" for this niche case.

Usage::

    $ util/grep.py -h
    usage: grep.py [-h] [-v] [-s SEPARATOR] [-f FIELD] [-q] ids target [out]

    This script is basically an efficient version of the GNU "grep -f" utility for
    table-like files, and functions with a similar sintax.

    positional arguments:
      ids                   The file of patterns to extract
      target                The file to filter
      out                   The output file

    optional arguments:
      -h, --help            show this help message and exit
      -v, --reverse         Equivalent to the "-v" grep option
      -s SEPARATOR, --separator SEPARATOR
                            The field separator. Default: consecutive
                            whitespace(s)
      -f FIELD, --field FIELD
                            The field to look in the target file.
      -q, --quiet           No logging.

merge_junction_bed12.py
~~~~~~~~~~~~~~~~~~~~~~~

This script will merge [Portcullis]_-like junctions into a single BED12, using the thick start/ends as unique keys.

Usage::

    $ merge_junction_bed12.py --help
    usage: Script to merge BED12 files *based on the thickStart/End features*.
        Necessary for merging junction files such as those produced by TopHat
           [-h] [--delim DELIM] [-t THREADS] [--tophat] [-o OUTPUT] bed [bed ...]

    positional arguments:
      bed                   Input BED files. Use "-" for stdin.

    optional arguments:
      -h, --help            show this help message and exit
      --delim DELIM         Delimiter for merged names. Default: ;
      -t THREADS, --threads THREADS
                            Number of threads to use for multiprocessing. Default:
                            1
      --tophat              Flag. If set, tophat-like junction style is assumed.
                            This means that junctions are defined using the
                            blockSizes rather than thickStart/End. The script will
                            convert the lines to this latter format. By default,
                            the script assumes that the intron start/end are
                            defined using thickStart/End like in portcullis.
                            Mixed-type input files are not supported.
      -o OUTPUT, --output OUTPUT
                            Output file. Default: stdout


remove_from_embl.py
~~~~~~~~~~~~~~~~~~~

Quick script to remove sequences from a given organism from SwissProt files, and print them out in FASTA format. Used to produce the BLAST datasets for the Mikado paper. Usage::

    $ remove_from_embl.py -h
    usage: Script to remove sequences specific of a given organism from a SwissProt file.
           [-h] -o ORGANISM [--format {fasta}] input [out]

    positional arguments:
      input
      out

    optional arguments:
      -h, --help            show this help message and exit
      -o ORGANISM, --organism ORGANISM
                            Organism to be excluded
      --format {fasta}      Output format. Choices: fasta. Default: fasta.

sanitize_blast_db.py
~~~~~~~~~~~~~~~~~~~~

Simple script to clean the header of FASTA files, so to avoid runtime errors and incrongruencies with BLAST and other tools which might be sensitive to long descriptions or the presence of special characters.

Usage::

    $ sanitize_blast_db.py --help
    usage: sanitize_blast_db.py [-h] [-o OUT] fasta [fasta ...]

    positional arguments:
      fasta

    optional arguments:
      -h, --help         show this help message and exit
      -o OUT, --out OUT


split_fasta.py
~~~~~~~~~~~~~~

This script is used to split a FASTA file in a fixed number of files, with an approximate equal number of sequences in each. If the number of sequences in the input file is lower than the number of requested splits, the script will create the necessary number of empty files. Used in :ref:`Daijin` for preparing the input data for the BLAST analysis. Usage::

    $ split_fasta.py --help
    usage: Script to split FASTA sequences in a fixed number of multiple files.
           [-h] [-m NUM_FILES] fasta [out]

    positional arguments:
      fasta                 Input FASTA file.
      out                   Output prefix. Default: filename+split

    optional arguments:
      -h, --help            show this help message and exit
      -m NUM_FILES, --num-files NUM_FILES
                            Number of files to create. Default: 1000

trim_long_introns.py
~~~~~~~~~~~~~~~~~~~~

This script parses an annotation file and truncates any transcript which has *UTR* introns over the provided threshold. In such cases, the UTR section after the long intron is simply removed. Usage::

    $ trim_long_introns.py --help
    usage: This script truncates transcript with UTR exons separated by long introns.
           [-h] [-mi MAX_INTRON] gff [out]

    positional arguments:
      gff
      out

    optional arguments:
      -h, --help            show this help message and exit
      -mi MAX_INTRON, --max-intron MAX_INTRON
                            Maximum intron length for UTR introns.

