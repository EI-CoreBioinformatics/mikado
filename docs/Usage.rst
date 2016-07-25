Usage
=====


Mikado is composed of four different programs (*configure, prepare, serialise, pick*) which have to be executed serially to go from an ensemble of different assemblies to the final dataset. In addition to these core programs, Mikado provides a utility to compare annotations, similarly to CuffCompare and ParsEval (*compare*), and various other minor utilities to perform operations such as extracting regions from a GFF, convert between different gene annotation formats, etc.

Mikado pipeline
---------------


Compare
-------

Mikado provides a utility to compare two different annotations. Please see the :ref:`relevant page <compare>` for details.


Utilities
---------

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

    $ mikado util metrics


.. _stat-command:

stat
~~~~

This command generates a statistics file for GFF3/GTF files. The output is a table including Average, Mode, and various quantiles for different features present in a typical GFF file (genes, introns, exons, cDNAs, etc.). The operation can be quite time consuming for large files, in which case it is advisable to ask for multiple processors.

.. warning:: GTF files have to have valid "transcript" features as top-level for CDS/exons.

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

A typical example statistics file can be found :download:`here, for the TAIR10 annotation <TAIR10.stats>`.

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

