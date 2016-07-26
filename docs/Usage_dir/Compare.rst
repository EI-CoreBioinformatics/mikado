.. _F1: https://en.wikipedia.org/wiki/F1_score
.. _Cuffcompare: http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html
.. _ParsEval: https://aegean.readthedocs.io/en/v0.16.0/parseval.html

.. _Compare:

Compare
=======

Overview
~~~~~~~~

This Mikado utility allows the user to compare the transcripts from any two annotations. Its output allows:
  - To understand which reference transcript each prediction is most similar to
  - To understand which prediction transcript best represent each reference model
  - To have a summary information about the similarity between the two annotations.

Mikado compare has been directly inspired by the popular `Cuffcompare`_ [Cufflinks]_ utility and by `ParsEval`_ [ParsEval]_. Please note that while superficially similar to Cuffcompare in the style of the output files, Mikado compare is more philosophically similar to ParsEval, as it will not try to aggregate transcripts in loci but will perform a pure comparison between the two annotation files. Both GTF and GFF files are accepted, in any combination.

Usage
~~~~~

Mikado compare is invoked by specifying the *reference* annotation and the desired mode of analysis. There are three possible options:

 #. In its default mode, compare will ask for a *prediction* annotation to compare the reference against.
 #. In the *"self"* mode, compare will do a self-comparison of the reference against itself, excluding as possible results the matches between a transcript and itself. It can be useful to glean the relationships between transcripts and genes in an annotation.
 #. In the *"internal"* mode of operations, compare will again perform a self-comparison, focussed on multi-isoform genes. For those, compare will perform and report all possible comparisons. It is useful to understand the relationships between the transcripts in a single locus.


Mikado stores the information of the reference in a specialised index, with a ".midx" suffix, which will be created by the program upon its first execution with a new reference. If the index file is already present, Mikado will try to use it rather than read again the annotation.

Command line
------------

.. code-block:: bash

    $ mikado compare --help
    usage: Mikado compare [-h] -r REFERENCE
                          (-p PREDICTION | --self | --internal | --index)
                          [--distance DISTANCE] [-pc] [-o OUT] [--lenient] [-eu]
                          [-n] [-l LOG] [-v]

    optional arguments:
      -h, --help            show this help message and exit
      --distance DISTANCE   Maximum distance for a transcript to be considered a
                            polymerase run-on. Default: 2000
      -pc, --protein-coding
                            Flag. If set, only transcripts with a CDS (both in
                            reference and prediction) will be considered.
      -o OUT, --out OUT     Prefix for the output files. Default: mikado_compare
      --lenient             If set, exonic statistics will be calculated leniently
                            in the TMAP as well - ie they will consider an exon as
                            match even if only the internal junction has been
                            recovered.
      -eu, --exclude-utr    Flag. If set, reference and prediction transcripts
                            will be stripped of their UTRs (if they are coding).
      -n, --no-index, --no-save-index
                            Unless this flag is set, compare will save an index of
                            the reference to quicken multiple calls.
      -l LOG, --log LOG
      -v, --verbose

    Prediction and annotation files.:
      -r REFERENCE, --reference REFERENCE
                            Reference annotation file. By default, an index will
                            be crated and saved with the suffix ".midx".
      -p PREDICTION, --prediction PREDICTION
                            Prediction annotation file.
      --self                Flag. If set, the reference will be compared with
                            itself. Useful for understanding how the reference
                            transcripts interact with each other.
      --internal            Flag. If set, for each gene with more than one
                            transcript isoform each will be compared to the
                            others. Useful for understanding the structural
                            relationships between the transcripts in each gene.
      --index               Flag. If set, compare will stop after having generated
                            the GFF index for the reference.

Output files
~~~~~~~~~~~~

Mikado compare produces two tabular files, tmap_ and refmap_, and one :ref:`statistics <stats>` file.

.. _tmap:

TMAP files
----------

TMAP are tabular files that store the information regarding the best match for each prediction in the reference. The columns are as follows:

#. **ref_id**: Transcript ID of the matched reference model(s).
#. **ref_gene**: Gene ID of the matched reference model(s).
#. **ccode**: class code of the match. See :ref:`the relevant section on Class codes <ccodes>`.
#. **tid**: Transcript ID of the prediction model.
#. **gid**: Gene ID of the prediction model.
#. **tid_num_exons**: Number of exons of the prediction model.
#. **ref_num_exons**: Number of exons of the reference model.
#. **n_prec**: Nucleotide precision of the prediction ( TP / (length of the prediction))
#. **n_recall**: Nucleotide recall of the reference (TP / (length of the reference))
#. **n_f1**: `F1`_ of recall and precision at the nucleotide level.
#. **j_prec**: Splice junction precision of the prediction model ( TP / (number of splice sites in the prediction))
#. **j_recall**: Splice junction recall of the reference model ( TP / (number of splice sites in the reference))
#. **j_f1**: `F1`_ of recall and precision at the splice junction level.
#. **e_prec**: Exon precision of the prediction model ( TP / (number of exons in the prediction)). **NB**: this value is calculated "leniently", ie terminal exons count as a match if the *internal* border is called correctly and the exon is terminal in both prediction and reference.
#. **e_recall**: Exon recall of the reference model ( TP / (number of exons in the reference))
#. **e_f1**: `F1`_ of recall and precision at the exon level.
#. **distance**: Distance of the model from its putative match.

An example of TMAP file is as follows::

    ref_id	ref_gene	ccode	tid	gid	tid_num_exons	ref_num_exons	n_prec	n_recall	n_f1	j_prec	j_recall	j_f1	e_prec	e_recall	e_f1	distance
    AT5G66610.1	AT5G66610	=	mikado.Chr5G1.2	mikado.Chr5G1	11	11	97.77	99.16	98.46	100.00	100.00	100.00	81.82	81.82	81.82	0
    AT5G66610.1	AT5G66610	j	mikado.Chr5G1.1	mikado.Chr5G1	11	11	92.93	94.74	93.82	95.00	95.00	95.00	81.82	81.82	81.82	0
    AT5G66620.1,AT5G66630.1,AT5G66631.1	AT5G66620,AT5G66630,AT5G66631	f,j,J,O	st_Stringtie_STAR.21710.6	Stringtie_STAR.21710	22	11,10,1	27.84,34.47,28.63	99.79,99.13,100.00	43.54,51.16,44.52	45.24,42.86,0.00	95.00,100.00,0.00	61.29,60.00,0.00	36.36,36.36,0.00	72.73,80.00,0.00	48.48,50.00,0.00	0

You can notice that the third example is particular as the prediction transcript matches not one but multiple reference transcripts. This is a fusion_ event.

.. _refmap:

RefMap files
------------

RefMap files are tabular files which store the information regarding the best match for each reference transcript, among all possible prediction models. The columns of the file are as follows:

#. **ref_id**: Transcript ID of the reference model.
#. **ccode**: class code of the match. See :ref:`the relevant section on Class codes <ccodes>`.
#. **tid**: Transcript ID of the prediction model.
#. **gid**: Gene ID of the prediction model.
#. **nF1**: `F1`_ of recall and precision at the nucleotide level.
#. **jF1**: `F1`_ of recall and precision at the splice junction level.
#. **eF1**: `F1`_ of recall and precision at the exon level. **NB**: this value is calculated "leniently", ie terminal exons count as a match if the *internal* border is called correctly and the exon is terminal in both prediction and reference.
#. **ref_gene**: Gene ID of the reference model.
#. **best_ccode**: Best possible class code found for any of the transcripts of the gene.
#. **best_tid**: Transcript ID of the prediction model which fit best one of the transcript models of the reference gene.
#. **best_gid**: Gene ID of the prediction model which fit best one of the transcript models of the reference gene.
#. **best_nF1**: `F1`_ of recall and precision at the nucleotide level, for the best possible comparison.
#. **best_jF1**: `F1`_ of recall and precision at the splice junction level, for the best possible comparison.
#. **best_eF1**: `F1`_ of recall and precision at the exon level, for the best possible comparison.

An example of a RefMap file is as follows::

    ref_id	ccode	tid	gid	nF1	jF1	eF1	ref_gene	best_ccode	best_tid	best_gid	best_nF1	best_jF1	best_eF1
    AT5G66610.1	=	mikado.Chr5G1.2	mikado.Chr5G1	98.46	100.0	81.82	AT5G66610	=	mikado.Chr5G1.2	mikado.Chr5G1	98.46	100.0	81.82
    AT5G66610.2	J	mikado.Chr5G1.2	mikado.Chr5G1	93.91	94.74	76.19	AT5G66610	=	mikado.Chr5G1.2	mikado.Chr5G1	98.46	100.0	81.82
    AT5G66630.1	f,n	tr_c58_g1_i4.mrna1.6	c58_g1_i4.path1.6	66.32	94.74	76.19	AT5G66630	f,n	tr_c58_g1_i4.mrna1.6	c58_g1_i4.path1.6	66.32	94.74	76.19

Please note that the third example (AT5G66630.1) has as best possible match a fusion_ event.

.. _stats:

Stats files
-----------

These files provide a summary of the comparison between the reference and the annotation. An example is as follows::

    Command line:
    /usr/users/ga002/venturil/py351/bin/mikado compare -r reference.gff3 -p mikado.loci.gff3 -o compare -l compare.log
    7 reference RNAs in 5 genes
    15 predicted RNAs in  8 genes
    --------------------------------- |   Sn |   Pr |   F1 |
                            Base level: 85.74  64.73  73.77
                Exon level (stringent): 63.83  42.86  51.28
                  Exon level (lenient): 80.00  52.94  63.72
                          Intron level: 89.47  59.65  71.58
                    Intron chain level: 33.33  14.29  20.00
          Transcript level (stringent): 0.00  0.00  0.00
      Transcript level (>=95% base F1): 28.57  13.33  18.18
      Transcript level (>=80% base F1): 42.86  20.00  27.27
             Gene level (100% base F1): 0.00  0.00  0.00
            Gene level (>=95% base F1): 40.00  25.00  30.77
            Gene level (>=80% base F1): 60.00  37.50  46.15

    #   Matching: in prediction; matched: in reference.

                Matching intron chains: 2
                 Matched intron chains: 2
       Matching monoexonic transcripts: 1
        Matched monoexonic transcripts: 1
            Total matching transcripts: 3
             Total matched transcripts: 3

              Missed exons (stringent): 17/47  (36.17%)
               Novel exons (stringent): 40/70  (57.14%)
                Missed exons (lenient): 9/45  (20.00%)
                 Novel exons (lenient): 32/68  (47.06%)
                        Missed introns: 4/38  (10.53%)
                         Novel introns: 23/57  (40.35%)

                    Missed transcripts: 0/7  (0.00%)
                     Novel transcripts: 6/15  (40.00%)
                          Missed genes: 0/5  (0.00%)
                           Novel genes: 2/8  (25.00%)

The first section of the file describes:

  #. Concordance of the two annotations at the base level (recall, precision, and F1)
  #. Concordance of the two annotation at the exonic level (recall, precision, and F1), in two ways:

     * *"stringent"*: only perfect exonic matches are considered.
     * *"lenient"*: in this mode, terminal exons are counted as a match if the **internal** border is matched. See the RGASP paper [RGASP]_ for details on the rationale.

  #. Concordance of the two annotations at the intron level.
  #. Concordance of the two annotations at the intron chain level - how many intron chains of the reference are found identical in the prediction. Only multiexonic models are considered for this level.
  #. Concordance of the two annotations at the transcript level, in three different modes:

     * *"stringent"*: in this mode, only perfect matches are considered.
     * *"95% base F1"*: in this mode, we only count instances where the nucleotide F1 is greater than *95%* and, for multiexonic transcripts, the intron chain is reconstructed perfectly.
     * *"80% base F1"*: in this mode, we only count instances where the nucleotide F1 is greater than *80%* and, for multiexonic transcripts, the intron chain is reconstructed perfectly.

  #. Concordance of the two annotations at the gene level, in three different modes:

     * *"stringent"*: in this mode, we consider reference genes for which it was possible to find at least one perfect match for one of its transcripts.
     * *"95% base F1"*: in this mode, we only count instances where the nucleotide F1 is greater than *95%* and, for multiexonic transcripts, the intron chain is reconstructed perfectly.
     * *"80% base F1"*: in this mode, we only count instances where the nucleotide F1 is greater than *80%* and, for multiexonic transcripts, the intron chain is reconstructed perfectly.

In the second section, the file reports how many of the intron chains, monoexonic transcripts and total transcripts in the **reference** were *matched* by at least one *matching* **prediction** transcript. Finally, in the third section the file reports the number of missed (present in the reference but not in the prediction) or novel (viceversa - present in the prediction but not in the reference) features.

.. note:: Please note that a gene might be considered as "found" even if its best match is intronic, on the opposite strand, or not directly overlapping it, or is in the opposite strand (see :ref:`next section <ccodes>`, in particular the *Intronic*, *Fragment* and *No overlap* categories).


.. _ccodes:

Class codes
~~~~~~~~~~~

In addition to recall, precision and F1 values, Mikado assign each comparison between two transcripts a *class code*, which summarises the relationship between the two transcripts. The idea is lifted from the popular tool `Cuffcompare`_, although Mikado greatly extends the catalogue of possible class codes.
All class codes fall within one of the following categories:

 - **Match**: class codes of this type indicate concordance between the two transcript models.
 - **Extension**: class codes of this type indicate that one of the two models extends the intron chain of the other, without internal interruptions. The extension can be from either perspective - either the prediction extends the reference, or it is instead *contained* within the reference (so that switching perspectives, the reference would "extend" the prediction).
 - **Alternative splicing**: the two exon chains overlap but differ in significant ways.
 - **Intronic**: either the prediction is completely contained within the introns of the reference, or viceversa.
 - **Fragment**: the prediction is a fragment of the reference, in most cases because they are on opposite strands.
 - **No overlap**: the prediction and the reference are near but do not directly overlap.

 .. _fusion:

 - **Fusion**: this special class code is a qualifier and it never appears on its own. When a transcript is defined as a fusion,  its class code in the *tmap* file will be an "f" followed by the class codes of the individual transcript matches, sperated by comma. So a prediction which matches two reference models, one with a "j" and another with a "o", will have a class code of **"f,j,o"**. In the *refmap* file, if the fusion is the best match, the class code will be "f" followed by the class code for the individual reference transcript; e.g., **"f,j"**



+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| Class code   | Definition                   | Is the       | Is the        | Nucleotide:       | Junction:         | Reverse           | Category          |
|              |                              | reference    | prediction    | Recall,           | Recall,           | class code        |                   |
|              |                              | transcript   | transcript    | Precision,        | Precision,        |                   |                   |
|              |                              | multiexonic? | multiexonic?  | F1                | F1                |                   |                   |
|              |                              |              |               |                   |                   |                   |                   |
+==============+==============================+==============+===============+===================+===================+===================+===================+
| **=**        | Complete intron chain match. | True         | True          | NA                | 100%, 100%, 100%  | **=**             | **Match**         |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **_**        | Complete match between two   | False        | False         | NA, NA, **>=80%** | NA                | **_**             | **Match**         |
| (underscore) | monoexonic transcripts.      |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **m**        | Generic match between two    | False        | False         | NA, NA, **< 80%** | NA                | **m**             | **Match**         |
|              | monoexonic transcripts.      |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **n**        | Intron chain extension, ie.  | True         | True          | **100%**, < 100%, | 100%, < 100%,     | **c**             | **Extension**     |
|              | both transcripts are         |              |               | < 100%            | < 100%            |                   |                   |
|              | multiexonic and the          |              |               |                   |                   |                   |                   |
|              | prediction has novel         |              |               |                   |                   |                   |                   |
|              | splice sites *outside* of    |              |               |                   |                   |                   |                   |
|              | the reference transcript     |              |               |                   |                   |                   |                   |
|              | boundaries.                  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **J**        | Intron chain extension,      | True         | True          | < 100%, <= 100%,  | **100%**, < 100%, | **C**             | **Extension**     |
|              | both transcripts are         |              |               | < 100%            | < 100%            |                   |                   |
|              | multiexonic and the          |              |               |                   |                   |                   |                   |
|              | prediction has novel         |              |               |                   |                   |                   |                   |
|              | splice sites *inside* of the |              |               |                   |                   |                   |                   |
|              | reference transcript         |              |               |                   |                   |                   |                   |
|              | boundaries.                  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **c**        | The prediction               | True         | NA            | < 100%, **100%**  | < 100%, **100%**  | **n**             | **Extension**     |
|              | is either multiexonic and    |              |               | NA                | NA                |                   |                   |
|              | with its intron chain        |              |               |                   |                   |                   |                   |
|              | completely contained within  |              |               |                   |                   |                   |                   |
|              | that of the reference, or    |              |               |                   |                   |                   |                   |
|              | monoexonic and contained     |              |               |                   |                   |                   |                   |
|              | within one of the reference  |              |               |                   |                   |                   |                   |
|              | exons.                       |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **C**        | The prediction intron chain  | True         | True          | <= 100%, < 100%,  | < 100%, **100%**, | **J**             | **Extension**     |
|              | is completely contained      |              |               | < 100%            | < 100%            |                   |                   |
|              | within that of the           |              |               |                   |                   |                   |                   |
|              | reference transcript, but    |              |               |                   |                   |                   |                   |
|              | it partially debords either  |              |               |                   |                   |                   |                   |
|              | into its introns or outside  |              |               |                   |                   |                   |                   |
|              | of the reference boundaries. |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **j**        | Alternative splicing event.  | True         | True          | NA                | <= 100%, < 100%,  | **j**             | **Alternative     |
|              |                              |              |               |                   | < 100%            |                   | splicing**        |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **h**        | Structural match between two | True         | True          | > 0%, > 0%, > 0%  | 0%, 0%, 0%        | **h**             | **Alternative     |
|              | models where no splice site  |              |               |                   |                   |                   | splicing**        |
|              | is conserved but **at least**|              |               |                   |                   |                   |                   |
|              | one intron of the reference  |              |               |                   |                   |                   |                   |
|              | and one intron of the        |              |               |                   |                   |                   |                   |
|              | prediction partially overlap.|              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **o**        | Generic overlap between two  | True         | True          | > 0%, > 0%, > 0%  | 0%, 0%, 0%        | **o**             | **Alternative     |
|              | multiexonic transcripts,     |              |               |                   |                   |                   | splicing**        |
|              | which do not share **any**   |              |               |                   |                   |                   |                   |
|              | overlap among their introns. |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **g**        | The monoexonic prediction    | True         | False         | > 0%, > 0%,       | 0%                | **G**             | **Alternative     |
| ("mo" before | overlaps one or more exons of|              |               | between 0 and 100%|                   |                   | splicing**        |
| release 1)   | the reference transcript; the|              |               |                   |                   |                   |                   |
|              | borders of the prediction    |              |               |                   |                   |                   |                   |
|              | cannot fall inside the       |              |               |                   |                   |                   |                   |
|              | introns of the reference.    |              |               |                   |                   |                   |                   |
|              | The prediction transcript    |              |               |                   |                   |                   |                   |
|              | can bridge multiple exons    |              |               |                   |                   |                   |                   |
|              | of the reference model.      |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **G**        | Generic match of a           | False        | True          | > 0%, > 0%, > 0%  | 0%                | **g**             | **Alternative     |
| ("O" before  | multiexonic prediction       |              |               |                   |                   |                   | splicing**        |
| release 1)   | transcript versus a          |              |               |                   |                   |                   |                   |
|              | monoexonic reference.        |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **i**        | Monoexonic prediction        | True         | False         | 0%                | 0%                | **ri**            | **Intronic**      |
|              | completely contained within  |              |               |                   |                   |                   |                   |
|              | one intron of the reference  |              |               |                   |                   |                   |                   |
|              | transcript.                  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **I**        | Prediction completely        | True         | True          | 0%                | 0%                | **rI**            | **Intronic**      |
|              | contained within the introns |              |               |                   |                   |                   |                   |
|              | of the reference transcript. |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **rI**       | Reference completely         | True         | True          | 0%                | 0%                | **I**             | **Intronic**      |
|              | contained within the introns |              |               |                   |                   |                   |                   |
|              | of the prediction transcript.|              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **ri**       | Reverse intron transcript -  | False        | True          | 0%                | 0%                | **i**             | **Intronic**      |
|              | the monoexonic reference is  |              |               |                   |                   |                   |                   |
|              | completely contained within  |              |               |                   |                   |                   |                   |
|              | one intron of the prediction |              |               |                   |                   |                   |                   |
|              | transcript.                  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **f**        | Fusion - this special code   | NA           | NA            | **> 10%**, NA, NA | **> 0%**, NA, NA  | NA                | **Fusion**        |
|              | is applied when a prediction |              |               |                   |                   |                   |                   |
|              | intersects more than one     |              |               |                   |                   |                   |                   |
|              | reference transcript. To be  |              |               |                   |                   |                   |                   |
|              | considered for fusions,      |              |               |                   |                   |                   |                   |
|              | candidate references must    |              |               |                   |                   |                   |                   |
|              | **either** share at least one|              |               |                   |                   |                   |                   |
|              | splice junction with the     |              |               |                   |                   |                   |                   |
|              | prediction, **or** have at   |              |               |                   |                   |                   |                   |
|              | least 10% of its bases       |              |               |                   |                   |                   |                   |
|              | recalled. If two or more     |              |               |                   |                   |                   |                   |
|              | reference transcripts fit    |              |               |                   |                   |                   |                   |
|              | these constraints, then the  |              |               |                   |                   |                   |                   |
|              | prediction model is          |              |               |                   |                   |                   |                   |
|              | classified as a **fusion**.  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **e**        | Single exon transcript       | True         | False         | > 0%, > 0%,       | 0%                | **G**             | **Fragment**      |
|              | overlapping *one* reference  |              |               | between 0 and 100%|                   |                   |                   |
|              | exon and at least 10 bps of a|              |               |                   |                   |                   |                   |
|              | reference intron, indicating |              |               |                   |                   |                   |                   |
|              | a possible pre-mRNA fragment.|              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **x**        | Monoexonic match on the      | NA           | False         | >= 0%             | 0%                | **x** or **X**    | **Fragment**      |
|              | *opposite* strand.           |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **X**        | Multiexonic match on the     | NA           | True          | >= 0%             | 0%                | **x** or **X**    | **Fragment**      |
|              | *opposite* strand.           |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **p**        | The prediction is on the same| NA           | NA            | 0%                | 0%                | **p**             | **No overlap**    |
|              | strand of a neighbouring but |              |               |                   |                   |                   |                   |
|              | non-overlapping transcript.  |              |               |                   |                   |                   |                   |
|              | Probable polymerase run-on.  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **P**        | The prediction is on the     | NA           | NA            | 0%                | 0%                | **P**             | **No overlap**    |
|              | *opposite* strand of a       |              |               |                   |                   |                   |                   |
|              | neighbouring but             |              |               |                   |                   |                   |                   |
|              | non-overlapping transcript.  |              |               |                   |                   |                   |                   |
|              | Probable polymerase run-on.  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+
| **u**        | Unknown - no suitable model  | NA           | NA            | 0%                | 0%                | NA                | **No overlap**    |
|              | has been found near enough   |              |               |                   |                   |                   |                   |
|              | the prediction to perform a  |              |               |                   |                   |                   |                   |
|              | comparison.                  |              |               |                   |                   |                   |                   |
+--------------+------------------------------+--------------+---------------+-------------------+-------------------+-------------------+-------------------+


Technical details
~~~~~~~~~~~~~~~~~

Mikado compare conceptualizes the reference annotation as a collection of interval trees, one per chromosome or scaffold, where each node corresponds to an array of genes at the location. The gene and transcript objects are stored separately. The location of each transcript model in the prediction is queried against the tree, with a padding (default 2kbps) to allow for neighouring but non-overlapping genes, and the transcript itself is subsequently compared with each reference transcript contained in the hits. Each comparison will yield precision, recall and F1 values for the nucleotide, splice junction and exonic levels, together with an associated class code. The best match for the prediction is selected for by choosing the comparison yielding the best splice junction F1 and the best nucleotide F1, in this order. If the prediction transcript overlaps two or more genes on the same strand, and for at least two it has one match each with either 10% nucleotide recall or junction recall over 0%, it is deemed as a fusion_ event, and its line in the tmap_ file will report the best match against each of the fused genes, separated by comma.

Each calculated match against a reference transcript is stored as a potential *best match* for the reference transcript. At the end of the run, the hits for each reference transcript will be ordered using the following function:

.. code-block:: python
    :linenos:

    @staticmethod
    def result_sorter(result):

        """
        Method to sort the results for the refmap. Order:
        - CCode does not contain "x", "P", "p" (i.e. fragments on opposite strand or
        polymerase run-on fragments)
        - Exonic F1 (e_f1)
        - Junction F1 (j_f1)
        - "f" in ccode (i.e. transcript is a fusion)
        - Nucleotide F1 (n_f1)

        :param result: a resultStorer object
        :type result: ResultStorer
        :return: (int, float, float, float)
        """

        bad_ccodes = ["x", "X", "P", "p"]
        bad_ccodes = set(bad_ccodes)

        orderer = (len(set.intersection(bad_ccodes, set(result.ccode))) == 0,
                   result.j_f1, result.e_f1,
                   result.n_f1,
                   "f" in result.ccode)

        return orderer

This function is used to select both for the best match *for the transcript*, as well as to select among these matches for the best match *for the gene*.

The interval tree data structure is created using Cython code originally part of the `bx-python <https://bitbucket.org/james_taylor/bx-python/overview>`_, kindly provided by `Dr. Taylor <mailto:james@taylorlab.org>`_ for modification and inclusion in Mikado. The code has been slightly modified for making it Python3 compliant.

The comparison code is written in Cython and is crucial during the :ref:`picking phase of Mikado <pick>`, not just for the functioning of the comparison utility.