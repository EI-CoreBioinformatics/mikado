.. _SQLAlchemy: http://www.sqlalchemy.org/
.. _Portcullis: https://github.com/maplesond/portcullis
.. _BED12: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

.. _configure-scoring-tutorial:

How to create a scoring configuration file
==========================================

The current version of Mikado relies upon the experimenter to determine how desirable transcripts should look like, and
how to prioritise across different possible isoforms. These instructions are provided through **scoring configuration files**,
whose detailed description can be found in the :ref:`general documentation <scoring_files>`. In this section, we will
expose a general way to create your own configuration file.

First step: defining transcript requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in the process is to define the minimum attributes of a transcript that should be retained in your annotation. For example, if the experimenter is interested only in coding transcripts and would like to ignore any non-coding RNA, this strict requirement would have to be encoded in this section. Transcripts that do not pass this preliminary filter are completely excluded from any subsequent analysis.

In general, this filter should be very gentle - being too stringent at this stage risks completely throwing out all the transcripts present in a given locus. It is rather preferable to strongly penalise dubious transcripts later at the scoring stage, so that they will be kept in the final annotation if no alternative is present, but discarded in most cases as there will be better alternatives.

The requirements block is composed of two sub-sections:
    - a *parameters* section, which details which :ref:`metrics <Metrics>` will have to be evaluated, and according to which :ref:`operators <operators>`.
    - an *expression* section, detailing how the various parameters have to be considered together.

As an example, let us imagine a quite stringent experiment in which we are interested only in transcripts that respect the following conditions:
    - minimum transcript length of 300 bps:
      -  *cdna_length: {operator: ge, value: 300}*
    - if they are multiexonic, at least one of their junctions must be validated by a junction checker such as Portcullis
      - *exon_num: {operator: ge, value: 2}*
      - *verified_introns_num: {operator: gt, value: 0}*
    - if they are multiexonic, their introns must be within 5 and 2000 bps:
      - *exon_num: {operator: ge, value: 2}*
      - *min_intron_length: {operator: ge, value: 5}*
      - *max_intron_length: {operator: le, value: 2000}*
    - if they are monoexonic, they must be coding:
      - *exon_num: {operator: eq, value: 1}*
      - *combined_cds_length: {operator: gt, value: 0}*

Having defined the parameters, we can now put them together in an *expression*:

    cdna_length and ((exon_num and verified_introns_num and min_intron_length and max_intron_length) or (exon_num and combined_cds_length))

Notice that we have used a property twice (*exon_num*). This would confuse Mikado and cause an error. In order to tell the program that these are actually two different values, we can prepend a suffix, starting with a mandatory dot sign:
    - exon_num.**multi**:  {operator: ge, value: 2}
    - exon_num.**mono**: {operator: eq, value: 1}

The expression now becomes:

    cdna_length and ((**exon_num.multi** and verified_introns_num and min_intron_length and max_intron_length) or (**exon_num.mono** and combined_cds_length))

Putting it all together, this is how the section in the configuration file would look like:

.. code-block:: yaml
  :emphasize-lines: 2,5
  :lineno-start: 5

  requirements:
      expression:
      - cdna_length and ((exon_num and verified_introns_num and min_intron_length
      - and max_intron_length) or (exon_num and combined_cds_length))
      parameters:
      - cdna_length: {operator: ge, value: 300}
      - exon_num.multi:  {operator: ge, value: 2}
      - verified_introns_num: {operator: gt, value: 0}
      - min_intron_length: {operator: ge, value: 5}
      - max_intron_length: {operator: le, value: 2000}
      - exon_num.mono: {operator: eq, value: 1}
      - combined_cds_length: {operator: gt, value: 0}

:warning: The example in this section is more stringent than the standard selection provided by the included scoring files.


Second step: prioritising transcripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After removing transcripts which are not good enough for our annotation, Mikado will analyse any remaining models and assign a score each. How to score models in Mikado is, explicitly, a procedure left to the experimenter, so to allow specific tailoring for each different species. In our own experiments, we have abided by these principles:

- Good transcripts should preferentially be protein coding, with homology to known proteins in other species, and sport a complete protein.
- Good coding transcripts should contain only one ORF, not multiple, and such an ORF should be complete. The total length of the CDS should be within 60 and 80% of the transcript length, ideally (with the value changing by species, on the basis of available data).
- All else equal, good coding transcripts should have the longest ORF among those present in the locus.
- Good coding transcripts should have a defined UTR, on both sides; however, if the UTR goes beyond a certain limit, the transcript should be penalised instead. For 5'UTR, we preferentially look at transcripts with at most four UTR exons, and preferentially **two**, for a total length of ideally 100 bps and maximally of 2500. For 3'UTR, based on literature and the phenomenon of nonsense mediated decay (NMD), we look for transcripts with at most **two** UTR exons and ideally **one**; the total length of this UTR should be ideally of 200 bps, and at most of 2,500.
- Good transcripts should be multiexonic, but if they are, they should have at least some of their junctions confirmed by Portcullis.
- The maximum distance between the stop codon and the last junction in the transcript should not exceed 55 bps (as discovered by studies on NMD).

These rules have to be encoded in a format that Mikado understands,



