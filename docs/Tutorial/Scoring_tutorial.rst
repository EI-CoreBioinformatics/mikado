.. _SQLAlchemy: http://www.sqlalchemy.org/
.. _Portcullis: https://github.com/maplesond/portcullis
.. _BED12: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

.. _configure-scoring-tutorial:

Tutorial: how to create a scoring configuration file
----------------------------------------------------

The current version of Mikado relies upon the experimenter to determine what desirable transcripts should look like, and
how to prioritise across different possible isoforms. These instructions are provided through **scoring configuration files**,
whose detailed description can be found in the :ref:`general documentation <scoring_files>`. In this section, we will
expose a general way to create your own configuration file.

The purpose of scoring files: introductory note
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mikado is an annotation pipeline, and as such, it does not have, at its core, the mission to represent every possible expressed transcript. Some sequenced transcripts may originate from transcriptional errors, exist only transiently and may have no functional role. When performing a sequencing experiment, we create a “snapshot” of the transient transcriptional state of the cell: inclusive of these erroneous events, along with immature transcripts and artifacts arising from misassemblies, fragmentation, genomic inclusion, etc.

In Mikado, we make a choice regarding what transcripts we want to retain in our annotation (as is the case with any genome annotation). For example, as you will see in this tutorial, the standard configuration penalises transcripts with retained CDS introns, and transcripts that are NMD targets. This choice does not mean that we think that these transcripts are artifacts; rather, it signifies that we prioritise those transcripts that are more likely to represent functionally active genes. Annotators will differ on the point: for example, the human reference annotation retains and marks these events, rather than discarding them.

Mikado allows the experimenter to make these choices simply and transparently. Our pre-defined configuration files strive to replicate the choices made by annotators over the years, and thus allow to replicate - as much as possible - the reference annotations of various species starting from partial sequencing data. However, if as an experimenter you are interested in a more stringent approach - say, only coding transcripts with a complete ORF and a UTR - or you would like instead to perform only a minimal cleaning up - say, just discarding obvious NMD targets - Mikado will allow you to do so. This tutorial will show you how.

Obtaining pre-defined scoring files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mikado comes with pre-packaged scoring files for multiple species. These can be found in the installation directory, under "configuration/scoring_files"; or, when launching mikado configure, you can request to copy a template file to the working directory ("--copy-scoring" command flag). In the rest of the tutorial, however, we will presume that no suitable scoring file exists and that a new one has to be created from scratch.

.. _scoring-tutorial-first-reqs:

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

.. code-block:: yaml

    expression:
        cdna_length and ((exon_num and verified_introns_num and min_intron_length and max_intron_length) or (exon_num and combined_cds_length))

Notice that we have used a property twice (*exon_num*). This would confuse Mikado and cause an error. In order to tell the program that these are actually two different values, we can prepend a suffix, starting with a mandatory dot sign:

    - exon_num.\ **multi**:  {operator: ge, value: 2}
    - exon_num.\ **mono**: {operator: eq, value: 1}

.. note::

    The suffix must be ASCII alphanumeric in character. Therefore, the following suffices are valid:

        - .mono
        - .mono1
        - .mono_first

.. code-block:: yaml

    The expression now becomes:
        cdna_length and ((exon_num.multi and verified_introns_num and min_intron_length and max_intron_length) or (exon_num.mono and combined_cds_length))

.. warning::

    if no expression is provided, Mikado will create a default one by connecting all the parameters with an and. This will make life simpler for simple cases (e.g. we only have a couple of parameters we want to check). In complex, conditional scenarios like this one, however, this might well lead to discarding all transcripts!

Putting it all together, this is how the section in the configuration file would look like:

.. code-block:: yaml
  :emphasize-lines: 2,5
  :lineno-start: 1

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

.. warning:: 

    The example in this section is more stringent than the standard selection provided by the included scoring files.

.. _scoring-tutorial-second-prior:

Second step: prioritising transcripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After removing transcripts which are not good enough for our annotation, Mikado will analyse any remaining models and assign each a score. How to score models in Mikado is, explicitly, a procedure left to the experimenter, to allow specific tailoring for each different species. In our own experiments, we have abided by the following principles:

1. Good transcripts should preferentially be protein coding, with homology to known proteins in other species, and sport both start and stop codon.
2. Good coding transcripts should contain only one ORF, not multiple; if they have multiple, most of the CDS should be within the primary.
3. The total length of the CDS should be within 60 and 80% of the transcript length, ideally (with the value changing by species, on the basis of available data).
4. All else equal, good coding transcripts should have a long ORF, contain most of the coding bases in the locus, and have that most of their introns are between coding exons.
5. All else equal, good transcripts should be longer and have more exons; however, there should be no preference between mono- and di-exonic transcripts.
6. Good coding transcripts should have a defined UTR, on both sides; however, if the UTR goes beyond a certain limit, the transcript should be penalised instead. For 5'UTR, we preferentially look at transcripts with at most four UTR exons, and preferentially **two**, for a total length of ideally 100 bps and maximally of 2500. For 3'UTR, based on literature and the phenomenon of nonsense mediated decay (NMD), we look for transcripts with at most **two** UTR exons and ideally **one**; the total length of this UTR should be ideally of 200 bps, and at most of 2500.
7. Multiexonic transcripts should have at least some of their junctions confirmed by Portcullis, ideally all of them. Ideally and all else equal, they should contain all of the verified junctions in the locus.
8. The distance between the stop codon and the last junction in the transcript should be the least possible, and in any case, not exceed 55 bps (as discovered by studies on NMD).

The first step is to associate each of these requirements with the proper :ref:`metric <Metrics>`. In order:

1. Good transcripts should preferentially be protein coding, with a good BLAST coverage of homologous proteins, and sport start and stop codon:

    - snowy_blast_score: look for the maximum value
    - is_complete: look for ``true``
    - has_start_codon: look for ``true``
    - has_stop_codon: look for ``true``

Looking at the documentation on :ref:`scoring files <scoring_files>`, we can write it down thus:

.. code-block:: yaml

        - snowy_blast_score: {rescaling: max}
        - is_complete: {rescaling: target, value: true}
        - has_start_codon: {rescaling: target, value: true}
        - has_stop_codon: {rescaling: target, value: true}

Applying the same procedure to the rest of the conditions:

2. Good coding transcripts should contain only one ORF, not multiple; if they have multiple, most of the CDS should be within the primary.

    - number_internal_orfs: look for a target of 1
    - cds_not_maximal: look for the **minimum** value
    - cds_not_maximal_fraction: look for the **minimum** value

.. code-block:: yaml

        - number_internal_orfs: {rescaling: target, value: 1}
        - cds_not_maximal: {rescaling: min}
        - cds_not_maximal_fraction: {rescaling: min}

3. The total length of the CDS should be within 60 and 80% of the transcript length, ideally (with the value changing by species, on the basis of available data).

    - selected_cds_fraction: look for a target of x *(where x depends on the species and is between 0 and 1)*, for example, let us set it to 0.7

.. code-block:: yaml

        - selected_cds_fraction: {rescaling: target, value: 0.7}

4. All else equal, good coding transcripts should have a long ORF, contain most of the coding bases in the locus, and have that most of their introns are between coding exons.

    - cdna_length: look for the maximum value
    - selected_cds_length: look for the maximum value
    - selected_cds_intron_fraction: look for the maximum value

.. code-block:: yaml

        - selected_cds_length: {rescaling: max}
        - selected_cds_intron_fraction: {rescaling: max}
        - selected_cds_intron_fraction: {rescaling: max}

5. All else equal, good transcripts should be longer and have more exons; however, there should be no preference between mono- and di-exonic transcripts.

    - cdna_length: look for the maximum value
    - exon_num: look for the maximum value, ignore for any transcript with one or two exons.

.. code-block:: yaml

        - cdna_length: {rescaling: max}
        - exon_num: {rescaling: max, filter: {operator: ge, value: 3}

6. Good coding transcripts should have a defined UTR, on both sides; however, if the UTR goes beyond a certain limit, the transcript should be penalised instead.

    - For 5'UTR, we preferentially look at transcripts with at most three UTR exons, and preferentially **two**, for a total length of ideally 100 bps and maximally of 2500.

        * five_utr_num: look for a target of 2, ignore anything with four or more 5' UTR exons
        * five_utr_length: look for a target of 100, ignore anything with 2500 or more bps
    - For 3'UTR, based on literature and the phenomenon of nonsense mediated decay (NMD), we look for transcripts with at most **two** UTR exons and ideally **one**; the total length of this UTR should be ideally of 200 bps, and at most of 2500.

        * three_utr_num: look for a target of 1, ignore anything with three or more 3'UTR exons
        * three_utr_length: look for a target of 200, ignore anything with 2500 bps or more

.. code-block:: yaml

        - five_utr_num: {rescaling: target, value: 2, filter: {operator: lt, value: 4}}
        - five_utr_length: {rescaling: target, value: 100, filter: {operator: le, value: 2500}}
        - three_utr_num: {rescaling: target, value: 1, filter: {operator: lt, value: 3}}
        - three_utr_length: {rescaling: target, value: 200, filter: {operator: lt, value: 2500}}

7. Multiexonic transcripts should have at least some of their junctions confirmed by Portcullis, ideally all of them. Ideally and all else equal, they should contain most of the verified junctions in the locus.

    - proportion_verified_introns_inlocus: look for the maximum value
    - non_verified_introns_num: look for the minimum value

.. code-block:: yaml

        - proportion_verified_introns_inlocus: {rescaling: max}
        - non_verified_introns_num: {rescaling: min}

8. The distance between the stop codon and the last junction in the transcript should be the least possible, and in any case, not exceed 55 bps (as discovered by studies on NMD).

    - end_distance_from_junction: look for the minimum value, discard anything over 55

.. code-block:: yaml

        - end_distance_from_junction: {rescaling: min, filter: {operator: lt, value: 55}}

Putting everything together:

.. code-block:: yaml

    scoring:
        - snowy_blast_score: {rescaling: max}
        - is_complete: {rescaling: target, value: true}
        - has_start_codon: {rescaling: target, value: true}
        - has_stop_codon: {rescaling: target, value: true}
        - number_internal_orfs: {rescaling: target, value: 1}
        - cds_not_maximal: {rescaling: min}
        - cds_not_maximal_fraction: {rescaling: min}
        - selected_cds_fraction: {rescaling: target, value: 0.7}
        - selected_cds_length: {rescaling: max}
        - selected_cds_intron_fraction: {rescaling: max}
        - selected_cds_intron_fraction: {rescaling: max}
        - cdna_length: {rescaling: max}
        - exon_num: {rescaling: max, filter: {operator: ge, value: 3}
        - five_utr_num: {rescaling: target, value: 2, filter: {operator: lt, value: 4}}
        - five_utr_length: {rescaling: target, value: 100, filter: {operator: le, value: 2500}}
        - three_utr_num: {rescaling: target, value: 1, filter: {operator: lt, value: 3}}
        - three_utr_length: {rescaling: target, value: 200, filter: {operator: lt, value: 2500}}
        - proportion_verified_introns_inlocus: {rescaling: max}
        - non_verified_introns_num: {rescaling: min}
        - end_distance_from_junction: {rescaling: min, filter: {operator: lt, value: 55}}

.. _scoring-tutorial-third-reqs:

Third step: defining acceptable alternative splicing events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After selecting a primary transcript for the locus, we have to define what would make a transcript **inherently** unacceptable as an alternative splicing event. This is done in a similar way to how we defined the :ref:`minimal requirements for all transcripts <scoring-tutorial-first-reqs>`.

.. warning::
    Keep in mind that this section defines the **inherent** requirements. **Relative** requirements, such as acceptable class codes, percentage of the score of the primary transcript, etc., are defined in the general configuration file, :ref:`under the "alternative_splicing" section <configure-alternative-splicing>`. By default, we also control whether to accept or refuse retained intron events there, rather than here.

Throughout our experiments, we have defined this section quite gently; potential candidates are discarded more due to their relationship to the primary transcript (:ref:`class code <ccodes>`, score percentage, etc.) rather than due to some inherent defect. This is how we generally selected:

- Minimum cDNA length of 200
- Combined UTR length less than 2500 bps
- No suspicious splicing event (ie junctions that would be canonical if ported on the opposite strand)

.. code-block:: yaml

    as_requirements:
      expression: [cdna_length and three_utr_length and five_utr_length and utr_length and suspicious_splicing]
      parameters:
        cdna_length: {operator: ge, value: 200}
        utr_length: {operator: le, value: 2500}
        five_utr_length: {operator: le, value: 2500}
        three_utr_length: {operator: le, value: 2500}
        suspicious_splicing: {operator: ne, value: true}

Fourth step: defining potential fragments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The final step in the selection process is to detect and discard potential transcript fragments present in the neighbourhood of good loci. Usually these originate by mismappings or polymerase run-ons, and can be easily identified "by eye" as short, non- or minimally coding transcripts near better looking loci. Mikado will use the requirements defined in this section to identify such spurious loci, and discard them.

.. note:: 
    The maximum distance between loci, for them to be considered for this step, is defined :ref:`in the general configuration file <clustering_specifics>` by the "flank" parameter. Any locus beyond this distance will **not** be evaluated as a potential fragment.

For our experiments, in general, this is how we defined potential fragments:

- If multiexonic:

    * Shorter than 300 bps
    * Or with an ORF shorter than 300 bps
    
- If monoexonic:

    * Non-coding and without any BLAST homology
    * Coding with an ORF lower than 600 bps

In the format understood by Mikado:

.. code-block:: yaml

    not_fragmentary:
        expression: [((exon_num.multi and (cdna_length.multi or selected_cds_length.multi)), or, (exon_num.mono and ((snowy_blast_score and selected_cds_length.zero)  or selected_cds_length.mono)))]
        parameters:
            selected_cds_length.zero: {operator: gt, value: 300} # 600
            exon_num.multi: {operator: gt, value: 2}
            cdna_length.multi: {operator: ge, value: 300}
            selected_cds_length.multi: {operator: gt, value: 250}
            exon_num.mono: {operator: eq, value: 1}
            snowy_blast_score: {operator: gt, value: 0}  # 0.3
            selected_cds_length.mono: {operator: gt, value: 600} # 900
            exon_num.mono: {operator: le, value: 2}

Pointing Mikado at the new configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the new scoring file is complete, we can point Mikado pick at it in two ways:

- Either transiently, with the "--scoring-file" switch, followed by the file name.
- Or in the configuration file for the project, by putting the file name under :ref:`the pick/scoring_file field <misc-settings>`.

When Mikado pick will be launched, it will validate - before starting - the validity of the scoring file. Common mistakes:

- Using a metric which does not exist.
- Using an invalid combination of "operator", "value" and "rescaling" parameters; for example using a value of "true" with "gt", ie "greater than" (see the :ref:`section on operators <operators>`).
- Using an invalid connector in the "expression" statements: only "and", "or", "xor", "not" and brackets are accepted (see :ref:`the requirements section <requirements-section>`)

Mikado should emit an error that will help you understand how to correct the issue.
