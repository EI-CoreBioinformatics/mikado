.. _scoring_files:

Scoring files
=============

Mikado employs user-defined configuration files to define the desirable features in genes. These files are in TOML, YAML or JSON format (default YAML) and are composed of five sections:

  #. a *requirements* section, specifying the minimum requirements that a transcript must satisfy to be considered as valid. **Any transcript failing these requirements will be scored at 0 and purged.**
  #. a *cds_requirements* section, specifying minimal conditions for a transcript to retain its CDS. If a transcript fails this check, it will be stripped of its coding component. If this section is not provided, the default will be to copy the *requirements* section above (in practice disabling it).
  #. a *not_fragmentary* section, specifying the minimum requirements that the primary transcript of a locus has to satisfy in order for the locus **not** to be considered as a putative fragment.
  #. an *as_requirements* section, which specifies the minimum requirements for transcripts for them to be considered as possible valid alternative splicing events.
  #. a *scoring* section, specifying which features Mikado should look for in transcripts, and how each of them will be weighted.

Conditions are specified using a strict set of :ref:`available operators <operators>` and the values they have to consider.

.. important:: Although at the moment Mikado does not offer any method to establish machine-learning based scoring configurations, it is a topic we plan to investigate in the future. Mikado already supports `Random Forest Regressors as scorers through Scikit-learn <http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html>`_, but we have yet to devise a proper way to create such regressors.

We provide a guide on how to write your own scoring files in a :ref:`separate tutorial <configure-scoring-tutorial>`.

.. _operators:

Operators
~~~~~~~~~

Mikado allows the following operators to express a relationship inside the scoring files:

* *eq*: equal to (:math:`=`). Valid for comparisons with numbers, boolean values, and strings.
* *ne*: different from (:math:`\neq`). Valid for comparisons with numbers, boolean values, and strings.
* *lt*: less than (:math:`<`). Valid for comparisons with numbers.
* *gt*: greater than (:math:`>`). Valid for comparisons with numbers.
* *le*: less or equal than (:math:`\le`). Valid for comparisons with numbers.
* *ge*: greater or equal than (:math:`\ge`). Valid for comparisons with numbers.
* *in*: member of (:math:`\in`). Valid for comparisons with arrays or sets.
* *not in*: not member of (:math:`\notin`). Valid for comparisons with arrays or sets.
* *within*: value comprised in the range of the two values, inclusive.
* *not within*: value *not* comprised in the range of the two values, inclusive.

Mikado will fail if an operator not present on this list is specified, or if the operator is assigned to compare against the wrong data type (eg. *eq* with an array).

.. _requirements-section:

The "requirements", "cds_requirements", "as_requirements" and "not_fragmentary" sections
----------------------------------------------------------------------------------------

These sections specifies the minimum requirements for a transcript at various stages.

* A transcript failing to pass the *requirements* check will be discarded outright (if "purge" is selected) or given a
  score of 0 otherwise.
* After passing the *rquirements* section, if the transcript is coding, Mikado will check whether its CDS passes the
  requirements specified in *cds_requirements*. If the check fails, the transcripts will be **stripped of its CDS**, and
  will only be considered further as a non-coding transcript. *This check can be used to properly consider transcripts
  that have a suspicious coding structure - e.g. a single coding exon in a transcript with five or more exons, or a low
  Coding Potential score coming from a third-party tool*.
* If a transcript has not been selected as the primary transcript of a locus, it has to pass the *as_requirements* check
  to be considered as a valid alternative splicing event.
* Finally, after loci have been defined, the primary transcript of each locus will be checked against the
  *not_fragmentary* expression. Any locus failing this check will be marked as a potential fragment and, if in the
  vicinity of other loci, might be purged out of the final output or clearly marked as a fragment (depending on whether
  the *purge* switch is set to true or false, respectively).

**It is strongly advised to use lenient parameters in the first requirements section**, as failing to do so might result
in discarding whole loci. Typically, transcripts filtered at this step should be obvious artefacts, eg monoexonic
transcripts produced by RNA-Seq with a total length lower than the *library* fragment length.

Each of these sections follows the same template, and they are composed by two parts:

* *parameters*: a list of the metrics to be considered. Each metric can be considered multiple times, by suffixing
  it with a ".<id>" construct (eg cdna_length.*mono* vs. cdna_length.*multi* to distinguish two uses of the cdna_length
  metric - once for monoexonic and once for multiexonic transcripts). Any parameter which is not a :ref:`valid metric
  name <Metrics>`, after removal of the suffix, **will cause an error**. Parameters have to specify the following:

    * a *value* that the metric has to be compared against
    * an *operator* that specifies the target operation. See :ref:`the operators section <operators>`.

* *expression*: a string *array* that will be compiled into a valid boolean expression. All the metrics present in the
  expression string **must be present in the parameters section**. If an unrecognized metric is present, Mikado will
  crash immediately, complaining that the scoring file is invalid. Apart from brackets, Mikado accepts only the
  following boolean operators to chain the metrics:

    * *and*
    * *or*
    * *not*
    * *xor*

.. hint:: if no *expression* is specified, Mikado will construct one by chaining all the provided parameters with and
   *and* operator. Most of the time, this would result in an unexpected behaviour - ie Mikado assigning a score of 0 to
   most transcripts. It is **strongly advised** to explicitly provide a valid expression.

As an example, the following snippet replicates a typical requirements section found in a scoring file:

.. code-block:: yaml

    requirements:
      expression: [((exon_num.multi and cdna_length.multi and max_intron_length and min_intron_length), or,
        (exon_num.mono and cdna_length.mono))]
      parameters:
        cdna_length.mono: {operator: gt, value: 50}
        cdna_length.multi: {operator: ge, value: 100}
        exon_num.mono: {operator: eq, value: 1}
        exon_num.multi: {operator: gt, value: 1}
        max_intron_length: {operator: le, value: 20000}
        min_intron_length: {operator: ge, value: 5}

In the parameters section, we ask for the following:

        * *exon_num.mono*: monoexonic transcripts must have one exon ("eq")
        * *exon_num.multi*: multiexonic transcripts must have more than one exon ("gt")
        * *cdna_length.mono*: monoexonic transcripts must have a length greater than 50 bps (the ".mono" suffix is
          arbitrary, as long as it is unique for all calls of *cdna_length*)
        * *cdna_length.multi*: multiexonic transcripts must have a length greater than or equal to 100 bps (the ".multi"
          suffix is arbitrary, as long as it is unique for all calls of *cdna_length*)
        * *max_intron_length*: multiexonic transcripts should not have any intron longer than 200,000 bps.
        * *min_intron_length*: multiexonic transcripts should not have any intron smaller than 5 bps.

The *expression* field will be compiled into the following expression::

        (exon_num > 1 and cdna_length >= 100 and max_intron_length <= 200000 and min_intron_length >= 5) or (exon_num == 1 and cdna_length > 50)


Any transcript for which the expression evaluates to ``false`` will be assigned a score of 0 outright and discarded,
unless the user has chosen to disable the purging of such transcripts.

.. _scoring-section:

The scoring section
~~~~~~~~~~~~~~~~~~~

This section specifies which metrics will be used by Mikado to score the transcripts. Each metric to be used is
specified as a subsection of the configuration, and will have the following attributes:

* *rescaling*: the only compulsory attribute. It specifies the kind of scoring that will be applied to the metric, and
  it has to be one of "max", "min", or "target". See :ref:`the explanation on the scoring algorithm <scoring_algorithm>`
  for details.
* *value*: compulsory if the chosen rescaling algorithm is "target". This should be either a number or a boolean value.
* *multiplier*: the weight assigned to the metric in terms of scoring. This parameter is optional; if absent, as it is
  in the majority of cases, Mikado will consider the multiplier to equal to 1. This is the :math:`w_{m}` element in the
  :ref:`equations above <scoring_algorithm>`.
* *filter*: It is possible to specify a filter which the metric has to fulfill to be considered for scoring, eg,
  "cdna_length >= 200". If the transcript fails to pass this filter, the score *for this metric only* will be set to 0.
  **The filter can apply to a different metric**, so it is possible to e.g. assign a score of 0 for, say, "combined_cds"
  to any transcript whose "cdna_length" is, say, below 150 bps.
  A "filter" subsection has to specify the following:

    * *operator*: the operator to apply for the boolean expression. See the :ref:`relative section <operators>`.
    * *value*: value that will be used to assess the metric.
    * *metric*: *optional*. The metric that this filter refers to. If omitted, this will be identical to the metric
      under examination.

.. hint:: the purpose of the *filter* section is to allow for fine-tuning of the scoring mechanism; ie it allows to
          penalise transcripts with undesirable qualities (eg a possible retained intron) without discarding them
          outright. As such, it is a less harsh version of the :ref:`requirements section <requirements-section>` and
          it is the preferred way of specifying which transcript features Mikado should be wary of.

For example, this is a snippet of a scoring section:

.. code-block:: yaml

    scoring:
        blast_score: {rescaling: max}
        cds_not_maximal: {rescaling: min}
        combined_cds_fraction: {rescaling: target, value: 0.8, multiplier: 2}
        five_utr_length:
            filter: {operator: le, value: 2500}
            rescaling: target
            value: 100
        end_distance_from_junction:
            filter: {operator: lt, value: 55}
            rescaling: min
        non_verified_introns_num:
            rescaling: max
            multiplier: -10
            filter:
                operator: gt
                value: 1
                metric: exons_num


Using this snippet as a guide, Mikado will score transcripts in each locus as follows:

* Assign a full score (one point, as no multiplier is specified) to transcripts which have the greatest *blast_score*
* Assign a full score (one point, as no multiplier is specified) to transcripts which have the lowest amount of CDS
  bases in secondary ORFs (*cds_not_maximal*)
* Assign a full score (**two points**, as a multiplier of 2 is specified) to transcripts that have a total amount of CDS
  bps approximating 80% of the transcript cDNA length (*combined_cds_fraction*)
* Assign a full score (one point, as no multiplier is specified) to transcripts that have a 5' UTR whose length is
  nearest to 100 bps (*five_utr_length*); if the 5' UTR is longer than 2,500 bps, this score will be 0
  (see the filter section)
* Assign a full score (one point, as no multiplier is specified) to transcripts which have the lowest distance between
  the CDS end and the most downstream exon-exon junction (*end_distance_from_junction*). If such a distance is greater
  than 55 bps, assign a score of 0, as it is a probable target for NMD (see the filter section).
* Assign a maximum penalty (**minus 10 points**, as a **negative** multiplier is specified) to the transcript with the
  highest number of non-verified introns in the locus.

  * Again, we are using a "filter" section to define which transcripts will be exempted from this scoring
    (in this case, a penalty)
  * However, please note that we are using the keyword **metric** in this section. This tells Mikado to check a
    *different* metric for evaluating the filter. Nominally, in this case we are excluding from the penalty any
    *monoexonic* transcript. This makes sense as a monoexonic transcript will never have an intron to be confirmed to
    start with.

.. note:: The possibility of using different metrics for the "filter" section is present from Mikado 1.3 onwards.

.. _Metrics:

Metrics
~~~~~~~

These are all the metrics available to quantify transcripts. The documentation for this section has been generated using
the :ref:`metrics utility <metrics-command>`.

Metrics belong to one of the following categories:

* **Descriptive**: these metrics merely provide a description of the transcript (eg its ID) and are not used for scoring.

* **cDNA**: these metrics refer to basic features of any transcript such as its number of exons, its cDNA length, etc.

* **Intron**: these metrics refer to features related to the number of introns and their lengths.

* **CDS**: these metrics refer to features related to the CDS assigned to the transcript.

* **UTR**: these metrics refer to features related to the UTR of the transcript. In the case in which a transcript has
  been assigned multiple ORFs, unless otherwise stated the UTR metrics will be derived only considering the *selected*
  ORF, not the combination of all of them.

* **Locus**: these metrics refer to features of the transcript in relationship to all other transcripts in its locus, eg
  how many of the introns present in the locus are present in the transcript. These metrics are calculated by Mikado during the picking phase, and as such their value can vary during the different stages as the transcripts are shifted to different groups.

* **External**: these metrics are derived from accessory data that is recovered for the transcript during the run time.
  Examples include data regarding the number of introns confirmed by external programs such as Portcullis, or the BLAST
  score of the best hits.

* **Attributes**: these metrics are extracted at runtime from attributes present in the input files. An example of this
  could be the TPM or FPKM values assigned to transcripts by rna expression analysis software.

.. hint:: Starting from version 1 beta8, Mikado allows to use externally defined metrics for the transcripts. These can
          be accessed using the keyword "external.<name of the metrics>" within the *configuration* file. See the
          :ref:`relevant section <external-metrics>` for details.

.. hint:: Starting from version 2, Mikado allows to use attribute defined metrics for the transcripts. These can be
          accessed using the keyword "attributes.<name of the metric>" within the *scoring* file. See the
          :ref:`relevant section <attributes-metrics>` for details.

.. important:: Starting from Mikado 1 beta 8, it is possible to use metrics with values between 0 and 1 directly as
               scores, without rescaling. This feature is available only for metrics whose values naturally lie between
               0 and 1, or that are boolean in nature.

.. topic:: Available metrics

+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| Metric name                         | Description                                               | Category    | Data type   | Usable raw   |
+=====================================+===========================================================+=============+=============+==============+
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| tid                                 | ID of the transcript - cannot be an undefined value.      | Descriptive | str         | False        |
|                                     | Alias of id.                                              |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| parent                              | Name of the parent feature of the transcript.             | Descriptive | str         | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| score                               | Numerical value which summarizes the reliability of the   | Descriptive | str         | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| external_scores                     | **SPECIAL** this Namespace contains all the information   | External    | Namespace   | True         |
|                                     | regarding external scores for the transcript. If an       |             |             |              |
|                                     | absent property is not defined in the Namespace, Mikado   |             |             |              |
|                                     | will set a default value of 0 into the Namespace and      |             |             |              |
|                                     | return it.                                                |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| alias                               | This property returns the alias of the transcript, if     | Descriptive | str         | False        |
|                                     | present, else its ID                                      |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| best_bits                           | Metric that returns the best BitS associated with the     | External    | float       | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| blast_identity                      | This metric will return the alignment identity for the    | External    | float       | True         |
|                                     | best BLAST hit according to the evalue. If no BLAST hits  |             |             |              |
|                                     | are available for the sequence, it will return 0.         |             |             |              |
|                                     | :return: :return:                                         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| blast_query_coverage                | This metric will return the **query** coverage for the    | External    | float       | True         |
|                                     | best BLAST hit according to the evalue. If no BLAST hits  |             |             |              |
|                                     | are available for the sequence, it will return 0.         |             |             |              |
|                                     | :return:                                                  |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| blast_score                         | Interchangeable alias for testing different blast-related | External    | float       | False        |
|                                     | scores. Current: best bit score.                          |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| blast_target_coverage               | This metric will return the **target** coverage for the   | External    | float       | True         |
|                                     | best BLAST hit according to the evalue. If no BLAST hits  |             |             |              |
|                                     | are available for the sequence, it will return 0.         |             |             |              |
|                                     | :return: :return:                                         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| canonical_intron_proportion         | This metric returns the proportion of canonical introns   | Intron      | float       | True         |
|                                     | of the transcript on its total number of introns.         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| cdna_length                         | This property returns the length of the transcript.       | cDNA        | int         | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| cds_disrupted_by_ri                 | This property describes whether the CDS is interrupted    | Locus       | bool        | True         |
|                                     | within a retained intron.                                 |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| cds_not_maximal                     | This property returns the length of the CDS excluded from | CDS         | int         | False        |
|                                     | the selected ORF.                                         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| cds_not_maximal_fraction            | This property returns the fraction of bases not in the    | CDS         | float       | True         |
|                                     | selected ORF compared to the total number of CDS bases in |             |             |              |
|                                     | the cDNA.                                                 |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_cds_fraction               | This property return the percentage of the CDS part of    | CDS         | float       | True         |
|                                     | the transcript vs. the cDNA length                        |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_cds_intron_fraction        | This property returns the fraction of CDS introns of the  | Locus       | float       | True         |
|                                     | transcript vs. the total number of CDS introns in the     |             |             |              |
|                                     | Locus. If the transcript is by itself, it returns 1.      |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_cds_length                 | This property return the length of the CDS part of the    | CDS         | int         | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_cds_locus_fraction         | This metric returns the fraction of CDS bases of the      | Locus       | float       | True         |
|                                     | transcript vs. the total of CDS bases in the locus.       |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_cds_num                    | This property returns the number of non-overlapping CDS   | CDS         | int         | False        |
|                                     | segments in the transcript.                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_cds_num_fraction           | This property returns the fraction of non-overlapping CDS | CDS         | float       | True         |
|                                     | segments in the transcript vs. the total number of exons  |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_utr_fraction               | This property returns the fraction of the cDNA which is   | UTR         | float       | True         |
|                                     | not coding according to any ORF. Complement of            |             |             |              |
|                                     | combined_cds_fraction                                     |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| combined_utr_length                 | This property return the length of the UTR part of the    | UTR         | int         | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| end_distance_from_junction          | This metric returns the cDNA distance between the stop    | CDS         | int         | False        |
|                                     | codon and the last junction in the transcript. In many    |             |             |              |
|                                     | eukaryotes, this distance cannot exceed 50-55 bps         |             |             |              |
|                                     | otherwise the transcript becomes a target of NMD. If the  |             |             |              |
|                                     | transcript is not coding or there is no junction          |             |             |              |
|                                     | downstream of the stop codon, the metric returns 0. This  |             |             |              |
|                                     | metric considers the combined CDS end.                    |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| end_distance_from_tes               | This property returns the distance of the end of the      | CDS         | int         | False        |
|                                     | combined CDS from the transcript end site. If no CDS is   |             |             |              |
|                                     | defined, it defaults to 0.                                |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| exon_fraction                       | This property returns the fraction of exons of the        | Locus       | float       | True         |
|                                     | transcript which are contained in the sublocus. If the    |             |             |              |
|                                     | transcript is by itself, it returns 1. Set from outside.  |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| exon_num                            | This property returns the number of exons of the          | cDNA        | int         | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| five_utr_length                     | Returns the length of the 5' UTR of the selected ORF.     | UTR         | float       | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| five_utr_num                        | This property returns the number of 5' UTR segments for   | UTR         | int         | False        |
|                                     | the selected ORF.                                         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| five_utr_num_complete               | This property returns the number of 5' UTR segments for   | UTR         | int         | False        |
|                                     | the selected ORF, considering only those which are        |             |             |              |
|                                     | complete exons.                                           |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| has_start_codon                     | Boolean. True if the selected ORF has a start codon.      | CDS         | bool        | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| has_stop_codon                      | Boolean. True if the selected ORF has a stop codon.       | CDS         | bool        | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| highest_cds_exon_number             | This property returns the maximum number of CDS segments  | CDS         | int         | False        |
|                                     | among the ORFs; this number can refer to an ORF           |             |             |              |
|                                     | *DIFFERENT* from the maximal ORF.                         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| highest_cds_exons_num               | Returns the number of CDS segments in the selected ORF    | CDS         | int         | False        |
|                                     | (irrespective of the number of exons involved)            |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| intron_fraction                     | This property returns the fraction of introns of the      | Locus       | float       | True         |
|                                     | transcript vs. the total number of introns in the Locus.  |             |             |              |
|                                     | If the transcript is by itself, it returns 1. Set from    |             |             |              |
|                                     | outside.                                                  |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| is_complete                         | Boolean. True if the selected ORF has both start and end. | CDS         | bool        | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| is_reference                        | Checks whether the transcript has been marked as          | External    | bool        | False        |
|                                     | reference by Mikado prepare                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| max_exon_length                     | This metric will return the length of the biggest exon in | cDNA        | int         | False        |
|                                     | the transcript.                                           |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| max_intron_length                   | This property returns the greatest intron length for the  | Intron      | int         | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| min_exon_length                     | This metric will return the length of the biggest exon in | cDNA        | int         | False        |
|                                     | the transcript.                                           |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| min_intron_length                   | This property returns the smallest intron length for the  | Intron      | int         | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| non_verified_introns_num            | This metric returns the number of introns of the          | External    | int         | False        |
|                                     | transcript which are not validated by external data.      |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| num_introns_greater_than_max        | This metric returns the number of introns greater than    | Intron      | int         | False        |
|                                     | the maximum acceptable intron size indicated in the       |             |             |              |
|                                     | constructor.                                              |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| num_introns_smaller_than_min        | This metric returns the number of introns smaller than    | Intron      | int         | False        |
|                                     | the mininum acceptable intron size indicated in the       |             |             |              |
|                                     | constructor.                                              |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| number_internal_orfs                | This property returns the number of ORFs inside a         | CDS         | int         | False        |
|                                     | transcript.                                               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| only_non_canonical_splicing         | This metric will return True if the canonical_number is 0 | Intron      | bool        | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| original_source                     | This property returns the original source assigned to the | Descriptive | str         | False        |
|                                     | transcript (before Mikado assigns its own final source    |             |             |              |
|                                     | value).                                                   |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| proportion_verified_introns         | This metric returns, as a fraction, how many of the       | External    | float       | True         |
|                                     | transcript introns are validated by external data.        |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| proportion_verified_introns_inlocus | This metric returns, as a fraction, how many of the       | Locus       | float       | True         |
|                                     | verified introns inside the Locus are contained inside    |             |             |              |
|                                     | the transcript. In loci without *any* verified introns,   |             |             |              |
|                                     | this metric will be set to 1.                             |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| retained_fraction                   | This property returns the fraction of the cDNA which is   | Locus       | float       | True         |
|                                     | contained in retained introns.                            |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| retained_intron_num                 | This property records the number of introns in the        | Locus       | int         | False        |
|                                     | transcripts which are marked as being retained. See the   |             |             |              |
|                                     | corresponding method in the sublocus class.               |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_cds_exons_fraction         | Returns the fraction of CDS segments in the selected ORF  | CDS         | float       | True         |
|                                     | (irrespective of the number of exons involved)            |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_cds_fraction               | This property calculates the fraction of the selected CDS | CDS         | float       | True         |
|                                     | vs. the cDNA length.                                      |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_cds_intron_fraction        | This property returns the fraction of CDS introns of the  | CDS         | float       | True         |
|                                     | selected ORF of the transcript vs. the total number of    |             |             |              |
|                                     | CDS introns in the Locus (considering only the selected   |             |             |              |
|                                     | ORF). If the transcript is by itself, it should return 1. |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_cds_length                 | This property calculates the length of the CDS selected   | CDS         | int         | False        |
|                                     | as best inside the cDNA.                                  |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_cds_locus_fraction         | This metric returns the fraction of CDS bases of the      | Locus       | float       | True         |
|                                     | transcript vs. the total of CDS bases in the locus.       |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_cds_num                    | This property calculates the number of CDS exons for the  | CDS         | int         | False        |
|                                     | selected ORF                                              |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_cds_number_fraction        | This property returns the proportion of best possible CDS | CDS         | float       | False        |
|                                     | segments vs. the number of exons. See                     |             |             |              |
|                                     | selected_cds_number.                                      |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_end_distance_from_junction | This metric returns the distance between the stop codon   | CDS         | int         | False        |
|                                     | and the last junction of the transcript. In many          |             |             |              |
|                                     | eukaryotes, this distance cannot exceed 50-55 bps,        |             |             |              |
|                                     | otherwise the transcript becomes a target of NMD. If the  |             |             |              |
|                                     | transcript is not coding or there is no junction          |             |             |              |
|                                     | downstream of the stop codon, the metric returns 0.       |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_end_distance_from_tes      | This property returns the distance of the end of the best | CDS         | int         | False        |
|                                     | CDS from the transcript end site. If no CDS is defined,   |             |             |              |
|                                     | it defaults to 0.                                         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| selected_start_distance_from_tss    | This property returns the distance of the start of the    | CDS         | int         | False        |
|                                     | best CDS from the transcript start site. If no CDS is     |             |             |              |
|                                     | defined, it defaults to 0.                                |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| snowy_blast_score                   | Metric that indicates how good a hit is compared to the   | External    | float       | False        |
|                                     | competition, in terms of BLAST similarities. As in        |             |             |              |
|                                     | SnowyOwl, the score for each hit is calculated by taking  |             |             |              |
|                                     | the coverage of the target and dividing it by (2 *        |             |             |              |
|                                     | len(self.blast_hits)). IMPORTANT: when splitting          |             |             |              |
|                                     | transcripts by ORF, a blast hit is added to the new       |             |             |              |
|                                     | transcript only if it is contained within the new         |             |             |              |
|                                     | transcript. This WILL screw up a bit the homology score.  |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| source_score                        | This metric returns a score that is assigned to the       | External    | float       | False        |
|                                     | transcript in virtue of its origin.                       |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| start_distance_from_tss             | This property returns the distance of the start of the    | CDS         | int         | False        |
|                                     | combined CDS from the transcript start site. If no CDS is |             |             |              |
|                                     | defined, it defaults to 0.                                |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| suspicious_splicing                 | This metric will return True if the transcript either has | Intron      | bool        | False        |
|                                     | canonical introns on both strands (probably a chimeric    |             |             |              |
|                                     | artifact between two neighbouring loci, or if it has no   |             |             |              |
|                                     | canonical splicing event but it would if it were assigned |             |             |              |
|                                     | to the opposite strand (probably a strand misassignment   |             |             |              |
|                                     | on the part of the assembler/predictor).                  |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| three_utr_length                    | Returns the length of the 5' UTR of the selected ORF.     |             | int         | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| three_utr_num                       | This property returns the number of 3' UTR segments       | UTR         | int         | False        |
|                                     | (referred to the selected ORF).                           |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| three_utr_num_complete              | This property returns the number of 3' UTR segments for   | UTR         | int         | False        |
|                                     | the selected ORF, considering only those which are        |             |             |              |
|                                     | complete exons.                                           |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| utr_fraction                        | This property calculates the length of the UTR of the     | UTR         | float       | True         |
|                                     | selected ORF vs. the cDNA length.                         |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| utr_length                          | Returns the sum of the 5'+3' UTR lengths                  | UTR         | int         | False        |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| utr_num                             | Returns the number of UTR segments (referred to the       | UTR         | int         | False        |
|                                     | selected ORF).                                            |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| utr_num_complete                    | Returns the number of UTR segments which are complete     | UTR         | int         | False        |
|                                     | exons (referred to the selected ORF).                     |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+
| verified_introns_num                | This metric returns the number of introns of the          | External    | int         | False        |
|                                     | transcript which are validated by external data.          |             |             |              |
+-------------------------------------+-----------------------------------------------------------+-------------+-------------+--------------+


.. _external-metrics:

External metrics
~~~~~~~~~~~~~~~~

Starting from version 1 beta 8, Mikado allows to load external metrics into the database, to be used for evaluating transcripts. Metrics of this kind **must** have a value comprised between 0 and 1.
The file can be provided either by specifying it in the :ref:`coonfiguration file <conf_anatomy>`, under "serialise/files/external_scores", or on the command line with the "--external-scores" parameters to mikado :ref:`serialise <serialise>`.
The external scores file should have the following format:

+--------------------------+------------------+------------------+-------+-----------------+
| TID                      | Metric_one       | Metric_two       | ...   |  Metric_N       |
+==========================+==================+==================+=======+=================+
| Transcript_one           | value            | value            |       | value           |
+--------------------------+------------------+------------------+-------+-----------------+
| Transcript_two           | value            | value            |       | value           |
+--------------------------+------------------+------------------+-------+-----------------+
| ...                      | ...              | ...              |       | ...             |
+--------------------------+------------------+------------------+-------+-----------------+
| Transcript_N             | value            | value            |       | value           |
+--------------------------+------------------+------------------+-------+-----------------+


Please note the following:

* the header is mandatory.
* the metric names at the head of the table should **not** contain any space or spcecial characters, apart from the underscore (_)
* the header provides the name of the metric as will be seen by Mikado. As such, it is advised to choose sensible and informative names (e.g. "fraction_covered") rather than uninformative ones (e.g. the "metric_one" from above)
* Column names **must be unique**.
* The transcript names present in the first column **must** be present in the FASTA file.
* The table should be tab-separated.
* Values can be of any numerical or boolean type. However, only values that are determined **at serialisation** to be comprised within 0 and 1 (inclusive) can be used as raw values.

A proper way of generating and using external scores would, therefore, be the following:

* Run Mikado prepare on the input dataset.
* Run all necessary supplementary analyses (ORF calling and/or homology analysis through DIAMOND or BLAST).
* Run supplementary analyses to assess the transcripts, e.g. expression analysis. Normalise results so that they can be expressed with values between 0 and 1.

  * Please note that boolean results (e.g. presence or absence) can be expressed with 0 and 1 intead of "False" and "True". Customarily, in Python 0 stands for False and 1 for True, but you can choose to switch the order if you so desire.
* Aggregate all results in a text table, like the one above, tab separated.
* Call mikado serialise, specifying the location of this table either through the configuration file or on the command line invocation.

Given the open ended nature of the external scores, the Daijin pipeline currently does not offer any system to generate these scores. This might change in the future.

Adding external scores to the scoring file
------------------------------------------

Once the external metrics have been properly loaded, it is necessary to tell Mikado how to use them. This requires :ref:`modifying the scoring file itself <configure-scoring-tutorial>`. The header that we used in the table above does provide the names of the metrics as they will be seen by Mikado.

Let us say that we have performed an expression analysis on our transcripts, and we have created and loaded the following three metrics:

* "fraction_covered", ie the percentage of the transcript covered by at least X reads (where X is decided by the experimenter)
* "samples_expressed", ie the percentage of samples where the transcript was expressed over a certain threshold (e.g. 1 TPM)
* "has_coverage_gaps", ie a boolean metrics that indicates whether there are windows *within* the transcript that lowly or not at all covered (e.g. a 100bp stretch with no coverage between two highly covered regions, indicating a possilble intron retention or chimera). For this example, a value of "0" indicates that there no coverage gaps (ie. it is *False* that there are gaps), "1" otherwise (it is *True* that there are coverage gaps).

We can now use these metrics as normal, by invoking them as "external." followed by the name of the metrics: e.g., "external.fraction_covered".
So for example, if we wished to prioritise transcripts that are expressed in the highest number of samples and are completely covered by RNASeq data upon reads realignment, under "scoring", we can add the following:

.. code-block:: yaml

    scoring:
        # [ ... other metrics ... ]
        - external.samples_expressed: {rescaling: max}
        - external.fraction_covered: {rescaling: max}

And if we wanted to consider any primary transcript with coverage gaps as a potential fragment, under the "fragmentary" section we could do so:

.. code-block:: yaml

    not_fragmentary:
      expression:
        # other metrics ..
        - and (external.has_coverage_gaps)
        # Finished expression
      parameters:
        # other metrics ...
	external.has_coverage_gaps: {operator: eq, value: 0}  # Please note, again, that here "0" means "no coverage gaps detected".
	# other metrics ...

As external metrics allow Mikado to accept any arbitrary metric for each transcript, they allow the program to assess transcripts in any way the experimenter desires. However, currently we do not provide any way of automating the process.

.. note:: also for external metrics, it is necessary to add a suffix to them if they are invoked more than once in an expression (see the :ref:`tutorial <scoring-tutorial-first-reqs>`). An invocation of e.g. "external.samples_expressed.mono" and "external.samples_expressed.multi", to distinguish between monoexonic and multiexonic transcripts, would be perfectly valid and actually *required* by Mikado. Notice the double use of the dot (".") as separator. Its usage as such is the reason that it cannot be present in the name of the metric itself (so, for example, "has.coverage.gaps" would be an invalid metric name).

.. _attributes-metrics:

Attributes metrics
~~~~~~~~~~~~~~~~~~
Starting from version 2, Mikado allows the usage of metrics defined in the attributes of the input files, these metrics
behave as the rest of the metrics but they are gathered at runtime from the input datasets. It is important to note that
these metrics must be equivalent in all the inputs and are by default initialised to "0" when a transcript does not have
an attribute defining the metric. The default initialisation value can be overridden in the scoring file.

Attribute metrics along with the required **rescaling** parameter, can define a *rtype* parameter as one of (float, int
or bool) which will be used to cast the value of the attribute internally, and a *percentage* boolean which indicates
that the values are in the 0-100 range and enables a transformation to the 0-1 range so that these can be used as 'raw'
scores (see the :ref:`scoring algorithm section <_scoring_algorithm>`).

An example for the usage of these metrics could be::

        Chr5	Cufflinks	transcript	26581218	26583874	1000	-	.	gene_id "cufflinks_star_at.23551";transcript_id "cufflinks_star_at.23551.1";exon_number "1";FPKM "0.4343609420";conf_hi "0.577851";frac "0.751684";cov "11.982854";conf_lo "0.293994";percentage_score "42.42"
        Chr5	Cufflinks	exon	26581218	26581528	.	-	.	gene_id "cufflinks_star_at.23551";transcript_id "cufflinks_star_at.23551.1";
        Chr5	Cufflinks	exon	26583335	26583874	.	-	.	gene_id "cufflinks_star_at.23551";transcript_id "cufflinks_star_at.23551.1";


If the scoring file defines:

.. code-block:: yaml

    scoring:
        # [ ... other metrics ... ]
        - attributes.FPKM: {rescaling: max}
        - attributes.frac: {rescaling: max, use_raw: true}
        - attributes.percentage_score: {rescaling: max, use_raw: true, percentage: true}

The same scoring rules defined previously will apply to metrics obtained from the transcript's attributes.