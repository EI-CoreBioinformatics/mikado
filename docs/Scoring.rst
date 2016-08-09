.. _scoring_files:

Transcript measurements and scoring
===================================

In order to determine the best transcript for each locus, Mikado measures each available candidate according to various different :ref:`metrics <Metrics>` and assigning a specific score for each of those. Similarly to `RAMPART <https://github.com/TGAC/RAMPART>`_ [Rampart]_, Mikado will assign a score to each transcript for each metric by assessing it relatively to the other transcripts in the locus. The particular feature rescaling equation used for a given metric depends on the type of feature it represents:

* metrics where higher values represent better transcript assemblies ("maximum").
* metrics where lower values represent better transcript assemblies ("minimum")
* metrics where values closer to a defined value represent better assemblies ("target")

To allow for this tripartite scoring system with disparate input values, we have to employ rescaling equations so that each metric for each transcript will be assigned a score between 0 and 1. Optionally, each metric might be assigned a greater weight so that its maximum possible value will be greater or smaller than 1. Formally, let metric :math:`m` be one of the available metrics :math:`M`, :math:`t` a transcript in locus :math:`L`, :math:`w_{m}` the weight assigned to metric :math:`m`, and :math:`r_{mt}` the raw value of metric :math:`m` for :math:`t`. Then, the score to metric :math:`m` for transcript :math:`t`, :math:`s_{mt}`, will be derived using one of the following three different rescaling equations:

* If higher values are best:
    :math:`s_{mt} = w_{m} * (\frac{r_{mt} - min(r_m)}{max(r_m)-min(r_m)})`
* If lower values are best:
    :math:`s_{mt} = w_{m} * (1 - \frac{r_{mt} - min(r_m)}{max(r_m)-min(r_m)})`
* If values closer to a target :math:`v_{m}` are best:
    :math:`s_{mt} = w_{m} * (1 - \frac{|r_{mt} - v_{m}|}{max(|r_{m} - v_{m}|)})`

Finally, the scores for each metric will be summed up to produce a final score for the transcript:
    :math:`s_{t} = \sum_{m \forall m \in M} s_{mt}`.

Not all the available metrics will be necessarily used for scoring; the choice of which to employ and how to score and weight each of them is left to the experimenter, although Mikado provides some pre-configured scoring files.

Scoring files
~~~~~~~~~~~~~

Mikado employs user-defined configuration files to define the desirable features in genes. These files are composed of two sections:

  #. a *requirements* section, specifying the minimum requirements that a transcript must satisfy to be considered as valid. **Any transcript failing these requirements will be scored at 0 and purged.**
  #. a *scoring* section, specifying which features Mikado should look for in transcripts, and how each of them will be weighted.

.. hint:: Although at the moment Mikado does not offer any method to establish machine-learning based scoring configurations, it is a topic we plan to investigate in the future. Mikado already supports `Random Forest Regressors as scorers through Scikit-learn <http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html>`_, but we have yet to devise a proper way to create such regressors.

The requirements section
------------------------

This section specifies the minimum requirements for a transcript. As transcripts failing to pass these checks will be discarded outright, **it is strongly advised to use lenient parameters in this section**. Being too stringent might end up removing valid models and potentially missing valid loci outright. Typically, transcripts filtered at this step should be obvious fragments, eg monoexonic transcripts produced by RNA-Seq with a total length lower than the *library* fragment length.
This section is composed by two parts:

* *parameters*: a list of the metrics to be considered. Each metric can be considered multiple times, by suffixing it with a ".<id>" construct (eg cdna_length.*mono* vs. cdna_length.*multi* to distinguish two uses of the cdna_length metric - once for monoexonic and once for multiexonic transcripts). Any parameter which is not a :ref:`valid metric name <Metrics>`, after removal of the suffix, **will cause an error**. Parameters have to specify the following:
    * a *value* that the metric has to be compared against
    * an *operator* that specifies the target operation. See :ref:`the operators section <operators>`.

* *expression*: a string array that will be compiled into a valid boolean expression. All the metrics present in the expression string **must be present in the parameters section**. If an unrecognized metric is present, Mikado will crash immediately, complaining that the scoring file is invalid. Apart from brackets, Mikado accepts only the following boolean operators to chain the metrics:
    * *and*
    * *or*
    * *not*
    * *xor*

.. hint:: if no *expression* is specified, Mikado will construct one by chaining all the provided parameters with and *and* operator. Most of the time, this would result in an unexpected behaviour - ie Mikado assigning a score of 0 to most transcripts. It is **strongly advised** to explicitly provide a valid expression.

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

In order:
    * In the parameters section, we ask for the following:
        * *exon_num.mono*: monoexonic transcripts must have one exon ("eq")
        * *exon_num.multi*: multiexonic transcripts must have more than one exon ("gt")
        * *cdna_length.mono*: monoexonic transcripts must have a length greater than 50 bps (the ".mono" suffix is arbitrary, as long as it is unique for all calls of *cdna_length*)
        * *cdna_length.multi*: multiexonic transcripts must have a length greater than or equal to 100 bps (the ".multi" suffix is arbitrary, as long as it is unique for all calls of *cdna_length*)
        * *max_intron_length*: multiexonic transcripts should not have any intron longer than 200,000 bps.
        * *min_intron_length*: multiexonic transcripts should not have any intron smaller than 5 bps.

    * the *expression* field will be compiled into the following expression::

        (exon_num > 1 and cdna_length >= 100 and max_intron_length <= 200000 and min_intron_length >= 5) or (exon_num == 1 and cdna_length > 50)


Any transcript for which the expression evaluates to :math:`False` will be assigned a score of 0 outright and therefore discarded.

The scoring section
-------------------


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

.. _operators:

Operators
---------

Mikado allows the following operators to express a relationship inside the scoring files:

* *eq*: equal to (:math:`=`). Valid for comparisons with numbers, boolean values, and strings.
* *ne*: different from (:math:`\neq`). Valid for comparisons with numbers, boolean values, and strings.
* *lt*: less than (:math:`<`). Valid for comparisons with numbers.
* *gt*: greater than (:math:`>`). Valid for comparisons with numbers.
* *le*: less or equal than (:math:`\le`). Valid for comparisons with numbers.
* *ge*: greater or equal than (:math:`\ge`). Valid for comparisons with numbers.
* *in*: member of (:math:`\in`). Valid for comparisons with arrays or sets.
* *not in*: not member of (:math:`\notin`). Valid for comparisons with arrays or sets.

Mikado will fail if an operator not present on this list is specified, or if the operator is assigned to compare against the wrong data type (eg. *eq* with an array).

.. _Metrics:

Metrics
~~~~~~~

These are all the metrics available to quantify transcripts. The documentation for this section has been generated using the :ref:`metrics utility <metrics-command>`.

Metrics belong to one of the following categories:

* **Descriptive**: these metrics merely provide a description of the transcript (eg its ID) and are not used for scoring.
* **cDNA**: these metrics refer to basic features of any transcript such as its number of exons, its cDNA length, etc.
* **Intron**: these metrics refer to features related to the number of introns and their lengths.
* **CDS**: these metrics refer to features related to the CDS assigned to the transcript.
* **UTR**: these metrics refer to features related to the UTR of the transcript. In the case in which a transcript has been assigned multiple ORFs, unless otherwise stated the UTR metrics will be derived only considering the *selected* ORF, not the combination of all of them.
* **Locus**: these metrics refer to features of the transcript in relationship to all other transcripts in its locus, eg how many of the introns present in the locus are present in the transcript. These metrics are calculated by Mikado during the picking phase, and as such their value can vary during the different stages as the transcripts are shifted to different groups.
* **External**: these metrics are derived from accessory data that is recovered for the transcript during the run time. Examples include data regarding the number of introns confirmed by external programs such as PortCullis, or the BLAST score of the best hits.


+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| Metric name                                    | Description                                               | Data type    | Category        |
|                                                |                                                           |              |                 |
+================================================+===========================================================+==============+=================+
| *tid*                                          | Name of the transcript. Not used for scoring.             | String       | **Descriptive** |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *parent*                                       | Name of the transcript parent. Not used for scoring.      | String       | **Descriptive** |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *score*                                        | Final score of the transcript.                            | Float        | **Descriptive** |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *best_bits*                                    | Best Bit Score associated with the transcript.            | Float        | **External**    |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *blast_score*                                  | Alias for either *best_bits* or *snowy_blast_score*. Set  | Float        | **External**    |
|                                                | currently to *best_bits*                                  |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *canonical_intron_proportion*                  | This metric returns the proportion of canonical introns of| Float        | **Intron**      |
|                                                | the transcript on its total number of introns.            | Float        |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *cdna_length*                                  | This property returns the length of the transcript.       | Int          | **cDNA**        |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *cds_not_maximal*                              | This property returns the length of the CDS excluding     | Int          | **CDS**         |
|                                                | that contained in the selected ORF. If the transcript only|              |                 |
|                                                | has one ORF, this metric returns a value of 0.            |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *cds_not_maximal_fraction*                     | This property returns the fraction of bases not in the    | Float        | **CDS**         |
|                                                | selected ORF compared to the total number of CDS bases    |              |                 |
|                                                | in the cDNA.                                              |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *combined_cds_fraction*                        | This property return the percentage of the CDS part of the| Float        | **CDS**         |
|                                                | transcript vs. the cDNA length.                           |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *combined_cds_intron_fraction*                 | This property returns the fraction of CDS introns of the  | Float        | **Locus**       |
|                                                | transcript vs. the total number of CDS introns in the     |              |                 |
|                                                | Locus. If the transcript is by itself, it returns 1.      |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *combined_cds_length*                          | This property returns the fraction of CDS introns of the  | Float        | **CDS**         |
|                                                | transcript, across all its ORFs.                          |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *combined_cds_num*                             | This property returns the number of non-overlapping CDS   | Int          | **CDS**         |
|                                                | segments in the transcript.                               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *combined_cds_num_fraction*                    | This property returns the fraction of non-overlapping CDS | Float        | **CDS**         |
|                                                | segments in the transcript vs. the total number of exons. |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *combined_utr_fraction*                        | This property returns the fraction of the cDNA which is   | Float        | **UTR**         |
|                                                | not coding according to any ORF. Complement of            |              |                 |
|                                                | *combined_cds_fraction*                                   |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *combined_utr_length*                          | This property return the length of the UTR part of the    | Int          | **UTR**         |
|                                                | transcript.                                               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *end_distance_from_junction*                   | This metric returns the cDNA distance between the stop    | Int          | **CDS**         |
|                                                | and the last junction of the transcript. In many          |              |                 |
|                                                | eukaryotes, this distance cannot exceed 50-55 bps,        |              |                 |
|                                                | otherwise the transcript becomes a target for NMD. If the |              |                 |
|                                                | transcript is not coding or there is no junction          |              |                 |
|                                                | downstream of the stop codon, the metric returns 0.       |              |                 |
|                                                | This metric considers the combined CDS end.               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *end_distance_from_tes*                        | This property returns the distance of the end of the      | Int          | **CDS**         |
|                                                | combined CDS from the transcript end site. If no CDS is   |              |                 |
|                                                | defined, it defaults to 0.                                |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *exon_fraction*                                | This property returns the fraction of exons of the        | Float        | **Locus**       |
|                                                | transcript which are contained in the sublocus. If the    |              |                 |
|                                                | transcript is by itself, it returns 1.                    |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *exon_num*                                     | This property returns the number of exons of the          | Int          | **cDNA**        |
|                                                | transcript.                                               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *five_utr_length*                              | Returns the length of the 5' UTR of the *selected* ORF.   | Int          | **UTR**         |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *five_utr_num*                                 | This property returns the number of 5' UTR segments for   | Int          | **UTR**         |
|                                                | the selected ORF.                                         |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *five_utr_num_complete*                        | This property returns the number of 5' UTR segments for   | Int          | **UTR**         |
|                                                | the selected ORF, considering only those which are        |              |                 |
|                                                | complete exons.                                           |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *has_start_codon*                              | True if the selected ORF has a start codon, False         | Bool         | **CDS**         |
|                                                | otherwise                                                 |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *has_stop_codon*                               | True if the selected ORF has a stop codon, False otherwise| Bool         | **CDS**         |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *highest_cds_exon_number*                      | This property returns the maximum number of CDS segments  | Int          | **CDS**         |
|                                                | among the ORFs; this number can refer to an ORF           |              |                 |
|                                                | *DIFFERENT* from the maximal ORF.                         |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *highest_cds_exons_num*                        | Returns the number of CDS segments in the selected ORF    | Int          | **CDS**         |
|                                                | (irrespective of the number of exons involved)            |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *intron_fraction*                              | This property returns the fraction of introns of the      | Float        | **Locus**       |
|                                                | transcript vs. the total number of introns in the Locus.  |              |                 |
|                                                | If the transcript is by itself, it returns 1.             |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *is_complete*                                  | Boolean. True if the selected ORF has both start and end. | Bool         | **CDS**         |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *max_intron_length*                            | This property returns the greatest intron length for the  | Int          | **Intron**      |
|                                                | transcript.                                               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *min_intron_length*                            | This property returns the smallest intron length for the  | Int          | **Intron**      |
|                                                | transcript.                                               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *non_verified_introns_num*                     | This metric returns the number of introns of the          | Int          | **External**    |
|                                                | transcript which are not validated by external data.      |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *num_introns_greater_than_max*                 | This metric returns the number of introns greater than the| Int          | **Intron**      |
|                                                | maximum acceptable intron size indicated in the           |              |                 |
|                                                | constructor.                                              |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *num_introns_smaller_than_min*                 | This metric returns the number of introns smaller than the| Int          | **Intron**      |
|                                                | mininum acceptable intron size indicated in the           |              |                 |
|                                                | constructor.                                              |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *number_internal_orfs*                         | This property returns the number of ORFs inside a         | Int          | **CDS**         |
|                                                | transcript.                                               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *proportion_verified_introns*                  | This metric returns, as a fraction, how many of the       | Float        | **External**    |
|                                                | transcript introns are validated by external data.        |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *proportion_verified_introns_inlocus*          | This metric returns, as a fraction, how many of the       | Float        | **Locus**       |
|                                                | verified introns inside the Locus are contained inside the|              |                 |
|                                                | transcript.                                               |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *retained_fraction*                            | This property returns the fraction of the cDNA which is   | Float        | **Locus**       |
|                                                | contained in retained introns.                            |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *retained_intron_num*                          | This property records the number of introns in the        | Int          | **Locus**       |
|                                                | transcripts which are marked as being retained.           |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_cds_exons_fraction*                  | Returns the fraction of CDS segments in the selected ORF  | Float        | **CDS**         |
|                                                | (irrespective of the number of exons involved)            |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_cds_fraction*                        | This property calculates the fraction of the selected CDS | Float        | **CDS**         |
|                                                | vs. the cDNA length.                                      |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_cds_intron_fraction*                 | This property returns the fraction of CDS introns of the  | Float        | **CDS**         |
|                                                | selected ORF of the transcript vs. the total number of    |              |                 |
|                                                | CDS introns in the Locus (considering only the selected   |              |                 |
|                                                | ORF). If the transcript is by itself, it should return 1. |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_cds_length*                          | This property calculates the length of the CDS selected   | Int          | **CDS**         |
|                                                | as best inside the cDNA.                                  |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_cds_num*                             | This property calculates the number of CDS exons for the  | Int          | **CDS**         |
|                                                | selected ORF.                                             |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_cds_number_fraction*                 | This property returns the proportion of best possible CDS | Float        | **CDS**         |
|                                                | segments vs. the number of exons. See selected_cds_number.|              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_end_distance_from_junction*          | This metric returns the distance between the stop codon   | Int          | **CDS**         |
|                                                | and the last junction of the transcript. In many          |              |                 |
|                                                | eukaryotes, this distance cannot exceed 50-55 bps,        |              |                 |
|                                                | otherwise the transcript becomes a target for NMD. If the |              |                 |
|                                                | transcript is not coding or there is no junction          |              |                 |
|                                                | downstream of the stop codon, the metric returns 0.       |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_end_distance_from_tes*               | This property returns the distance of the end of the best | Int          | **CDS**         |
|                                                | CDS from the transcript end site. If no CDS is defined,   |              |                 |
|                                                | it defaults to 0.                                         |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *selected_start_distance_from_tss*             | This property returns the distance of the start of the    | Int          | **CDS**         |
|                                                | best CDS from the transcript start site. If no CDS is     |              |                 |
|                                                | defined, it defaults to 0.                                |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *snowy_blast_score*                            | Metric that indicates how good a hit is compared to the   | Float        | **External**    |
|                                                | competition, in terms of BLAST similarities. As in        |              |                 |
|                                                | SnowyOwl [SnowyOwl]_, the score for each hit is calculated|              |                 |
|                                                | by taking the percentage of positive matches and dividing |              |                 |
|                                                | it by (2 * len(self.blast_hits)). IMPORTANT: when         |              |                 |
|                                                | splitting transcripts by ORF, a blast hit is added to the |              |                 |
|                                                | new transcript only if it is contained within it. This    |              |                 |
|                                                | will influnce directly this metric.                       |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *source_score*                                 | This metric returns a score that is assigned to the       | Float        | **External**    |
|                                                | transcript solely in virtue of its origin.                |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *start_distance_from_tss*                      | This property returns the distance of the start of the    | Int          | **CDS**         |
|                                                | combined CDS from the transcript start site.              |              |                 |
|                                                | If no CDS is defined, it defaults to 0.                   |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *three_utr_length*                             | Returns the length of the 5' UTR of the selected ORF.     | Int          | **UTR**         |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *three_utr_num*                                | This property returns the number of 3' UTR segments       | Int          | **UTR**         |
|                                                | (referred to the selected ORF).                           |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *three_utr_num_complete*                       | This property returns the number of 3' UTR segments for   | Int          | **UTR**         |
|                                                | the selected ORF, considering only those which are        |              |                 |
|                                                | complete exons.                                           |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *utr_fraction*                                 | This property calculates the length of the UTR of the     | Float        | **UTR**         |
|                                                | selected ORF vs. the cDNA length.                         |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *utr_length*                                   | Returns the sum of the 5'+3' UTR lengths.                 | Int          | **UTR**         |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *utr_num*                                      | Returns the number of UTR segments.                       | Int          | **UTR**         |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *utr_num_complete*                             | Returns the number of UTR segments which are complete     | Int          | **UTR**         |
|                                                | exons.                                                    |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+
| *verified_introns_num*                         | This metric returns the number of introns of the          | Int          | **External**    |
|                                                | transcript which are validated by external data.          |              |                 |
+------------------------------------------------+-----------------------------------------------------------+--------------+-----------------+

Technical details
~~~~~~~~~~~~~~~~~

Metrics are an extension of the ``property`` construct in Python3. Compared to normal properties, they are distinguished only by the optional ``category`` attribute, whose function is only descriptive. The main reason to subclass ``property`` is to allow Mikado to be self-aware of which properties will be used for scoring transcripts, and which will not. So, for example, in the following snippet from the :ref:`Transcript class definition <transcript-class>`:

.. code-block:: python

    @property
    def combined_cds(self):
        """This is a list which contains all the non-overlapping CDS
        segments inside the cDNA. The list comprises the segments
        as duples (start,end)."""
        return self.__combined_cds

    @combined_cds.setter
    def combined_cds(self, combined):
        """
        Setter for combined_cds. It performs some basic checks,
        e.g. that all the members of the list are integer duplexes.

        :param combined: list
        :type combined: list[(int,int)]
        """

        if ((not isinstance(combined, list)) or
                any(self.__wrong_combined_entry(comb) for comb in combined)):
            raise TypeError("Invalid value for combined CDS: {0}".format(combined))

    @Metric
    def combined_cds_length(self):
        """This property return the length of the CDS part of the transcript."""
        c_length = sum([c[1] - c[0] + 1 for c in self.combined_cds])
        if len(self.combined_cds) > 0:
            assert c_length > 0
        return c_length

    combined_cds_length.category = "CDS"

    @Metric
    def combined_cds_num(self):
        """This property returns the number of non-overlapping CDS segments
        in the transcript."""
        return len(self.combined_cds)

    combined_cds_num.category = "CDS"

    @Metric
    def has_start_codon(self):
        """Boolean. True if the selected ORF has a start codon.
        :rtype: bool"""
        return self.__has_start_codon

    @has_start_codon.setter
    def has_start_codon(self, value):
        """Setter. Checks that the argument is boolean.
        :param value: boolean flag
        :type value: bool
        """

        if value not in (None, False, True):
            raise TypeError(
                "Invalid value for has_start_codon: {0}".format(type(value)))
        self.__has_start_codon = value

    has_start_codon.category = "CDS"

Mikado will recognize that "derived_children" is a normal property, while "combined_cds_length", "combined_cds_num" and "has_start_codon" are Metrics (and as such, we assign them a "category" - by default, that attribute will be ``None``.). Please note that Metrics behave and are coded like normal properties in any other regard - including docstrings and setters/deleters.

The requirements expression is evaluated using ``eval``.

.. warning:: While we took pains to ensure that the expression is properly sanitised and inspected **before** ``eval``, Mikado might prove itself to be permeable to clever code injection attacks. Do **not** execute Mikado with super user privileges if you do not want to risk from such attacks, and always inspect third-party YAML scoring files before execution!