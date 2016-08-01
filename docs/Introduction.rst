.. _Introduction:

.. |python_badge| image:: https://img.shields.io/pypi/pyversions/snakemake.svg?style=flat-square
   :target: https://www.python.org/
.. |snake_badge| image:: https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)]
   :target: http://snakemake.bitbucket.org

Introduction
============

|python_badge| |snake_badge|


Introduction test


.. _Metrics:

Metrics
-------

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
|                                                | codon and the last junction in the transcript. In many    |              |                 |
|                                                | eukaryotes, this distance cannot exceed 50-55 bps         |              |                 |
|                                                | otherwise the transcript becomes a target of NMD. If the  |              |                 |
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
|                                                | and the nearest downstream junction. In many eukaryotes,  |              |                 |
|                                                | this distance cannot exceed 50-55 bps, otherwise the      |              |                 |
|                                                | transcript becomes a target of NMD. If the transcript is  |              |                 |
|                                                | not coding or there is no junction downstream of the stop |              |                 |
|                                                | codon, the metric returns 0.                              |              |                 |
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
