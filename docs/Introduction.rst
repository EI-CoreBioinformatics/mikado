.. _Introduction:

Introduction
============

Introduction test


.. _Metrics:

Metrics
-------

These are all the metrics available to quantify transcripts. The documentation for this section has been generated using the :ref:`metrics utility <metrics-command>`.

- *tid*:	ID of the transcript - cannot be an undefined value. Alias of id. :rtype str
- *parent*:	Name of the parent feature of the transcript.
- *score*:	Numerical value which summarizes the reliability of the transcript.
- *best_bits*:	Metric that returns the best BitS associated with the transcript.
- *blast_score*:	 Interchangeable alias for testing different blast-related scores. Current: best bit score. :return:
- *canonical_intron_proportion*:	 This metric returns the proportion of canonical introns of the transcript on its total number of introns. :return:
- *cdna_length*:	This property returns the length of the transcript.
- *cds_not_maximal*:	This property returns the length of the CDS excluded from the selected ORF.
- *cds_not_maximal_fraction*:	This property returns the fraction of bases not in the selected ORF compared to the total number of CDS bases in the cDNA.
- *combined_cds_fraction*:	This property return the percentage of the CDS part of the transcript vs. the cDNA length
- *combined_cds_intron_fraction*:	This property returns the fraction of CDS introns of the transcript vs. the total number of CDS introns in the Locus. If the transcript is by itself, it returns 1.
- *combined_cds_length*:	This property return the length of the CDS part of the transcript.
- *combined_cds_num*:	This property returns the number of non-overlapping CDS segments in the transcript.
- *combined_cds_num_fraction*:	This property returns the fraction of non-overlapping CDS segments in the transcript vs. the total number of exons
- *combined_utr_fraction*:	This property returns the fraction of the cDNA which is not coding according to any ORF. Complement of combined_cds_fraction
- *combined_utr_length*:	This property return the length of the UTR part of the transcript.
- *end_distance_from_junction*:	This metric returns the cDNA distance between the stop codon and the last junction in the transcript. In many eukaryotes, this distance cannot exceed 50-55 bps otherwise the transcript becomes a target of NMD. If the transcript is not coding or there is no junction downstream of the stop codon, the metric returns 0. This metric considers the combined CDS end.
- *end_distance_from_tes*:	This property returns the distance of the end of the combined CDS from the transcript end site. If no CDS is defined, it defaults to 0.
- *exon_fraction*:	This property returns the fraction of exons of the transcript which are contained in the sublocus. If the transcript is by itself, it returns 1. Set from outside.
- *exon_num*:	This property returns the number of exons of the transcript.
- *five_utr_length*:	Returns the length of the 5' UTR of the selected ORF.
- *five_utr_num*:	This property returns the number of 5' UTR segments for the selected ORF.
- *five_utr_num_complete*:	This property returns the number of 5' UTR segments for the selected ORF, considering only those which are complete exons.
- *has_start_codon*:	Boolean. True if the selected ORF has a start codon. :rtype: bool
- *has_stop_codon*:	Boolean. True if the selected ORF has a stop codon. :rtype bool
- *highest_cds_exon_number*:	This property returns the maximum number of CDS segments among the ORFs; this number can refer to an ORF *DIFFERENT* from the maximal ORF.
- *highest_cds_exons_num*:	Returns the number of CDS segments in the selected ORF (irrespective of the number of exons involved)
- *intron_fraction*:	This property returns the fraction of introns of the transcript vs. the total number of introns in the Locus. If the transcript is by itself, it returns 1. Set from outside.
- *is_complete*:	Boolean. True if the selected ORF has both start and end.
- *max_intron_length*:	This property returns the greatest intron length for the transcript.
- *min_intron_length*:
- *non_verified_introns_num*:	 This metric returns the number of introns of the transcript which are not validated by external data. :rtype : int
- *num_introns_greater_than_max*:	 This metric returns the number of introns greater than the maximum acceptable intron size indicated in the constructor. :rtype : int
- *num_introns_smaller_than_min*:	 This metric returns the number of introns smaller than the mininum acceptable intron size indicated in the constructor. :rtype : int
- *number_internal_orfs*:	This property returns the number of ORFs inside a transcript.
- *proportion_verified_introns*:	This metric returns, as a fraction, how many of the transcript introns are validated by external data.
- *proportion_verified_introns_inlocus*:	This metric returns, as a fraction, how many of the verified introns inside the Locus are contained inside the transcript.
- *retained_fraction*:	This property returns the fraction of the cDNA which is contained in retained introns.
- *retained_intron_num*:	This property records the number of introns in the transcripts which are marked as being retained. See the corresponding method in the sublocus class.
- *selected_cds_exons_fraction*:	Returns the fraction of CDS segments in the selected ORF (irrespective of the number of exons involved)
- *selected_cds_fraction*:	This property calculates the fraction of the selected CDS vs. the cDNA length.
- *selected_cds_intron_fraction*:	This property returns the fraction of CDS introns of the selected ORF of the transcript vs. the total number of CDS introns in the Locus (considering only the selected ORF). If the transcript is by itself, it should return 1.
- *selected_cds_length*:	This property calculates the length of the CDS selected as best inside the cDNA.
- *selected_cds_num*:	This property calculates the number of CDS exons for the selected ORF
- *selected_cds_number_fraction*:	This property returns the proportion of best possible CDS segments vs. the number of exons. See selected_cds_number.
- *selected_end_distance_from_junction*:	This metric returns the distance between the stop codon and the nearest downstream junction. In many eukaryotes, this distance cannot exceed 50-55 bps, otherwise the transcript becomes a target of NMD. If the transcript is not coding or there is no junction downstream of the stop codon, the metric returns 0.
- *selected_end_distance_from_tes*:	This property returns the distance of the end of the best CDS from the transcript end site. If no CDS is defined, it defaults to 0.
- *selected_start_distance_from_tss*:	This property returns the distance of the start of the best CDS from the transcript start site. If no CDS is defined, it defaults to 0.
- *snowy_blast_score*:	 Metric that indicates how good a hit is compared to the competition, in terms of BLAST similarities. As in SnowyOwl, the score for each hit is calculated by taking the percentage of positive matches and dividing it by (2 * len(self.blast_hits)). IMPORTANT: when splitting transcripts by ORF, a blast hit is added to the new transcript only if it is contained within the new transcript. This WILL screw up a bit the homology score. :return
- *source_score*:	This metric returns a score that is assigned to the transcript in virtue of its origin.
- *start_distance_from_tss*:	This property returns the distance of the start of the combined CDS from the transcript start site. If no CDS is defined, it defaults to 0.
- *three_utr_length*:	Returns the length of the 5' UTR of the selected ORF.
- *three_utr_num*:	This property returns the number of 3' UTR segments (referred to the selected ORF).
- *three_utr_num_complete*:	This property returns the number of 3' UTR segments for the selected ORF, considering only those which are complete exons.
- *utr_fraction*:	This property calculates the length of the UTR of the selected ORF vs. the cDNA length.
- *utr_length*:	Returns the sum of the 5'+3' UTR lengths
- *utr_num*:	Returns the number of UTR segments (referred to the selected ORF).
- *utr_num_complete*:	Returns the number of UTR segments which are complete exons (referred to the selected ORF).
- *verified_introns_num*:	 This metric returns the number of introns of the transcript which are validated by external data. :rtype : int