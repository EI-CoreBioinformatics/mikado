#Locus brewer: a pipeline to determine and select the best RNA-Seq prediction

##Overview

The locus brewer is a lightweight Python3 pipeline whose purpose is to facilitate the identification
of expressed loci from RNA-Seq data * and to select the best model in each locus.

The logic of the pipeline is as follows:

1. In a first step, the annotation (provided in GTF/GFF3 format) is parsed to locate *superloci* of overlapping features on the **same strand**.
2. The superloci are divided into different *subloci*, each of which is defined as follows:

    * For multiexonic transcripts, to belong to the same sublocus they must share at least a splicing junction (i.e. an intron)
    * For monoexonic transcripts, they must overlap for at least one base pair
    * All subloci must contain either only multiexonic or only monoexonic transcripts
3. In each sublocus, the pipeline selects the best transcript according to a user-defined prioritization scheme.
4. The resulting *monosubloci* are merged together, if applicable, into *monosubloci_holders*
5. The best non-overlapping transcripts are selected, in order to define the *loci* contained inside the superlocus.

    * At this stage, monoexonic and multiexonic transcript are checked for overlaps
    * Moreover, two multiexonic transcripts are considered to belong to the same locus if they share a splice *site* (not junction)

The resulting, cleaned file does contain therefore only one transcript per locus.
The criteria used to select the "*best*" transcript are left to the user's discretion.

##Implementation
The pipeline is implemented in Python3, and does not require any non-standard library.
It has been tested on UNIX and Windows systems, using Python 3.4.3.
Beware that Cygwin at the time of this writing does not provide a version of Python3 more recent than 3.2.3,
which might generate problems. 

###Implementation details

In a first step, transcripts are recovered from the input files and inserted into a *superlocus* object, irrespective
of strand.
Superloci objects are then split into superloci stranded objects, according to their strand (method "split_by_strand",
which returns an iterator).
In each stranded superlocus, transcripts are grouped into subloci according to the logic outlined in the "Overview" paragraph,
using the "define_subloci" method.
At this point, the superlocus object uses the "define_monosubloci" method to choose the best transcript for each sublocus.
In each sublocus, transcripts are analyzed to give them a score according to the specifications in the JSON configuration
(see "Scoring transcripts" section), using the "calculate_scores" method.
The JSON configuration might also contain a "requirements" field. Any transcript not passing these requirements will
be given a score of 0.
The "best" transcripts are used to create "monosubloci" objects. All the monosubloci are gathered into
monosubloci_holder objects, and the process of scoring and selection begins again, using the "define_loci" method.
Surviving transcripts are used to define the loci which will be the final output.

All "loci" objects are defined as implementations of the "abstractlocus" class.
As such, each of them must implement the following methods:
* *is_intersecting*		This method is used to define when two transcripts are intersecting. It varies from locus to locus, so it must be specified by each of the objects.
* __init__
* __str__

The abstractlocus class provides many methods (both class and instance ones) that are used throughout the program. 

##Input files
The pipeline can use as input files any properly GTF or GFF3 file. The only requirements are:
* GTF files must have a "transcript/mRNA" line for each transcript present in the annotation
  * Standard GTF files, which have only CDS/exon/start/stop features, are not interpreted correctly
* All files must be sorted by chromosome, start and stop.
  * A utility called "sort_gtf.py", present in the util folder, can be used to sort GTFs.
  For GFF3s, our recommendation is to use [GenomeTools](http://genometools.org/), with e.g. the following command line:
  "gt gff3 -sort -retainids -addids no -tidy -force -o <output> <input>"


Additionally, the pipeline can also use the output of coding potential programs such as [TransDecoder](https://transdecoder.github.io/ "TransDecoder (Find Coding Regions Within Transcripts)").
In order to give this type of data to the pipeline, the coding regions should be encoded into a BED12
file, with the "thickStart/thickEnd" fields (7 and 8) indicating the start and stop of the CDS.
As in transdecoder, it is expected that the coding potential has been defined using oriented transcript sequences.
Any multiexonic CDS found on the negative transcript strand will be therefore ignored, as it happens in TD.

##Output files
The pipeline can produce output files at the following stages:
* Definition of subloci
* Definition of monosubloci
* Definition of loci

At each stage it is possible to obtain both a a GFF of the transcripts currently considered, and all the statistics
of the remaining transcripts.
GFFs can be printed using a simple "print" call, while the metrics have to be retrieved as dictionaries using the
"print_metrics" file. These dictionaries can be written out using the DictWriter method from the standard csv library.

If the "purge" option is selected, all transcripts which do not pass the minimum requirements will be excluded.
This option can factively lead to the loss of whole loci, if they contain only low-quality transcripts.

GTF output is not supported, as the output is more hierarchical than what is supportable in a GTF file.

##Scoring transcripts

###Available metrics
Metrics are defined into a text file, "metrics.txt". Each of them is defined as a property inside the "transcript" class.
The documentation for each of them can be generated with the utility "generate_metric_docs.py" inside the util folder.  

- *tid*:	ID of the transcript - cannot be an undefined value. Alias of id.
- *parent*:	Name of the parent feature of the transcript.
- *score*:	Numerical value which summarizes the reliability of the transcript.
- *cdna_length*:	This property returns the length of the transcript.
- *cds_not_maximal*:	This property returns the length of the CDS excluded from the selected ORF.
- *cds_not_maximal_fraction*:	This property returns the fraction of bases not in the selected ORF compared to the total number of CDS bases in the cDNA.
- *combined_cds_fraction*:	This property return the percentage of the CDS part of the transcript vs. the cDNA length
- *combined_cds_intron_fraction*:	This property returns the fraction of CDS introns of the transcript vs. the total number of CDS introns in the locus. If the transcript is by itself, it returns 1.
- *combined_cds_length*:	This property return the length of the CDS part of the transcript.
- *combined_cds_num*:	This property returns the number of non-overlapping CDS segments in the transcript.
- *combined_cds_num_fraction*:	This property returns the fraction of non-overlapping CDS segments in the transcript vs. the total number of exons
- *combined_utr_fraction*:	This property returns the fraction of the cDNA which is not coding according to any ORF. Complement of combined_cds_fraction
- *combined_utr_length*:	This property return the length of the UTR part of the transcript.
- *end_distance_from_tes*:	This property returns the distance of the end of the combined CDS from the transcript end site. If no CDS is defined, it defaults to 0.
- *exon_fraction*:	This property returns the fraction of exons of the transcript which are contained in the sublocus. If the transcript is by itself, it returns 1. Set from outside.
- *exon_num*:	This property returns the number of exons of the transcript.
- *five_utr_length*:	Returns the length of the 5' UTR of the selected ORF.
- *five_utr_num*:	This property returns the number of 5' UTR segments for the selected ORF.
- *five_utr_num_complete*:	This property returns the number of 5' UTR segments for the selected ORF, considering only those which are complete exons.
- *has_start_codon*:	Boolean. True if the selected ORF has a start codon.
- *has_stop_codon*:	Boolean. True if the selected ORF has a stop codon.
- *highest_cds_exon_number*:	This property returns the maximum number of CDS segments among the ORFs; this number can refer to an ORF *DIFFERENT* from the maximal ORF.
- *intron_fraction*:	This property returns the fraction of introns of the transcript vs. the total number of introns in the locus. If the transcript is by itself, it returns 1. Set from outside.
- *is_complete*:	Boolean. True if the selected ORF has both start and end.
- *number_internal_orfs*:	This property returns the number of ORFs inside a transcript.
- *retained_fraction*:	This property returns the fraction of the cDNA which is contained in retained introns.
- *retained_intron_num*:	This property records the number of introns in the transcripts which are marked as being retained. See the corresponding method in the sublocus class.
- *selected_cds_exons_fraction*:	Returns the fraction of CDS segments in the selected ORF (irrespective of the number of exons involved)
- *selected_cds_fraction*:	This property calculates the fraction of the selected CDS vs. the cDNA length.
- *selected_cds_intron_fraction*:	This property returns the fraction of CDS introns of the selected ORF of the transcript vs. the total number of CDS introns in the locus (considering only the selected ORF). If the transcript is by itself, it should return 1.
- *selected_cds_length*:	This property calculates the length of the CDS selected as best inside the cDNA.
- *selected_cds_num*:	This property calculates the number of CDS exons for the selected ORF
- *selected_end_distance_from_tes*:	This property returns the distance of the end of the best CDS from the transcript end site. If no CDS is defined, it defaults to 0.
- *selected_start_distance_from_tss*:	This property returns the distance of the start of the best CDS from the transcript start site. If no CDS is defined, it defaults to 0.
- *start_distance_from_tss*:	This property returns the distance of the start of the combined CDS from the transcript start site. If no CDS is defined, it defaults to 0.
- *three_utr_length*:	Returns the length of the 5' UTR of the selected ORF.
- *three_utr_num*:	This property returns the number of 3' UTR segments (referred to the selected ORF).
- *three_utr_num_complete*:	This property returns the number of 3' UTR segments for the selected ORF, considering only those which are complete exons.
- *utr_fraction*:	This property calculates the length of the UTR of the selected ORF vs. the cDNA length.
- *utr_length*:	Returns the sum of the 5'+3' UTR lengths
- *utr_num*:	Returns the number of UTR segments (referred to the selected ORF).
- *utr_num_complete*:	Returns the number of UTR segments which are complete exons (referred to the selected ORF).


###Defining new metrics
For the program to be able to use a novel metric, it must be implemented as a *metric* inside the
transcript class. The "metric" class is an alias of "property", and therefore metrics are coded as properties would:

```
@metric
def my_metric(self):
	return self.__my_metric

@my_metric.setter
def my_metric(self,value):
	self.__my_metric = value
``` 

If the metric is relative to other transcripts in the locus (e.g. retained introns),
it must still be implemented as above, but its value must be set by the "calculate_metrics" method in the
sublocus/monosublocus_holder class.

##Scoring configuration
At each selection stage, the best transcript for the span is selected according to user-defined criteria, supplied
through a JSON file.Each metric must be specified in the "parameters" head field of the JSON configuration.
For each parameter, it is possible to specify the following:

* "multiplier"		A number by which the metrics will be multiplied to get the final score.
* "rescaling"		This key controls the rescaling performed to calculate the score:
**"max"
**"min"
**"target"
    
Each parameter will have a score assigned which varies from 0 to 1. If the "target" rescaling is selected, it is *mandatory* to specify a "value" keyword. For details, see [RAMPART supporting material (section 2, page3)](http://bioinformatics.oxfordjournals.org/content/suppl/2015/01/29/btv056.DC1/supplementary.pdf)

Moreover, for each parameter it is possible to configure a "filter", i.e. boundaries after which the score for this parameter is set automatically to 0 (e.g. a 3'UTR total length over 2.5 kbps). Each "filter" subfield must contain the following:

* "operator": one of
** "eq": equal to
** "ne": different from
** "gt": greater than
** "ge": greater than or equal to
** "lt": less than
** "le": less than or equal to
** "in": value in array of valid values
** "not in": value not in array of invalid values

The comparisons are always made against the reference value.

  

