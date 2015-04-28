#Locus brewer: a pipeline to determine and select the best RNA-Seq prediction

##Overview

The locus brewer is a lightweight Python3 pipeline whose purpose is to facilitate the identification
of expressed loci from RNA-Seq data - and to select the best model in each locus.

The logic of the pipeline is as follows:

1. In a first step, the annotation (provided in GTF/GFF3 format) is parsed to locate *superloci* of overlapping features on the **same strand**.
2. The superloci are divided into different *subloci*, each of which is defined as follows:
   * For multiexonic transcripts, to belong to the same sublocus they must share at least a splicing junction (i.e. an intron)
   * For monoexonic transcripts, they must overlap for at least one base pair
   * All subloci must contain either only multiexonic or only monoexonic transcripts
3. In each sublocus, the pipeline select the best transcript (__WORK IN PROGRESS__)
4. The resulting *monosubloci* are merged together, if applicable
5. The best non-overlapping transcripts are selected, in order to define the *loci* contained inside the superlocus.
   * At this stage, monoexonic and multiexonic transcript are checked for overlaps
   * Moreover, two multiexonic transcripts are considered to belong to the same locus if they share a splice *site* (not junction)

The resulting, cleaned file does contain therefore only one transcript per locus.
The criteria used to select the "*best*" transcript are left to the user's discretion.

##Implementation
The pipeline is implemented in Python3, and does not require any non-standard library.
It should be usable on any system, but so far it has been tested only on Ubuntu 14.04.

##Input files
The pipeline can use as input files any properly GTF or GFF3 file. The only requirements are:
* GTF files must have a "transcript/mRNA" line for each transcript present in the annotation
  * Standard GTF files, which have only CDS/exon/start/stop features, are not interpreted correctly
* All files must be sorted by chromosome, start and stop.

Additionally, the pipeline can also use the output of coding potential programs such as [TransDecoder](https://transdecoder.github.io/ "TransDecoder (Find Coding Regions Within Transcripts)").
In order to give this type of data to the pipeline, the coding regions should be encoded into a BED12
file, with the "thickStart/thickEnd" fields (7 and 8) indicating the start and stop of the CDS.
As in transdecoder, it is expected that the coding potential has been defined using oriented transcript sequences.
Any multiexonic CDS found on the negative transcript strand will be therefore ignored, as it happens in TD.
