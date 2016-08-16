.. _Introduction:

.. |python_badge| image:: https://img.shields.io/pypi/pyversions/snakemake.svg?style=flat-square
   :target: https://www.python.org/
.. |snake_badge| image:: https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)]
   :target: http://snakemake.bitbucket.org

Introduction
============

|python_badge| |snake_badge|


Numerous algorithms have been proposed to analyse RNA-Seq data, both in terms of aligning the reads to a reference genome ([TopHat2]_, [STAR]_, [Hisat]_) or to assemble them to infer the sequence and structure of the original molecules present in the sample. The latter phase can be performed either by using alignment data against the reference genome ([Cufflinks]_, [StringTie]_, [Class2]_, [Trinity]_) or in the absence of such information ([Trinity]_, [Oases]_, [Bridger]_, ). Each of these methods has to contend with numerous sources of variability in the RNA-Seq data:
  * alternative splicing events at the same locus
  * extremely variable expression and therefore sequencing depth at different loci
  * presence of orthologous genes with very similar cDNAs.
  * tandem duplicated genes which can be artefactually reconstructed as a single, fused entity

Multiple assessments of these methods on real and simulated data. In particular, the RGASP consortium promoted a competition among some of the most popular methods in 2012 [RGASP]_; the results showed that the accuracy of the methods was dependent on the input data, the expression levels of the genes, and on the species under analysis. A different concern regards the completeness of each of the assembly methods, ie whether methods that are not necessarily the best across all gene loci might nonetheless be capable of outperforming other competitors for specific contigs. Recently, approaches have been proposed to integrate multiple RNA-Seq assemblies by inferring their protein-coding content and collapsing them accordingly [EviGeneTobacco]_ or to determine the best proposed assembly using multiple measures related to how the transcript compares to the input RNA-Seq reads ([Transrate]_, [RSEMeval]_).

A different but related line of research has focused on how to integrate data coming from multiple samples. While some researchers advocate for merging together the input data and assembling it [Trinity]_, others have developed methods to integrate multiple assemblies into a single coherent annotation ([CuffMerge]_, Stringtie-merge). Another common approach relies on a third party tool to integrate multiple assemblies together with other data - typically protein alignments and *ab initio* predictions - to select for the best model. The latter approach has been chosen by some of the most popular pipelines for genome annotation in the last years ([EVM]_, [Maker]_).

Our tool, Mikado, contributes to this area of research by proposing a novel algorithm to integrate multiple transcript assemblies, leveraging additional data regarding the position of ORFs (generally derived with Transdecoder [Trinity]_), sequence similarity (derived with BLAST [Blastplus]_) and reliable splicing junctions (generally derived using Portcullis [Portcullis]_). Through the combination of these input data sources, we are capable of identifying artefacts - especially gene fusions and small fragments - and retrieve the original transcripts in each locus with high accuracy and sensitivity.

An important caveat is that our approach does not look for the *best* candidate in terms of the input data, as Transrate [Transrate]_ or RSEM-Eval [RSEMeval]_ do. Rather, our approach is more similar to EvidentialGene [EviGeneTobacco]_ as Mikado will try to find the candidates most likely to represent the isoforms which would be annotated as canonical by a human annotator. Biological events such as intron retentions or trans-splicing that might be present in the sample are explicitly selected against, in order to provide a clean annotation of the genome.


The core algorithms
~~~~~~~~~~~~~~~~~~~

.. sidebar:: The Mikado algorithm

    .. figure:: Mikado_algorithm.jpeg
        :align: center
        :scale: 50%

    Schematic representation of the steps in the Mikado pipeline.

Mikado analyses an ensemble of RNA-Seq assemblies by looking for overlapping transcripts based on their genomic position. Such clusters of genes, or *superloci*, are then analysed to locate the gene loci that will form the final annotation. The algorithm presupposes that the transcripts have been ordered on the genome based on their location (see the second step in the pipeline, :ref:`mikado prepare <prepare>`.). Additionally, transcript data is integrated with external data sources in real time, while traversing the annotation file (see the third step in the pipeline, :ref:`mikado serialise <serialise>`, for details on which data can be integrated). The selection algorithm proper is implemented in the :ref:`last step of the pipeline <pick>`.

Transcripts are scored and selected according to user-defined rules, based on many different features of the transcripts themselves (cDNA length, CDS length, UTR fraction, number of reliable junctions, etc.; please see the :ref:`dedicated section on scoring <scoring_files>` for details on the scoring algorithm).

The detection and analysis of a locus proceeds as follows:

#. When the first transcript is detected, Mikado will create a *superlocus* - a container of transcripts sharing the same genomic location - and assign the transcript to it.
#. While traversing the genome, as long as any new transcript is within the maximum allowed flanking distance, it will be added to the superlocus.
#. When the last transcript is added, Mikado performs the following preliminary operations:
    #. Integrate all the data from the database (including ORFs, reliable junctions in the region, and BLAST homology).
    #. If a transcript is monoexonic, assign or reverse its strand if the ORF data supports the decision
    #. If requested and the ORF data supports the operation, split chimeric transcripts - ie those that contain two or more non-overlapping ORFs on the same strand.
    #. Split the superlocus into groups of transcripts that:
        * share the same strand
        * have at least 1bp overlap
    #. Analyse each of these novel "stranded" superloci separately.
#. Create *subloci*, ie group transcripts so to minimize the probability of mistakenly merging multiple gene loci due to chimeras. These groups are defined as follows:
    * if the transcripts are multiexonic, they must share at least one intron, inclusive of the borders
    * if the transcripts are monoexonic, they must overlap by at least 1bp.
    * Monoexonic and multiexonic transcripts *cannot* be part of the same sublocus.
#. Select the best transcript inside each sublocus:
    #. Score the transcripts (see the :ref:`section on scoring <scoring_files>`)
    #. Select as winner the transcript with the highest score and assign it to a *monosublocus*
    #. Discard any transcript which is overlapping with it, according to the definitions in the point above
    #. Repeat the procedure from point 2 until no transcript remains in the sublocus
#. *Monosubloci* are gathered together into *monosubloci holders*, ie the seeds for the gene loci. Monosubloci holder have more lenient parameters to group transcripts, as the first phase should have already discarded most chimeras. Once a holder is created by a single *monosublocus*, any subsequent candidate *monosublocus* will be integrated only if the following conditions are satisfied:
    * if the candidate is monoexonic, its exon must overlap at least one exon of a transcript already present in the holder
    * if the candidate is multiexonic and the holder contains only monoexonic transcripts, apply the same criterion, ie check whether its exons overlap the exons of at least one of the transcripts already present
    * if the candidate is multiexonic and the holder contains multiexonic transcripts, check whether its introns overlap the introns of at least one of the transcripts in the holder.
#. Once the holders are created, apply the same scoring and selection procedure of the sublocus selection step. The winning transcripts are assigned to the final *loci*. These are called the *primary transcripts of the loci*.
#. Once the loci are created, track back to the original transcripts of the superlocus:
    #. discard any transcript overlapping more than one locus, as these are probably chimeras.
    #. For those transcripts that are overlapping to a single locus, verify that they are valid alternative splicing events using the :ref:`class code <ccode>` of the comparison against the primary transcript.
#. Finally detect and either tag or discard fragments inside the initial *superlocus* (irrespective of strand):
    #. Check whether the primary transcript of any locus meets the criteria to be defined as a fragment (by default, maximum ORF of 30AA and maximum 2 exons - any transcript exceeding either criterion will be considered as non-fragment by default)
    #. If so, verify whether they are near enough any valid locus to be considered as a fragment (in general, class codes which constitute the "Intronic", "Fragmentary" and "No overlap" categories).
    #. If these conditions are met, tag the locus as a fragment. If requested, Mikado will just discard these transcripts (advised).

These steps help Mikado identify and solve fusions, detect correctly the gene loci, and define valid alternative splicing events.

Technical details
~~~~~~~~~~~~~~~~~

Most of the selection (ie "pick") stage of the pipeline relies on the implementation of the objects in the :ref:`loci submodule <sub-loci>`. In particular, the library defines an abstract class, "Abstractlocus", which requires all its children to implement a version of the "is_intersecting" method. Each implementation of the method is specific to the stage. So the *superlocus* class will require in the "is_intersecting" method only overlap between the transcripts, optionally with a flanking and optionally restricting the groups to transcripts that share the same strand. The *sublocus* class will implement a different algorithm, and so on.
The scoring is effectuated by first asking to recalculate the metrics (.calculate_metrics) and subsequently
to calculate the scores (.calculate_scores). Mikado will try to cache and avoid recalculation of metrics and scores as much as possible, to make the program faster.
