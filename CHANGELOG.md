#Version 1.0

Changes in this release:

- **MAJOR**: solved a bug which caused a failure of clustering into loci in rare occasions. Through the graph clustering, now Mikado is guaranteed to group monoloci correctly.
- **MAJOR**: When looking for fragments, now Mikado will consider transcripts without a strand as being on the **opposite** strand of neighbouring transcripts. This prevents many monoexonic, non-coding fragments from being retained in the final output.
- **MAJOR**: now Mikado serialise also stores the ***frame*** information of transcripts. Hits on the opposite strand will be **ignored**. This requires to **regenerate all Mikado databases**.
- **MAJOR**: Added the final configuration files used for the article.
- Added three new metrics, "blast_target_coverage", "blast_query_coverage", "blast_identity"
- Changed the *default* repertoire of valid AS events to J, j, G, h (removed C and g).          
- **Bug fix**: now Mikado will consider the cDNA/CDS overlap also for monoexonic transcripts, even when the "simple_overlap_for_monoexonic_loci" flag is set to true.
- Solved some issues with the Daijin schemas, which prevented correct referencing.
- Bug fix for finding retained introns - Mikado was not accounting for cases where an exon started within an intron and crossed multiple subsequent junctions.
- BF: Loci will never purge transcripts
- After creating the final loci, now Mikado will check for, and remove, any AS event transcript which would cross into the AS event.

#Version 1.0.0beta10

Changes in this release:

- **MAJOR**: re-written the clustering algorithm for the MonosublocusHolder stage. Now a holder will accept another monosublocus if:
    - the cDNA and CDS overlap is over a user-specified threshold
    *OR* 
    - there is some intronic overlap
    *OR*
    - one intron of either transcript is completely contained within an exon of the other.
    *OR*
    - at least one of the transcripts is monoexonic and there is some overlap of any kind. This behaviour (which was the default until this release) can be switched off through pick/clustering/simple_overlap_for_monoexonic (default true).
- **MAJOR**: changed slightly the anatomy of the configuration files. Now "pick" has two new subsections, "clustering" and "fragments".
    - Clustering: dedicated to how to cluster the transcripts in the different steps. Currently it contains the keys:
        - "flank"
        - "min_cdna_overlap" and "min_cds_overlap" (for the second clustering during the monosublocusHolder phase)
        - "cds_only": to indicate whether we should only consider the CDS for clustering after the initial merging in the Superlocus.
        - "simple_overlap_for_monoexonic": to switch on/off the old default behaviour with monoexonic transcripts
        - "purge": whether to completely exclude failed loci, previously under "run_options"
    - Fragments: dedicated to how to identify and treat putative fragments. Currently it contains the keys:
        - "remove": whether to exclude fragments, previously under "run_options"
        - "valid_class_codes": which class codes constitute a fragment match. Only class codes in the "Intronic", "Overlap" (inclusive of _) and "Fragment" categories are allowed.
        - max_distance: for non-overlapping fragments (ie p and P), maximum distance from the gene.
- Solved a long-standing bug which caused Mikado compare to consider as fusion also hits on the opposite strand only.
- Mikado compare now also provides the location of the matches in TMAP and REFMAP files.
- Introduced a new utility, "class_codes", to print out the information of the class codes. The definition of class codes is now contained in a subpackage of "scales".
- The "metrics" utility now allows for interactive querying based on category or metric name.
- The class code repertoire for putative fragments has been expanded, and made configurable through the "fragments" section.
- When printing out putative fragments, now Mikado will indicate the class code of the fragment, the match against which it was deemed a fragment of, and the distance of said fragment (if they are not overlapping). 
- Deprecated the "discard_definition" flag in Mikado serialise. Now Mikado will infer on its own whether to use the definition or the ID for serialising BLAST results.
- Now AbstractLocus implementations have a private method to check the correctness of the json_conf. As a corollary, Transcript and children have been moved to their own subpackage ("transcripts") in order to break the circular dependency Mikado.loci.Abstractlocus <- Mikado.configurator <- Mikado.loci.Transcript. *Technical note*: checking the consinstency of the configuration is an expensive operation, so it will be executed on demand rather than automatically.
- The methods to calculate scores and metrics have been moved to the AbstractLocus class, so to minimize the incidence of bugs due to code duplication and diversion.
- Made the checks for the scoring files more robust.
- Re-written the "find_retained_introns" method of AbstractLocus, to solve some bugs found during the utilisation of last beta. As a corollary, expanded the intervaltree module to allow searches for "tagged" intervals.
- Now the "monoloci_out" files contain the Monosublocus**Holder** step, not the Monosublocus step. This should help during fine-tuning. 
- Minimal requirements for alternative splicing events are now specified with a syntax analogous to that of minimal requirements, and that for not considering a locus as a putative fragment, under the tag "as_requirements".
- Fixed a bug which caused transcript requirements to be ignored if pick/clustering/purge was set to False.
- Mikado now supports also Python3.6.


#Version 1.0.0beta9 - "External scores"

Changes in this release:

- Mikado prepare, stats and compare are capable of using GFFs with "match/match_part" or "cDNA_match" features as input. This allows eg. to obtain sensible statistics for the alignments of Trinity transfrags or long reads when using GMAP inside Daijin.
- When the option "subloci_from_cds_only" is set to true (or the equivalent command line switch is invoked), AS events class codes will be calculated on the coding portion only of the transcript.
- Mikado now allows to perform unittest on the installed module, through the function "Mikado.test". This is simply a wrapper to NoseTest, through a call to the numpy tester.	
- **IMPORTANT**: in strand-specific mode, now Mikado prepare will not flip transcripts which have been misassigned to the opposite strand. Those transcripts will be kept on the original strand, and tagged.
- **IMPORTANT**: now Mikado will use all the verified introns found in a **Superlocus** to determine the fraction of verified introns per locus in each transcript. At the stage of Locus, ie creation of genes, this will revert to check only the verified introns in the locus. Also, solved a bug related to this metric.
- **IMPORTANT**: at the stage of the creation of monosubloci-holders, now Mikado groups together also transcripts **for which at least one intron is completely contained within the exon of another**. This should solve spurious cases where we called multiple loci instead of one, while probably avoiding trans-splicing.
- **IMPORTANT**: otionally now Mikado can perform the second clustering of transcripts based on simple overlap (either of the CDS or of all exonic features). The option can be invoked from the command line. 
- Two new metrics:
  - "suspicious_splicing" will allow to identify transcripts with mixed_splices, or which would have at least one canonical intron if moved to the opposite strand.
  - "only_non_canonical_splicing" will allow to identify transcripts whose splicing sites are all non-canonical.
- It is now possible to give Mikado a tab-delimited file of pre-calculated metrics (which must be numeric), during serialise. The file should have the transcript ids in the first column and have a header as first line; this header must have "TID" as first field, and no repeated fields afterwards. External metrics can be specified in the scoring configuration using the syntax "external.{name of the score}". If an inexistent metric is asked for, Mikado will assign a default value of 0 to it.
- It is now possible to use metrics with values between 0 and 1, inclusive directly as scoring, by specifying the parameter "use_raw: True". This is available only for metrics which have been tagged as being "usable raw", or with externally provided metrics. The option is valid only when looking for the maximum or minimum value for a metric, not when looking for a target. If an incorrect configuration is specified, Mikado will crash.
- Mikado prepare in "lenient" mode will keep also transcripts with a mixture of strands for the splicing junctions. Such transcripts are marked with the "suspicious_splicing" GTF attribute.
- Mikado prepare can be asked to keep all transcripts, even if they are redundant. The new behaviour (disabled by default) is switched on by the boolean parameter "prepare/keep_redundant".
- Mikado pick can consider transcripts with CDS ending within a CDS intron as truncated due to a retained intron event. This potentially allows Mikado to detect retained introns even when only CDSs are provided. The behaviour is disabled by default, and can be switched on using the boolean configuration parameter "pick/run_options/consider_truncated_for_retained".
- Some bugs have been detected and solved thanks to the collaboration with Hugo Darras.
- Many tests, including system ones, have been added to the test suite. Tests have been fixed for Python3.4 as well.

#Version 1.0.0beta7.2

Mainly fixes for Daijin, in order to be able to handle different versions of Trinity, GMAP and PortCullis.

#Version 1.0.0beta7.1

BugFix release:

- Corrected a bug identified by Hugo Darras, whereby Mikado crashed when asked not to report alternative splicing events
- Mikado compare now supports compressing the outputs
- Bug fix for Mikado util stats - now it functions also when exonic features contain "derived_from" tags
- Bug fix for bam2gtf.py

#Version 1.0.0beta7 - "Storing the stacks"

Changes:

- Now we use ujson and SQLite to store out-of-memory the data loaded during the "prepare" step, massively reducing the memory needed for this step.
- TransDecoder is now optional in the Daijin pipeline
- Mikado can be instructed to pad transcripts at their ends, to make all transcripts at a locus compatible in terms of their ends (eventually extended the CDS, if possible). The protocol was inspired by the AtRTD2 release: http://biorxiv.org/content/early/2016/05/06/051938

#Version 1.0.0beta6

Beta 6, Hopefully final release for the article. Highlights:

- Compatibility with DIAMOND
- Essential bugfix for handling correctly verified introns
- Updated scoring files

#Version 1.0.0beta5

First public release of Mikado.

Changelog:

- Added a flank option in Daijin
- Documentation ready, inclusive of a tutorial for Daijin itself
- Bug fixes, especially to Daijin

#Version 1.0.0beta4

Changelog:

- Daijin now has a proper schema and a proper configure command. Still to implement: schema check on starting Daijin itself
- Documentation is now more complete (see eg sections on Compare and Configuration)
- Now mo and O have been renamed to g and G respectively

#Version 1.0.0beta3

Changelog:

- Solved a bug which prevented correct loading of ORFs in monoexonic transcripts with an ORF on the negative strand
- Dagger is now daijin
- Added controls in Daijin for:
    - Minimum length of TD ORFs
    - rerun from a particular target
- Minor polishes to command line interfaces and the like.

#Version 1.0.0beta2

Changes:

- Small bug fixes for DAGGER
- Introduced a very basic configuration utility for DAGGER
- Now Mikado programs have a sensible help message, including description.

#Version 1.0.0beta1

Finally almost ready! Major changes:

- mikado is now completely integrated with DAGGER, h/t Daniel Mapleson
- Both mikado and dagger are called as "mikado" and "dagger", without any trailing ".py"
- We are now in feature freeze, with the exception of dagger.

#Version 0.24.0 "Blade sharpening"

This release is mainly focused on two aims: integration with Dagger, and the possibility of performing a good self-comparison with mikado.py compare.

Detailed changes for the release:

- scoring files are now inside a proper subfolder, and can be copied during configuration, for user modification
- the configuration process has now been widely rehauled, together with the configuration files, to make Mikado compatible with DAGGER
- Introduced a functioning deep copy method for gene objects
- Now redundant adding of exons into a transcript object will not throw an exception but will be silently ignored
- Mikado parsers now support GZIPped and BZIPped files, through python-magic
- added the possibility of specifying which subset of input assemblies are strand-specific during the preparation step
- configure and prepare now accept a tab-separated text file (format: , ) as input
- compare now can perform two types of self-analysis:
    - "internal", which outputs only a tmap file and performs a comparison all-vs-all within each multi-isoform gene locus
    - "self", in which each transcript is analysed as it would during a normal run, but removing itself from the reference.
- Now serialise can accept FASTA references, pyFaidx is used to recover the correct FAI
- pick and configure can now specify the mode of action (nosplit, split, stringent, lenient, permissive) by CL instead of having to modify the configuration file
- cleaned up the sample_data directory, plus now the configuration file is created on the spot by Snakemake (ensuring compatibility with the latest version)

#Version 0.23.2 "Stats bugfix"

Changes for this micro release:

- Bug fixes for some debug log messages
- Bug fixes for the stats utility
- added the flank option for mikado

#Version 0.23.1 "Faster stats"

Changes:

- Stats is now much faster and lighter, as objects are discarded upon analysis, not kept around for the whole run duration
- Pick now preferentially outputs only one ORF per transcript, behaviour regulated by pick/output_format/report_all_orfs
- BugFix for the detection of fragments during pick

#Version 0.23.0 "Clap-o-meter"

New features:

- Added a converter for (GTF <=> GFF3) > BED12 conversions in mikado.py util convert
- Now it is possible to provide a predefined malus/bonus to transcripts according to their sources using the following syntax in the configuration file:

        pick:
          source_score:
            source1: score1
            source2: score
            ...

It is generally *not* advised to use such an option unless the user really wants to emphasize or penalise a particular set of transcripts, eg. those coming from PacBio data.

#Version 0.22.0 "FastXML"

Transitioned the XML parser to the experimental Bio.SearchIO module, which uses internally C-level APIs and is therefore twice as fast than the canonical implementation. This required some rewriting at the class level in Mikado and a cloning of the original class to make it a proper iterative class.
Moreover, now requirements for the scoring can accept also eg boolean values.

#Version 0.21.0 "Trim galore"

Major changes for this release:

- Trim is now functioning
- ORFs from transdecoder are treated correctly - even truncated and internal ones. No more faux UTRs!
- Changes in the DB structure to accomodate the ORF information - note, this means that the databases have to be regenerated

