
.. _adapting-mikado:

Adapting Mikado to specific user-cases
======================================

Although Mikado provides generally sane defaults for most projects and species, one of its key advantages is its flexibility, which allows it to be tailored for the needs of various projects. In this section we provide an overview on how to adapt the workflow to specific user-cases.

In general, this tailoring can be performed in these sections of the workflow:

.. add links to these bullet points

- when launching ``mikado configure``, to set up influential parameters such as the flank clustering distance or the expected intron size range.
- by modifying directly the resulting configuration file.
- by modifying the :ref:`scoring file <configure-scoring-tutorial>`.

.. _adapting-case-one:

Case study 1: adapting Mikado to your genome of interest
--------------------------------------------------------

When adapting Mikado to a new species, some of the most important factors to be considered are:

- how compact is the species' genome - ie, what is the expected genomic distance between two neighbouring genes?
- what is the expected intron size range?
- what is the expected UTR/CDS ratio, ie, will transcripts generally have short UTR sections (e.g. *Arabidopsis thaliana*, with its average ~80% coding section) or will they instead have long, multiexonic UTRs (as is the case in e.g. *Homo sapiens*)?
- will the organism mostly possess multi-exonic, long genes with many splicing variants (as is the case for mammals, e.g. our own species) or does it instead harbour mostly short transcripts with a low number of exons - potentially, even, mostly monoexonic (as is the case for many fungi, e.g., *Saccharomyces cerevisiae*)?

.. we are using pombe rather than cerevisiae because the annotation for cerevisiae is very simplistic: only monoexonic genes without UTR or splicing events.

A good starting point for understanding how to answer these questions is to head over to `EnsEMBL <https://www.fungi.ensembl.org/index.html>`_, and look for an annotated species sufficiently similar to the one under examination. For example, let us say that we want to annotate a yeast in the same family of *S. pombe*. We can fetch its `GTF annotation <ftp://ftp.ensemblgenomes.org/pub/fungi/release-38/gff3/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.38.gff3.gz>`_ from the `FTP download page of EnsEMBL <http://fungi.ensembl.org/info/website/ftp/index.html>`_ and start analysing it. We can use to this end mikado util stats:

.. code-block:: bash

    mikado util stats Schizosaccharomyces_pombe.ASM294v2.38.gff3.gz Schizosaccharomyces_pombe.ASM294v2.38.stats

which returns the following:

================================  ========  =========  ============  ======  ======  ======  ======  =====  ========  =====  =====  =====  =====  ======
Stat                              Total     Average    Mode          Min     1%      5%      10%     25%    Median    75%    90%    95%    99%    Max
================================  ========  =========  ============  ======  ======  ======  ======  =====  ========  =====  =====  =====  =====  ======
Number of genes                   7268      NA         NA            NA      NA      NA      NA      NA     NA        NA     NA     NA     NA     NA
Number of genes (coding)          5145      NA         NA            NA      NA      NA      NA      NA     NA        NA     NA     NA     NA     NA
Number of monoexonic genes        4715      NA         NA            NA      NA      NA      NA      NA     NA        NA     NA     NA     NA     NA
Transcripts per gene              7269      1.00       1             1       1       1       1       1      1         1      1      1      1      2
Coding transcripts per gene       5146      0.71       1             0       0       0       0       0      1         1      1      1      1      2
CDNA lengths                      12610347  1,734.81   72            47      72      83      306     846    1,504     2,322  3,370  4,103  5,807  15,022
CDNA lengths (mRNAs)              10592917  2,058.48   1315          75      283     546     753     1,201  1,792     2,613  3,675  4,425  6,404  15,022
CDS lengths                       7178717   987.58     0             0       0       0       0       0      750       1,461  2,308  3,008  4,944  14,775
CDS lengths (mRNAs)               NA        1,395.01   354;750       75      201     307     387     669    1,137     1,752  2,661  3,420  5,502  14,775
CDS/cDNA ratio                    NA        67.48      100.0         5       15      28      37      53     70        84     93     100    100    100
Monoexonic transcripts            4715      1,656.09   72            47      71      74      119     695    1,413     2,276  3,386  4,124  5,950  14,362
MonoCDS transcripts               2753      1,485.99   375;432;4002  93      218     332     418     732    1,191     1,806  2,849  3,784  5,875  14,154
Exons per transcript              12633     1.74       1             1       1       1       1       1      1         2      3      4      7      16
Exons per transcript (mRNAs)      3081      2.03       1             1       1       1       1       1      1         3      4      5      7      16
Exon lengths                      NA        998.21     72            2       25      62      79      171    538       1,474  2,522  3,298  5,002  14,362
Exon lengths (mRNAs)              NA        1,013.39   106           3       24      59      86      175    499       1,516  2,605  3,414  5,117  14,362
Intron lengths                    NA        83.69      49            1       17      38      41      46     56        85     162    226    411    2,526
Intron lengths (mRNAs)            NA        84.16      46            1       36      39      41      46     56        86     162    227    412    2,526
CDS exons per transcript          2091      1.41       1             0       0       0       0       0      1         2      3      4      7      16
CDS exons per transcript (mRNAs)  2091      1.99       1             1       1       1       1       1      1         3      4      5      7      16
CDS exon lengths                  7178717   702.76     69;99         1       8       29      49      106    307       990    1,797  2,474  4,370  14,154
CDS Intron lengths                414501    81.77      44;45         0       35      38      40      45     55        84     157    220    395    2,525
5'UTR exon number                 5146      0.95       1             0       0       0       1       1      1         1      1      1      2      3
3'UTR exon number                 5146      0.93       1             0       0       0       1       1      1         1      1      1      2      3
5'UTR length                      1372304   266.67     0             0       0       0       17      71     154       309    586    932    1,935  4,397
3'UTR length                      2041896   396.79     0             0       0       0       46      126    243       441    865    1,386  2,644  5,911
Stop distance from junction       NA        7.58       0             0       0       0       0       0      0         0      0      0      23     3,385
Intergenic distances              NA        -60.33     -66           -9,461  -3,753  -2,180  -1,382  -161   64        365    842    1,279  2,698  31,961
Intergenic distances (coding)     NA        297.72     -66           -7,815  -3,477  -1,598  -440    -32    176       600    1,302  1,913  3,924  78,421
================================  ========  =========  ============  ======  ======  ======  ======  =====  ========  =====  =====  =====  =====  ======


From this table we can already see the following:

    - Most genes (5145 out of 7268, or 70.9%) are monoexonic
    - The average and modal intergenic distance between genes are very small, with almost half of the recorded distances being negative - indicating that most genes are actually *overlapping*.
    - Only a very small handful of genes (less than 1%) is annotated as having any splicing event
    - On average, UTRs occupy 33% of the length of coding transcripts (CDS/cDNA ratio is at 67%, on average) but most often transcripts actually lack an UTR at all (mode of 100%).
    - 98% of the introns have a length between 36 and 412 bps.

On the basis of this information, we can now start to customize the behaviour of Mikado for the species.

Creating the scoring file
^^^^^^^^^^^^^^^^^^^^^^^^^

The first step in the process is for us to create a scoring file, following the :ref:`tutorial on the subject <configure-scoring-tutorial>`. We will call it "spombe.yaml"; as detailed in the link before, we will write it in the textual `YAML format <http://yaml.org/spec/1.2/spec.html>`_.

Following the indications above and those in the tutorial, we should make the following changes in terms of priority for transcripts:

- we want mostly monoexonic transcripts
- transcripts with a UTR ratio under 33%
- we should look at most to 1 UTR exon, each way, targeting 0 (most transcripts are monoexonic and have their UTR contained in the same exon as the ORF).
- the distance of the stop codon from the nearest junction should be 0 (again this follows from having mostly monoexonic transcripts).

The scoring section would therefore end up looking like this:

.. code-block:: yaml
    :emphasize-lines: 6-9,17,21,22,26,29-30,32,40

    scoring:
          snowy_blast_score: {rescaling: max}
          cdna_length: {rescaling: max}
          cds_not_maximal: {rescaling: min}
          cds_not_maximal_fraction: {rescaling: min}
          exon_num: {
            rescaling: target,
            value: 1
          }
          five_utr_length:
            filter: {operator: le, value: 2500}
            rescaling: target
            value: 100
          five_utr_num:
            filter: {operator: lt, value: 4}
            rescaling: target
            value: 0
          end_distance_from_junction:
            filter: {operator: lt, value: 23}
            rescaling: min
          highest_cds_exon_number: {rescaling: target, value: 1}
          intron_fraction: {rescaling: min}
          is_complete: {rescaling: target, value: true}
          number_internal_orfs: {rescaling: target, value: 1}
          non_verified_introns_num: {rescaling: min}
          # proportion_verified_introns_inlocus: {rescaling: max}
          retained_fraction: {rescaling: min}
          retained_intron_num: {rescaling: min}
          selected_cds_fraction: {rescaling: target, value: 1, filter: {operator: gt, value: 0.7 }}
          # selected_cds_intron_fraction: {rescaling: max}
          selected_cds_length: {rescaling: max}
          selected_cds_num: {rescaling: target, value: 1}
          three_utr_length:
            filter: {operator: le, value: 2500}
            rescaling: target
            value: 200
          three_utr_num:
            filter: {operator: lt, value: 2}
            rescaling: target
            value: 0
          combined_cds_locus_fraction: {rescaling: max}

Now that we have codified the scoring part, the next step is to determine the :ref:`requirements <scoring-tutorial-first-reqs>` regarding the transcripts that should be accepted into our annotation. Given the simplicity of the organism, we can satisfy ourselves with the following two requirements:

    - No transcript should be shorter than 75 bps (minimum length for coding transcripts)
    - No transcript should have an intron longer than ~2600 bps (in the annotation the maximum is 2,526); we can be slightly more permissive here and set the limit at 3,000 bps.

This will yield the following, very simple requirements section:

.. code-block:: yaml

    requirements:
        expression:
            - cdna_length and max_intron_length
        parameters:
            cdna_length: {operator: ge, value: 75}
            max_intron_length: {operator: lt, value: 3000}


Modifying the general configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The second step in the customization process is to personalize the general configuration. On the basis of what we know of *S. pombe*, we have to intervene here in the following way:

- set the intron range: per above, a reasonable setting should be 36-412.
- set the clustering flank: given the very compact size of the genome, we should aim for something very small - probably 50bps is plenty.
- given the very compact size of the genome and the general lack of splicing, it is also advised to set Mikado to split any chimeric transcripts - the chances are very, very high that any such occurrence is artifactual.
- make alternative splicing calling a very rare occurrence

First of all, we will download `our genome <ftp://ftp.ensemblgenomes.org/pub/fungi/release-38/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna_sm.toplevel.fa.gz>`_ in a single file (genome.fasta). We will use the pretty boilerplate *A. thaliana* scoring configuration as our starting block, and we will ask Daijin to copy it to the current location.

.. code-block:: bash

    daijin configure --scoring spombe.yaml \
        --flank 50 \
        -i 36 412 \
        -m split \
        -o configuration.yaml \
        --genome genome.fasta

Once the configuration file has been created, we have to perform another couple of modifications, to make Mikado more stringent in terms of alternative splicing events. Look for the section :ref:`mikado/pick <pick>`. Here we can do the following:

1. If you are completely uninterested in alternative splicing events, you can just set the "report" flag to false. This will disable AS calling completely.
2. If you want to still report AS events but at a far lower rate, you can:

    - reduce the number of maximum isoforms reported: from 5 to 2, for example. **Note**: reducing this number to 1 will have the same effect as disabling AS calling completely.
    - restrict the types of AS events we call (see :ref:`the class code section <class-codes>` for more details). We can for example restrict the calling to "j" and "G", and potentially add "g" (i.e. consider as a valid alternative splicing event for a multiexonic transcript a monoexonic one).
    - increase the minimum score percentage of an AS event for it to be reported, to extremely high values (such as 0.9 to 0.99). This will ensure that only a small amount of isoforms will be called.
    - increase the minimum cDNA/CDS overlap between the AS events and the primary transcript. This cannot go up to 100% for both, otherwise no AS event will ever be reported. However, you could for example set the CDS overlap to 100%, if you are only interested in alternative UTR splicing.
    - leave the "keep_retained_introns" field as false, and "only_confirmed_introns" field as "true".

Once these modifications have been made, Mikado is ready to be run.

.. _adapting-case-two:

Case study 2: noisy RNA-Seq data
--------------------------------

With RNA-Seq, a relatively common happenstance is the presence of noise in the data - either experimentally, through the presence of pre-mRNA, genomic contamination, or otherwise erroneous transcripts; or from computational artifacts, e.g. an explicit choice on the part of the experimenter to retrieve from the data even isoforms and loci with little coverage support, in an attempt to boost the sensitivity of the analysis at the cost of decreased precision.

In such instances, it might make sense to make Mikado more stringent than usual. In this tutorial we will focus on the following:

- Making Mikado more aggressive in filtering out putative fragments
- Making Mikado more aggressive in splitting chimeric transcripts
- Making Mikado more aggressive in filtering out incorrect alternative splicing events such as retained introns

For ease of discussion, we will suppose that we are working in a species similar in features to *D. melanogaster*. We will, therefore, be using a copy of the dmelanogaster_scoring.yaml file included in the distribution of Mikado.

.. _adapting-case-two-general:

Modifying the general configuration file and obtaining a copy of the original template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before touching the scoring file, this time we will call the Daijin configurator in order to obtain a copy of the original *D. melanogaster* scoring file.
We will suppose to have relevant proteins in "proteins.fasta" (e.g. a dataset assembled from SwissProt), and that - like for *D. melanogaster* - the acceptable intron size range is between 50 and 26000 bps. As the data is quite noisy, we have to expect that there will be fragments derived from mis-alignments or genomic contamination; we will, therefore, enlarge the normal flanking area to 2000 bps. This will allow to catch more of these events, when we check for potential fragments in the neighbourhood of good loci. Regarding probable chimeric events, we will be quite aggressive - we will split any chimeric event which is not supported by a good blast hit against the database ("-m permissive").

.. code-block:: bash

    daijin configure \
        --scoring dmelanogaster_scoring.yaml --copy-scoring noisy.yaml  \
        --flank 2000 \
        -i 50 26000 \
        -m permissive \
        -o configuration.yaml \
        --genome genome.fasta \
        --prot-db proteins.fasta

Once created, the configuration file should be modified as follows:

    - in the pick/alternative_splicing section:

        - increase the stringency for calling an alternative splicing event:

            - min_score_percentage: from 0.5 to 0.75
            - max_isoforms: from 5 to 3
    
    - in the pick/fragments section:

        - add "I" (multi-exonic and within an intron of the reference locus) to the list of valid_class_codes

Please note that by default Mikado will look for alternative splicing events that have all introns not shared with the primary transcript to be confirmed externally. Also, it will exclude any transcript with retained introns. We should keep these options on their default value, as they will already contribute a significantly to reducing the number of spurious splicing events.

Customising the scoring file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Looking at the scoring section of the file, we do not need to apply anything particular here - the predefined definitions will already reward coding, homology-supported transcripts.

.. admonition:: Default scoring for *D. melanogaster*
    :class: toggle

    .. code-block:: yaml
       :linenos:

        scoring:
          snowy_blast_score: {rescaling: max}
          cdna_length: {rescaling: max}
          cds_not_maximal: {rescaling: min}
          cds_not_maximal_fraction: {rescaling: min}
          # exon_fraction: {rescaling: max}
          exon_num: {
            rescaling: max,
            filter: {
            operator: ge,
            value: 3}
          }
          five_utr_length:
            filter: {operator: le, value: 2500}
            rescaling: target
            value: 100
          five_utr_num:
            filter: {operator: lt, value: 4}
            rescaling: target
            value: 2
          end_distance_from_junction:
            filter: {operator: lt, value: 55}
            rescaling: min
          highest_cds_exon_number: {rescaling: max}
          intron_fraction: {rescaling: max}
          is_complete: {rescaling: target, value: true}
          number_internal_orfs: {rescaling: target, value: 1}
          # proportion_verified_introns: {rescaling: max}
          non_verified_introns_num: {rescaling: min}
          proportion_verified_introns_inlocus: {rescaling: max}
          retained_fraction: {rescaling: min}
          retained_intron_num: {rescaling: min}
          selected_cds_fraction: {rescaling: target, value: 0.8}
          selected_cds_intron_fraction: {rescaling: max}
          selected_cds_length: {rescaling: max}
          selected_cds_num: {rescaling: max}
          three_utr_length:
            filter: {operator: le, value: 2500}
            rescaling: target
            value: 200
          three_utr_num:
            filter: {operator: lt, value: 3}
            rescaling: target
            value: 1
          combined_cds_locus_fraction: {rescaling: max}

We can and should, however, modify the minimum requirements for transcripts in general, for alternative splicing events, and for not considering a given locus as a putative fragment.

First off, for the minimum requirements, we will tweak the requirements in this way:

    - discard any multiexonic transcript without verified introns. Normally we would discard such transcripts only if there are verified introns in the region. In this case, we would like to get rid of these transcripts altogether:

        - verified_introns_num: {operator: gt, value: 0}
        - If we would like to be really stringent, we could instead exclude any transcript with any amount of non-verified introns:

            - non_verified_introns_num: {operator: eq, value: 0}
    - discard any transcript with suspicious splicing events (ie splicing events that would be canonical if transferred on the opposite strand):

        - suspicious_splicing: {operator: eq, value: false}
    - let us also be more stringent on the maximum intron length, and decrease it from the permissive 150,000 to a much more stringent 30,000 (slightly higher than the 26,000 used for the "acceptable" intron range, above).

        - max_intron_length: {operator: le, value: 30000}
    - discard any monoexonic transcript without a CDS. This is more stringent than the default setting (where we keep non-coding monoexonic transcripts that have a some homology to a protein in the supplied database).

        - selected_cds_length.mono: {operator: gt, value: 0}

Altogether, this becomes:

.. code-block:: yaml

    requirements:
        expression:
            - ((exon_num.multi and (cdna_length.multi or selected_cds_length.multi)
            - and
            - max_intron_length and min_intron_length and verified_introns_num and suspicious_splicing)
            - or
            - (exon_num.mono and selected_cds_length.mono)))
        parameters:
            selected_cds_length.mono: {operator: gt, value: 300}
            cdna_length.multi: {operator: ge, value: 400}
            selected_cds_length.multi: {operator: gt, value: 200}
            exon_num.mono: {operator: eq, value: 1}
            exon_num.multi: {operator: gt, value: 1}
            max_intron_length: {operator: le, value: 30000}
            min_intron_length: {operator: ge, value: 20}
            verified_introns_num: {operator: gt, value: 1}
            suspicious_splicing: {operator: eq, value: false}

We should also adapt the requirements for alternative splicing events. Compared with the default settings, we can now remove the "suspicious_splicing" requirement - it is already present in the general requirements for a transcript, so it will never be invoked. However, we will make certain that no transcript with more than one ORF will ever be selected as an alternative splicing event: these transcripts are often generated by retained intron events or trans-splicing. It should be a rare event, but by putting a requirement here, we will ensure that no transcript of this kind will be brought back as ASE.

.. code-block:: yaml
    :emphasize: 2,8

    as_requirements:
      expression: [cdna_length and three_utr_length and five_utr_length and utr_length and number_internal_orfs]
      parameters:
        cdna_length: {operator: ge, value: 200}
        utr_length: {operator: le, value: 2500}
        five_utr_length: {operator: le, value: 2500}
        three_utr_length: {operator: le, value: 2500}
        number_internal_orfs: {operator: le, value: 1}

Finally, we will consider as fragmentary any non-coding transcript in the neighbourhood of a coding locus. We will consider as potentially fragmentary also any coding transcript with a short ORF (<100 aa, or 300 bps). The expression will be, in this case, very simple:

.. code-block:: yaml

    expression: [selected_cds_length]
    parameters:
        selected_cds_length: {operator: gt, value: 300}


Case study 3: comprehensive splicing catalogue
----------------------------------------------

There are cases in which we would like our annotation to be as comprehensive as possible, ie. to include transcripts that we would normally exclude from consideration. For example, we might want to study the prevalence of retained intron events in a sample, or keep events that do not have a good read coverage and whose introns might, therefore, be recognised as invalid by Portcullis. It is possible to tweak Mikado's behaviour to this end.

Modifying the general configuration file and obtaining a copy of the original template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Like in :ref:`the second case <adapting-case-two-general>`, we will presume to be working in a similar species to *D. melanogaster*. Again, we will create the configuration file thus:

.. code-block:: bash

    daijin configure \
        --scoring dmelanogaster_scoring.yaml --copy-scoring comprehensive.yaml  \
        --flank 200 \
        -i 50 26000 \
        -m permissive \
        -o configuration.yaml \
        --genome genome.fasta \
        --prot-db proteins.fasta

Notice that compared to the previous example we reduced the flanking distance to the standard value (200 bps instead of 2000 bps) as we are less worried of fragmentary loci.

In the configuration file, we will change the following:

    - under pick/alternative_splicing:

        - switch "keep_retained_introns" to true
        - switch "only_confirmed_introns" to false
        - potentially, increase the number of isoforms from 5 to 10 or higher
        - consult the documentation on :ref:`class codes <class-codes>` to verify which additional AS events you would like to keep; by default, Mikado will include cases where the transcript has at least a different splicing site (j), no splicing site in common with the original transcript but introns roughly coincident (h), novel introns in the terminal exons (J) or within the primary mono-exonic transcript (G).

            - For a comprehensive catalogue, we would recommend to include at least "C" (transcript roughly contained, but with "spilling" within the intron(s) of the primary transcript).
        - To include transcripts quite dissimilar from the primary, potentially lower the percentages for:

            - min_cds_overlap
            - min_cdna_overlap
            - min_score_perc

.. warning::
  The heuristics we are touching in this section are core to the precision of Mikado. For example, allowing Mikado to bring back retained intron events will, by definition, bring into the annotation transcripts that are normally ignored. Please consider this when configuring the run and later, when reviewing the results.

Customising the scoring file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case, as we are interested in retaining a greater variety of splicing events, we will concentrate only on one section of the file, ie the "as_requirements" section. Compared to the default settings, we are going to remove the UTR requirements, and bring back transcripts with long UTRs. This is because many transcripts with retained intron events will have, by default, longer UTRs than usual. We will rely on the general prioritisation instead to penalise these transcripts in general (and thus avoid bringing them back if they are or really poor quality). We will still exclude cases with "suspicious_splicing", ie cases most probably generated by trans-splicing.


.. code-block:: yaml

    as_requirements:
      expression: [cdna_length and suspicious_splicing]
      parameters:
        cdna_length: {operator: ge, value: 200}
        suspicious_splicing: {operator: ne, value: true}
