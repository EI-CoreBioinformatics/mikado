.. _pick:

Mikado pick
===========

This is the final stage of the pipeline, in which Mikado identifies gene loci and selects the best transcripts.

Input files
~~~~~~~~~~~

``mikado pick`` requires as input files the following:

#. A sorted GTF files with unique transcript names, derived through the :ref:`prepare stage <prepare>`.
#. A database containing all the external data to be integrated with the transcript structure, derived through the :ref:`serialisation stage <serialise>`.
#. A scoring file specifying the minimum requirements for transcripts and the relevant metrics for scoring. See the :ref:`section on scoring files <scoring_files>` for details.

.. _pick-output:

Output files
~~~~~~~~~~~~

``mikado pick`` will produce three kinds of output files: a GFF3 file, a *metrics* file, and a *scores* file. This triad will be produced for the *loci* level, and optionally also for the *subloci* and *monoloci* level.

GFF3 files
----------

This output file is a standard-compliant GFF, with the addition of the *superloci* to indicate the original spans. An example with two superloci, one on the negative and one on the positive strand, follows::

    ##gff-version 3
    ##sequence-region Chr5 1 26975502
    Chr5    Mikado_loci     superlocus      26584796        26601707        .       +       .       ID=Mikado_superlocus:Chr5+:26584796-26601707;Name=superlocus:Chr5+:26584796-26601707
    Chr5    Mikado_loci     gene    26584796        26587912        23      +       .       ID=mikado.Chr5G1;Name=mikado.Chr5G1;multiexonic=True;superlocus=Mikado_superlocus:Chr5+:26584796-26601707
    Chr5    Mikado_loci     mRNA    26584796        26587912        24      +       .       ID=mikado.Chr5G1.2;Parent=mikado.Chr5G1;Name=mikado.Chr5G1.2;alias=st_Stringtie_STAR.21710.1;canonical_junctions=1,2,3,4,5,6,7,8,9,10;canonical_number=10;canonical_proportion=1.0;ccode=j;cov=25.165945;primary=False
    Chr5    Mikado_loci     exon    26584796        26584879        .       +       .       ID=mikado.Chr5G1.2.exon1;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     five_prime_UTR  26584796        26584879        .       +       .       ID=mikado.Chr5G1.2.five_prime_UTR1;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26585220        26585273        .       +       .       ID=mikado.Chr5G1.2.exon2;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     five_prime_UTR  26585220        26585222        .       +       .       ID=mikado.Chr5G1.2.five_prime_UTR2;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26585223        26585273        .       +       0       ID=mikado.Chr5G1.2.CDS1;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26585345        26585889        .       +       0       ID=mikado.Chr5G1.2.CDS2;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26585345        26585889        .       +       .       ID=mikado.Chr5G1.2.exon3;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26585982        26586102        .       +       1       ID=mikado.Chr5G1.2.CDS3;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26585982        26586102        .       +       .       ID=mikado.Chr5G1.2.exon4;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26586217        26586294        .       +       0       ID=mikado.Chr5G1.2.CDS4;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26586217        26586294        .       +       .       ID=mikado.Chr5G1.2.exon5;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26586420        26586524        .       +       0       ID=mikado.Chr5G1.2.CDS5;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26586420        26586524        .       +       .       ID=mikado.Chr5G1.2.exon6;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26586638        26586850        .       +       0       ID=mikado.Chr5G1.2.CDS6;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26586638        26586850        .       +       .       ID=mikado.Chr5G1.2.exon7;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26586934        26586996        .       +       0       ID=mikado.Chr5G1.2.CDS7;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26586934        26586996        .       +       .       ID=mikado.Chr5G1.2.exon8;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26587084        26587202        .       +       0       ID=mikado.Chr5G1.2.CDS8;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26587084        26587202        .       +       .       ID=mikado.Chr5G1.2.exon9;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26587287        26587345        .       +       1       ID=mikado.Chr5G1.2.CDS9;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26587287        26587345        .       +       .       ID=mikado.Chr5G1.2.exon10;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     CDS     26587427        26587755        .       +       2       ID=mikado.Chr5G1.2.CDS10;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     exon    26587427        26587912        .       +       .       ID=mikado.Chr5G1.2.exon11;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     three_prime_UTR 26587756        26587912        .       +       .       ID=mikado.Chr5G1.2.three_prime_UTR1;Parent=mikado.Chr5G1.2
    Chr5    Mikado_loci     mRNA    26584930        26587912        23      +       .       ID=mikado.Chr5G1.1;Parent=mikado.Chr5G1;Name=mikado.Chr5G1.1;alias=st_Stringtie_STAR.21710.3;canonical_junctions=1,2,3,4,5,6,7,8,9,10;canonical_number=10;canonical_proportion=1.0;cov=2.207630;primary=True
    Chr5    Mikado_loci     exon    26584930        26585023        .       +       .       ID=mikado.Chr5G1.1.exon1;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     five_prime_UTR  26584930        26585023        .       +       .       ID=mikado.Chr5G1.1.five_prime_UTR1;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26585220        26585273        .       +       .       ID=mikado.Chr5G1.1.exon2;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     five_prime_UTR  26585220        26585222        .       +       .       ID=mikado.Chr5G1.1.five_prime_UTR2;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26585223        26585273        .       +       0       ID=mikado.Chr5G1.1.CDS1;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26585345        26585889        .       +       0       ID=mikado.Chr5G1.1.CDS2;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26585345        26585889        .       +       .       ID=mikado.Chr5G1.1.exon3;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26585982        26586102        .       +       1       ID=mikado.Chr5G1.1.CDS3;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26585982        26586102        .       +       .       ID=mikado.Chr5G1.1.exon4;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26586217        26586294        .       +       0       ID=mikado.Chr5G1.1.CDS4;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26586217        26586294        .       +       .       ID=mikado.Chr5G1.1.exon5;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26586420        26586524        .       +       0       ID=mikado.Chr5G1.1.CDS5;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26586420        26586524        .       +       .       ID=mikado.Chr5G1.1.exon6;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26586638        26586850        .       +       0       ID=mikado.Chr5G1.1.CDS6;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26586638        26586850        .       +       .       ID=mikado.Chr5G1.1.exon7;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26586934        26586996        .       +       0       ID=mikado.Chr5G1.1.CDS7;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26586934        26586996        .       +       .       ID=mikado.Chr5G1.1.exon8;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26587084        26587202        .       +       0       ID=mikado.Chr5G1.1.CDS8;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26587084        26587202        .       +       .       ID=mikado.Chr5G1.1.exon9;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26587287        26587345        .       +       1       ID=mikado.Chr5G1.1.CDS9;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26587287        26587345        .       +       .       ID=mikado.Chr5G1.1.exon10;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     CDS     26587427        26587755        .       +       2       ID=mikado.Chr5G1.1.CDS10;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     exon    26587427        26587912        .       +       .       ID=mikado.Chr5G1.1.exon11;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     three_prime_UTR 26587756        26587912        .       +       .       ID=mikado.Chr5G1.1.three_prime_UTR1;Parent=mikado.Chr5G1.1
    Chr5    Mikado_loci     gene    26588402        26592561        20      +       .       ID=mikado.Chr5G2;Name=mikado.Chr5G2;multiexonic=True;superlocus=Mikado_superlocus:Chr5+:26584796-26601707
    Chr5    Mikado_loci     mRNA    26588402        26592561        24      +       .       ID=mikado.Chr5G2.2;Parent=mikado.Chr5G2;Name=mikado.Chr5G2.2;alias=st_Stringtie_STAR.21710.9.split1;canonical_junctions=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21;canonical_number=21;canonical_proportion=1.0;ccode=j;cov=0.000000;primary=False
    Chr5    Mikado_loci     exon    26588402        26588625        .       +       .       ID=mikado.Chr5G2.2.exon1;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     five_prime_UTR  26588402        26588625        .       +       .       ID=mikado.Chr5G2.2.five_prime_UTR1;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26589203        26589279        .       +       .       ID=mikado.Chr5G2.2.exon2;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     five_prime_UTR  26589203        26589237        .       +       .       ID=mikado.Chr5G2.2.five_prime_UTR2;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26589238        26589279        .       +       0       ID=mikado.Chr5G2.2.CDS1;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26589386        26590167        .       +       0       ID=mikado.Chr5G2.2.CDS2;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26589386        26590167        .       +       .       ID=mikado.Chr5G2.2.exon3;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26590261        26590393        .       +       1       ID=mikado.Chr5G2.2.CDS3;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26590261        26590393        .       +       .       ID=mikado.Chr5G2.2.exon4;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26590495        26590566        .       +       0       ID=mikado.Chr5G2.2.CDS4;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26590495        26590566        .       +       .       ID=mikado.Chr5G2.2.exon5;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26590641        26590739        .       +       0       ID=mikado.Chr5G2.2.CDS5;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26590641        26590739        .       +       .       ID=mikado.Chr5G2.2.exon6;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26590880        26591092        .       +       0       ID=mikado.Chr5G2.2.CDS6;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26590880        26591092        .       +       .       ID=mikado.Chr5G2.2.exon7;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26591174        26591236        .       +       0       ID=mikado.Chr5G2.2.CDS7;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26591174        26591236        .       +       .       ID=mikado.Chr5G2.2.exon8;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26591324        26591442        .       +       0       ID=mikado.Chr5G2.2.CDS8;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26591324        26591442        .       +       .       ID=mikado.Chr5G2.2.exon9;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26591520        26591578        .       +       1       ID=mikado.Chr5G2.2.CDS9;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26591520        26591578        .       +       .       ID=mikado.Chr5G2.2.exon10;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26591681        26592002        .       +       2       ID=mikado.Chr5G2.2.CDS10;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26591681        26592002        .       +       .       ID=mikado.Chr5G2.2.exon11;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     CDS     26592528        26592561        .       +       1       ID=mikado.Chr5G2.2.CDS11;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     exon    26592528        26592561        .       +       .       ID=mikado.Chr5G2.2.exon12;Parent=mikado.Chr5G2.2
    Chr5    Mikado_loci     mRNA    26588402        26592561        20      +       .       ID=mikado.Chr5G2.1;Parent=mikado.Chr5G2;Name=mikado.Chr5G2.1;alias=st_Stringtie_STAR.21710.6.split3;canonical_junctions=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21;canonical_number=21;canonical_proportion=1.0;cov=0.000000;primary=True
    Chr5    Mikado_loci     exon    26588402        26588625        .       +       .       ID=mikado.Chr5G2.1.exon1;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     five_prime_UTR  26588402        26588625        .       +       .       ID=mikado.Chr5G2.1.five_prime_UTR1;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26589196        26589279        .       +       .       ID=mikado.Chr5G2.1.exon2;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     five_prime_UTR  26589196        26589237        .       +       .       ID=mikado.Chr5G2.1.five_prime_UTR2;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26589238        26589279        .       +       0       ID=mikado.Chr5G2.1.CDS1;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26589386        26590167        .       +       0       ID=mikado.Chr5G2.1.CDS2;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26589386        26590167        .       +       .       ID=mikado.Chr5G2.1.exon3;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26590261        26590393        .       +       1       ID=mikado.Chr5G2.1.CDS3;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26590261        26590393        .       +       .       ID=mikado.Chr5G2.1.exon4;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26590495        26590566        .       +       0       ID=mikado.Chr5G2.1.CDS4;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26590495        26590566        .       +       .       ID=mikado.Chr5G2.1.exon5;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26590641        26590739        .       +       0       ID=mikado.Chr5G2.1.CDS5;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26590641        26590739        .       +       .       ID=mikado.Chr5G2.1.exon6;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26590880        26591092        .       +       0       ID=mikado.Chr5G2.1.CDS6;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26590880        26591092        .       +       .       ID=mikado.Chr5G2.1.exon7;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26591174        26591236        .       +       0       ID=mikado.Chr5G2.1.CDS7;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26591174        26591236        .       +       .       ID=mikado.Chr5G2.1.exon8;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26591324        26591442        .       +       0       ID=mikado.Chr5G2.1.CDS8;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26591324        26591442        .       +       .       ID=mikado.Chr5G2.1.exon9;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26591520        26591578        .       +       1       ID=mikado.Chr5G2.1.CDS9;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26591520        26591578        .       +       .       ID=mikado.Chr5G2.1.exon10;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26591681        26592002        .       +       2       ID=mikado.Chr5G2.1.CDS10;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26591681        26592002        .       +       .       ID=mikado.Chr5G2.1.exon11;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     CDS     26592528        26592561        .       +       1       ID=mikado.Chr5G2.1.CDS11;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     exon    26592528        26592561        .       +       .       ID=mikado.Chr5G2.1.exon12;Parent=mikado.Chr5G2.1
    Chr5    Mikado_loci     gene    26592649        26595691        19      +       .       ID=mikado.Chr5G3;Name=mikado.Chr5G3;multiexonic=True;superlocus=Mikado_superlocus:Chr5+:26584796-26601707
    Chr5    Mikado_loci     mRNA    26592720        26595691        19      +       .       ID=mikado.Chr5G3.1;Parent=mikado.Chr5G3;Name=mikado.Chr5G3.1;alias=st_Stringtie_STAR.21710.7.split2;canonical_junctions=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20;canonical_number=20;canonical_proportion=1.0;cov=0.000000;primary=True
    Chr5    Mikado_loci     CDS     26592720        26593365        .       +       0       ID=mikado.Chr5G3.1.CDS1;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26592720        26593365        .       +       .       ID=mikado.Chr5G3.1.exon1;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26593449        26593836        .       +       2       ID=mikado.Chr5G3.1.CDS2;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26593449        26593836        .       +       .       ID=mikado.Chr5G3.1.exon2;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26593930        26594062        .       +       1       ID=mikado.Chr5G3.1.CDS3;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26593930        26594062        .       +       .       ID=mikado.Chr5G3.1.exon3;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26594172        26594243        .       +       0       ID=mikado.Chr5G3.1.CDS4;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26594172        26594243        .       +       .       ID=mikado.Chr5G3.1.exon4;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26594318        26594416        .       +       0       ID=mikado.Chr5G3.1.CDS5;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26594318        26594416        .       +       .       ID=mikado.Chr5G3.1.exon5;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26594569        26594772        .       +       0       ID=mikado.Chr5G3.1.CDS6;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26594569        26594772        .       +       .       ID=mikado.Chr5G3.1.exon6;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26594860        26594922        .       +       0       ID=mikado.Chr5G3.1.CDS7;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26594860        26594922        .       +       .       ID=mikado.Chr5G3.1.exon7;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26595003        26595121        .       +       0       ID=mikado.Chr5G3.1.CDS8;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26595003        26595121        .       +       .       ID=mikado.Chr5G3.1.exon8;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26595210        26595268        .       +       1       ID=mikado.Chr5G3.1.CDS9;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26595210        26595268        .       +       .       ID=mikado.Chr5G3.1.exon9;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     CDS     26595366        26595691        .       +       2       ID=mikado.Chr5G3.1.CDS10;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     exon    26595366        26595691        .       +       .       ID=mikado.Chr5G3.1.exon10;Parent=mikado.Chr5G3.1
    Chr5    Mikado_loci     mRNA    26592649        26595268        21      +       .       ID=mikado.Chr5G3.2;Parent=mikado.Chr5G3;Name=mikado.Chr5G3.2;abundance=2.390309;alias=cl_Chr5.6283;canonical_junctions=1,2,3,4,5,6,7,8;canonical_number=8;canonical_proportion=1.0;ccode=j;primary=False
    Chr5    Mikado_loci     exon    26592649        26593365        .       +       .       ID=mikado.Chr5G3.2.exon1;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     five_prime_UTR  26592649        26592719        .       +       .       ID=mikado.Chr5G3.2.five_prime_UTR1;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26592720        26593365        .       +       0       ID=mikado.Chr5G3.2.CDS1;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26593449        26593836        .       +       2       ID=mikado.Chr5G3.2.CDS2;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26593449        26593836        .       +       .       ID=mikado.Chr5G3.2.exon2;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26593930        26594095        .       +       1       ID=mikado.Chr5G3.2.CDS3;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26593930        26594095        .       +       .       ID=mikado.Chr5G3.2.exon3;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26594172        26594243        .       +       0       ID=mikado.Chr5G3.2.CDS4;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26594172        26594243        .       +       .       ID=mikado.Chr5G3.2.exon4;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26594318        26594416        .       +       0       ID=mikado.Chr5G3.2.CDS5;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26594318        26594416        .       +       .       ID=mikado.Chr5G3.2.exon5;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26594569        26594772        .       +       0       ID=mikado.Chr5G3.2.CDS6;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26594569        26594772        .       +       .       ID=mikado.Chr5G3.2.exon6;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26594860        26594922        .       +       0       ID=mikado.Chr5G3.2.CDS7;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26594860        26594922        .       +       .       ID=mikado.Chr5G3.2.exon7;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26595003        26595121        .       +       0       ID=mikado.Chr5G3.2.CDS8;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26595003        26595121        .       +       .       ID=mikado.Chr5G3.2.exon8;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     CDS     26595210        26595268        .       +       1       ID=mikado.Chr5G3.2.CDS9;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     exon    26595210        26595268        .       +       .       ID=mikado.Chr5G3.2.exon9;Parent=mikado.Chr5G3.2
    Chr5    Mikado_loci     gene    26596207        26598231        20      +       .       ID=mikado.Chr5G4;Name=mikado.Chr5G4;multiexonic=False;superlocus=Mikado_superlocus:Chr5+:26584796-26601707
    Chr5    Mikado_loci     mRNA    26596207        26598231        20      +       .       ID=mikado.Chr5G4.1;Parent=mikado.Chr5G4;Name=mikado.Chr5G4.1;alias=st_Stringtie_STAR.21710.6.split3;canonical_junctions=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21;canonical_number=21;canonical_proportion=1.0;cov=0.000000;primary=True
    Chr5    Mikado_loci     CDS     26596207        26598192        .       +       0       ID=mikado.Chr5G4.1.CDS1;Parent=mikado.Chr5G4.1
    Chr5    Mikado_loci     exon    26596207        26598231        .       +       .       ID=mikado.Chr5G4.1.exon1;Parent=mikado.Chr5G4.1
    Chr5    Mikado_loci     three_prime_UTR 26598193        26598231        .       +       .       ID=mikado.Chr5G4.1.three_prime_UTR1;Parent=mikado.Chr5G4.1
    Chr5    Mikado_loci     gene    26599417        26601137        20      +       .       ID=mikado.Chr5G5;Name=mikado.Chr5G5;multiexonic=True;superlocus=Mikado_superlocus:Chr5+:26584796-26601707
    Chr5    Mikado_loci     mRNA    26599417        26601137        20      +       .       ID=mikado.Chr5G5.1;Parent=mikado.Chr5G5;Name=mikado.Chr5G5.1;abundance=0.371780;alias=cl_Chr5.6286;canonical_junctions=1,2,3,4,5,6;canonical_number=6;canonical_proportion=1.0;primary=True
    Chr5    Mikado_loci     exon    26599417        26599654        .       +       .       ID=mikado.Chr5G5.1.exon1;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     five_prime_UTR  26599417        26599612        .       +       .       ID=mikado.Chr5G5.1.five_prime_UTR1;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     CDS     26599613        26599654        .       +       0       ID=mikado.Chr5G5.1.CDS1;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     CDS     26599767        26600053        .       +       0       ID=mikado.Chr5G5.1.CDS2;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     exon    26599767        26600053        .       +       .       ID=mikado.Chr5G5.1.exon2;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     CDS     26600151        26600244        .       +       1       ID=mikado.Chr5G5.1.CDS3;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     exon    26600151        26600244        .       +       .       ID=mikado.Chr5G5.1.exon3;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     CDS     26600314        26600394        .       +       0       ID=mikado.Chr5G5.1.CDS4;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     exon    26600314        26600394        .       +       .       ID=mikado.Chr5G5.1.exon4;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     CDS     26600497        26600616        .       +       0       ID=mikado.Chr5G5.1.CDS5;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     exon    26600497        26600616        .       +       .       ID=mikado.Chr5G5.1.exon5;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     CDS     26600696        26600908        .       +       0       ID=mikado.Chr5G5.1.CDS6;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     exon    26600696        26600908        .       +       .       ID=mikado.Chr5G5.1.exon6;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     CDS     26600987        26601085        .       +       0       ID=mikado.Chr5G5.1.CDS7;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     exon    26600987        26601137        .       +       .       ID=mikado.Chr5G5.1.exon7;Parent=mikado.Chr5G5.1
    Chr5    Mikado_loci     three_prime_UTR 26601086        26601137        .       +       .       ID=mikado.Chr5G5.1.three_prime_UTR1;Parent=mikado.Chr5G5.1
    ###
    Chr5    Mikado_loci     superlocus      26575364        26579730        .       -       .       ID=Mikado_superlocus:Chr5-:26575364-26579730;Name=superlocus:Chr5-:26575364-26579730
    Chr5    Mikado_loci     ncRNA_gene      26575364        26579730        18      -       .       ID=mikado.Chr5G6;Name=mikado.Chr5G6;multiexonic=True;superlocus=Mikado_superlocus:Chr5-:26575364-26579730
    Chr5    Mikado_loci     ncRNA   26575711        26579730        18      -       .       ID=mikado.Chr5G6.1;Parent=mikado.Chr5G6;Name=cl_Chr5.6271;abundance=1.141582;alias=cl_Chr5.6271;canonical_junctions=1,2,3,4,5,6,7,8,9,10;canonical_number=10;canonical_proportion=1.0;primary=True
    Chr5    Mikado_loci     exon    26575711        26575797        .       -       .       ID=mikado.Chr5G6.1.exon1;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26575885        26575944        .       -       .       ID=mikado.Chr5G6.1.exon2;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26576035        26576134        .       -       .       ID=mikado.Chr5G6.1.exon3;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26576261        26577069        .       -       .       ID=mikado.Chr5G6.1.exon4;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26577163        26577288        .       -       .       ID=mikado.Chr5G6.1.exon5;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26577378        26577449        .       -       .       ID=mikado.Chr5G6.1.exon6;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26577856        26577937        .       -       .       ID=mikado.Chr5G6.1.exon7;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26578239        26578792        .       -       .       ID=mikado.Chr5G6.1.exon8;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26579079        26579161        .       -       .       ID=mikado.Chr5G6.1.exon9;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26579301        26579395        .       -       .       ID=mikado.Chr5G6.1.exon10;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     exon    26579602        26579730        .       -       .       ID=mikado.Chr5G6.1.exon11;Parent=mikado.Chr5G6.1
    Chr5    Mikado_loci     ncRNA   26578496        26579563        13      -       .       ID=mikado.Chr5G6.3;Parent=mikado.Chr5G6;Name=tr_c73_g1_i1.mrna1.160;alias=tr_c73_g1_i1.mrna1.160;canonical_junctions=1;canonical_number=1;canonical_proportion=1.0;ccode=j;gene_name=c73_g1_i1;primary=False
    Chr5    Mikado_loci     exon    26578496        26578518        .       -       .       ID=mikado.Chr5G6.3.exon1;Parent=mikado.Chr5G6.3
    Chr5    Mikado_loci     exon    26579301        26579563        .       -       .       ID=mikado.Chr5G6.3.exon2;Parent=mikado.Chr5G6.3
    Chr5    Mikado_loci     ncRNA   26575364        26578163        16      -       .       ID=mikado.Chr5G6.2;Parent=mikado.Chr5G6;Name=cuff_cufflinks_star_at.23553.1;alias=cuff_cufflinks_star_at.23553.1;fpkm=2.9700103727;canonical_junctions=1,2,3,4,5,6,7,8;canonical_number=8;canonical_proportion=1.0;ccode=j;conf_hi=3.260618;conf_lo=2.679403;cov=81.895309;frac=0.732092;primary=False
    Chr5    Mikado_loci     exon    26575364        26575410        .       -       .       ID=mikado.Chr5G6.2.exon1;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26575495        26575620        .       -       .       ID=mikado.Chr5G6.2.exon2;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26575711        26575797        .       -       .       ID=mikado.Chr5G6.2.exon3;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26575885        26575944        .       -       .       ID=mikado.Chr5G6.2.exon4;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26576035        26576134        .       -       .       ID=mikado.Chr5G6.2.exon5;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26576261        26577069        .       -       .       ID=mikado.Chr5G6.2.exon6;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26577163        26577288        .       -       .       ID=mikado.Chr5G6.2.exon7;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26577378        26577449        .       -       .       ID=mikado.Chr5G6.2.exon8;Parent=mikado.Chr5G6.2
    Chr5    Mikado_loci     exon    26577856        26578163        .       -       .       ID=mikado.Chr5G6.2.exon9;Parent=mikado.Chr5G6.2
    ###

Things to note:
    * multiple RNAs for the same gene are identified by progressive enumeration after a "." (eg. mikado.Chr5G5.1, mikado.Chr5G5.2, etc.).
    * All RNAs retain their old name under the attribute "alias". If a transcript was split due to the presence of multiple ORFs, its alias will end with ".split<progressive ID>".
    * RNAs have the boolean attribute "primary", which identifies them as the primary transcript of the gene or as an alternative splicing isoform.
    * Non-primary RNAs have the additional "ccode" field, which identifies the :ref:`class code <ccodes>` assigned to them when they were compared to the primary transcript.
    * multiexonic RNAs have the attributes "canonical_junctions", "canonical_number", and "canonical_proportion" assigned to them. These properties are calculated by Mikado during the :ref:`prepare stage <prepare>`.

Metrics files
-------------

These are tabular files that enumerate all the :ref:`metrics raw values <Metrics>` for each transcript. This is the section of the metrics file corresponding to the GFF3 file above::

    tid	parent	score	best_bits	blast_score	canonical_intron_proportion	cdna_length	cds_not_maximal	cds_not_maximal_fraction	combined_cds_fraction	combined_cds_intron_fractioncombined_cds_length	combined_cds_num	combined_cds_num_fraction	combined_utr_fraction	combined_utr_length	end_distance_from_junction	end_distance_from_tes	exon_fraction	exon_num	five_utr_length	five_utr_num	five_utr_num_complete	has_start_codon	has_stop_codon	highest_cds_exon_number	highest_cds_exons_num	intron_fraction	is_complete	max_intron_length	min_intron_length	non_verified_introns_num	num_introns_greater_than_max	num_introns_smaller_than_min	number_internal_orfs	proportion_verified_introns	proportion_verified_introns_inlocusretained_fraction	retained_intron_num	selected_cds_exons_fraction	selected_cds_fraction	selected_cds_intron_fraction	selected_cds_length	selected_cds_num	selected_cds_number_fraction	selected_end_distance_from_junction	selected_end_distance_from_tes	selected_start_distance_from_tss	snowy_blast_score	source_score	start_distance_from_tss	three_utr_length	three_utr_num	three_utr_num_complete	utr_fraction	utr_length	utr_num	utr_num_complete	verified_introns_num
    mikado.Chr5G1.2	mikado.Chr5G1	19.0	1086.25	1086.25	1.0	1927	0	0.0	0.87	1.0	1683	10	0.91	0.13	244	0	157	0.92	11	87	2	1	TrueTrue	10	10	0.91	True	340	71	10	0	0	1	0.0	0	0.0	0	0.91	0.87	1.0	1683	10	0.91	0	157	87	13.78	0	87	157	1	0	0.13	244	3	1	0
    mikado.Chr5G1.1	mikado.Chr5G1	21.89	1086.63	1086.63	1.0	1937	0	0.0	0.87	1.0	1683	10	0.91	0.13	254	0	157	0.92	11	97	2	1	TrueTrue	10	10	0.91	True	196	71	10	0	0	1	0.0	0	0.0	0	0.91	0.87	1.0	1683	10	0.91	0	157	97	13.78	0	97	157	1	0	0.13	254	3	1	0
    mikado.Chr5G2.2	mikado.Chr5G2	19.04	1140.95	1140.95	1.0	2197	0	0.0	0.88	1.0	1938	11	0.92	0.12	259	0	0	0.92	12	259	2	1	TrueTrue	11	11	0.92	True	577	74	11	0	0	1	0.0	0	0.0	0	0.92	0.88	1.0	1938	11	0.92	0	0	259	16.66	0	259	0	0	0	0.12	259	2	1	0
    mikado.Chr5G2.1	mikado.Chr5G2	20.06	1140.95	1140.95	1.0	2204	0	0.0	0.88	1.0	1938	11	0.92	0.12	266	0	0	0.92	12	266	2	1	TrueTrue	11	11	0.92	True	570	74	11	0	0	1	0.0	0	0.0	0	0.92	0.88	1.0	1938	11	0.92	0	0	266	16.66	0	266	0	0	0	0.12	266	2	1	0
    mikado.Chr5G3.2	mikado.Chr5G3	8.59	1193.72	1193.72	1.0	1887	0	0.0	0.96	0.8	1816	9	1.0	0.04	71	0	0	0.75	9	71	1	0	TrueFalse	9	9	0.8	False	152	74	8	0	0	1	0.0	0	0.0	0	1.0	0.96	0.8	1816	9	1.0	0	0	71	14.16	0	71	0	0	0	0.04	71	1	0	0
    mikado.Chr5G3.1	mikado.Chr5G3	19.0	1353.19	1353.19	1.0	2109	0	0.0	1.0	0.9	2109	10	1.0	0.0	0	0	0	0.83	10	0	0	0	TrueTrue	10	10	0.9	True	152	74	9	0	0	1	0.0	0	0.0	0	1.0	1.0	0.9	2109	10	1.0	0	0	0	16.66	0	00	0	0	0.0	0	0	0	0
    mikado.Chr5G4.1	mikado.Chr5G4	20.0	1258.43	1258.43	1.0	2025	0	0.0	0.98	0	1986	1	1.0	0.02	39	0	39	1.0	1	0	0	0	TrueTrue	1	1	0	True	0	0	0	0	0	1	0	0	0.0	0	1.0	0.98	0	1986	1	1.0	0	39	0	16.66	0	039	1	0	0.02	39	1	0	0
    mikado.Chr5G5.1	mikado.Chr5G5	20.0	565.46	565.46	1.0	1184	0	0.0	0.79	1.0	936	7	1.0	0.21	248	0	52	1.0	7	196	1	0	TrueTrue	7	7	1.0	True	112	69	6	0	0	1	0.0	0	0.0	0	1.0	0.79	1.0	936	7	1.0	0	52	196	13.67	0	196	52	1	0	0.21	248	2	0	0
    mikado.Chr5G6.2	mikado.Chr5G6	17.1	0	0	1.0	1735	0	0	0.0	0	0	0	0.0	1.0	0	0	0	0.56	9	0	0	0	False	False	0	0	0.62	False	406	84	0	0	0	0	1.0	0.89	0.0	0	0.0	0.0	0	0	0	0.0	0	0	0	0	00	0	0	0	1.0	0	0	0	8
    mikado.Chr5G6.1	mikado.Chr5G6	17.9	0	0	1.0	2197	0	0	0.0	0	0	0	0.0	1.0	0	0	0	0.69	11	0	0	0	False	False	0	0	0.77	False	406	87	3	0	0	0	0.9	1.0	0.0	0	0.0	0.0	0	0	0	0.0	0	0	0	0	00	0	0	0	1.0	0	0	0	9
    mikado.Chr5G6.3	mikado.Chr5G6	13.0	0	0	1.0	286	0	0	0.0	0	0	0	0.0	1.0	0	0	0	0.12	2	0	0	0	False	False	0	0	0.08	False	782	782	1	0	0	0	0.0	0.0	0.0	0	0.0	0.0	0	0	0	0.0	0	0	0	0	00	0	0	0	1.0	0	0	0	0

As it can be noted, metrics can assume values in a very wide range. We direct you to the :ref:`metrics section of the documentation <Metrics>` for further details.

Scoring files
-------------

This file contains the scores assigned to each metric for each transcript. Only metrics which have been used for the scoring will be present. This is the section of the metrics file corresponding to the above GFF3 file::

    tid	parent	score	blast_score	cdna_length	cds_not_maximal	cds_not_maximal_fraction	combined_cds_fraction	combined_cds_intron_fraction	combined_cds_length	combined_cds_num	end_distance_from_junction	exon_fraction	exon_num	five_utr_length	five_utr_num	highest_cds_exon_number	intron_fraction	number_internal_orfs	proportion_verified_introns	retained_fraction	retained_intron_num	selected_cds_fraction	selected_cds_intron_fraction	selected_cds_length	selected_cds_num	source_score	three_utr_length	three_utr_num
    mikado.Chr5G1.2	mikado.Chr5G1	19.0	0.0	0.0	1	1	0.0	1	1	1	1	1	1	0.0	1.0	1	1	1.0	1	1	1	0.0	1	11	0	0.0	1.0
    mikado.Chr5G1.1	mikado.Chr5G1	21.89	1.0	1.0	1	1	0.06	1	1	1	1	1	1	0.77	1.0	1	1	1.0	1	1	1	0.06	1	11	0	0.0	1.0
    mikado.Chr5G2.2	mikado.Chr5G2	19.04	1	0.0	1	1	0.0	1	1	1	1	1	1	0.04	1.0	1	1	1.0	1	1	1	0.0	1	11	0	0.0	0.0
    mikado.Chr5G2.1	mikado.Chr5G2	20.06	1	1.0	1	1	0.03	1	1	1	1	1	1	0.0	1.0	1	1	1.0	1	1	1	0.03	1	11	0	0.0	0.0
    mikado.Chr5G3.1	mikado.Chr5G3	19.0	1.0	1.0	1	1	0.0	1.0	1.0	1.0	1	1.0	1.0	0.0	0.0	1.0	1.0	1.0	1	1	1	0.0	1.01.0	1.0	0	0.0	0.0
    mikado.Chr5G3.2	mikado.Chr5G3	8.59	0.0	0.0	1	1	0.19	0.0	0.0	0.0	1	0.0	0.0	0.71	0.5	0.0	0.0	1.0	1	1	1	0.19	0.00.0	0.0	0	0.0	0.0
    mikado.Chr5G4.1	mikado.Chr5G4	20.0	1	1	1	1	0.0	1	1	1	1	1	1	0.0	0.0	1	1	1.0	1	1	1	0.0	1	11	0	0.0	1.0
    mikado.Chr5G5.1	mikado.Chr5G5	20.0	1	1	1	1	0.0	1	1	1	1	1	1	0.0	0.0	1	1	1.0	1	1	1	0.0	1	11	0	0.0	1.0
    mikado.Chr5G6.3	mikado.Chr5G6	13.0	1	0.0	1	1	0.0	1	1	1	1	0.0	0.0	0.0	0.0	1	0.0	0.0	0.0	1	1	0.0	1	11	0	0.0	0.0
    mikado.Chr5G6.2	mikado.Chr5G6	17.1	1	0.76	1	1	0.0	1	1	1	1	0.78	0.78	0.0	0.0	1	0.78	0.0	1.0	1	1	0.0	1	11	0	0.0	0.0
    mikado.Chr5G6.1	mikado.Chr5G6	17.9	1	1.0	1	1	0.0	1	1	1	1	1.0	1.0	0.0	0.0	1	1.0	0.0	0.9	1	1	0.0	1	11	0	0.0	0.0


The final score value is obtained by summing all the individual metrics.

.. important:: If you compare the scores assigned to transcripts at the loci level with those assigned at the subloci level, you will notice that the scores are different and that even some of the raw metrics values are. The former phenomenon is due to the fact that :ref:`the Mikado scoring system is not absolute but relative <scoring_algorithm>`; the latter, to the fact that :ref:`some metrics are locus-dependent <Metrics>`, ie their values change due the presence or absence of other transcripts. A typical example is given by the "retained_intron" metrics; retained introns are identified by looking for non-coding regions of transcript which fall inside the intron of another transcript. Changing the transcripts in the locus will change the value associated to this metric, as non-coding sections will or will not be classified as "retained introns", and therefore the score associated with both the metric and the transcript.


Transcript padding
__________________

After calculating the final loci, Mikado can try to uniform the ends of transcripts present in the locus, by extending
the shorter ones so that their ends coincide with those of longer transcripts in the locus. The procedure is explained more
in detail in the :ref:`dedicated section in the Algorithms page <padding>`. The approach has been inspired by the consolidation
approach taken by the Araport annotation for *Arabidopsis thaliana* [AraPort]_.

Usage
~~~~~

``mikado pick`` allows to modify some of the parameters regarding the run at runtime. However, some sections - such as most of the settings regarding alternative splicing - are left untouched by the utility, and are best modified by editing the :ref:`configuration file itself <configure>`. The available parameters are as follows:

* *json-conf*: required. This is the configuration file created in the :ref:`first step <configure>` of the pipeline.
* *gff*; optionally, it is possible to point Mikado prepare to the GTF it should use here on the command line. This file should be the output of the :ref:`preparation step <prepare>`. Please note that this file should be in GTF format, sorted by chromosome and position; if that is not the case, Mikado will fail.
* *db*: Optionally, it is possible to specify the database to Mikado on the command line, rather than on the configuration file. Currently, this option *supports SQLite databases only*.
* Options related to how Mikado will treat the data:

    * *intron_range*: this option expects a couple of positive integers, in ascending order, indicating the 98% CI where most intron lengths should fall into. Gene models with introns whose lengths fall outside of this range might be penalized, depending on the scoring system used. If uncertain, it is possible to use the :ref:`included stats utility <stat-command>` on the gene annotation of a closely related species.
    * *no-purge*: flag. If set, Mikado will not not exclude putative fragments from the output, but will report them (appropriately flagged).
    * *flank*: for the purposes of identifying fragments, it is useful to consider together loci which are not necessarily overlapping but which are lying relatively near on the genome sequence. This parameter (a positive integer) specifies the maximum distance for Mikado for gathering data together for this purpose.
    * *mode*: how Mikado will treat BLAST and ORF data in the presence of putative chimeras. Please refer to the :ref:`algorithms section <chimera_splitting_algorithm>` for details.
* Options regarding the output files:

    * *output-dir*: Output directory. By default, Mikado will write all files and the log on the current directory.
    * *loci_out*: required. This it the main output file, in GFF format.
    * *prefix*: this ID will be prefixed to all gene and transcript models. IN general, IDs will be of the form "<prefix>.<chromosome><progressive ID>". Default: Mikado.
    * *source*: source field prefix for the output files. Useful for eg loading Mikado runs into `WebApollo <http://genomearchitect.org/>`_ [Apollo]_.
    * *no_cds*: if present, this flg will indicate to Mikado not to print out the CDS of selected models but only their transcript structures.
    * *subloci_out*: If requested, Mikado can output the data regarding the first intermediate step, ie the *subloci*. See the :ref:`introduction <Introduction>` for details.
    * *monoloci_out*: If requested, Mikado can output the data regarding the second intermediate step, ie the *monosubloci*. See the :ref:`introduction <Introduction>` for details.
* Options regarding the resources to be used:

    * *procs*: number of processors to use.
    * *start-method*: multiprocessing start method. See the :ref:`explanation on Python multiprocessing <scheduler-multiprocessing>`
    * *single*: flag. If present, multiprocessing will be disabled.
* Options regarding logging:

    * *log*: name of the log file. By default, "pick.log"
    * *verbose*: sets the log level to DEBUG. Please be advised that the debug mode is **extremely** verbose and is bestly invoked only for real, targeted debugging sessions.
    * *noverbose*: sets the log level to ERROR. If set, in most cases, the log file will be practically empty.
    * *log-level*: this flag directly sets the log level. Available values: DEBUG, INFO, WARNING, ERROR.
* Options related to padding:

    * *pad*: if set, this option will enforce transcript padding. The default is inferred from the configuration (on by default).
    * *no-pad*: if set, this option will disable transcript padding. The default is inferred from the configuration (on by default).
    * *pad-max-splices*: maximum amount of splicing sites that an expanded exon can cross. Default is inferred from the configuration file (currently default is 1)
    * *pad-max-distance*: Maximum amount of basepairs that transcripts can be padded with (per side). Default is inferred from the configuration file (default 300 bps)
    * *fasta*: genome FASTA file. **Required if the padding is switched on**. Default: inferred from the configuration file.

Usage::

    $ mikado pick --help
    usage: Mikado pick [-h] [--start-method {fork,spawn,forkserver}] [-p PROCS]
                       --json-conf JSON_CONF [--scoring-file SCORING_FILE]
                       [-i INTRON_RANGE INTRON_RANGE] [--fasta FASTA]
                       [--no-pad | --pad] [--pad-max-splices PAD_MAX_SPLICES]
                       [--pad-max-distance PAD_MAX_DISTANCE]
                       [--subloci_out SUBLOCI_OUT] [--monoloci_out MONOLOCI_OUT]
                       [--loci_out LOCI_OUT] [--prefix PREFIX] [--source SOURCE]
                       [--no_cds] [--flank FLANK] [--no-purge] [--cds-only]
                       [--only-reference-update]
                       [--consider-truncated-for-retained] [-db SQLITE_DB]
                       [-od OUTPUT_DIR] [--single] [-l LOG] [-v | -nv]
                       [-lv {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                       [--mode {nosplit,stringent,lenient,permissive,split}]
                       [gff]

    positional arguments:
      gff

    optional arguments:
      -h, --help            show this help message and exit
      --start-method {fork,spawn,forkserver}
                            Multiprocessing start method. (default: None)
      -p PROCS, --procs PROCS
                            Number of processors to use. Default: look in the
                            configuration file (1 if undefined) (default: None)
      --json-conf JSON_CONF
                            JSON/YAML configuration file for Mikado. (default:
                            None)
      --scoring-file SCORING_FILE
                            Optional scoring file for the run. It will override
                            (default: None)
      -i INTRON_RANGE INTRON_RANGE, --intron-range INTRON_RANGE INTRON_RANGE
                            Range into which intron lengths should fall, as a
                            couple of integers. Transcripts with intron lengths
                            outside of this range will be penalised. Default: (60,
                            900) (default: None)
      --fasta FASTA         Genome FASTA file. Required if pad is enabled
                            (default). (default: None)
      --no-pad              Disable transcript padding. (default: None)
      --pad                 Whether to pad transcripts in loci. (default: None)
      --pad-max-splices PAD_MAX_SPLICES
                            Maximum splice sites that can be crossed during
                            transcript padding. (default: None)
      --pad-max-distance PAD_MAX_DISTANCE
                            Maximum amount of bps that transcripts can be padded
                            with (per side). (default: None)
      --no_cds              Flag. If set, not CDS information will be printed out
                            in the GFF output files. (default: False)
      --flank FLANK         Flanking distance (in bps) to group non-overlapping
                            transcripts into a single superlocus. Default:
                            determined by the configuration file. (default: None)
      --no-purge            Flag. If set, the pipeline will NOT suppress any loci
                            whose transcripts do not pass the requirements set in
                            the JSON file. (default: False)
      --cds-only            "Flag. If set, Mikado will only look for overlap in
                            the coding features when clustering transcripts
                            (unless one transcript is non-coding, in which case
                            the whole transcript will be considered). Default:
                            False, Mikado will consider transcripts in their
                            entirety. (default: False)
      --only-reference-update
                            Flag. If switched on, Mikado will only keep loci where
                            at least one of the transcripts is marked as
                            "reference". CAUTION: new and experimental. If no
                            transcript has been marked as reference, the output
                            will be completely empty! (default: None)
      --consider-truncated-for-retained
                            Flag. If set, Mikado will consider as retained intron
                            events also transcripts which lack UTR but whose CDS
                            ends within a CDS intron of another model. (default:
                            False)
      -db SQLITE_DB, --sqlite-db SQLITE_DB
                            Location of an SQLite database to overwrite what is
                            specified in the configuration file. (default: None)
      -od OUTPUT_DIR, --output-dir OUTPUT_DIR
                            Output directory. Default: current working directory
                            (default: None)
      --single              Flag. If set, Creator will be launched with a single
                            process. Useful for debugging purposes only. (default:
                            False)
      --mode {nosplit,stringent,lenient,permissive,split}
                            Mode in which Mikado will treat transcripts with
                            multiple ORFs. - nosplit: keep the transcripts whole.
                            - stringent: split multi-orf transcripts if two
                            consecutive ORFs have both BLAST hits and none of
                            those hits is against the same target. - lenient:
                            split multi-orf transcripts as in stringent, and
                            additionally, also when either of the ORFs lacks a
                            BLAST hit (but not both). - permissive: like lenient,
                            but also split when both ORFs lack BLAST hits - split:
                            split multi-orf transcripts regardless of what BLAST
                            data is available. (default: None)

    Options related to the output files.:
      --subloci_out SUBLOCI_OUT
      --monoloci_out MONOLOCI_OUT
      --loci_out LOCI_OUT   This output file is mandatory. If it is not specified
                            in the configuration file, it must be provided here.
                            (default: None)
      --prefix PREFIX       Prefix for the genes. Default: Mikado (default: None)
      --source SOURCE       Source field to use for the output files. (default:
                            None)

    Log options:
      -l LOG, --log LOG     File to write the log to. Default: decided by the
                            configuration file. (default: None)
      -v, --verbose         Flag. If set, the debug mode will be activated.
                            (default: False)
      -nv, --noverbose      Flag. If set, the log will report only errors and
                            critical events. (default: False)
      -lv {DEBUG,INFO,WARNING,ERROR,CRITICAL}, --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                            Logging level. Default: retrieved by the configuration
                            file. (default: None)



.. block end


Technical details
~~~~~~~~~~~~~~~~~

``mikado pick`` uses a divide-et-impera algorithm to find and analyse loci separately. As the data to be integrated with the transcripts is stored on the database rather than be calculated on the fly, rerunning ``pick`` with different options takes little time and resources.
To keep the data sorted, Mikado will write out temporary files during the operation and merge them at the end of the run (see function ``merge_loci_gff`` in :ref:`the picking module <sub-picking-loci>`.