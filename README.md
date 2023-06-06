[![GitHub Downloads](https://img.shields.io/github/downloads/EI-CoreBioinformatics/mikado/total.svg?style=social&logo=github&label=download)](https://github.com/EI-CoreBioinformatics/mikado/releases)
[![Release](https://img.shields.io/github/release/EI-CoreBioinformatics/mikado.svg)](https://github.com/EI-CoreBioinformatics/mikado/releases)
[![PyPI](https://img.shields.io/pypi/v/mikado.svg?style=flat)](https://pypi.python.org/pypi/mikado)
[![Build Status](https://github.com/EI-CoreBioinformatics/mikado/workflows/Mikado/badge.svg)](https://github.com/EI-CoreBioinformatics/mikado/actions/workflows/python-package.yml)

# Mikado - pick your transcript: a pipeline to determine and select the best RNA-Seq prediction

Mikado is a lightweight Python3 pipeline to identify the most useful or “best” set of transcripts from multiple transcript assemblies. Our approach leverages transcript assemblies generated by multiple methods to define expressed loci, assign a representative transcript and return a set of gene models that selects against transcripts that are chimeric, fragmented or with short or disrupted CDS. Loci are first defined based on overlap criteria and each transcript therein is scored based on up to 50 available metrics relating to ORF and cDNA size, relative position of the ORF within the transcript,  UTR length and presence of multiple ORFs. Mikado can also utilize blast data to score transcripts based on proteins similarity and to identify and split chimeric transcripts. Optionally, junction confidence data as provided by [Portcullis][Portcullis] can be used to improve the assessment. The best-scoring transcripts are selected as the primary transcripts of their respective gene loci; additionally, Mikado can bring back other valid splice variants that are compatible with the primary isoform.

Mikado uses GTF or GFF files as mandatory input. Non-mandatory but highly recommended input data can be generated by obtaining a set of reliable splicing junctions with Portcullis_, by locating coding ORFs on the transcripts using either [Transdecoder][Transdecoder] or [Prodigal][Prodigal], and by obtaining homology information through either [BLASTX][Blast+] or [DIAMOND][DIAMOND].

Our approach is amenable to include sequences generated by *de novo* Illumina assemblers or reads generated from long read technologies such as Pacbio.

Extended documentation is hosted on ReadTheDocs: http://mikado.readthedocs.org/

## Installation

Using mamba

download mamba using pip

```bash
pip install mamba=0.27.0
```

Create a mamba environment using the environment.yml file

```bash
mamba env create -f environment.yml
conda activate mikado2
```

Check and run mikado

```bash
mikado --help
```



Mikado can also be be installed from PyPI with pip (**deprecated**):

``pip3 install mikado``

Alternatively, you can clone the repository from source and install with:

    pip wheel -w dist .
    pip install dist/*whl

You can verify the correctness of the installation with the unit tests (*outside of the source folder*, as otherwise Python will get confused and try to use the `Mikado` source folder instead of the system installation):

    python -c "import Mikado; Mikado.test(); Mikado.test(label='slow')"

An alternative way of installing using `setuptools`:

    pip install -r requirements.txt
    pip install Cython
    python setup.py bdist_wheel
    pip install dist/*whl

The steps above will ensure that any additional python dependencies will be installed correctly. A full list of library dependencies can be found in the file ``requirements.txt``

### Additional dependencies

Mikado by itself does require only the presence of a database solution, such as SQLite (although we do support MySQL and PostGRESQL as well).
However, the Daijin pipeline requires additional programs to run.

For driving Mikado through Daijin, the following programs are required:

- [DIAMOND][DIAMOND] or [Blast+][Blast+] to provide protein homology. DIAMOND is preferred for its speed.
- [Prodigal][Prodigal] or [Transdecoder][Transdecoder] to calculate ORFs. The versions of Transdecoder that we tested scale poorly in terms of runtime and disk usage, depending on the size of the input dataset. Prodigal is much faster and lighter, however, the data on our paper has been generated through Transdecoder - not Prodigal. Currently we set Prodigal as default.
- Mikado also makes use of a dataset of RNA-Seq high-quality junctions. We are using [Portcullis][Portcullis] to calculate this data alongside the alignments and assemblies.

If you plan to generate the alignment and assembly part as well through Daijin, the pipeline requires the following:

- SAMTools
- If you have short-read RNA-Seq data:
  - At least one short-read RNA-Seq aligner, choice between [GSNAP], [GMAP][GMAP], [STAR][STAR], [TopHat2][TopHat2], [HISAT2][HISAT2]
  - At least one RNA-Seq assembler, choice between [StringTie][StringTie], [Trinity][Trinity], [Cufflinks], [CLASS2][CLASS2]. Trinity additionally requires [GMAP][GMAP].
  - [Portcullis][Portcullis] is optional, but highly recommended to retrieve high-quality junctions from the data
- If you have long-read RNA-Seq data:
  - At least one long-read RNA-Seq aligner, current choice between [STAR][STAR] and [GMAP][GMAP]

## Development guide

We provide source trail files ([https://www.sourcetrail.com/](https://www.sourcetrail.com/)) to aid in development.
As required by the SourceTrail application, these files are present in the master directory, as "Mikado.srctrl*".

## Citing Mikado

If you use Mikado in your work, please consider to cite:

> Venturini L., Caim S., Kaithakottil G., Mapleson D.L., Swarbreck D. Leveraging multiple transcriptome assembly methods for improved gene structure annotation. GigaScience, Volume 7, Issue 8, 1 August 2018, giy093, [doi:10.1093/gigascience/giy093](https://doi.org/10.1093/gigascience/giy093)

If you also use Portcullis to provide reliable junctions to Mikado, either independently or as part of the Daijin pipeline, please consider to cite:

> Mapleson D.L., Venturini L., Kaithakottil G., Swarbreck D. Efficient and accurate detection of splice junctions from RNAseq with Portcullis. GigaScience, Volume 7, Issue 12, 12 December 2018, giy131, [doi:10.1093/gigascience/giy131](https://doi.org/10.1093/gigascience/giy131)

[Portcullis]: https://github.com/maplesond/portcullis
[Blast+]: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
[Transdecoder]: https://github.com/TransDecoder/TransDecoder/
[Prodigal]: https://github.com/hyattpd/Prodigal
[DIAMOND]: https://github.com/bbuchfink/diamond/
[GMAP]: http://research-pub.gene.com/gmap/
[STAR]: https://github.com/alexdobin/STAR
[TopHat2]: https://ccb.jhu.edu/software/tophat/index.shtml
[HISAT2]: http://ccb.jhu.edu/software/hisat2
[StringTie]: https://ccb.jhu.edu/software/stringtie/
[Trinity]: https://github.com/trinityrnaseq/trinityrnaseq
[CLASS2]: http://ccb.jhu.edu/people/florea/research/CLASS2/
