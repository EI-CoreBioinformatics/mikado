"""Setup file for PyPI"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    long_description = description.read()

setup(

    name = "shanghai",
    version = "0.2",
    
    description="A Python3 annotation program to select the best gene model in each locus",
    long_description=long_description,

    url="http://stash.tgac.ac.uk/users/venturil/repos/locus_pipeline/",

    author="Luca Venturini",
    author_email="luca.venturini@tgac.ac.uk",

    license="GPL3",

    classifiers = [

        "Development Status :: 3 - Alpha",
        "Topic :: Gene Annotation ",
        "License :: OSI Approved :: GPL3",
        'Programming Language :: Python :: 3.4',
    ],

    keywords="rna-seq annotation genomics transcriptomics",

    packages = find_packages(),

    install_requires=["pyyaml", "networkx", "sqlalchemy", "biopython" ],

    data_files = [ ("sample_data", ["sample_data/*"])],

)
