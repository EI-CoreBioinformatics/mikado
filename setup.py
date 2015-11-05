# coding: utf-8

"""Setup file for PyPI"""

from setuptools import setup, find_packages
from codecs import open
from os import path
import glob

here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    long_description = description.read()

setup(

    name="Mikado",
    version="0.9.5",

    description="A Python3 annotation program to select the best gene model in each locus",
    long_description=long_description,

    url="https://github.com/lucventurini/mikado_lib.git",

    author="Luca Venturini",
    author_email="luca.venturini@tgac.ac.uk",

    license="GPL3",

    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Gene Annotation",
        "License :: OSI Approved :: GPL3",
        'Programming Language :: Python :: 3.4',
        "Operating System :: Linux"
    ],

    zip_safe=False,

    keywords="rna-seq annotation genomics transcriptomics",

    packages=find_packages(),

    scripts=glob.glob("bin/*.py"),

    install_requires=["pyyaml",
                      "jsonschema",
                      "numpy",
                      "networkx",
                      "sqlalchemy>=1",
                      "sqlalchemy_utils",
                      "biopython>=1.6",
                      "intervaltree"
                      ],

    extras_require={
        "postgresql": ["psycopg2"],
        "mysql": ["mysqlclient>=1.3.6"],
    },

    data_files=[("mikado_lib/configuration",
                 glob.glob("mikado_lib/configuration/*json") + glob.glob("mikado_lib/configuration/*yaml") )],

)
