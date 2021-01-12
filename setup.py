# coding: utf-8

"""Setup file for PyPI"""

from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from Cython.Build import cythonize
from codecs import open
from os import path
import glob
import re
import sys
import numpy as np
from scipy._build_utils import numpy_nodepr_api


here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    long_description = description.read()

version = {}
with open(path.join(here, "Mikado", "version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]

if version is None:
    print("No version found, exiting", file=sys.stderr)
    sys.exit(1)

if sys.version_info.major != 3:
    raise EnvironmentError("""Mikado is a pipeline specifically programmed for python3,
    and is not compatible with Python2. Please upgrade your python before proceeding!""")

extensions = [Extension("Mikado.utilities.overlap",
                        sources=[path.join("Mikado", "utilities", "overlap.pyx")],
                        **numpy_nodepr_api),
              Extension("Mikado.utilities.f1",
                        sources=[path.join("Mikado", "utilities", "f1.pyx")],
                        **numpy_nodepr_api),
              Extension("Mikado.scales.contrast",
                        sources=[path.join("Mikado", "scales", "contrast.pyx")],
                        **numpy_nodepr_api),
              Extension("Mikado.utilities.intervaltree",
                        sources=[path.join("Mikado", "utilities", "intervaltree.pyx")],
                        **numpy_nodepr_api),
              Extension("Mikado.serializers.blast_serializer.btop_parser",
                        include_dirs=[np.get_include()],
                        language="c++",
                        sources=[path.join("Mikado", "serializers", "blast_serializer", "btop_parser.pyx")],
                        **numpy_nodepr_api),
              Extension("Mikado.serializers.blast_serializer.aln_string_parser",
                        include_dirs=[np.get_include()],
                        language="c++",
                        sources=[path.join("Mikado", "serializers", "blast_serializer", "aln_string_parser.pyx")],
                        **numpy_nodepr_api)
              ]

setup(
    name="Mikado",
    version=version,
    description="A Python3 annotation program to select the best gene model in each locus",
    long_description=long_description,
    url="https://github.com/EI-CoreBioinformatics/mikado",
    author="Luca Venturini",
    author_email="lucventurini@gmail.com",
    license="LGPL3",
    tests_require=["pytest"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: POSIX :: Linux",
        "Framework :: Pytest",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        'Programming Language :: Python :: 3.7'
    ],
    ext_modules=cythonize(extensions, compiler_directives = {"language_level": "3"}),
    zip_safe=False,
    keywords="rna-seq annotation genomics transcriptomics",
    packages=find_packages(),
    scripts=glob.glob("util/*.py"),
    entry_points={"console_scripts": ["mikado = Mikado.__main__:main",
                                      "daijin = Mikado.daijin:main",
                                      ]},
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    extras_require={
        "postgresql": ["psycopg2"],
        "mysql": ["mysqlclient>=1.3.6"],
        "bam": ["pysam>=0.8"]
    },
    # test_suite="nose2.collector.collector",
    package_data={
        "Mikado.configuration":
            glob.glob("Mikado/configuration/*json") + glob.glob("Mikado/configuration/*yaml"),
        "Mikado.configuration.scoring_files":
            glob.glob("Mikado/configuration/scoring_files/*"),
        "Mikado":
            glob.glob(path.join("Mikado", "daijin", "*smk")) +
            glob.glob(path.join("Mikado", "daijin", "*yaml")) + glob.glob("Mikado/daijin/*json"),
        "Mikado.utilities.overlap": [path.join("Mikado", "utilities", "overlap.pxd")],
        "Mikado.utilities.intervaltree": [path.join("Mikado", "utilities", "intervaltree.pxd")],
        "Mikado.serializers.blast_serializers": glob.glob(path.join("Mikado", "serializers", "blast_serializers",
                                                                    "*pxd")),
        "Mikado.tests.blast_data": glob.glob(path.join("Mikado", "tests", "blast_data", "*"))
        },
    include_package_data=True
)
