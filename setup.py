# coding: utf-8

"""Setup file for PyPI"""

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from codecs import open
from os import path
import glob
import re
import sys

here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    long_description = description.read()

version = None
with open(path.join(here, "Mikado", "__init__.py")) as main:
    for line in main:
        if "__version__" in line:
            version = re.sub("'", "", re.sub('"', "", line.rstrip().split()[-1]))
assert version is not None

if sys.version_info.major != 3:
    raise EnvironmentError("""Mikado is a pipeline specifically programmed for python3,
    and is not compatible with Python2. Please upgrade your python before proceeding!""")

setup(
    name="Mikado",
    version=version,
    description="A Python3 annotation program to select the best gene model in each locus",
    long_description=long_description,
    url="https://github.com/lucventurini/mikado",
    author="Luca Venturini",
    author_email="luca.venturini@tgac.ac.uk",
    license="GPL3",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: POSIX :: Linux",
        'Programming Language :: Python :: 3.4',
        "Programming Language :: Python :: 3.5",
    ],
    ext_modules=cythonize([Extension(path.join("Mikado.utilities.overlap"),
                                     [path.join("Mikado", "utilities", "overlap.pyx")]),
                           Extension(path.join("Mikado.scales.f1"),
                                     [path.join("Mikado", "scales", "f1.pyx")]),
                           Extension(path.join("Mikado.scales.contrast"),
                                     [path.join("Mikado", "scales", "contrast.pyx")]),
                           Extension(path.join("Mikado.scales.intervaltree"),
                                     [path.join("Mikado", "scales", "intervaltree.pyx")]),
                           ]),
    zip_safe=False,
    keywords="rna-seq annotation genomics transcriptomics",
    packages=find_packages(),
    # scripts=glob.glob("bin/*.py") + glob.glob("util/*.py"),
    scripts=glob.glob("util/*.py"),
    entry_points={"console_scripts": ["mikado = Mikado:main",
                                      "daijin = Mikado.daijin:main"]},
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    extras_require={
        "postgresql": ["psycopg2"],
        "mysql": ["mysqlclient>=1.3.6"],
    },
    test_suite="Mikado.test",
    package_data={
        "Mikado.configuration":
            glob.glob("Mikado/configuration/*json") + glob.glob("Mikado/configuration/*yaml"),
        "Mikado.configuration.scoring_files":
            glob.glob("Mikado/configuration/scoring_files/*"),
        "Mikado":
            glob.glob(path.join("Mikado", "daijin", "*yaml")) +
            glob.glob("Mikado/daijin/*json") + \
            glob.glob("Mikado/daijin/*snakefile")
        },
    include_package_data=True
    # data_files=[
    #     ("Mikado/configuration",
    #      glob.glob("Mikado/configuration/*json") + glob.glob("Mikado/configuration/*yaml"))],

)
