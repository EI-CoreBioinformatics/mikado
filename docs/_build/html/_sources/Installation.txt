:: Installation_

Installation
============


Download
--------

Mikado is available on PyPI, so it is possible to install it with

``pip3 install mikado``

The source for the latest release can be obtained with

``pip3 download mikado``

As the package contains some core Cython components, it might be necessary to download and compile the source code instead of relying on the wheel.

Alternatively, Mikado can be installed from source by obtaining it from our `GitHub`_ repository. Either download the `latest release <https://github.com/lucventurini/mikado/releases/latest>`_  or download the latest stable development snapshot with

``git clone https://github.com/lucventurini/mikado.git; cd mikado``

If you desire, the unstable development version can be obtained with the command

``git checkout development``

in the Git directory. Please note that development can proceed quite rapidly.

Building and installing from source
-----------------------------------

If you desire to install Mikado from source, this can be achieved with

``python3 setup.py bdist_wheel``

Followed by

``pip3 install dist/*whl``

We advise to test whether the distribution has been built successfully by executing the unit test suite with

``python3 setup.py test``

Although code coverage is not perfect yet, it is over 50% for the whole package and considerably higher for the core components.

Python Dependencies
-------------------

Mikado has been written for Python 3.4 and 3.5. It is dependent on the following Python3 modules:

* wheel>=0.28.0
* pyyaml [PyYaml]_
* jsonschema
* Cython [Cython]_
* cython [Cython]_
* numpy [Numpy]_
* networkx>=1.10 [NetworkX]_
* sqlalchemy>=1
* sqlalchemy_utils
* biopython>=1.66 [BioPython]_
* intervaltree
* nose
* pyfaidx
* scikit-learn>=0.17.0 [SciKit]_
* scipy>=0.15.0 [Scipy]_
* frozendict
* python-magic
* drmaa [DRMAA]_
* snakemake [Snake]_
* docutils

These dependencies will be installed automatically by PIP.

.. _GitHub: https://github.com/lucventurini/mikado

Additional dependencies
-----------------------

Mikado relies on relational databases for its functioning, so one of SQLite, PosGRESql or MySQL have to present for if to function properly.
Additionally, the Daijin pipeline, which drives Mikado, requires BLAST+ and TransDecoder for the Mikado stage, and at least one RNA-Seq aligner and one assembler, to be installed. If you are planning to execute it on a cluster, we do support job submission either with DRMAA or without on SLURM, LSF and PBS clusters.