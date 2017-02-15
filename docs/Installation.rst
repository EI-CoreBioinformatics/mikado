.. _Installation:

Installation
============

System requirements
-------------------

Mikado requires CPython 3.4 or later to function (Python2 is not supported). Additionally, it requires a functioning installation of one among SQLite, PostgreSQL and MySQL. Mikado has additionally the following python dependencies:

.. literalinclude:: ../requirements.txt

Mikado can run with little system requirements, being capable of analysing human data with less than 4GB of RAM in all stages. It benefits from the availability of multiple processors, as many of the steps of the pipeline are parallelised.

Mikado is at its core a data integrator. The :ref:`Daijin pipeline <Daijin>` has been written to facilitate the creation of a functioning workflow. To function, it requires Snakemake [Snake]_ and the presence of at least one supported RNA-Seq aligner and one supported transcript assembler. If you are planning to execute it on a cluster, we do support job submission on SLURM, LSF and PBS clusters, either in the presence or absence of DRMAA.


Download
--------

Mikado is available on PyPI, so it is possible to install it with

``pip3 install mikado``

The source for the latest release on PyPI can be obtained with

``pip3 download mikado``

As the package contains some core Cython components, it might be necessary to download and compile the source code instead of relying on the wheel.

Alternatively, Mikado can be installed from source by obtaining it from our `GitHub`_ repository. Either download the `latest release <https://github.com/lucventurini/mikado/releases/latest>`_  or download the latest stable development snapshot with

``git clone https://github.com/lucventurini/mikado.git; cd mikado``

If you desire, the unstable development version can be obtained with the command

``git checkout development``

in the Git directory. Please note that at the time of this writing development proceeds quite rapidly.

Building and installing from source
-----------------------------------

If you desire to install Mikado from source, this can be achieved with

``python3 setup.py bdist_wheel``

Followed by

``pip3 install dist/*whl``

..note:
    If you want to update your installation of Mikado, the command to be executed is
    ``pip install -U dist/*whl``

Testing the installed module
----------------------------

It is possible to test whether Mikado has been built successfully by opening a python3 interactive session and digiting::

    >>  import Mikado
    >>  Mikado.test()

.. block-end

This will run all the tests included in the suite. Although code coverage is not perfect yet, it is at 70% for the whole package and over 80% for the core components.

.. _GitHub: https://github.com/lucventurini/mikado
