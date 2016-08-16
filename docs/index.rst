.. Mikado documentation master file, created by
   sphinx-quickstart on Mon Jul 18 14:33:33 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |python_badge| image:: https://img.shields.io/pypi/pyversions/snakemake.svg?style=flat-square
   :target: https://www.python.org/
.. |snake_badge| image:: https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)]
   :target: https://bitbucket.org/snakemake/snakemake/wiki/Home

============================
Mikado: pick your transcript
============================

|python_badge| |snake_badge|

:Authors:
    Venturini Luca,
    Caim Shabhonam,
    Mapleson Daniel,
    Kaithakottil Gemy George,
    Swarbreck David
:Version: 1.0 (July 2016)

Mikado is a lightweight Python3 pipeline to extract the best transcript models from multiple transcript assemblies.

Contents:

.. toctree::
   :maxdepth: 2
   :numbered:

   Introduction
   Installation
   Scoring
   Usage/index
   Tutorial/index
   References
   Library/index


Citing
~~~~~~

We are currently working on our paper, and we will be releasing a pre-print shortly.
In the meantime, if you use Mikado please reference our github page: `https://github.com/lucventurini/mikado <https://github.com/lucventurini/mikado>`_


Availability and License
~~~~~~~~~~~~~~~~~~~~~~~~

Open source code available on github: `https://github.com/lucventurini/mikado <https://github.com/lucventurini/mikado>`_

This documentation is hosted publicly on read the docs: `https://mikado.readthedocs.org/en/latest/ <https://mikado.readthedocs.org/en/latest/>`_

Mikado is available under `GNU GLP V3 <http://www.gnu.org/licenses/gpl.txt>`_.

Acknowledgements
~~~~~~~~~~~~~~~~

Mikado has greatly benefitted from the public libraries, in particular the [NetworkX]_ library, Scipy and Numpy ([Scipy]_, [Numpy]_), Intervaltree [PYinterval]_, and the BX library for a Cython implementation of interval trees. Mikado has also been constantly optimised using Snakeviz, a tool which proved invaluable during the development process.


Credits
~~~~~~~

 - Luca Venturini (The software architect and developer)
 - Shabhonam Caim (Primary tester and analytic consultancy)
 - Daniel Mapleson (Developer of PortCullis and of the Daijin pipeline)
 - Gemy Kaithakottil (Tester and analytic consultancy)
 - David Swarbreck (Annotation guru and ideator of the pipeline)