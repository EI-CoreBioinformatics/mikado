.. _SQLAlchemy: http://www.sqlalchemy.org/
.. _Portcullis: https://github.com/maplesond/portcullis
.. _BED12: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

.. _configure-scoring-tutorial:

How to create a scoring configuration file
==========================================

The current version of Mikado relies upon the experimenter to determine how desirable transcripts should look like, and
how to prioritise across different possible isoforms. These instructions are provided through **scoring configuration files**,
whose detailed description can be found in the :ref:`general documentation <scoring_files>`. In this section, we will
expose a general way to create your own configuration file.

First step: defining transcript requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When defining which transcripts are desirable for Mikado, the first step is to define which are the minimum attributes that a model should have to be even considered. This information is defined in the "requirements" section