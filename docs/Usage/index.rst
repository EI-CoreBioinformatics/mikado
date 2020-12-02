.. _Portcullis: https://github.com/maplesond/portcullis
.. _Usage:

Usage
=====


Mikado is composed of four different programs (``configure``, ``prepare``, ``serialise``, ``pick``) which have to be executed serially to go from an ensemble of different assemblies to the final dataset. In addition to these core programs, Mikado provides a utility to compare annotations, similarly to CuffCompare and ParsEval (*compare*), and various other minor utilities to perform operations such as extracting regions from a GFF, convert between different gene annotation formats, etc.


Mikado pipeline stages
~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   Configure
   Prepare
   Serialise
   Pick

The Mikado pipeline is composed of four different stages, that have to be executed serially:

#. :ref:`configure`, for creating the configuration file that will be used throughout the run.
#. :ref:`prepare`, for collapsing the input assemblies into a single file. After this step, it is possible to perform additional analyses on the data such as TransDecoder (highly recommended), Portcullis_, or BLAST.
#. :ref:`serialise`, to gather all external data into a single database.
#. :ref:`pick`, to perform the actual selection of the best transcripts in each locus.

Compare
~~~~~~~

.. toctree::
   :maxdepth: 1

   Compare


Mikado also comprises a dedicated utility, :ref:`Mikado compare <Compare>`, to assess the similarity between two annotations.

Daijin
~~~~~~

.. toctree::
   :maxdepth: 1

   Daijin

Mikado provides a pipeline manager, :ref:`Daijin <Daijin>`, to align and assemble transcripts with multiple methods and subsequently choose the best assemblies among the options. The pipeline is implemented using Snakemake [Snake]_.


Miscellaneous utilities
~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   Utilities

Finally, Mikado provides some dedicated utilities to perform common tasks.

* Some of the utilities are :ref:`integral to the Mikado suite <utils>` and can be accessed as subcommands of Mikado. These utilities comprise programs to calculate annotation statistics, retrieve or exclude specific loci from a file, etc.
* Other utilities are provided as :ref:`stand-alone scripts <included_scripts>`; while some of them directly depend on the Mikado library, this is not necessarily the case for them all.
