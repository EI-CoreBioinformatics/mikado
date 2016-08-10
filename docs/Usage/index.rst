.. _Portcullis: https://github.com/maplesond/portcullis

Usage
=====


Mikado is composed of four different programs (*configure, prepare, serialise, pick*) which have to be executed serially to go from an ensemble of different assemblies to the final dataset. In addition to these core programs, Mikado provides a utility to compare annotations, similarly to CuffCompare and ParsEval (*compare*), and various other minor utilities to perform operations such as extracting regions from a GFF, convert between different gene annotation formats, etc.


Mikado pipeline stages
~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   :numbered:

   Configure
   Prepare
   Serialise
   Pick

The Mikado pipeline is composed of four different stages, that have to be executed serially:

#. :ref:`configure`, for creating the configuration file that will be used throughout the run.
#. :ref:`prepare`, for collapsing the input assemblies into a single file. After this step, it is possible to perform additional analyses on the data such as TransDecoder (highly recommended), Portcullis_, or BLAST.
#. :ref:`serialise`, to gather all external data into a single database.
#. :ref:`pick`, to perform the actual selection of the best transcripts in each locus.

Mikado utilities
~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   :numbered:

   Compare
   Daijin
   Utilities

In addition with the pipeline proper, Mikado includes the following:

#. A dedicated utility, :ref:`Mikado compare <Compare>`, to assess the similarity between two annotations.

#. A pipeline manager, :ref:`Daijin`, to align and assemble transcripts with multiple methods and subsequently drive Mikado on the assemblies.

#. Assorted utilities, some of them :ref:`part of the Mikado suite <utils>` and others provided as :ref:`accessory scripts <included_scripts>`, to perform standard operations such as calculating statistics from a GFF file or recover a specific suite of transcripts from an annotation file.
