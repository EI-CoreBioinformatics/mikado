Usage
=====


Mikado is composed of four different programs (*configure, prepare, serialise, pick*) which have to be executed serially to go from an ensemble of different assemblies to the final dataset. In addition to these core programs, Mikado provides a utility to compare annotations, similarly to CuffCompare and ParsEval (*compare*), and various other minor utilities to perform operations such as extracting regions from a GFF, convert between different gene annotation formats, etc.

Mikado pipeline
---------------

The Mikado pipeline is composed of four different stages, that have to be executed serially:

#. :ref:`configure`, for creating the configuration file that will be used throughout the run.
#. :ref:`prepare`, for collapsing the input assemblies into a single file. After this step, it is possible to perform additional analyses on the data such as TransDecoder (highly recommended), Portcullis, or BLAST.
#. :ref:`serialise`, to gather all external data into a single database.
#. :ref:`pick`, to perform the actual selection of the best transcripts in each locus.

Compare
-------

Mikado provides a :ref:`dedicated utility <compare>` to compare two different annotations.

Utilities
---------

Mikado provides :ref:`some utilities <utils>` directly accessible from the command line, to perform standard operations such as obtaining basic statistics from a GTF/GFF file.


Mikado stages
=============

.. toctree::
   :maxdepth: 1
   :numbered:

   Usage_dir/Configure
   Usage_dir/Prepare
   Usage_dir/Serialise
   Usage_dir/Pick

Mikado utilities
================

.. toctree::
   :maxdepth: 1
   :numbered:

   Usage_dir/Compare
   Usage_dir/Daijin
   Usage_dir/Utilities
