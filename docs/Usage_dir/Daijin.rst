.. _Snakemake: https://bitbucket.org/snakemake/snakemake/wiki/Home
.. _YAML: http://www.yaml.org/spec/1.2/spec.html

.. _Daijin:

The Daijin pipeline for driving Mikado
======================================

.. image:: https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)]
   :target: http://snakemake.bitbucket.org

No emperor or empress can lead its nation without a trusty chancellor to help him or her in organising the bureaucracy. Daijin, the Japanese minister, has the same role in Mikado - it smooths the path to go from a collection of read inputs (both RNA-Seq or long reads) to a polished assembly. The pipeline is based on Snakemake_ [Snake]_ and supports natively three cluster management systems (SLURM, PBS and LSF), plus any DRMAA-compliant batch submission system.

.. _daijin-configure:

Configure
~~~~~~~~~

This utility creates the configuration file that will drive Daijin, in YAML_ format. The file will need to be edited


.. _daijin-assemble:

Assemble
~~~~~~~~

.. image:: daijin_assemble.svg


.. _daijin-mikado:

Mikado
~~~~~~

.. image:: daijin_mikado.svg