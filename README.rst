========
Overview
========
Elastic network model using Structure Mechanics Statistical Learning.

* Free software: BSD license

Introduction
============

SMSL is a different approach to parameterize the multiscale models for protein structural analysis.
Typically, an elastic network model (ENM) creates springs between coarse-grain
sites and then uses normal mode analysis (NMA) to determine the vibrational
information that occurs within the structure; most ENMs employ study proteins
representing each residue at the C-alpha position. SMSL tool has been
programmed to incorporate different multiscale models for more molecules of interest, including nucleic acids,
solvent, and ions. Residues can either be represented on the alpha-carbon; an
alpha-carbon and the sidechain; or the amino N, carboxyl O, and
chemically specific atoms of sidechain.

Previous versions of fluctuation matching strictly employed CHARMM for all
calculations (average structure, initial bond statistics, and NMA). Furthermore,
the directory structure was such that all average structures were used to
determine the springs available within the system. The current code base as been
completely retooled allowing for easier definitions of additional coarse-grain
models, additions to analysis, and implementation of other MD packages. Because
MDAnalysis 0.16.2+ has also been employed (compared with 0.10.0 for fluctmatch
2.0), greater improvements have been made in the efficiency of the code.

PLease Cite:
1. N. Raj, T. Click, H. Yang and J.-W. Chu, Comput. Struct.Biotechnol. J., 2021, 19, 5309–5320

2. N. Raj, T. Click, H. Yang and J.-W. Chu, Chem. Sci. RSC,Chem. Sci., 2022,13, 3688-3696.
https://doi.org/10.1039/D1SC06184D

========
Overview
========

.. start-badges

 * - docs
   - |docs|
 * - tests
   - |requires|

.. |docs| image:: https://readthedocs.org/projects/python-fluctmatch/badge/?style=flat
    :target: https://readthedocs.org/projects/python-fluctmatch
    :alt: Documentation Status

.. |requires| image:: https://requires.io/github/nixnmtm/SMSL/requirements.txt?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/nixnmtm/SMSL/requirements.txt?branch=master

.. end-badges

Installation
============

::

1. Create a conda environment:
    (requirements.txt will install all the require packages)
conda create --name [env name] --file requirements.txt

2. conda activate [env name]

3. Download the Master branch in .zip format
4. pip install [zip file]

Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox

