========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |appveyor| |requires|
        | |coveralls| |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|

.. |docs| image:: https://readthedocs.org/projects/python-fluctmatch/badge/?style=flat
    :target: https://readthedocs.org/projects/python-fluctmatch
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/tclick/python-fluctmatch.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/tclick/python-fluctmatch

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/github/tclick/python-fluctmatch?branch=master&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/tclick/python-fluctmatch

.. |requires| image:: https://requires.io/github/tclick/python-fluctmatch/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/tclick/python-fluctmatch/requirements/?branch=master

.. |coveralls| image:: https://coveralls.io/repos/tclick/python-fluctmatch/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/tclick/python-fluctmatch

.. |codecov| image:: https://codecov.io/github/tclick/python-fluctmatch/coverage.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/tclick/python-fluctmatch

.. |version| image:: https://img.shields.io/pypi/v/fluctmatch.svg
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/fluctmatch

.. |commits-since| image:: https://img.shields.io/github/commits-since/tclick/python-fluctmatch/v3.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/tclick/python-fluctmatch/compare/v3.0.0...master

.. |wheel| image:: https://img.shields.io/pypi/wheel/fluctmatch.svg
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/fluctmatch

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/fluctmatch.svg
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/fluctmatch

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/fluctmatch.svg
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/fluctmatch


.. end-badges

Elastic network model using fluctuation matching.

* Free software: BSD license

Introduction
============

Fluctuation matching is a different approach to protein structural analysis.
Typically, an elastic network model (ENM) creates springs between coarse-grain
sites and then uses normal mode analysis (NMA) to determine the vibrational
information that occurs within the structure; most ENMs emply study proteins
representing each residue at the C-alpha position. Fluctuation matching has been
programmed to incorporate more molecules of interest, including nucleic acids,
solvent, and ions. Residues can either be represented on the alpha-carbon; an
alpha-carbon and the sidechain; or the amino group, carboxyl group, and
sidechain.

Previous versions of fluctuation matching strictly employed CHARMM for all
calculations (average structure, initial bond statistics, and NMA). Furthermore,
the directory structure was such that all average structures were used to
determine the springs available within the system. The current code base as been
completely retooled allowing for easier definitions of additional coarse-grain
models, additions to analysis, and implementation of other MD packages. Because
MDAnalysis 0.16.2+ has also been employed (compared with 0.10.0 for fluctmatch
2.0), greater improvements have been made in the efficiency of the code.

Installation
============

::

    pip install fluctmatch

Documentation
=============

https://python-fluctmatch.readthedocs.io/

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
