.. image:: https://raw.githubusercontent.com/jeetsukumaran/DendroPy/DendroPy4/doc/source/_static/dendropy_logo.png
   :align: right
   :alt: DendroPy

DendroPy is a Python library for phylogenetic computing.
It provides classes and functions for the simulation, processing, and
manipulation of phylogenetic trees and character matrices, and supports the
reading and writing of phylogenetic data in a range of formats, such as NEXUS,
NEWICK, NeXML, Phylip, FASTA, etc.  Application scripts for performing some
useful phylogenetic operations, such as data conversion and tree posterior
distribution summarization, are also distributed and installed as part of the
libary.  DendroPy can thus function as a stand-alone library for phylogenetics,
a component of more complex multi-library phyloinformatic pipelines, or as a
scripting "glue" that assembles and drives such pipelines.

The primary home page for DendroPy, with detailed tutorials and documentation, is at:

    http://dendropy.org/

DendroPy is also hosted in the official Python repository:

    http://packages.python.org/DendroPy/

Requirements and Installation
=============================

DendroPy 4.x runs under Python 3 (all versions > 3.1) and Python 2 (Python 2.7 only).

You can install DendroPy by running::

    $ sudo pip install dendropy

More information is available here:

    http://dendropy.org/downloading.html

Documentation
=============

Full documentation is available here:

    http://dendropy.org/

This includes:

    -   `A comprehensive "getting started" primer <http://dendropy.org/primer/index.html>`_ .
    -   `API documentation <http://dendropy.org/library/index.html>`_ .
    -   `Descriptions of data formats supported for reading/writing <http://dendropy.org/schemas/index.html>`_ .

and more.

Testing
=======

.. note::

    Note that some tests rely on PAUP* being available on your system.
    You will need to set the environmental variable ``DENDROPY_PAUP_EXECUTABLE_PATH`` to the path
    of the PAUP* executable for these tests to be run, e.g.::

        DENDROPY_PAUP_EXECUTABLE_PATH=/usr/local/bin/paup python setup.py test

    or::

        DENDROPY_PAUP_EXECUTABLE_PATH=/usr/local/bin/paup python -m dendropy.test

    If this variable is not set or set to "NONE", then any tests that rely on
    PAUP* will NOT be run.

Tests can be run by typing::

    $ python -m dendropy.test

By default, all tests are run. You can run specific by providing the
fully-qualified name of the modules, test cases, or specific test methods to
run, e.g.::

    $ python -m dendropy.test test_tokenizer
    $ python -m dendropy.test test_tokenizer.TestCase
    $ python -m dendropy.test test_tokenizer.TestCase.test1
    $ python -m dendropy.test test_tokenizer test_datamodel_taxon

Or special pre-defined sub-groups of tests, e.g.::

    $ python -m dendropy.test @datamodel
    $ python -m dendropy.test @dataio
    $ python -m dendropy.test @datamodel @newick

A list of all available sub-groups can be seen by::

    $ python -m dendropy.test --help-testgroups

For any tests run, you can set the level at which the test progress is logged
by::

    $ python -m dendropy.test -l DEBUG all

For all options, type::

    $ python -m dendropy.test --help

License and Warranty
====================

Please see the file "LICENSE.rst" for details.
