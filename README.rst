Introduction
============

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

    http://packages.python.org/DendroPy/

Requirements and Installation
=============================

DendroPy 4.x runs under Python 3 (all versions > 3.1) and Python 2 (Python 2.7 only).

DendroPy 4.x is under development, and is not suitable yet for public or production use.

But it will soon be!

Testing
=======

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