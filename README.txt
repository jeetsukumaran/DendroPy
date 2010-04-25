Introduction
============

DendroPy is a Python library for phylogenetic computing, supporting object-oriented reading, writing, simulation, processing, and manipulation of phylogenetic data.
It is designed both as a stand-alone library for phylogenetics, as well as to support cross-library phyloinformatic pipelines.
DendroPy provides classes and functions for working with phylogenetic data such as trees and character matrices, and supports reading and writing of the data in a range of standard phylogenetic data formats, such as NEXUS, NEWICK, NeXML, Phylip, FASTA, etc.
In addition, scripts for performing some useful phylogenetic computations are distributed as part of the libary, such as SumTrees, which summarizes the support for splits or clades given by a posterior sample of phylogenetic trees.

The primary home page for DendroPy, with detailed tutorials and documentation, is at:

    http://packages.python.org/DendroPy/

Requirements and Installation
=============================

DendroPy runs under any version of Python 2 from 2.4 upwards (i.e., Python 2.4, 2.5, 2.6, 2.7, etc.). At present, it does not run under Python 3.

If you have `setuptools <http://pypi.python.org/pypi/setuptools>`_ installed, you can install the latest public release of DendroPy directly from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_ by running::

    $ sudo easy_install -U dendropy

Alternatively, if you have `pip <http://pypi.python.org/pypi/pip>`_ installed, you can run::

    $ sudo pip install dendropy

If you download the source code archive from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_, you can unarchive it and install it from the local source by running::

    $ tar -xvzf DendroPy-3.x.x.tar.gz
    $ cd DendroPy-3.x.x
    $ sudo python setup.py install

Documentation
=============

A detailed example-rich "cookbook"-style tutorial on using DendroPy can be found here:

    http://packages.python.org/DendroPy/tutorial/index.html

While the API reference can be found here:

    http://packages.python.org/DendroPy/library/index.html

Source Code Repository
======================

The DendroPy source code is version-controlled using `Git <http://git-scm.com/>`_, and the `DendroPy Git repository <http://github.com/jeetsukumaran/DendroPy>`_ can be cloned by running::

    $ git clone git://github.com/jeetsukumaran/DendroPy.git

More Information
=================
More information, including documentation, tutorials, citation, license, etc. can be found in the `DendroPy home page <http://packages.python.org/DendroPy/>`_:

    http://packages.python.org/DendroPy/

