Introduction
============

DendroPy is a Python library for phylogenetic scripting, simulation, data processing and data manipulation, with its primary home page at:

    http://packages.python.org/DendroPy/

DendroPy provides classes and functions for working with phylogenetic data such as trees and character matrices, and supports reading and writing of the data in a range of standard phylogenetic data formats, such as NEXUS, NEWICK, NeXML, Phylip, FASTA, etc.

In addition, scripts for performing some useful phylogenetic computations are distributed as part of the libary, such as SumTrees, which summarizes the support for splits or clades given by a posterior sample of phylogenetic trees.

Requirements and Installation
=============================

DendroPy runs under any version of **Python 2.x greater than Python 2.4 up to, but not including, Python 3.0**.

If you have `setuptools <http://pypi.python.org/pypi/setuptools>`_ installed, you can install the latest public release of DendroPy directly from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_ by running::

    $ sudo easy_install -U dendropy

Alternatively, if you have `pip <http://pypi.python.org/pypi/pip>`_ installed, you can run::

    $ sudo pip install dendropy

If you have downloaded the source code distribution archive from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_, you can unarchive it and install it from the local source by running::

    $ tar -xvzf DendroPy-3.0.0.tar.gz
    $ cd DendroPy-3.0.0
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
