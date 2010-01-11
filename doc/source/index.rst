#######################################
DendroPy Phylogenetic Computing Library
#######################################

.. toctree::
    :hidden:
    :maxdepth: 2

    tutorial/index.rst
    library/index.rst
    scripts/sumtrees.rst

Introduction
============

DendroPy is a |Python|_ library for phylogenetic scripting, simulation, data processing and manipulation by |js|_ and |mth|_.
DendroPy provides classes and functions for working with phylogenetic data such as trees and character matrices, and supports reading and writing of the data in a range of standard phylogenetic data formats, such as NEXUS, NEWICK, NeXML, Phylip, FASTA, etc.
In addition, scripts for performing some useful phylogenetic computations are distributed as part of the libary, such as SumTrees, which summarizes the support for splits or clades given by a posterior sample of phylogenetic trees.

Prerequisites
=============

DendroPy runs under any version of **Python 2.x greater than Python 2.4 up to, but not including, Python 3.0**.

Installation
============

If you have `setuptools <http://pypi.python.org/pypi/setuptools>`_ installed, you can install the latest release of DendroPy directly from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_ by running::

    $ sudo easy_install -U dendropy

Alternatively, if you have `pip <http://pypi.python.org/pypi/pip>`_ installed, you can run::

    $ sudo pip install dendropy

If you have downloaded the source code distribution archive from the `Python Package Index <http://pypi.python.org/pypi/DendroPy/>`_, you can unarchive it and install it from the local source by running::

    $ tar -xvzf DendroPy-3.0.0.tar.gz
    $ cd DendroPy-3.0.0
    $ sudo python setup.py install

Documentation
==============

DendroPy Tutorial
-----------------

The :doc:`DendroPy Tutorial </tutorial/index>` is a "cookbook"-style tutorial, consisting of annotated practical examples and code walk-throughs, and can be found here:

    :doc:`DendroPy Tutorial </tutorial/index>`

DendroPy Library Reference
--------------------------
The :doc:`DendroPy Library Reference </library/index>` documents the modules, classes and methods of the DendroPy library, and can be found here:

    :doc:`DendroPy Library Reference </library/index>`

Much of this information is also readily available from the |Python|_ interpreter by invoking ``help()`` on an object or name.

..  Indices and Tables
    ------------------

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`

Included Scripts
================
* :doc:`SumTrees </scripts/sumtrees>` is a script that faciliates the summarization of non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.

Repository Access
=================
#. If you do not already have it, get yourself |Git|_:
    - `Source <http://www.kernel.org/pub/software/scm/git/git-1.6.6.tar.gz>`_
    - `OS X binaries <http://code.google.com/p/git-osx-installer/downloads/list?can=3>`_
    - `Microsoft Windows <http://code.google.com/p/msysgit/downloads/list>`_.

#. Clone the `DendroPy Git repository <http://github.com/jeetsukumaran/DendroPy>`_::

    $ git clone git://github.com/jeetsukumaran/DendroPy.git

#. Install DendroPy in development mode::

    $ cd DendroPy
    $ sudo python setup.py develop

Citation
=========

If you use this library either in whole or in part in your analysis, or use any code derived from it, please cite it as (replacing the version numbers with ones corresponding to the version that you used):

    |dendropy_citation|

Bugs, Suggestions, Comments, etc.
=================================

If you encounter any problems, errors, crashes etc. while using this program, please let me know at jeet@ku.edu.
Please feel free to contact me if you have any other questions, suggestions or comments as well.


.. include:: license.inc
.. include:: acknowledgements.inc




