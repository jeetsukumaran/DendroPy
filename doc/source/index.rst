%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DendroPy Phylogenetic Computing Library
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 3
    :hidden:

    Introduction <self>
    downloading.rst
    tutorial/index.rst
    library/index.rst
    scripts/sumtrees.rst
    changes.rst

DendroPy is a |Python|_ library for phylogenetic scripting, simulation, data processing and manipulation by |js|_ and |mth|_.
DendroPy provides classes and functions for working with phylogenetic data such as trees and character matrices, and supports reading and writing of the data in a range of standard phylogenetic data formats, such as NEXUS, Newick, NeXML, PHYLIP, FASTA, etc.
In addition, scripts for performing some useful phylogenetic computations are distributed as part of the libary, such as SumTrees, which summarizes the support for splits or clades given by a posterior sample of phylogenetic trees.

Prerequisites
=============

DendroPy is a pure-Python library with no dependencies, and runs under any version of Python 2 from 2.4 upwards (i.e., Python 2.4, 2.5, 2.6, 2.7, etc.). At present, it does not run under Python 3.

.. versionchanged:: 3.2
   Python 2.4 support added.

Installing
==========

DendroPy is fully easy-installable and can be installed using |pip|_::

    $ sudo pip install dendropy

or |setuptools|_::

    $ sudo easy_intall -U dendropy

if these are available on your system.

This, and other ways of obtaining and installing DendroPy (e.g., by downloading the |dendropy_source_archive|_, or by cloning the `DendroPy Git repository <http://github.com/jeetsukumaran/DendroPy>`_), are discussed in detail in the ":doc:`downloading`" section.

Documentation
==============

    :doc:`Downloading and Installing DendroPy </downloading>`

        The many ways to get DendroPy up and running on your system.

    :doc:`DendroPy Tutorial </tutorial/index>`

        A detailed tutorial on how to use the DendroPy library, with lots of annotated practical examples and code walk-throughs.

    :doc:`DendroPy Library API Reference </library/index>`

        The technical details of the modules, classes and methods of the DendroPy library. Almost all of this information is also readily available from the |Python|_ interpreter by invoking ``help()`` on an object or name.

    :doc:`SumTrees User Manual </scripts/sumtrees>`

        How to use SumTrees, an application script bundled with DendroPy that faciliates the summarization of non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.

    :doc:`Change History <changes>`

        A summary of major changes (new features, bug fixes, bug creations, etc.) of each release of DendroPy.

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




