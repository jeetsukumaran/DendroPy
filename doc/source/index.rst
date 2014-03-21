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
    glossary.rst
    migration.rst
    changes.rst

DendroPy is a |Python|_ library for phylogenetic computing.
It provides classes and functions for the simulation, processing, and manipulation of phylogenetic trees and character matrices, and supports the reading and writing of phylogenetic data in a range of formats, such as NEXUS, NEWICK, NeXML, Phylip, FASTA, etc.
Application scripts for performing some useful phylogenetic operations, such as data conversion and tree posterior distribution summarization, are also distributed and installed as part of the libary.
DendroPy can thus function as a stand-alone library for phylogenetics, a component of more complex multi-library phyloinformatic pipelines, or as a scripting "glue" that assembles and drives such pipelines.

Prerequisites
=============

DendroPy is a pure-Python library with no dependencies, and runs under any version of Python 3 and Python 2.7. If you want to run DendroPy under Python version of less than 2.4, you should download the `DendroPy legacy release <https://github.com/jeetsukumaran/DendroPy/releases/tag/v3.12.1>`_.

.. versionchanged:: 4.0.0
   Python 2.4, 2.5, and 2.6 support removed.
   Python 3.x support added.

Installing
==========

DendroPy is fully easy-installable and can be installed using |pip|_::

    $ pip install dendropy

or |setuptools|_::

    $ easy_install -U dendropy

if these are available on your system.

These, and other ways of obtaining and installing DendroPy (e.g., by downloading the |dendropy_source_archive|_, or by cloning the `DendroPy Git repository <http://github.com/jeetsukumaran/DendroPy>`_), are discussed in detail in the ":doc:`downloading`" section.

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

    :doc:`Migration Guide </migration>`

        DendroPy 4 is improved in *many* ways over DendroPy 3. However, some of the changes are significant enough to break code written for DendroPy 3. This reference provides an overview of these changes, as well as notes on how to fix, handle or otherwise work-around issues that might result from trying to run code written for DendroPy 3 under DendroPy 4.


    :doc:`Change History </changes>`

        A summary of major changes (new features, bug fixes, bug creations, etc.) of each release of DendroPy.

Citation
=========

If you use this library either in whole or in part in your analysis, or use any code derived from it, please cite it as (replacing the version numbers with ones corresponding to the version that you used):

    |dendropy_citation|

Suggestions, Comments, Help, or Bug Reports
===========================================

The |dendropy_users|_ discussion group is available for general discussion on the usage of the DendroPy library.
If you want more information on how to use a particular feature or carry out specific task, would like to request a new feature or change an existing one, or even just want to make some general comments on the DendroPy library, you should post a message to the |dendropy_users|_ discussion group.
I follow and respond to messages on this group fairly regularly, so you should get a response within 24 hours.
If you prefer, you can also contact me directly via e-mail at jeet@ku.edu, but by posting to the |dendropy_users|_ group forum, you allow other users, current as well as future, to benefit from the discussion and results.

The |dendropy_users|_ discussion group is also a suitable forum for the reporting of bugs or errors encountered with DendroPy, as most other users would also appreciate learning of problems and issues with library.
However, again, if you prefer you can contact me directly via e-mail at jeet@ku.edu.


.. include:: license.inc
.. include:: acknowledgements.inc

