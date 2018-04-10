%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DendroPy Phylogenetic Computing Library
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. toctree::
    :maxdepth: 3
    :hidden:

    Introduction <self>
    downloading.rst
    primer/index.rst
    library/index.rst
    schemas/index.rst
    programs/index.rst
    glossary.rst
    migration.rst
    changes.rst

DendroPy is a |Python|_ library for phylogenetic computing.
It provides classes and functions for the simulation, processing, and manipulation of phylogenetic trees and character matrices, and supports the reading and writing of phylogenetic data in a range of formats, such as NEXUS, NEWICK, NeXML, Phylip, FASTA, etc.
Application scripts for performing some useful phylogenetic operations, such as data conversion and tree posterior distribution summarization, are also distributed and installed as part of the libary.
DendroPy can thus function as a stand-alone library for phylogenetics, a component of more complex multi-library phyloinformatic pipelines, or as a scripting "glue" that assembles and drives such pipelines.

Prerequisites
=============

DendroPy is a pure-Python library with no dependencies, and runs under any version of Python 3 and Python 2.7. If you want to run DendroPy under Python version of less than 2.7, you should download the `DendroPy legacy release <https://github.com/jeetsukumaran/DendroPy/releases/tag/v3.12.1>`_.

.. versionchanged:: 4.0.0
   Python 2.4, 2.5, and 2.6 support removed.
   Python 3.x support added.

Installing
==========

DendroPy is fully easy-installable and can be installed using |pip|_::

    $ sudo pip install -U dendropy

Note: the "``sudo``" command should only be used if installing system-wide on a machine on which you have administrative privileges. Otherwise, you would use the "``--user``" flag for a local user install::

    $ pip install --user -U dendropy

You can install directly from the main GitHub repository using::

    $ pip install git+https://github.com/jeetsukumaran/DendroPy.git
    $ pip install git+git://github.com/jeetsukumaran/DendroPy.git

Alternatively, if you want to install from a particular branch, e.g., the latest development branch, ``development-master``::

    $ pip install git+https://github.com/jeetsukumaran/DendroPy.git@development-master

Or::

    $ pip install git+git://github.com/jeetsukumaran/DendroPy.git@development-master

If you do not have |pip|_ installed, you should *definitely* `install it <https://pip.pypa.io/en/latest/installing.html>`_ !
Other ways of obtaining and installing DendroPy (e.g., by downloading the |dendropy_source_archive|_, or by cloning the `DendroPy Git repository <http://github.com/jeetsukumaran/DendroPy>`_), are discussed in detail in the ":doc:`/downloading`" section.

Documentation
==============

    :doc:`Downloading and Installing DendroPy </downloading>`

        The many ways to get DendroPy up and running on your system.

    :doc:`DendroPy Primer </primer/index>`

        A detailed primer on how to use the DendroPy library, with lots of annotated practical examples and code walk-throughs.

    :doc:`DendroPy Library API Reference </library/index>`

        The technical details of the modules, classes and methods of the DendroPy library. Almost all of this information is also readily available from the |Python|_ interpreter by invoking ``help()`` on an object or name.

    :doc:`SumTrees User Manual </programs/sumtrees>`

        How to use SumTrees, an application script bundled with DendroPy that faciliates the summarization of non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.

    :doc:`Migration Guide </migration>`

        DendroPy 4 is improved in *many* ways over DendroPy 3. However, some of the changes are significant enough to break code written for DendroPy 3. This reference provides an overview of these changes, as well as notes on how to fix, handle or otherwise work-around issues that might result from trying to run code written for DendroPy 3 under DendroPy 4.


    :doc:`Change History </changes>`

        A summary of major changes (new features, bug fixes, bug creations, etc.) of each release of DendroPy.

Feature Requests
================

We are constantly expanding and improving the DendroPy library.
And we are constantly looking for ideas on *how* to expand and improve DendroPy.
If you have an idea that you would like to see implemented in DendroPy, or a need/requirement for some sort of functionality that DendroPy does not provide or does not provide easily, or, for that matter, a suggestion/idea for how something DendroPy already does that can be done better, *please* let us know by posting it as an issue on the |dendropy_issues|_ page.
We take these very seriously, and look forward to ideas/suggestions/requests/critiques such as these, not only because they improve the library for the entire community of users, but also because they sometimes provide for interesting side-projects for us.

Bug Reports and Other Issues
============================

Everyone has issues.
And sometimes, not pointing them out is enabling rather than helping. So, if you encounter or think you have encountered a problem in the way DendroPy works, please report it on the |dendropy_issues|_ page.
The DendroPy library has a really large suite of tests, but it is also a really, really, really, *really* large and complex library, and there are many nooks and crannies in which can scurry many crooks and nannies, so you will be doing us *and* the community a really generous favor by reporting bugs and problems.
Even if you are not sure if you are actually dealing with a bug, please do not hesitate to report it as an issue: false positives are better than false negatives in this context.

Having said that, *please* do take the trouble to write-up a decent bug report.
It really is quite frustrating to receive vague and lackadaisical "reports" that are little more than stream-of-consciousness responses to stimuli from the monitor that were almost accidentally blurted into the keyboard (e.g., "doesn't work").
At the very least, you should provide:

    -   A brief explication of what occurred, describing the operation that you were trying to do, and the result (or non-result) that you observed that led you to think there is an error.
    -   The environment in which this error occurred. This includes, at the very least, the operating system, as well as things such as the python version, DendroPy version, the installation locations of DendroPy and the Python libraries, etc. Apart from the operating system, you can retrieve all this information by running the following command::

            $ python -m dendropy

        The version of SumTrees that ships with DendroPy 4 onwards includes a special flag, "``--describe``" (i.e., you would type "``sumtrees.py --describe``"), that also provides this information.
        Including the details of either of these commands in their entirety along with the operating system is not only useful, but essential.
    -   We have gone to great lengths to write sensible and meaningful error messages, and that chunk of text that usually appears on the screen on error (i.e., the "stack trace") is packed with useful information, and should be included in its entirety should it appear.
    -   We need to be able to reproduce the error on our side to be able to fix. Thus, providing a *self-contained* and *minimum* example of the error is crucial. This might involve a little bit of work on your side, to extract the essential bits of code from their context and to ensure that it can run, at least, up to the point where it does not due to the error. The "steps to reproduce" section of a bug report is so important and useful that anyone who makes sure to include a good one into their bug report would be perfectly justified in expecting that strawberry cheesecake and a orange mocha frappuccino get delivered to their doorstep the moment the report is posted.

Help, Discussion, and Comments
==============================

The |dendropy_users|_ discussion group is available for general discussion on the usage of the DendroPy library.
If you want more information on how to use a particular feature or carry out specific task or want to make some general comments on the DendroPy library, you should post a message to the |dendropy_users|_ discussion group.
I follow and respond to messages on this group fairly regularly, so you should get a response within 24 hours.


Of course, we accept all bug reports, bad or good, and I guess we would prefer to have a bad bug report than no report at all, so if all the above seem enough of a hassle to discourage you from posting an issue at all, feel free to go ahead and write it up any way you see fit.

.. include:: citation.inc
.. include:: license.inc
.. include:: acknowledgements.inc

