.. include:: common.inc

**********************************************************
|dendropymasthead| DendroPy Phylogenetic Computing Library
**********************************************************

Contents:

.. toctree::
   :maxdepth: 2
   
   tutorial.rst
   sumtrees.rst

Introduction
============

|DendroPy|_ is a |Python|_ library for phylogenetic scripting, simulation, data processing and manipulation by |js|_ and |mth|_.

Installation
============

The current public release of DendroPy is available from the Python Package Index:

    |dendropy_download_url|
    
If you have `setuptools <http://pypi.python.org/pypi/setuptools>`_ installed on your system, then you can install DendroPy by simply running the following command::

    $ sudo easy_install -U dendropy
    
Otherwise you will have to download the distribution archive from `the Python Package Index <http://pypi.python.org/pypi/DendroPy>`_, unpack it, and run the setup yourself::

    $ tar -xvzf DendroPy-2.5.1.tar.gz
    $ cd DendroPy-2.5.1
    $ sudo python setup.py install
            
Documentation
==============

DendroPy Cookbook Tutorial
--------------------------
A `"cookbook"-style tutorial <tutorial.html>`_, consisting of annotated practical examples and code walk-throughs can be found here:

        |dendropy_tutorial_url|
        

Library Reference
-----------------
Most of the modules, classes and methods in the |DendroPy|_ library are richly documented, and this information is readily available from the |Python|_ interpreter by invoking ``help()`` on an object or name.

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`  

Included Scripts
================
* |SumTrees|_ is a script that faciliates the summarization of non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.

Repository Access
=================
The |DendroPy|_ public-access |Git|_ repository can be cloned from:
    
        |dendropy_public_repo_url|
            
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




