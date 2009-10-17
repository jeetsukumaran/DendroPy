.. include:: common.rst

******************************************************
|dendropylogo| DendroPy Phylogenetic Computing Library
******************************************************

Introduction
============

|DendroPy|_ is a |Python|_ library for phylogenetic scripting, simulation, data processing and manipulation by |js| and |mth|.

Installation
============

The current public release of DendroPy is available from the Python Package Index:

    |dendropy_download_page|
    
If you have `setuptools <http://pypi.python.org/pypi/setuptools>`_ installed on your system, then you can install DendroPy by simply running the following command::

    $ sudo easy_install -U dendropy
    
Otherwise you will have to download the distribution archive from `the Python Package Index <http://pypi.python.org/pypi/DendroPy>`_, unpack it, and run the setup yourself::

    $ tar -xvzf DendroPy-2.5.1.tar.gz
    $ cd DendroPy-2.5.1
    $ sudo python setup.py install
    
Repository Access
=================
The DendroPy public-access `Git <http://git-scm.com/>`_ repository can be cloned from:
    
        |dendropy_public_repo|
        
Documentation
==============
The `DendroPy Cookbook <cookbook.html>`_, consisting of tutorials, walk-throughs, and annotated examples, can be found here:

        |dendropy_cookbook|
    
Included Scripts
================
* `SumTrees <sumtrees.html>`_ is a script that faciliates the summarization of non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.

Citation
=========

If you use this library either in whole or in part in your analysis, or use any code derived from it, please cite it as:

    Sukumaran, J. and Mark T. Holder. 2009. *DendroPy Phylogenetic Computing Library Version 2.5.1* (|dendropy_homepage|).
    
Bugs, Suggestions, Comments, etc.
=================================

If you encounter any problems, errors, crashes etc. while using this program, please let me know at jeet@ku.edu.
Please feel free to contact me if you have any other questions, suggestions or comments as well.

Copyright, License and Warranty
===============================

Copyright 2009 Jeet Sukumaran and Mark T. Holder.
  
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the `GNU General
Public License <http://www.gnu.org/licenses/gpl.html>`_ for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Acknowledgments
================

Portions of `DendroPy
<http://packages.python.org/DendroPy/>`_ were developed under `CIPRES <http://www.phylo.org>`_, a multi-site collaboration funded by the `NSF <http://www.nsf.gov/>`_ Information Technology Research (ITR) program grant entitled "`BUILDING THE TREE OF LIFE: A National Resource for Phyloinformatics and Computational Phylogenetics <http://www.phylo.org/about/acknowledgements>`_".

.. image:: logo_cipres.gif
    :height: 40   
    :target: http://www.phylo.org/

.. image:: nsf.gif
    :width: 40
    :target: http://www.nsf.gov/
