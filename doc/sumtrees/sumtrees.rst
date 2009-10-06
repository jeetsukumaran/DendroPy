Introduction
============

SumTrees is a program to summarize non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.

The basis of the support assessment is typically given by a set of non-parametric bootstrap replicate tree samples produced by programs such as GARLI or RAxML, or by a set of MCMC tree samples produced by programs such as Mr. Bayes or BEAST.
The proportion of trees out of the samples in which a particular split is found is taken to be the degree of support for that split as indicated by the samples.
The samples that are the basis of the support can be distributed across multiple files, and a burn-in option allows for an initial number of trees in each file to be excluded from the analysis if they are not considered to be drawn from the true support distribution.

The support for the splits will be mapped onto one or more target trees either in terms of node labels or branch lengths.
The target trees can be supplied by yourself, or, if no target trees are given, then a majority-rule clade consensus tree will be constructed based on the samples given.
In the latter case, you have the option of specifying the minimum posterior probability or proportional frequency threshold for a clade to be included on the consensus tree.

Where to Find the Package
=========================

SumTrees is distributed and installed as part of the `DendroPy
<http://pypi.python.org/pypi/DendroPy>`_ phylogenetic computation library, which  can be found here: 

    http://pypi.python.org/pypi/DendroPy
    
How to Install the Package
==========================

The Easy Way
------------

Simply type the following command in your shell/terminal::

    sudo easy_install -U dendropy
    
And that is all it takes!    

This requires you have `setuptools <http://pypi.python.org/pypi/setuptools>`_ already installed. If you get an error message about "``easy_install``" not being found, then you do *not* have `setuptools <http://pypi.python.org/pypi/setuptools>`_ installed, and you need to `get it <http://pypi.python.org/pypi/setuptools>`_ and install it first.
It is *highly recommended* that you use this method, even if it means that you have to pre-install `setuptools <http://pypi.python.org/pypi/setuptools>`_.

The Not So Easy Way
-------------------

Alternatively, you can install DendroPy yourself, by following these steps.

1.  Download the latest release of DendroPy from `here <http://pypi.python.org/pypi/DendroPy>`_:

    http://pypi.python.org/pypi/DendroPy

    Note that you might end up downloading a newer version of DendroPy, in which case the version numbers in the file and directory names may not correspond exactly to those given in the examples. For instance, you may end up with "``DendroPy-2.2.1rc3.tar.gz``" rather than "``DendroPy-2.1.3.tar.gz``". However, as long as you substitute the correct archive and directory name in the examples and discussion below, everything else should remain the same. 

#.  Expand the downloaded archive

    This step varies depending on the operating system and the particular programs that you have installed. 
    In most cases, simply double-clicking on the file that you have downloaded should kick off the process.
    Otherwise, open a terminal shell window and go to the directory in which you have downloaded the archive, and type "``tar -xvzf DendroPy-2.1.3.tar.gz``".
    For example, say you saved the downloaded file on your desktop.
    Then, opening up the terminal and entering the following commands will take you to your Desktop and expand the archive::
    
        $ cd ~/Desktop
        $ tar -xvzf DendroPy-2.1.3.tar.gz
    
    One way or another, you should end up with a directory called "``DendroPy-2.1.3``" or something similar, which contains the entire Dendropy package.
    
#.  Install the library

    In the terminal shell, go to the directory of the package that you have just archived and type "``sudo python setup.py install``".
    Continuing the example from above::

        $ cd DendroPy-2.1.3
        $ sudo python setup.py install
        
    The library installation will automatically create an executable script called "``sumtrees.py``" and place it on your system path for you, so that you can call it from anywhere.

#.  If the installation was successful, you should be able to type     "``sumtrees.py``" in the terminal window and see something like the following (with possibly a different date or version number)::

        ====================================================================
        SumTrees - Phylogenetic Tree Split Support Summarization
        Version 1.1.1 (Mar 03 2009)
        By Jeet Sukumaran and Mark T. Holder
        (using the DendroPy Phylogenetic Computation Library Version 2.1.3)
        ====================================================================
        
        No sources of support specified or could be found. Please provide the
        path to at least one (valid and existing) file containing non-
        parametric or MCMC tree samples to summarize.
        
    You can now delete the original downloaded archive and unpacked directory if you want.         

How to Use the Program
======================

SumTrees is typically invoked by providing it a list of one or more tree files to be summarized::

    $ sumtrees.py [OPTIONS] <TREEFILE> [<TREEFILE> [<TREEFILE> ...]]]

Common options include specification of a target topology onto which to map the support ("``-t``" or "``--target``"), an output file ("``-o``" or "``--output``"), and a burn-in value ("``-b``" or "``--burnin``").

Full help on program usage and options is given by using the "``--help``" option::
    
    $ sumtrees.py --help
    
    
Quick Recipes
=============
                        
Non-parametric Bootstrap Support of a Model Tree
------------------------------------------------
Calculate support for nodes on a specific tree, "``best.tre``" as given by a set of tree files, with support reported as percentages rounded to integers, and saving the result to "``results.sumtree``"::

    $ sumtrees.py --decimals=0 --percentages --target=best.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -d0 -p -t best.tre treefile1.tre treefile2.tre treefile3.tre

Summarization of Posterior Probabilities of Clades with a Consensus Tree
------------------------------------------------------------------------
Summarize a set of tree files using a 95% rule consensus tree, with support for clades indicated as proportions (posterior probabilities) and branch lengths the mean across all trees, dropping the first 200 trees in each file as a burn-in, and saving the result to "``results.sumtree``"::

    $ sumtrees.py --min-clade-freq=0.95 --burn-in=200 --support-as-labels --output=results.sumtrees treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -f0.95 -b200 -l -o results.sumtrees treefile1.tre treefile2.tre treefile3.tre
 

Tutorials and Examples
======================

At its most basic, you will need to supply SumTrees with the path to one or more tree files in NEWICK or NEXUS format that you want to summarize::
    
    $ sumtrees.py phylo.tre

The above command will construct a 50% majority-rule consensus tree of the all trees found in the file "``phylo.tre``", with the internal node labels of the resulting consensus tree indicating the proportion of trees in "``phylo.tre``" in which that clade was found, while the branch lengths of the resulting consensus tree being set to the mean of the branch lengths of that clade across all the trees in "``phylo.tre``".

If you have split searches across multiple runs (across, for example, multiple computers, so as to speed up the search time), such that you have multiple tree files ("``phylo.run1.tre``", "``phylo.run2.tre``", "``phylo.run3.tre``", ...), you can instruct SumTrees to consider all these files together when summarizing the support by simply listing them one after another separated by spaces::
    
    $ sumtrees.py phylo.run1.tre phylo.run2.tre phylo.run3.tre

As before, the above command will construct a 50% majority-rule consensus tree with clade supported indicated by internal node labels and branch lengths being the mean across all trees, but this time it will use all the trees found across all the files listed: "``phylo.run1.tre``", "``phylo.run2.tre``", and "``phylo.run3.tre``".

You will notice that the final resulting tree is displayed to the terminal and not saved anywhere.
It will probably be more useful if we can save it to a file for visualization for further analysis.
This can be done in one of two ways, either by redirecting the screen output to a file, using the standard (at least on UNIX and UNIX-like systems) redirection operator, ``>``::
    
    $ sumtrees.py phylo.tre > phylo.consensus.sumtrees

or by using the or "``--output``" option::
    
    $ sumtrees.py --output=phylo.consensus.sumtrees phylo.tre 

If the files are in different directories, or you are not in the same directory as the files, you should use the full directory path specification::
    
    $ sumtrees.py --output=/Users/myself/MyProjects/phylo1/final/phylo.consensus.sumtrees /Users/myself/MyProjects/phylo1/phylo.tre 
 
More extended options specify things like: where to save the output (by default it goes to the screen), the topology or tree to which to map the support (user-supplied or consensus tree), the output format (NEXUS by default, but it can also be NEWICK), whether support is indicated in terms of proportions or percentages etc. 
All of these options are specified on the command line when invoking the program, with multiple options separated by spaces.
Many of the options have two alternate forms, a long form (a word or phrase preceded by two dashes, e.g., "``--option``") and a short form (a single letter preceded by a single dash, "``-o``").
The long form of the options needs an equals sign before setting the paramater (e.g., "``--option=1``"), while the short one does not (e.g., "``-o1``" or "``-o 1``").
Most of the options have default values that will be used if not explicitly set when the program is invoked.
The order that the options are given does *not* matter, i.e., "``sumtrees.py --option1=something --option2=something``" is the same as "``sumtrees.py --option2=something --option1=something``".
As mentioned above, full details on these options, their long and short forms, as well as their default values will be given by invoking the program with the "``--help``" or "``-h``" option: "``sumtrees.py --help``".

Summarizing Non-Parametric Bootstrap Support with a Consensus Tree
------------------------------------------------------------------

Say you have completed a 1000-replicate non-parametric analysis of your dataset using a program such as GARLI or RAxML.
You want to construct a 70% majority-rule consensus tree of the replicates, with support indicated as percentages on the node labels.
If the bootstrap replicates are in the file "``phylo-boots.tre``", you would then enter the following command::
    
    $ sumtrees.py --min-clade-freq=0.7 --percentages --decimals=0 phylo-boots.tre 

Or, using the short option syntax::
    
    $ sumtrees.py -f0.7 -p -d0 phylo-boots.tre 

Here, the "``--min-clade-freq=0.7``" or "``-f0.7``" option lowers the minimum threshold for clade inclusion to 70%.
If you want a 95% majority-rule consensus tree instead, you would use "``--min-clade-freq=0.95``" or "``-f0.95``".
The default threshold if none is specified is 0.5 (50%).
The "``--percentages``" or "``-p``" option instructs SumTrees to report the support in terms of percentages, while the "``--decimals=0``" or "``-d 0``" option instructs SumTrees not to bother reporting any decimals. 
Note that even if you instruct SumTrees to report the support in terms of percentages, the minimum clade inclusion threshold is still given in terms of proportions.

Again, if we want to actually save the results to the file, we should use the "``--output``" option::
    
    $ sumtrees.py --output=phylo-mle-support.sumtrees --min-clade-freq=0.7 --percentages --decimals=0 phylo-boots.tre
    $ sumtrees.py -o phylo-mle-support.sumtrees -f0.7 --p --d0 phylo-boots.tre

Summarizing Non-Parametric Bootstrap Support of an Estimated Tree
-----------------------------------------------------------------

Say you also have a maximum likelihood estimate of the phylogeny, and want to annotate the nodes of the maximum likelihood tree with the proportion of trees out of the bootstrap replicates in which the node is found.
Then, assuming your maximum likelihood tree is in the file, "``phylo-mle.tre``", and the bootstrap tree file is "``phylo-boots.tre``", you would use the "``--target``" options, as in the following command::
    
    $ sumtrees.py --target=phylo-mle.tre phylo-boots.tre

Here, "``--target``" specifies the target topology onto which the support will be mapped, while the remaining (unprefixed) argument specifies the tree file that is the source of the support. 
An equivalent form of the same command, using the short option syntax is::
    
    $ sumtrees.py -t phylo-mle.tre phylo-boots.tre

If you want the support expressed in percentages instead of proportions, and the final tree saved to a file, you would enter::
    
    $ sumtrees.py --output phylo-mle-support.sumtrees --target phylo-mle.tre --proportions --decimals=0 phylo-boots.tre
    $ sumtrees.py -o phylo-mle-support.sumtrees -t phylo-mle.tre -p -d0 phylo-boots.tre

Summarizing MCMC Trees
----------------------

Say you have just completed a BEAST analysis resulting in a file of MCMC tree samples called "``phylo.trees``". 
While the program TreeAnnotator that is distributed along with BEAST does construct a tree summarizing the split support for you, it produces a MCCT topology as the summary tree.
This is not the same summarization strategy as used by Mr. Bayes using its "``sumt``" command, and thus the two summary trees are not truly directly comparable.
You can use SumTrees to construct a majority-rule clade consensus tree out of your BEAST MCMC samples, which you can then use to compare with your Mr. Bayes tree::
    
    $ sumtrees.py phylo.trees

This command will construct a 50% majority rule clade consensus tree out of the all the trees found in "``phylo.trees``", label each node with its posterior probability and output the resulting tree in NEXUS format to the terminal.

Of course, we want to discard the first few samples of trees, as these were probably not drawn in frequencies in proportion to the stationary distribution of the chain.
To do this::
    
    $ sumtrees.py --burnin=200 phylo.trees

The above command will cause SumTrees to ignore the first 200 trees it finds in the file for all its calculations.

Again, instead of displaying the tree to the screen we can save it directly to a file, either by redirecting the screen output to a file::
    
    $ sumtrees.py --burnin=200 phylo.trees > phylo.trees.sumtrees

or by using the "``-o``" or "``--output``" option::
    
    $ sumtrees.py --output=phylo.trees.sumtrees --burnin=200 phylo.trees

We might also have split up our analysis into multiple independent runs, resulting in multiple MCMC tree sample files (e.g., "``phylo1.trees``", "``phylo2.trees``" and "``phylo3.trees``").
We can ask SumTrees to summarize posterior probability from across all these runs, treating the first 200 trees in *each* sample file as a burn-in by typing the following::
    
    $ sumtrees.py --output=phylo.trees.sumtrees --burnin=200 phylo1.trees phylo2.trees phylo3.trees

Alternatively, we might be quite happy with the MCCT tree produced by BEAST, and in fact we want to see how the MCMC samples produced by Mr. Bayes map onto this tree (i.e., the posterior probability of the splits on the MCCT as given by the Mr. Bayes samples).
To do this, we would supply the Mr. Bayes ``.run.t``" files as the tree samples to be summarized, and use the "``-t``" or "``--target``" option to instruct SumTrees to map the posterior probabilities onto the BEAST MCMCT tree.
Thus, assuming that our Mr. Bayes runs are is in the files "``phylo.nex.run1.t``" and "``phylo.nex.run2.t``", and the BEAST summarized MCCT tree is in the file "``phylo.beast.tree``" we could type the following::
    
    $ sumtrees.py --target=phylo.beast.tree --output=phylo.mb-beast.sumtrees --burnin=200 phylo.nex.run1.t phylo2.nex.run2.t
    
Troubleshooting
===============

Prerequisites
-------------   

DendroPy is a `Python <http://www.python.org/>`_ library.
It requires and presupposes not only the existence of a Python installation on your system, but also that this Python installation is available on the system path.

The biggest problem faced by most users is not so much not having Python installed, but not having the correct version of Python installed. You can check which version of Python you have running by typing::

    $ python -V
    
SumTrees, and the DendroPy library that it is part of, works out-of-the-box with Python version 2.4 or greater, up to and including Python 2.6. 

SumTrees will not work with versions of Python prior to 2.4, such as Python 2.3. It can probably be made to work pretty easily, and if you have strong enough motiviation to use Python 2.3, it might be worth the effort for you.
It is not for me.

SumTrees (and DendroPy, and, for that matter, most existing Python code) is flat-out broken under Python 3.0.

All this can be summarized as the follows:

.. pull-quote::

    Then, shalt thou count to **2.5**.
    
    No more.     
    
    No less.     
    
    **2.5** shalt be the number thou shalt count, and the number of the counting shall be **2.5**.     
    
    **3.0** shalt thou not count, nor either count thou **2.3**, excepting that thou then proceed to **2.5**.     
    
    **4.0** is right out.

My Computer Does Not Know What a Python Is
-------------------------------------------

If you get a message like::

    python: command not found
    
it is either because Python is not installed on your system, or is not found on the system path.

SumTrees is a Python script, and, as such, you will need to have a Python interpreter installed on your system.

Otherwise, you must download and install Python 2.6 from: http://www.python.org/download/releases/2.6/.
For your convenience, the clicking on the following links should lead you directly to the appropriate pre-compiled download:

* `Mac OS X <http://www.python.org/ftp/python/2.6/python-2.6-macosx.dmg>`_
* `Microsoft Windows <http://www.python.org/ftp/python/2.6/python-2.6.msi>`_

For other platforms, the usual "``./configure``", "``make``", and "``sudo make install``" dance should get you up and running with the following:

* `Cross-platform Source <http://www.python.org/ftp/python/2.6/Python-2.6.tgz>`_

Microsoft Windows users should also refer to the `"Python Windows FAQ" <http://www.python.org/doc/faq/windows.html>`_
(http://www.python.org/doc/faq/windows.html)
after installing Python, and pay particular attention to the
`"How do I run a Python program under Windows?" <http://www.python.org/doc/faq/windows.html#id2>`_ section, as it will
help them greatly in getting Python up and running on the system path.

Manual Installation
===================

The DendroPy library is actually quite straightforward to install manually, especially if you have any familiarity with Python and how Python files are organized.
There are a couple of different things you could do:

* Add the current location of the "``dendropy``" subdirectory to your Python path environmental variable, "``$PYTHONPATH``", and place the file "``scripts\sumtrees.py``" on your system path. 

* Copy (or symlink) the "``dendropy``" directory to the "``site-packages``" directory of your Python installation, and place the file "``scripts\sumtrees.py``" on your system path. 

Repository Access
=================
The DendroPy public-access `Git <http://git-scm.com/>`_ repository can be cloned from:
    
        git://dendropy.git.sourceforge.net/gitroot/dendropy

Bugs, Suggestions, Comments, etc.
=================================

If you encounter any problems, errors, crashes etc. while using this program, please let me know at jeet@ku.edu. If you include the term "sumtrees" anywhere on the subject line (e.g. "Problem such-and-such with bootscore), it would help greatly with getting through the spam filter. Please include all the datafiles involved, as 
well the complete command used (with all the options and parameters) and the complete error message returned (simply cutting-and-pasting the terminal text should work fine).
Please feel free to contact me if you have any other questions, suggestions or comments as well.

How to Cite this Program
=========================

If you use this program in your analysis, please cite it as:

    Sukumaran, J. and Mark T. Holder. 2008. *SumTrees: Summarization of Split Support on Phylogenetic Trees. Version 1.0.2*. Part of the *DendroPy Phylogenetic Computation Library Version 2.1.3* (http://sourceforge.net/projects/dendropy).
    
Copyright, License and Warranty
===============================

SumTrees and DendroPy are: Copyright 2008 Jeet Sukumaran and Mark T. Holder.
  
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
SumTrees is part of the `DendroPy
<http://pypi.python.org/pypi/DendroPy>`_ library, which is authored by myself (`Jeet Sukumaran <http://people.ku.edu/~jeet>`_) and `Mark T. Holder <http://people.ku.edu/~mtholder>`_.

We would like to thank all the people who have contributed suggestions, bug reports and critiqes, and especially our *de facto* beta testers who contributed valuable time and trusted our program with their valuable data: `Charles W. Linkem <http://people.ku.edu/~cwlinkem>`_ and `Jamie Oaks <http://people.ku.edu/~joaks1>`_.

Portions of `DendroPy
<http://pypi.python.org/pypi/DendroPy>`_ were developed under `CIPRES <http://www.phylo.org>`_, a multi-site collaboration funded by the `NSF <http://www.nsf.gov/>`_ Information Technology Research (ITR) program grant entitled "`BUILDING THE TREE OF LIFE: A National Resource for Phyloinformatics and Computational Phylogenetics <http://www.phylo.org/about/acknowledgements>`_".

.. image:: logo_cipres.gif
    :height: 40   
    :target: http://www.phylo.org/

.. image:: nsf.gif
    :width: 40
    :target: http://www.nsf.gov/



