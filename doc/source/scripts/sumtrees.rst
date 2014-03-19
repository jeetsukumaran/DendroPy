########################################################
SumTrees: Phylogenetic Tree Summarization and Annotation
########################################################

Introduction
============

`SumTrees <sumtrees.html>`_ is a program by |js|_ and |mth|_ to summarize non-parameteric bootstrap or Bayesian posterior probability support for splits or clades on phylogenetic trees.

The basis of the support assessment is typically given by a set of non-parametric bootstrap replicate tree samples produced by programs such as GARLI or RAxML, or by a set of MCMC tree samples produced by programs such as Mr. Bayes or BEAST.
The proportion of trees out of the samples in which a particular split is found is taken to be the degree of support for that split as indicated by the samples.
The samples that are the basis of the support can be distributed across multiple files, and a burn-in option allows for an initial number of trees in each file to be excluded from the analysis if they are not considered to be drawn from the true support distribution.

The support for the splits will be mapped onto one or more target trees either in terms of node labels or branch lengths.
The target trees can be supplied by yourself, or, if no target trees are given, then a majority-rule clade consensus tree will be constructed based on the samples given.
In the latter case, you have the option of specifying the minimum posterior probability or proportional frequency threshold for a clade to be included on the consensus tree.

.. versionadded:: 3.3.1 (with DendroPy 3.8.0)

By default SumTrees will now provide summaries of edge lengths (i.e., mean, median, standard deviation, range, 95% HPD, 5% and 95% quantiles, etc.) as special node comments. These can be visualized in `FigTree <http://tree.bio.ed.ac.uk/software/figtree/>`_ by, for example, checking "Node Labels", then selecting one of "length_mean", "length_median", "length_sd", "length_hpd95", etc. If the trees are ultrametric and the "``--ultrametric``" flag is used, or edge lengths are set so that node ages on the output trees correspond to be mean or median of the node ages of the input trees ("``--edges=mean-age``" or "``--edges=median-age``"), then node ages will be summarized as well. In all cases, the flag "``--no-summary-metadata``" will suppress calculation and output of these summaries.

Where to Find the Package
=========================

SumTrees is distributed and installed as part of the |DendroPy|_ phylogenetic computing library, which has its homepage here:

        |dendropy_homepage_url|

with the main download page here:

        |dendropy_download_url|

How to Install the Package
==========================

The Easy Way
------------

Simply type the following command in your shell/terminal::

    sudo easy_install -U dendropy

And that is all it takes!

This requires you have |setuptools|_ already installed. If you get an error message about "``easy_install``" not being found, then you do *not* have |setuptools|_ installed, and you need to install it first.
It is *highly recommended* that you use this method, even if it means that you have to pre-install |setuptools|_.

The Not So Easy Way
-------------------

Alternatively, you can install |DendroPy|_ yourself, by following these steps.

1.  Download the latest release of |DendroPy|_ from:

    |dendropy_source_archive_url|

#.  Expand the downloaded archive

    This step varies depending on the operating system and the particular programs that you have installed.
    In most cases, simply double-clicking on the file that you have downloaded should kick off the process.
    Otherwise, open a terminal shell window and go to the directory in which you have downloaded the archive, and type:

        .. parsed-literal::

            tar -xvzf DendroPy-|version|.tar.gz

    For example, say you saved the downloaded file on your desktop.
    Then, opening up the terminal and entering the following commands will take you to your Desktop and expand the archive:

        .. parsed-literal::

            $ cd ~/Desktop
            $ tar -xvzf DendroPy-|version|.tar.gz

    One way or another, you should end up with a directory called "DendroPy-|version|" or something similar, which contains the entire Dendropy package.

#.  Install the library

    In the terminal shell, go to the directory of the package that you have just archived and type "``sudo python setup.py install``".
    Continuing the example from above:

        .. parsed-literal::

            $ cd DendroPy-|version|
            $ sudo python setup.py install

    The library installation will automatically create an executable script called "``sumtrees.py``" and place it on your system path for you, so that you can call it from anywhere.

#.  If the installation was successful, you should be able to type "``sumtrees.py``" in the terminal window and see something like the following (with possibly a different date or version number):

        .. parsed-literal::

            ======================================================================
            SumTrees - Phylogenetic Tree Split Support Summarization
            Version 3.0.0 (Aug 22 2010)
            By Jeet Sukumaran and Mark T. Holder
            (using the DendroPy Phylogenetic Computing Library Version |version|)
            ======================================================================

            SumTrees: No sources of input trees specified. Please provide the path to at
                      least one (valid and existing) file containing tree samples to
                      summarize. See '--help' for other options.

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
Calculate support for nodes on a specific tree, "``best.tre``" as given by a set of tree files, with support reported as percentages rounded to integers, and saving the result to "``result.tre``"::

    $ sumtrees.py --decimals=0 --percentages --output=result.tre --target=best.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -d0 -p -o result.tre -t best.tre treefile1.tre treefile2.tre treefile3.tre

Summarization of Posterior Probabilities of Clades with a Consensus Tree
------------------------------------------------------------------------
Summarize a set of tree files using a 95% rule consensus tree, with support for clades indicated as proportions (posterior probabilities) and branch lengths the mean across all trees, dropping the first 200 trees in each file as a burn-in, and saving the result to "``result.tre``"::

    $ sumtrees.py --min-clade-freq=0.95 --burn-in=200 --support-as-labels --output=result.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -f0.95 -b200 -l -o result.tre treefile1.tre treefile2.tre treefile3.tre

Set Node Ages of Consensus or Target Tree(s) to Mean/Median Node Age of Input Tree
----------------------------------------------------------------------------------

.. versionadded:: DendroPy 3.8.0 / SumTrees 3.3.0

Summarize a set of ultrametric tree files using a 95% rule consensus tree, with support for clades indicated as proportions (posterior probabilities) and branch lengths the mean across all trees, dropping the first 200 trees in each file as a burn-in, with node ages of the consensus tree set to the mean node ages of the input trees::

    $ sumtrees.py --min-clade-freq=0.95 --burn-in=200 --support-as-labels --edges=mean-age --output=result.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -f0.95 -b200 -o result.tre -l -e mean-age treefile1.tre treefile2.tre treefile3.tre

To use the median age instead::

    $ sumtrees.py --min-clade-freq=0.95 --burn-in=200 --support-as-labels --edges=median-age --output=result.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -f0.95 -b200 -o result.tre -l -e median-age treefile1.tre treefile2.tre treefile3.tre

Running in Parallel Mode
------------------------

.. versionadded:: DendroPy 3.6.0 / SumTrees 3.0.0

.. note::

    This feature is only available when running under Python 2.6 of greater.

Starting with DendroPy version 3.6 (SumTrees version 3.0), and when running Python 2.6 or greater, you can run SumTrees in parallel mode.
This will analyze each input source in its own independent process, with multiple processes running in parallel.
Multiprocessing analysis is invoked by adding the "``-m``" or "``--multiprocessing``"  flag to the SumTrees command, and passing in the maximum number of processes to run in parallel.
For example, if your machine has two cores, and you want to run the previous analyses using both of them, you would specify that SumTrees run in parallel mode with two processes by adding "``-m2``" or "``--multiprocessing=2``" to the SumTrees command invocation::

    $ sumtrees.py --multiprocessing=2 --decimals=0 --percentages --output=result.tre --target=best.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -m2 -d0 -p -o result.tre -t best.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py --multiprocessing=2 --min-clade-freq=0.95 --burn-in=200 --support-as-labels --output=result.tre treefile1.tre treefile2.tre treefile3.tre
    $ sumtrees.py -m2 -f0.95 -b200 -l -o result.tre treefile1.tre treefile2.tre treefile3.tre

You can specify as many processes as you want, up to the total number of tree support files passed as input sources.
If you specify fewer processes than input sources, then the files will be cycled through the processes.

Tutorials and Examples
======================

At its most basic, you will need to supply SumTrees with the path to one or more tree files in Newick or NEXUS format that you want to summarize::

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

More extended options specify things like: where to save the output (by default it goes to the screen), the topology or tree to which to map the support (user-supplied or consensus tree), the output format (NEXUS by default, but it can also be Newick), whether support is indicated in terms of proportions or percentages etc.
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

However as, this is a BEAST analysis, the trees are going to be rooted and ultrametric. We can tell SumTrees this by passing it the "``--ultrametric``" flag::

    $ sumtrees.py --ultrametric phylo.trees

This will result in node ages being summarized as well as edge lengths.

Of course, we want to discard the first few samples of trees, as these were probably not drawn in frequencies in proportion to the stationary distribution of the chain.
To do this::

    $ sumtrees.py --ultrametric --burnin=200 phylo.trees

The above command will cause SumTrees to ignore the first 200 trees it finds in the file for all its calculations.

Again, instead of displaying the tree to the screen we can save it directly to a file, either by redirecting the screen output to a file::

    $ sumtrees.py --ultrametric --burnin=200 phylo.trees > phylo.trees.sumtrees

or by using the "``-o``" or "``--output``" option::

    $ sumtrees.py --ultrametric --output=phylo.trees.sumtrees --burnin=200 phylo.trees

We might also have split up our analysis into multiple independent runs, resulting in multiple MCMC tree sample files (e.g., "``phylo1.trees``", "``phylo2.trees``" and "``phylo3.trees``").
We can ask SumTrees to summarize posterior probability from across all these runs, treating the first 200 trees in *each* sample file as a burn-in by typing the following::

    $ sumtrees.py --ultrametric --output=phylo.trees.sumtrees --burnin=200 phylo1.trees phylo2.trees phylo3.trees

Alternatively, we might be quite happy with the MCCT tree produced by BEAST, and in fact we want to see how the MCMC samples produced by Mr. Bayes map onto this tree (i.e., the posterior probability of the splits on the MCCT as given by the Mr. Bayes samples).
To do this, we would supply the Mr. Bayes ``.run.t``" files as the tree samples to be summarized, and use the "``-t``" or "``--target``" option to instruct SumTrees to map the posterior probabilities onto the BEAST MCMCT tree.
Thus, assuming that our Mr. Bayes runs are is in the files "``phylo.nex.run1.t``" and "``phylo.nex.run2.t``", and the BEAST summarized MCCT tree is in the file "``phylo.beast.tree``" we could type the following::

    $ sumtrees.py --target=phylo.beast.tree --output=phylo.mb-beast.sumtrees --burnin=200 phylo.nex.run1.t phylo2.nex.run2.

Summarizing Rooted and Ultrametric Trees
----------------------------------------

.. versionadded:: DendroPy 3.8.0 / SumTrees 3.3.0

SumTrees treats all trees as unrooted unless specified otherwise. You can force SumTrees to treat all trees as rooted by passing it the "``--rooted``" flag::

    $ sumtrees.py --rooted phylo.trees

If the trees are rooted **and** ultrametric, the "``--ultrametric``" flag will result in SumTrees summarizing node age information as well::

    $ sumtrees.py --ultrametric phylo.trees

Summarizing Edge Lengths and Node Ages
--------------------------------------

.. versionadded:: DendroPy 3.8.0 / SumTrees 3.3.0

When constructing a consensus tree onto which to map support, by default SumTrees sets the edge lengths of the consensus tree to the mean of the lengths of corresponding edges of the input trees.
However, if the "``--ultrametric``" flag is given (which requires the input trees to be rooted and ultrametric), then by default SumTrees adjusts the edge lengths of the consensus tree such that the ages of the subtended nodes are equal to the median of the ages of the corresponding nodes of the input trees.
If target trees are given, then SumTrees will *not* change the original target tree edges unless instructed otherwise.

You can explicitly request an alternate edge summarization strategy using the "``-e``"/"``--edges``" flag. This will result in the edge lengths of all target trees being adjusted, whether the target tree is the consensus tree constructed by SumTrees, or one or more trees specified by the "``-t``"/"``--target``" flag.

The "``-e``"/"``--edges``" flag can take one of the following values:

        - ``mean-length``: sets the edge lengths of the target/consensus tree(s) to the mean of the lengths of the corresponding edges of the input trees.
        - ``median-length``: sets the edge lengths of the target/consensus tree(s) to the median of the lengths of the corresponding edges of the input trees.
        - ``median-age``: adjusts the edge lengths of the target/consensus tree(s) such that the node ages correspond to the median age of corresponding nodes of the input trees [requires rooted ultrametric trees].
        - ``mean-age``: adjusts the edge lengths of the target/consensus tree(s) such that the node ages correspond to the mean age of corresponding nodes of the input trees [requires rooted ultrametric trees].

So, for example, to construct a consensus tree of a post-burnin set of ultrametric trees, with the node ages set to the *mean* instead of the median node age::

    $ sumtrees.py --edges=mean-age --burnin=200 beast1.trees beast2.trees beast3.trees
    $ sumtrees.py --e mean-age --b 200 beast1.trees beast2.trees beast3.trees

Or to set the edges of a user-specifed tree to the median edge length of the input trees::

    $ sumtrees.py --edges=median-length --target=mle.tre boots1.tre boots2.tre
    $ sumtrees.py --e median-length -t mle.tre boots1.tre boots2.tre

Parallelizing SumTrees
----------------------

.. versionadded:: 3.6

.. note::

    This feature is only available when running under Python 2.6 of greater.

Basics
^^^^^^

DendroPy version 3.6 (SumTrees version 3.0) added support for running multiple parallel processes when running under Python 2.6 or greater.

If you have multiple input (support) files, you can greatly increase the performance of SumTrees by running it in parallel mode.
In parallel mode, each input source will be handled in a separate process, resulting in a speed-up linearly proportional to the total number of processes running in parallel.
At its most basic, running in parallel mode involves nothing more than adding the "``-m``" or "``--multiprocessing`` option to the SumTrees invocation, and passing in the number of parallel processes to run.
So, for example, if you have four tree files that you want to summarize, and you want to run these using two processes in parallel::

    $ sumtrees.py --multiprocessing=2 -o result.tre phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre
    $ sumtrees.py -m 2 -o result.tre phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

Or, to run in four processes simultaneously::

    $ sumtrees.py --multiprocessing=4 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre
    $ sumtrees.py -m 4 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

Other options as described above, can, of course be added as needed::

    $ sumtrees.py --multiprocessing=4 --burnin=200 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre
    $ sumtrees.py -m 4 -b 200 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

    $ sumtrees.py --multiprocessing=4 --burnin=200 --min-clade-freq=0.75 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre
    $ sumtrees.py -m 4 -b 200 -f0.75 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

    $ sumtrees.py --multiprocessing=4 --burnin=200 --min-clade-freq=0.75 --output=con.tre phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre
    $ sumtrees.py -m 4 -b 200 -f0.75 -o con.tre phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

Parallelization Strategy: Deciding on the Number of Processes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In parallel mode, SumTrees parallelizes by input files, which means that the maximum number of processes that can be run is limited to the number of input files.
If you have four support tree files as input sources, as in the above example ("``phylo.run1.tre``", "``phylo.run2.tre``","``phylo.run3.tre``", and "``phylo.run4.tre``"), then you can run from 2 up to 4 processes in parallel.
If you have 40 tree source or support files, then you can run from 2 up to 40 processes in parallel.

Note that there is a difference between the number of *processes* that SumTrees runs, and the number of *processors* or cores available on your machine or in your given hardware context.
When running SumTrees in parallel mode, you can specify *any* number of parallel *processes*, up to the number of support files, even if the number of processes greatly exceeds the number *processors* available.
For example, you might have an octo-core machine available, which means that you have 8 processors available.
If your analysis has 40 independent support files, you can invoke SumTrees with 40 parallel processes::

    $ sumtrees.py --multiprocessing=40 -o result.tre t1.tre ... t40.tre

SumTrees will indeed launch 40 parallel processes to process the files, and it will seem like 40 processes are being executed in parallel.
In reality, however, the operating system is actually cycling the processes through the available processors in rapid succession, such that, on a nanosecond time-scale, only 8 processes are actually executing simultaneously: one on each of the available cores.
While your run should complete without problems when oversubscribing the hardware in this way, there is going to be some degree of performance hit due to the overhead involved in managing the cycling of the processes through the processors.
On some operating systems and hardware contexts, depending on the magnitude of oversubscription, this can be considerable.
Thus it probably is a good idea to match the number of processes to the number of processors available.

Another issue to consider is an even distribution of workload.
Assuming that each of your input support files have the same number of trees, then it makes sense to specify a number of processes that is a factor of the number of input files.
So, for example, if you have 8 input files to be summarized, you will get the best performance out of SumTrees by specifying 2, 4, or 8 processes, with the actual number given by the maximum number of processors available or that you want to dedicate to this task.

Running Parallel-Mode SumTrees in a Parallel Environment on a High-Performance Computing (HPC) Cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unfortunately, the diversity and idiosyncracies in various HPC configurations and set-ups are so great that it is very difficult to provide general recipes on how to run parallel-mode SumTrees in a parallel environment on an HPC cluster.
However, in all cases, all you really need to do is to set up an appropriate parallel execution environment on the cluster, requesting a specific number of processors from the cluster scheduler software, and then tell SumTrees to run the same number of processes.

For example, if your cluster uses the `Sun Grid Engine <http://gridengine.sunsource.net/>`_ as its scheduler, then you might use a job script such as the following:

    .. parsed-literal::

        # /bin/sh
        #$ -cwd
        #$ -V
        #$ -pe mpi 4
        sumtrees.py -m4 -o result.con.tre -burnin 100 mcmc1.tre mcmc2.tre mcmc3.tre mcmc4.tre

The "``#$ -pe mpi 4``" line tells the SGE scheduler to allocate four processors to the job in the "``mpi``" parallel environment, while the "``-m4``" part of the SumTrees command tells SumTrees to run 4 processes in parallel.

If you are using PBS/Torque as a scheduler, the equivalent job script might be:

    .. parsed-literal::

        # /bin/sh
        #PBS -l nodes=2:ncpus=2
        cd $PBS_O_WORKDIR
        sumtrees.py -m4 -o result.con.tre -burnin 100 mcmc1.tre mcmc2.tre mcmc3.tre mcmc4.tre

Here, the "``#PBS -l nodes=2:ncpus=2``" requests two processors on two nodes, while the "``-m4``" part of the SumTrees command, as before, tells SumTrees to run 4 processes in parallel.

The particular job scripts you use will almost certainly be different, varying with the cluster job/load management software, scheduler and computing resources.
Apart from the parallel environment name, number of processors and/or machine configuration, you might also need to provide a queue name, a maximum run time limit, and a soft or hard memory limit.
None of it should make any difference to how SumTrees is actually invoked: you would still just use the "``-m``" or "``--multiprocessing``" flags to specify that SumTrees runs in parallel mode with a particular number of parallel processes).
You just need to check with your cluster administrators or documentation to make sure your job script or execution context provides sufficient processors to match the number of processes that you want run (as well as other resources, e.g. wall time limit and memory limit, that SumTrees needs to finish its job).

Improving Performance
=====================

    * Run SumTrees in parallel mode (see "`Running in Parallel Mode`_" or "`Parallelizing SumTrees`_")::

        $ sumtrees.py -m 4 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

    * Reduce the tree-processing logging frequency.
      By default, SumTrees will report back every 500th tree in a tree file, just to let you know where it is and to give you a sense of how long it will take. You can use the "``-g``"  or "``--log-frequency``" flag to control this behavior. If you have very large files, you may be content to have SumTrees report back every 5000th or even every 10000th tree instead::

        $ sumtrees.py --log-frequency=10000 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre
        $ sumtrees.py -g10000 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

      If you are content to let SumTrees run without reporting its progress *within* each file (SumTrees will still report back whenever it begins or ends working on a file), then you can switch off tree processing logging altogether by specifying a logging frequency of 0::

        $ sumtrees.py --log-frequency=0 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre
        $ sumtrees.py -g0 phylo.run1.tre phylo.run2.tre phylo.run3.tre phylo.run4.tre

Troubleshooting
===============

Prerequisites
-------------

|DendroPy|_ is a |Python|_ library.
It requires and presupposes not only the existence of a |Python|_ installation on your system, but also that this Python installation is available on the system path.

The biggest problem faced by most users is not so much not having |Python|_ installed, but not having the correct version of Python installed. You can check which version of Python you have running by typing::

    $ python -V

SumTrees, and the |DendroPy|_ library that it is part of, works out-of-the-box with any 2.x Python version 2.4 or greater.

SumTrees will not work with versions of |Python|_ prior to 2.4, such as |Python|_ 2.3. It can probably be made to work pretty easily, and if you have strong enough motiviation to use Python 2.3, it might be worth the effort for you.
It is not for me.

SumTrees (and |DendroPy|) is currently not compatible with Python 3.

My Computer Does Not Know What a Python Is
-------------------------------------------

If you get a message like::

    python: command not found

it is either because |Python|_  is not installed on your system, or is not found on the system path.

SumTrees is a |Python|_ script, and, as such, you will need to have a |Python|_  interpreter installed on your system.

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

The |DendroPy|_ library is actually quite straightforward to install manually, especially if you have any familiarity with Python and how Python files are organized.
There are a couple of different things you could do:

* Add the current location of the "``dendropy``" subdirectory to your Python path environmental variable, "``$PYTHONPATH``", and place the file "``scripts/sumtrees.py``" on your system path.

* Copy (or symlink) the "``dendropy``" directory to the "``site-packages``" directory of your Python installation, and place the file "``scripts/sumtrees.py``" on your system path.

Repository Access
=================
The |DendroPy|_ public-access |Git|_ repository can be cloned from:

        |dendropy_public_repo_url|

Bugs, Suggestions, Comments, etc.
=================================

If you encounter any problems, errors, crashes etc. while using this program, please let me (|js|_) know at jeet@ku.edu. If you include the term "sumtrees" anywhere on the subject line (e.g. "Problem such-and-such with bootscore), it would help greatly with getting through the spam filter. Please include all the datafiles involved, as
well the complete command used (with all the options and parameters) and the complete error message returned (simply cutting-and-pasting the terminal text should work fine).
Please feel free to contact me if you have any other questions, suggestions or comments as well.

How to Cite this Program
=========================

If you use this program in your analysis, please cite it as:

    |dendropy_citation|

In the text of your paper, if you want to look like you know what you are doing, you should probably also mention explicitly that you specifically used the SumTrees program of the |DendroPy|_ package, as well as the particular version numbers of SumTrees and |DendroPy|_ that you used.

.. include:: ../license.inc
.. include:: ../acknowledgements.inc

