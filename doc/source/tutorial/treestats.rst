******************************************
Tree Statistics, Metrics, and Calculations
******************************************

Tree Length
===========

The :meth:`~dendropy.dataobject.tree.Tree.length()` method returns the sum of edge lengths of a |Tree| object, with edges that do not have any length assigned being treated as edges with length 0.
The following example shows how to identify the "critical" value for an `Archie-Faith-Cranston or PTP test <http://hymenoptera.tamu.edu/courses/ento606/Suggested%20Readings/Slowinksi_Crother_1998.pdf>`_ from a sample of |Tree| objects, i.e. a tree length equal to or greater than 95% of the trees in the sample:

.. literalinclude:: /examples/tree_length_crit.py
    :linenos:

Node Ages
=========

The :meth:`~dendropy.dataobject.tree.Tree.add_ages_to_nodes()` method calculates the age of a node (i.e., the sum of edge lengths from the node to a tip) and assigns it to a new attribute of the node: :attr:`~dendropy.dataobject.tree.Node.age`. The following example iterates through the post-burn-in of an MCMC sample of ultrametric trees, calculating the age of the MRCA of two taxa, and reports the mean age of the node.

.. literalinclude:: /examples/node_ages1.py
    :linenos:

Pybus-Harvey Gamma
==================

The Pybus-Harvey Gamma statistic is given by the :meth:`~dendropy.dataobject.tree.Tree.pybus_harvey_gamma()` instance method. The following example iterates through the post-burn-in of an MCMC sample of trees, reporting the mean Pybus-Harvey Gamma statistic:

.. literalinclude:: /examples/pbhg.py
    :linenos:

Patristic Distances
===================

The :class:`~dendropy.treecalc.PatristicDistanceMatrix` is the most efficient way to calculate the patristic distances between any pair of taxa on a tree.
Its constructor takes a |Tree| object as an argument, and the object return is callable, taking two |Taxon| objects as arguments and returning the sum of edge lengths between the two. The following example reports the pairwise distances between all taxa on the input tree:

.. literalinclude:: /examples/pdm.py
    :linenos:

Probability Under the Coalescent and Counting of Deep Coalescences
==================================================================

The :mod:`~dendropy.coalescent` module provides a range of methods for simulations and calculations under Kingman's coalescent framework and related models:

    :func:`~dendropy.coalescent.log_probability_of_coalescent_tree`
        Given a |Tree| object as the first argument, and the haploid population size as the second, returns the log probability of the |Tree| under the neutral coalescent.

    :func:`~dendropy.coalescent.num_deep_coalescences_with_fitted_tree`
        Given two |Tree| objects, a gene tree and a species tree, sharing the same leaf-set, this returns the number of deep coalescences resulting from fitting the gene tree to the species tree.

    :func:`~dendropy.coalescent.num_deep_coalescences_with_grouping`
        Given a |Tree| object as the first argument, and a list of lists of
        |Taxon| objects representing the expected monophyletic partitioning of the |TaxonSet| of the |Tree| as the second argument, this returns the number of deep coalescences found in the relationships implied by the |Tree| object, conditional on the taxon groupings given by the second argument.

Majority-Rule Consensus Tree from a Collection of Trees
=======================================================

To get the majority-rule consensus tree of a |TreeList| object, you can call the :meth:`~dendropy.dataobject.tree.TreeList.consensus()` instance method.
You can specify the frequency threshold for the consensus tree by the ``min_freq`` argument, which default to 0.5 (i.e., a 50% majority rule tree).
The following example aggregates the post-burn-in trees from four MCMC samples into a single |TreeList| object, and prints the 95% majority-rule consensus as a NEWICK string:

.. literalinclude:: /examples/majrule.py
    :linenos:

Frequency of a Split in a Collection of Trees
=============================================

The :meth:`~dendropy.dataobject.tree.TreeList.frequency_of_split()` method of a |TreeList| object returns the frequency of occurrence of a single split across all the |Tree| objects in the |TreeList|.
The split can be specified by passing a split bitmask directly using the ``split_bitmask`` keyword argument, as a list of |Taxon| objects using the ``taxa`` keyword argument, or as a list of taxon labels using the ``labels`` keyword argument.
The following example shows how to calculate the frequency of a split defined by two taxa, "Morelia amethistina" and "Morelia tracyae", from the post-burn-in trees aggregated across four MCMC samples:

.. literalinclude:: /examples/splitfreq.py
    :linenos:

Tree Distances
==============

The :mod:`~dendropy.treecalc` module provides a number of functions to calculate the distance between two trees passed as arguments:

    :func:`~dendropy.treecalc.symmetric_distance`
        This function returns the symmetric distance between two trees. The symmetric distance between two trees is the sum of splits found in one of the trees but not the other. It is common to see this statistic called the "Robinson-Foulds distance", but in DendroPy we reserve this term to apply to the Robinson-Foulds distance in the strict sense, i.e., the weighted symmetric distance (see below).

    :func:`~dendropy.treecalc.euclidean_distance`
        This function returns the "branch length distance" of Felsenstein (2004), i.e. the sum of absolute differences in branch lengths for equivalent splits between two trees, with the branch length for a missing split taken to be 0.0.

    :func:`~dendropy.treecalc.robinson_foulds_distance`
        This function returns the Robinsons-Foulds distance between two trees, i.e., the sum of the square of differences in branch lengths for equivalent splits between two trees, with the branch length for a missing split taken to be 0.0.
