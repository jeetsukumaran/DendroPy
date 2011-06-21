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

The :meth:`~dendropy.dataobject.tree.Tree.calc_node_ages()` method calculates the age of a node (i.e., the sum of edge lengths from the node to a tip) and assigns it to the :attr:`~dendropy.dataobject.tree.Node.age` attribute. The following example iterates through the post-burn-in of an MCMC sample of ultrametric trees, calculating the age of the MRCA of two taxa, and reports the mean age of the node.

.. literalinclude:: /examples/node_ages1.py
    :linenos:

Pybus-Harvey Gamma
==================

The Pybus-Harvey Gamma statistic is given by the :meth:`~dendropy.dataobject.tree.Tree.pybus_harvey_gamma()` instance method. The following example iterates through the post-burn-in of an MCMC sample of trees, reporting the mean Pybus-Harvey Gamma statistic:

.. literalinclude:: /examples/pbhg.py
    :linenos:

Patristic Distances
===================

The :class:`~dendropy.treecalc.PatristicDistanceMatrix` is the most efficient way to calculate the patristic distances between taxa or leaves on a tree, when doing multiple such calculations.
Its constructor takes a |Tree| object as an argument, and the object return is callable, taking two |Taxon| objects as arguments and returning the sum of edge lengths between the two. The following example reports the pairwise distances between all taxa on the input tree:

.. literalinclude:: /examples/pdm.py
    :linenos:

Probability Under the Coalescent Model
=======================================

The :mod:`~dendropy.coalescent` module provides a range of methods for simulations and calculations under Kingman's coalescent framework and related models. For example:

    :func:`~dendropy.coalescent.log_probability_of_coalescent_tree`
        Given a |Tree| object as the first argument, and the haploid population size as the second, returns the log probability of the |Tree| under the neutral coalescent.

    :func:`~dendropy.coalescent.kl_divergence_coalescent_trees`
        Reports the Kullback-Leilber divergence of a list of trees from the theoretical distribution of neutral coalescent trees. Requires the `de Hoon statistics package <http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/python/Statistics>`_ package to be installed.

Numbers of Deep Coalescences
============================

    :func:`~dendropy.reconcile.reconciliation_discordance`
        Given two |Tree| objects *sharing the same leaf-set*, this returns the number of deep coalescences resulting from fitting the first tree (e.g., a gene tree) to the second (e.g., a species tree). This is based on the algorithm described `Goodman, et al. <bioinformatics.oxfordjournals.org/cgi/reprint/14/9/819.pdf>`_ (Goodman, et al., 1979. Fitting the gene lineage into its species lineage,a parsimony strategy illustrated by cladograms constructed from globin sequences. Syst. Zool. 19: 99-113).

    .. versionchanged:: 3.3.0
        Renamed and moved to :mod:`~dendropy.reconcile` module.

    :func:`~dendropy.reconcile.monophyletic_partition_discordance`
        Given a |Tree| object as the first argument, and a list of lists of
        |Taxon| objects representing the expected monophyletic partitioning of the |TaxonSet| of the |Tree| as the second argument, this returns the number of deep coalescences found in the relationships implied by the |Tree| object, conditional on the taxon groupings given by the second argument. This statistic corresponds to the Slatkin and Maddison (1989) **s** statistic, as described `here <http://mesquiteproject.org/Mesquite_Folder/docs/mesquite/popGen/popGen.html#s>`_.

    .. versionchanged:: 3.3.0
        Renamed and moved to :mod:`~dendropy.reconcile` module.

Number of Deep Coalescences when Embedding One Tree in Another (e.g. Gene/Species Trees)
========================================================================================

Imagine we wanted to generate the distribution of the number of deep coalescences under two scenarios: one in which a population underwent sequential or step-wise vicariance, and another when there was simultaneous fragmentation.
In this case, the containing tree and the embedded trees have different leaf sets, and there is a many-to-one mapping of embedded tree taxa to containing tree taxa.

The :class:`~dendropy.reconcile.ContainingTree` class is designed to allow for counting deep coalescences in cases like this.
It requires a |TaxonSetMapping| object, which provides an association between the embedded taxa and the containing taxa.
The easiest way to get a |TaxonSetMapping| object is to call the special factory function :meth:`~dendropy.dataobject.taxon.TaxonSetMapping.create_contained_taxon_mapping()`.
This will create a new |TaxonSet| to manage the gene taxa, and create the associations between the gene taxa and the containing tree taxa for you.
It takes two arguments: the |TaxonSet| of the containing tree, and the number of genes you want sampled from each species.
If the gene-species associations are more complex, e.g., different numbers of genes per species, we can pass in a list of values as the second argument to `~dendropy.dataobject.taxon.TaxonSetMapping.create_contained_taxon_mapping()`.
This approach should be used with caution if we cannot be certain of the order of taxa (as is the case with data read in Newick formats). In these case, and in more complex cases, we might need to directly instantiate the :class:`~dendropy.dataobject.taxon.TaxonSetMapping` object. The API to describe the associations when constructing this object is very similar to that of the :class:`~dendropy.dataobject.taxon.TaxonSetPartition` object: you can use a function, attribute or dictionary.

The :class:`~dendropy.reconcile.ContainingTree` class has its own native contained coalescent simulator, :meth:`~dendropy.reconcile.ContainingTree.embed_contained_kingman()`, which simulates *and* embeds a contained coalescent tree at the same time.

.. literalinclude:: /examples/sim_and_count_deepcoal2.py

If you have used some other method to simulate your trees, you can use :meth:`~dendropy.reconcile.ContainingTree.embed_tree()` to embed the trees and count then number of deep coalescences.

.. literalinclude:: /examples/sim_and_count_deepcoal1.py

For more details on simulating contained coalescent trees and counting numbers of deep coalescences on them, see ":ref:`Simulating_Contained_Coalescent_Trees`" or ":ref:`Simulating_and_Counting_Deep_Coalescences`".

Majority-Rule Consensus Tree from a Collection of Trees
=======================================================

To get the majority-rule consensus tree of a |TreeList| object, you can call the :meth:`~dendropy.dataobject.tree.TreeList.consensus()` instance method.
You can specify the frequency threshold for the consensus tree by the ``min_freq`` argument, which default to 0.5 (i.e., a 50% majority rule tree).
The following example aggregates the post-burn-in trees from four MCMC samples into a single |TreeList| object, and prints the 95% majority-rule consensus as a Newick string:

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

Native |Tree| Methods
---------------------

.. versionadded:: 3.2

The |Tree| class provides methods for calculating distances between two trees:


    :meth:`~dendropy.dataobject.tree.Tree.symmetric_difference`
        This method returns the symmetric distance between two trees. The symmetric distance between two trees is the sum of the number of splits found in one of the trees but not the other. It is common to see this statistic called the "Robinson-Foulds distance", but in DendroPy we reserve this term to apply to the Robinson-Foulds distance in the strict sense, i.e., the weighted symmetric distance (see below).

    :meth:`~dendropy.dataobject.tree.Tree.false_positives_and_negatives`
        This method returns a tuple pair, with the first element the number of splits in the current |Tree| object *not* found in the |Tree| object to which it is being compared, while the second element is the number of splits in the second |Tree| object that are not in the current |Tree|. The sum of these two elements is exactly equal to the value reported by :meth:`~dendropy.dataobject.tree.Tree.symmetric_distance`.

    :meth:`~dendropy.dataobject.tree.Tree.euclidean_distance`
        This method returns the "branch length distance" of Felsenstein (2004), i.e. the sum of absolute differences in branch lengths for equivalent splits between two trees, with the branch length for a missing split taken to be 0.0.

    :meth:`~dendropy.dataobject.tree.Tree.robinson_foulds_distance`
        This method returns the Robinsons-Foulds distance between two trees, i.e., the sum of the square of differences in branch lengths for equivalent splits between two trees, with the branch length for a missing split taken to be 0.0.

For example::

    >>> import dendropy
    >>> s1 = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247)"
    >>> s2 = "((t5:2.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247)"
    >>> tree1 = dendropy.Tree.get_from_string(s1, 'newick')
    >>> tree2 = dendropy.Tree.get_from_string(s2, 'newick')
    >>> tree1.symmetric_difference(tree2)
    0
    >>> tree1.false_positives_and_negatives(tree2)
    (0, 0)
    >>> tree1.euclidean_distance(tree2)
    2.2232636377544162
    >>> tree1.robinson_foulds_distance(tree2)
    2.971031

Using the :mod:`~dendropy.treecalc` Module
------------------------------------------

The :mod:`~dendropy.treecalc` module provides for these operations as independent functions that take two |Tree| objects as arguments.
These independent functions require that both trees have the same |TaxonSet| reference, otherwise an exception is raised::

        >>> import dendropy
        >>> from dendropy import treecalc
        >>> s1 = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247)"
        >>> s2 = "((t5:2.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247)"
        >>> tree1 = dendropy.Tree.get_from_string(s1, 'newick')
        >>> tree2 = dendropy.Tree.get_from_string(s2, 'newick')
        >>> treecalc.symmetric_difference(tree1, tree2)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "treecalc.py", line 240, in symmetric_difference
            t = false_positives_and_negatives(tree1, tree2)
          File "treecalc.py", line 254, in false_positives_and_negatives
            % (hex(id(reference_tree.taxon_set)), hex(id(test_tree.taxon_set))))
        TypeError: Trees have different TaxonSet objects: 0x10111ec00 vs. 0x10111eaa0
        >>> treecalc.euclidean_distance(tree1, tree2)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "treecalc.py", line 236, in euclidean_distance
            value_type=value_type)
          File "treecalc.py", line 160, in splits_distance
            % (hex(id(tree1.taxon_set)), hex(id(tree2.taxon_set))))
        TypeError: Trees have different TaxonSet objects: 0x10111ec00 vs. 0x10111eaa0
        >>> treecalc.robinson_foulds_distance(tree1, tree2)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "treecalc.py", line 223, in robinson_foulds_distance
            value_type=float)
          File "treecalc.py", line 160, in splits_distance
            % (hex(id(tree1.taxon_set)), hex(id(tree2.taxon_set))))
        TypeError: Trees have different TaxonSet objects: 0x10111ec00 vs. 0x10111eaa0
        >>> tree3 = dendropy.Tree.get_from_string(s2, 'newick', taxon_set=tree1.taxon_set)
        >>> treecalc.symmetric_difference(tree1, tree3)
        0
        >>> treecalc.euclidean_distance(tree1, tree3)
        2.2232636377544162
        >>> treecalc.robinson_foulds_distance(tree1, tree3)
        2.971031
