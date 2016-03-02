****************************************************************
Tree Statistics, Metrics, Summarizations, and Other Calculations
****************************************************************

Some general tree metrics that are calculated without reference to any particular model or data and general report some tree metadata (e.g., tree length, node ages, etc.) are available as instance methods.
More specialized tree statistics, however, are available through functions in various other modules:

-   The :mod:`~dendropy.calculate.treemeasure` module provides for calculation of statistics that are typically calculated on a single tree.
-   The :mod:`~dendropy.calculate.treecompare` module provides for calculation of statistics that are typically calculated between trees
-   The :mod:`~dendropy.calculate.treescore` module provides for statistics that typically score a tree under a model and with reference to some sort of data.
-   The :mod:`~dendropy.model.coalescent` module provides for calcuations on trees under the coalescent model.

In addition, see the :doc:`PhylogeneticDistanceMatrix class </primer/phylogenetic_distances>` for statistics, operations, and inferences based on phylogenetic (taxon-to-taxon) distances, including rapid calculation of MRCA's, patristics distances, neighbor-joining (NJ) and Unweighted Pair Group with Mathematical Average (UPGMA) trees, phylogenetic community statistics (such as the Mean Pairwise Distance, MPD, or Mean Nearest Taxon Distance, MNTD), and more.

Native Tree Statistic and Metric Methods
========================================

Basic meta-information about tree structure are available as native |Tree| methods.

Tree Length
-----------

The :meth:`~dendropy.datamodel.treemodel.Tree.length()` method returns the sum of edge lengths of a |Tree| object, with edges that do not have any length assigned being treated as edges with length 0.
The following example shows how to identify the "critical" value for an `Archie-Faith-Cranston or PTP test <http://hymenoptera.tamu.edu/courses/ento606/Suggested%20Readings/Slowinksi_Crother_1998.pdf>`_ from a sample of |Tree| objects, i.e. a tree length equal to or greater than 95% of the trees in the sample:

.. literalinclude:: /examples/tree_length_crit.py

Node Ages
---------

The :meth:`~dendropy.datamodel.treemodel.Tree.calc_node_ages()` method calculates the age of a node (i.e., the sum of edge lengths from the node to a tip) and assigns it to the :attr:`~dendropy.datamodel.treemodel.Node.age` attribute. The following example iterates through the post-burnin of an MCMC sample of ultrametric trees, calculating the age of the MRCA of two taxa, and reports the mean age of the node.

.. literalinclude:: /examples/node_ages1.py

Number of Lineages at a Particular Time and Lineage Through Time Plots
----------------------------------------------------------------------

The :meth:`~dendropy.datamodel.treemodel.Tree.num_lineages_at()` method of the |Tree| class returns the number of lineages at a particular time given in terms of distance from the root.
The following example extracts the number of lineages at fixed intervals along the length of the tree to use in an Lineage Through Time (LTT) plot:

.. literalinclude:: /examples/ltt.py


Unary Tree Statistics and Metrics
=================================

Numerous specialized statistics and indexes of tree shape and structure (B1, Colless' imbalance, Pybus-Harvey-Gamma, etc.) are available through the :mod:`~dendropy.calculate.treemeasure` module:


.. literalinclude:: /examples/treemeasures1.py

Pybus-Harvey Gamma
------------------

The Pybus-Harvey Gamma statistic is given by the :meth:`~dendropy.datamodel.treemodel.Tree.pybus_harvey_gamma()` instance method. The following example iterates through the post-burn-in of an MCMC sample of trees, reporting the mean Pybus-Harvey Gamma statistic:

.. literalinclude:: /examples/pbhg.py

Patristic Distances
-------------------

The :class:`~dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix` is the most efficient way to calculate the patristic distances between taxa or leaves on a tree, when doing multiple such calculations.
The easiest way to get an object of this class for a particular tree is to call :meth:`~dendropy.datamodel.treemodel.Tree.phylogenetic_distance_matrix()`.
The object is callable, taking two |Taxon| objects as arguments and returning the sum of edge lengths between the two. The following example reports the pairwise distances between all taxa on the input tree:

.. literalinclude:: /examples/pdm.py

Note that the |PhylogeneticDistanceMatrix| object does not automatically update if the original |Tree| changes: it is essentially a snapshot of |Tree| at the point in which it is instantiated.
If the original |Tree| changes, you should create a new instance of the corresponding |PhylogeneticDistanceMatrix| object.

Comparing and Summarizing Trees
===============================

Distances Between Trees
-----------------------

Unweighted Robinson-Foulds Distance
...................................

The *unweighted* Robinson-Foulds distance (often referred to as just the Robinson-Foulds distance) is given by the :func:`dendropy.calculate.treecompare.symmetric_difference` function:


.. literalinclude:: /examples/symdiff1.py

Note that the two trees *must* share the same |TaxonNamespace| reference, otherwise an error will be raised::

    >> import dendropy
    >> from dendropy.calculate import treecompare
    >> s1 = "(a,(b,(c,d)));"
    >> s2 = "(a,(d,(b,c)));"
    >> tree1 = dendropy.Tree.get(data=s1, schema='newick')
    >> tree2 = dendropy.Tree.get(data=s2, schema='newick')
    >> print(treecompare.symmetric_difference(tree1, tree2))
    Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
            print(treecompare.symmetric_difference(tree1, tree2))
        File "/Users/jeet/Documents/Projects/Phyloinformatics/DendroPy/dendropy/dendropy/calculate/treecompare.py", line 85, in symmetric_difference
            is_bipartitions_updated=is_bipartitions_updated)
        File "/Users/jeet/Documents/Projects/Phyloinformatics/DendroPy/dendropy/dendropy/calculate/treecompare.py", line 221, in false_positives_and_negatives
            raise error.TaxonNamespaceIdentityError(reference_tree, comparison_tree)
        dendropy.utility.error.TaxonNamespaceIdentityError: Non-identical taxon namespace references: <TaxonNamespace object at 0x10052d310> is not <TaxonNamespace object at 0x101572210>

Note, too, that results very much depend on the rooting states of the tree:

.. literalinclude:: /examples/symdiff2.py

Weighted Robinson-Foulds Distance
.................................

The *weighted* Robinson-Foulds distance takes edge lengths into account, and is given by the :func:`dendropy.calculate.treecompare.weighted_robinson_foulds_distance`:

.. literalinclude:: /examples/weightedrf1.py

Euclidean Distance
..................

The Euclidean distance, like the weighted Robinson-Foulds distance takes edge lengths into account, but squares the edge lengths instead of taking the absolute distance, and is given by the :func:`dendropy.calculate.treecompare.euclidean_distance`:

.. literalinclude:: /examples/euctree.py

Majority-Rule Consensus Tree from a Collection of Trees
-------------------------------------------------------

To get the majority-rule consensus tree of a |TreeList| object, you can call the :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.consensus()` instance method.
You can specify the frequency threshold for the consensus tree by the ``min_freq`` argument, which default to 0.5 (i.e., a 50% majority rule tree).
The following example aggregates the post-burn-in trees from four MCMC samples into a single |TreeList| object, and prints the 95% majority-rule consensus as a Newick string:

.. literalinclude:: /examples/majrule.py

Frequency of a Split in a Collection of Trees
---------------------------------------------

The :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.frequency_of_split()` method of a |TreeList| object returns the frequency of occurrence of a single split across all the |Tree| objects in the |TreeList|.
The split can be specified by passing a split bitmask directly using the ``split_bitmask`` keyword argument, as a list of |Taxon| objects using the ``taxa`` keyword argument, or as a list of taxon labels using the ``labels`` keyword argument.
The following example shows how to calculate the frequency of a split defined by two taxa, "Morelia amethistina" and "Morelia tracyae", from the post-burn-in trees aggregated across four MCMC samples:

.. literalinclude:: /examples/splitfreq.py

The Maximum Clade Credibility Tree: The Tree that Maximizes the Product of Split Support
----------------------------------------------------------------------------------------

The Maximum Clade Credibility Tree (MCCT) is one that maximize the *product* of split support, and is returned for a collection of trees managed in a |TreeList| instance by the :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.maximum_product_of_split_support_tree` method:

.. literalinclude:: /examples/mcct.py

Unfortunately, terminology in usage and literature regarding this type of summary is *very* confusing, and sometimes the term "MCCT" is used to refer to the tree that maximizes the *sum* of split support and "MCT" to the tree that maximizes the product of split support.
If the tree that maximizes the *sum* of split support is the criteria required, then the :meth:`~dendropy.datamodel.treecollectionmodel.TreeList.maximum_sum_of_split_support_tree` method of the |TreeList| object should be used.

Scoring Trees Under the Coalescent
==================================

Probability Under the Coalescent Model
---------------------------------------

The :mod:`~dendropy.model.coalescent` module provides a range of methods for simulations and calculations under Kingman's coalescent framework and related models. For example:

    :func:`~dendropy.model.coalescent.log_probability_of_coalescent_tree`
        Given a |Tree| object as the first argument, and the haploid population size as the second, returns the log probability of the |Tree| under the neutral coalescent.

Numbers of Deep Coalescences
----------------------------

    :func:`~dendropy.model.reconcile.reconciliation_discordance`
        Given two |Tree| objects *sharing the same leaf-set*, this returns the number of deep coalescences resulting from fitting the first tree (e.g., a gene tree) to the second (e.g., a species tree). This is based on the algorithm described `Goodman, et al. <bioinformatics.oxfordjournals.org/cgi/reprint/14/9/819.pdf>`_ (Goodman, et al., 1979. Fitting the gene lineage into its species lineage,a parsimony strategy illustrated by cladograms constructed from globin sequences. Syst. Zool. 19: 99-113).

    :func:`~dendropy.model.reconcile.monophyletic_partition_discordance`
        Given a |Tree| object as the first argument, and a list of lists of
        |Taxon| objects representing the expected monophyletic partitioning of the |TaxonNamespace| of the |Tree| as the second argument, this returns the number of deep coalescences found in the relationships implied by the |Tree| object, conditional on the taxon groupings given by the second argument. This statistic corresponds to the Slatkin and Maddison (1989) **s** statistic, as described `here <http://mesquiteproject.org/Mesquite_Folder/docs/mesquite/popGen/popGen.html#s>`_.

Number of Deep Coalescences when Embedding One Tree in Another (e.g. Gene/Species Trees)
----------------------------------------------------------------------------------------

Imagine we wanted to generate the distribution of the number of deep coalescences under two scenarios: one in which a population underwent sequential or step-wise vicariance, and another when there was simultaneous fragmentation.
In this case, the containing tree and the embedded trees have different leaf sets, and there is a many-to-one mapping of embedded tree taxa to containing tree taxa.

The :class:`~dendropy.model.reconcile.ContainingTree` class is designed to allow for counting deep coalescences in cases like this.
It requires a |TaxonNamespaceMapping| object, which provides an association between the embedded taxa and the containing taxa.
The easiest way to get a |TaxonNamespaceMapping| object is to call the special factory function :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespaceMapping.create_contained_taxon_mapping()`.
This will create a new |TaxonNamespace| to manage the gene taxa, and create the associations between the gene taxa and the containing tree taxa for you.
It takes two arguments: the |TaxonNamespace| of the containing tree, and the number of genes you want sampled from each species.
If the gene-species associations are more complex, e.g., different numbers of genes per species, we can pass in a list of values as the second argument to :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespaceMapping.create_contained_taxon_mapping()`.
This approach should be used with caution if we cannot be certain of the order of taxa (as is the case with data read in Newick formats). In these case, and in more complex cases, we might need to directly instantiate the :class:`~dendropy.datamodel.taxonmodel.TaxonNamespaceMapping` object. The API to describe the associations when constructing this object is very similar to that of the :class:`~dendropy.datamodel.taxonmodel.TaxonNamespacePartition` object: you can use a function, attribute or dictionary.

The :class:`~dendropy.model.reconcile.ContainingTree` class has its own native contained coalescent simulator, :meth:`~dendropy.model.reconcile.ContainingTree.embed_contained_kingman()`, which simulates *and* embeds a contained coalescent tree at the same time.

.. literalinclude:: /examples/sim_and_count_deepcoal2.py

If you have used some other method to simulate your trees, you can use :meth:`~dendropy.model.reconcile.ContainingTree.embed_tree()` to embed the trees and count then number of deep coalescences.

.. literalinclude:: /examples/sim_and_count_deepcoal1.py

For more details on simulating contained coalescent trees and counting numbers of deep coalescences on them, see ":ref:`Simulating_Contained_Coalescent_Trees`" or ":ref:`Simulating_and_Counting_Deep_Coalescences`".



