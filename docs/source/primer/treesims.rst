******************************
Tree Simulation and Generation
******************************

The :mod:`~dendropy.simulate.treesim` module provides functions for the simulation of trees under a variety of theoretical models.
This module is actually just a namespace that aggregates functions and classes for tree simulation routines in one convenient place.
For example the :func:`~dendropy.model.birthdeath.birth_death_tree()` function
is actually defined in the :mod:`~dendropy.model.birthdeath` module, but is
exposed in the :mod:`~dendropy.simulate.treesim` for ease of access.

Birth-Death Process Trees
=========================

There are two different birth-death process tree simulation routines in DendroPy:

    :func:`~dendropy.simulate.treesim.birth_death_tree()`
        Returns a tree generated under a continuous-time birth-death process, with branch lengths in arbitrary time units.

    :func:`~dendropy.simulate.treesim.discrete_birth_death_tree()`
        Returns a tree generated under discrete-time birth-death process, with branch length in generation units.

Both of these functions have identical interfaces, and will grow a tree under a branching process with the specified birth-date and death-rate until the termination condition (pre-specified number of leaves or maximum amount of time) is met.

For example, to get a continuous-time tree with 10 leaves, generated under a birth rate of 1.0 and death rate of 0.5::

    >>> from dendropy.simulate import treesim
    >>> t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=10)
    >>> t.print_plot()
                  /-------------------------------------------- T1
                  |
    /-------------+                             /-------------- T2
    |             |              /--------------+
    |             \--------------+              \-------------- T3
    |                            |
    |                            \----------------------------- T4
    +
    |                            /----------------------------- T5
    |             /--------------+
    |             |              |              /-------------- T6
    |             |              \--------------+
    \-------------+                             \-------------- T7
                  |
                  |                             /-------------- T8
                  |              /--------------+
                  \--------------+              \-------------- T9
                                 |
                                 \----------------------------- T10

While to get a continuous time tree generated under the same rates after 6 time units::

    >>> t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, max_time=6.0)

If both conditions are given simultaneously, then tree growth will terminate when
*any* of the termination conditions (i.e., number of tips == ``ntax``, or number
of tips == len(taxon_namespace) or maximum time == ``max_time``) are met.

Specifying a |TaxonNamespace|
-----------------------------

By default, a new |Taxon| object will be created and associated with each leaf (labeled "T1", "T2", etc.),  all belonging to a new |TaxonNamespace| object associated with the resulting tree.

You can pass in an explicit |TaxonNamespace| object using the "``taxon_namespace``" keyword::

    >>> import dendropy
    >>> from dendropy.simulate import treesim
    >>> taxa = dendropy.TaxonNamespace(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
    >>> t = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa)
    >>> t.print_plot()
                /-------------------------------------- h
                |
    /-----------+                         /------------ c
    |           |            /------------+
    |           \------------+            \------------ a
    |                        |
    +                        \------------------------- g
    |
    |                                     /------------ e
    |                        /------------+
    |                        |            \------------ f
    \------------------------+
                             |            /------------ d
                             \------------+
                                          \------------ b


In this case, the branching process underlying the tree generation will terminate when the number of leaves in the tree equals the number of taxa in the |TaxonNamespace| "``taxa``", and the |Taxon| objects in "``taxa``" will be randomly assigned to the leaves.

The "``taxon_namespace``" keyword can be combined with the "``ntax``" keyword.
If the size of the |TaxonNamespace| object given by the ``taxon_namespace`` argument is greater than the specified target tree taxon number, then a random subset of |Taxon| object in the |TaxonNamespace| will be assigned to the leaves::

    >>> import dendropy
    >>> from dendropy.simulate import treesim
    >>> taxa = dendropy.TaxonNamespace(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
    >>> t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=5, taxon_namespace=taxa)
    >>> t.print_plot()
    /-------------------------------------------------- g
    |
    +                        /------------------------- a
    |           /------------+
    |           |            |            /------------ d
    \-----------+            \------------+
                |                         \------------ c
                |
                \-------------------------------------- f

If the size of the |TaxonNamespace| object is less than the target taxon number, then new |Taxon| objects will be created as needed and added to the |TaxonNamespace| object as well as associated with the leaves::

    >>> import dendropy
    >>> from dendropy.simulate import treesim
    >>> taxa = dendropy.TaxonNamespace(['a', 'b'])
    >>> t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=5, taxon_namespace=taxa)
    >>> t.print_plot()
                                     /---------------- a
    /--------------------------------+
    |                                \---------------- b
    +
    |               /--------------------------------- T3
    \---------------+
                    |                /---------------- T4
                    \----------------+
                                     \---------------- T5


Repeating Failed Branching Processes
------------------------------------

With a non-zero death rate, it is possible for all lineages of a tree to go extinct before the termination conditions are reached.
In this case, by default a :class:`~dendropy.simulate.treesim.TreeSimTotalExtinctionException` will be raised::

    >>> t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.9, ntax=10)

produces::

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/Users/jeet/Projects/DendroPy/dendropy/treesim.py", line 188, in birth_death
        raise TreeSimTotalExtinctionException()
    dendropy.simulate.treesim.TreeSimTotalExtinctionException

If the keyword argument "``repeat_until_success``" is given, then instead of raising an exception the process starts again and repeats until the termination condition is met::

    >>> t = treesim.birth_death_tree(birth_rate=1.0,
    ...                         death_rate=0.9,
    ...                         ntax=10,
    ...                         repeat_until_success=True)
    >>> t.print_plot()
                                           /------------------- T1
    /--------------------------------------+
    |                                      |         /--------- T2
    |                                      \---------+
    |                                                \--------- T3
    +
    |                                                /--------- T4
    |        /---------------------------------------+
    |        |                                       \--------- T5
    |        |
    \--------+                   /----------------------------- T6
             |         /---------+
             |         |         |         /------------------- T7
             |         |         \---------+
             \---------+                   |         /--------- T8
                       |                   \---------+
                       |                             \--------- T9
                       |
                       \--------------------------------------- T10

Suppressing Taxon Assignment
----------------------------
You can specify "``assign_taxa``" to be ``False``  to avoid taxa from being automatically assigned to a tree (for example, when you want to build a tree in stages -- see below).

Extending an Existing Tree
--------------------------

Both these functions also accept a Tree object (with valid branch lengths) as an argument passed using the keyword ``tree``.
If given, then this tree will be used as the starting point; otherwise a new one will be created.

Evolving Birth and Death Rates
------------------------------

The same functions can also produce trees generated under variable birth and death rates.
The "``birth_rate_sd``" keyword argument specifies the standard deviation of the normally-distributed error of birth rates as they evolve from parent to child node, while the "``death_rate_sd``" keyword argument specifies the same of the the death rates.
For example, to get a 10-taxon tree generated under a birth- and death-rate that evolves with a standard deviation of 0.1::

    >>> t = treesim.birth_death_tree(birth_rate=1.0,
                death_rate=0.5,
                birth_rate_sd=0.1,
                death_rate_sd=0.1,
                ntax=10)

Building a Tree in Multiple Stages under Different Conditions
-------------------------------------------------------------

You might want to generate a tree under different condition in different stages.
To do this, you would start with an empty tree and passing it to the birth-death function as an argument using the "``tree``" keyword argument, and at the same time suppress the automatic taxon assignment using the "``assign_taxa=False``" keyword argument to avoid taxa being assigned to what will eventually become internal nodes.
When the tree is ready, you will call the :meth:`~dendropy.datamodel.treemodel.Tree.randomly_assign_taxa()` function to assign taxa at random to the leaves.

For example, the following generates a birth-death tree with equal birth and death rates, but both rates shifting for a short while to a temporarily higher (though equal) rates:

.. literalinclude:: /examples/bdtree_multi1.py
    :linenos:

Another example draws birth and death rates from a normal distribution with the same mean and standard deviation in multiple stages:

.. literalinclude:: /examples/bdtree_multi2.py
    :linenos:

Star Trees
==========

The :func:`~dendropy.simulate.treesim.star_tree()` generates a simple polytomy tree, with a single node as the immediate ancestor to a set of leaves, with one leaf per |Taxon| in the |TaxonNamespace| object given by the ``taxon_namespace`` argument.
For example::

    >>> from dendropy.simulate import treesim
    >>> taxa = dendropy.TaxonNamespace(['a', 'b', 'c', 'd', 'e'])
    >>> tree = treesim.star_tree(taxa)
    >>> print(tree.as_ascii_plot())
    /-------------------------------------- a
    |
    |-------------------------------------- b
    |
    +-------------------------------------- c
    |
    |-------------------------------------- d
    |
    \-------------------------------------- e


(Pure Neutral) Coalescent Trees
===============================

The :func:`~dendropy.simulate.treesim.pure_kingman()` function returns a tree generated under an unconstrained neutral coalescent model. The first argument to this function, ``taxon_namespace``, is a |TaxonNamespace| object, where each member |Taxon| object represents a gene to be coalesced. The second argument, ``pop_size``, specifies the population size in terms of the number of gene copies in the population. This means that for a diploid population of size ``N``, ``pop_size`` should be ``N*2``, while for a haploid population of size ``N``, ``pop_size`` should be ``N``. If ``pop_size`` is |None|, 1, or 0, then the edge lengths of the returned gene tree will be in population units (i.e., 1 unit of edge length == to 2N generations if a diploid population or 1N generations if a haploid population). Otherwise, the edge lengths will be in generation units. For example:

.. literalinclude:: /examples/pure_kingman1.py

.. _Simulating_Contained_Coalescent_Trees:

Multispecies Coalescent ("Contained Coalescent" or "Censored Coalescent") Trees
===============================================================================

The :func:`~dendropy.simulate.treesim.contained_coalescent()` function returns a tree generated under a neutral coalescent model conditioned on population splitting times or events given by a containing species or population tree.
Such a tree is often referred to as a contained, embedded, censored, truncated, or constrained genealogy/tree.
At a minimum, this function takes two arguments: a |Tree| object representing the containing (species or population) tree, and a |TaxonNamespaceMapping| object describing how the sampled gene taxa map or are associated with the species/population |Tree| taxa.

The |Tree| object representing the containing species or population tree should be rooted and ultrametric.
If edge lengths are given in generations, then a meaningful population size needs to be communicated to the :func:`~dendropy.simulate.treesim.contained_coalescent()` function.
In general, for coalescent operations in DendroPy, unless otherwise specified, population sizes are the *haploid* population size, i.e. the number of genes in the population.
This is 2N for a diploid population with N individuals, or N for haploid population of N individuals.
If edge lengths are given in population units (e.g., N), then the appropriate population size to use is 1.

If the population size is fixed throughout the containing species/population tree, then simply passing in the appropriate value using the ``default_pop_size`` argument to the :func:`~dendropy.simulate.treesim.contained_coalescent()` function is sufficient.
If, on the other hand, the population size varies, the a special attribute must be added to each edge, "``pop_size``", that specifies the population size for that edge.
For example::

    tree = dendropy.Tree.get_from_path("sp.tre", "newick")
    for edge in tree.postorder_edge_iter():
            edge.pop_size = 100000

The easiest way to get a |TaxonNamespaceMapping| object is to call the special factory function :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespaceMapping.create_contained_taxon_mapping()`.
This will create a new |TaxonNamespace| to manage the gene taxa, and create the associations between the gene taxa and the containing tree taxa for you.
It takes two arguments: the |TaxonNamespace| of the containing tree, and the number of genes you want sampled from each species.

The following example shows how to create a |TaxonNamespaceMapping| using :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespaceMapping.create_contained_taxon_mapping()`, and then calls :meth:`~dendropy.simulate.treesim.contained_coalescent()` to produce a contained coalescent tree:

.. literalinclude:: /examples/contained_coalescent1.py

In the above example, the branch lengths were in haploid population units, so we did not specify a population size.
If the gene-species associations are more complex, e.g., different numbers of genes per species, we can pass in a list of values as the second argument to :meth:`~dendropy.datamodel.taxonmodel.TaxonNamespaceMapping.create_contained_taxon_mapping()`:


.. literalinclude:: /examples/contained_coalescent2.py

This approach should be used with caution if we cannot be certain of the order of taxa (as is the case with data read in Newick formats). In these case, and in more complex cases, we might need to directly instantiate the :class:`~dendropy.datamodel.taxonmodel.TaxonNamespaceMapping` object. The API to describe the associations when constructing this object is very similar to that of the :class:`~dendropy.datamodel.taxonmodel.TaxonNamespacePartition` object: you can use a function, attribute or dictionary.

.. _Simulating_and_Counting_Deep_Coalescences:

Simulating the Distribution of Number Deep Coalescences Under Different Phylogeographic History Scenarios
=========================================================================================================

A typical application for simulating censored coalescent trees is to produce a distribution of trees under different hypotheses of demographic or phylogeographic histories.

For example, imagine we wanted to generate the distribution of the number of deep coalescences under two scenarios: one in which a population underwent sequential or step-wise vicariance, and another when there was simultaneous fragmentation.
This can be achieved by generating trees under :meth:`~dendropy.simulate.treesim.contained_coalescent()`, and then using a :class:`~dendropy.reconcile.ContainingTree` object to embed the trees and count the number of deep coalescences.

.. literalinclude:: /examples/sim_and_count_deepcoal1.py

Actually, the  :class:`~dendropy.reconcile.ContainingTree` class has its own native contained coalescent simulator, :meth:`~dendropy.reconcile.ContainingTree.embed_contained_kingman()`, which simulates *and* embeds a contained coalescent tree at the same time. So a more practical approach might be:

.. literalinclude:: /examples/sim_and_count_deepcoal2.py
