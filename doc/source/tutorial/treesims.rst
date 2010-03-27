******************************
Tree Simulation and Generation
******************************

The :mod:`~dendropy.treesim` module provides functions for the simulation of trees under a variety of theoretical models.

Birth-Death Process Trees
=========================

There are two different birth-death process tree simulation routines in DendroPy:

    :func:`~dendropy.treesim.birth_death()`
        Returns a tree generated under a continuous-time birth-death process, with branch lengths in arbitrary time units.

    :func:`~dendropy.treesim.discrete_birth_death()`
        Returns a tree generated under discrete-time birth-death process, with branch length in generation units.

Both of these functions have identical interfaces, and will grow a tree under a branching process with the specified birth-date and death-rate until the termination condition (pre-specified number of leaves or maximum amount of time) is met.

For example, to get a continuous-time tree with 10 leaves, generated under a birth rate of 1.0 and death rate of 0.5::

    >>> from dendropy import treesim
    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, ntax=10)
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

    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, max_time=6.0)

If both conditions are given simultaneously, then tree growth will terminate when
*any* of the termination conditions (i.e., number of tips == `ntax`, or number
of tips == len(taxon_set) or maximum time == `max_time`) are met.

Specifying a |TaxonSet|
-----------------------

By default, a new |Taxon| object will be created and associated with each leaf (labeled "T1", "T2", etc.),  all belonging to a new |TaxonSet| object associated with the resulting tree.

You can pass in an explicit |TaxonSet| object using the "``taxon_set``" keyword::

    >>> import dendropy
    >>> from dendropy import treesim
    >>> taxa = dendropy.TaxonSet(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
    >>> t = treesim.birth_death(0.4, 0.1, taxon_set=taxa)
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


In this case, the branching process underlying the tree generation will terminate when the number of leaves in the tree equals the number of taxa in the |TaxonSet| "``taxa``", and the |Taxon| objects in "``taxa``" will be randomly assigned to the leaves.

The "``taxon_set``" keyword can be combined with the "``ntax``" keyword.
If the size of the |TaxonSet| object given by the ``taxon_set`` argument is greater than the specified target tree taxon number, then a random subset of |Taxon| object in the |TaxonSet| will be assigned to the leaves::

    >>> import dendropy
    >>> from dendropy import treesim
    >>> taxa = dendropy.TaxonSet(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, ntax=5, taxon_set=taxa)
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

If the size of the |TaxonSet| object is less than the target taxon number, then new |Taxon| objects will be created as needed and added to the |TaxonSet| object as well as associated with the leaves::

    >>> import dendropy
    >>> from dendropy import treesim
    >>> taxa = dendropy.TaxonSet(['a', 'b'])
    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, ntax=5, taxon_set=taxa)
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
In this case, by default a :class:`~dendropy.treesim.TreeSimTotalExtinctionException` will be raised::

    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.9, ntax=10)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/Users/jeet/Projects/DendroPy/dendropy/treesim.py", line 188, in birth_death
        raise TreeSimTotalExtinctionException()
    dendropy.treesim.TreeSimTotalExtinctionException

If the keyword argument "``repeat_until_success``" is given, then instead of raising an exception the process starts again and repeats until the termination condition is met::

    >>> t = treesim.birth_death(birth_rate=1.0,
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

Both these functions also accept a Tree object (with valid branch lengths) as an argument passed using the keyword `tree`.
If given, then this tree will be used as the starting point; otherwise a new one will be created.

Evolving Birth and Death Rates
------------------------------

The same functions can also produce trees generated under variable birth and death rates.
The "``birth_rate_sd``" keyword argument specifies the standard deviation of the normally-distributed error of birth rates as they evolve from parent to child node, while the "``death_rate_sd``" keyword argument specifies the same of the the death rates.
For example, to get a 10-taxon tree generated under a birth- and death-rate that evolves with a standard deviation of 0.1::

    >>> t = treesim.birth_death(birth_rate=1.0,
                death_rate=0.5,
                birth_rate_sd=0.1,
                death_rate_sd=0.1,
                ntax=10)

Building a Tree in Multiple Stages under Different Conditions
-------------------------------------------------------------

You might want to generate a tree under different condition in different stages.
To do this, you would start with an empty tree and passing it to the birth-death function as an argument using the "``tree``" keyword argument, and at the same time suppress the automatic taxon assignment using the "``assign_taxa=False``" keyword argument to avoid taxa being assigned to what will eventually become internal nodes.
When the tree is ready, you will call the :meth:`~dendropy.dataobject.Tree.randomly_assign_taxa()` function to assign taxa at random to the leaves.

For example, the following generates a birth-death tree with equal birth and death rates, but both rates shifting for a short while to a temporarily higher (though equal) rates:

.. literalinclude:: /examples/bdtree_multi1.py
    :linenos:

Another example draws birth and death rates from a normal distribution with the same mean and standard deviation in multiple stages:

.. literalinclude:: /examples/bdtree_multi2.py
    :linenos:

Star Trees
==========

The :func:`~dendropy.treesim.star_tree()` generates a simple polytomy tree, with a single node as the immediate ancestor to a set of leaves, with one leaf per |Taxon| in the |TaxonSet| object given by the ``taxon_set`` argument.
For example::

    >>> from dendropy import treesim
    >>> taxa = dendropy.TaxonSet(['a', 'b', 'c', 'd', 'e'])
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


Population Genetic Trees
========================

Coming soon: :func:`~dendropy.treesim.pop_gen_tree()`.

(Pure Neutral) Coalescent Trees
===============================

Coming soon: :func:`~dendropy.treesim.pure_kingman()`.


Censored/Constrained Coalescent Trees
=====================================

Coming soon: :func:`~dendropy.treesim.constrained_kingman()`.

