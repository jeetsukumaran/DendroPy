******************************
Tree Simulation and Generation
******************************

The :mod:`~dendropy.treesim` module provides functions for the simulation of trees under a variety of theoretical models.

Birth-Death Process Trees
=========================

Constant Birth and Death Rates
------------------------------

There are two different birth-death process tree simulation routines in DendroPy:

    :func:`~dendropy.treesim.birth_death()`
        Returns a tree generated under a continuous-time birth-death process, with branch lengths in arbitrary time units.

    :func:`~dendropy.treesim.discrete_birth_death()`
        Returns a tree generated under discrete-time birth-death process, with branch length in generation units.

Both of these functions have identical interfaces, and will grow a tree under a branching process with the specified birth-date and death-rate until the termination condition (pre-specified number of leaves or maximum amount of time) is met.

For example, return a continuous time tree with 10 leaves, generated under a birth rate of 1.0 and death rate of 0.5::

    >>> from dendropy import treesim
    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, ntax=10)

While to return a continuous time tree generated under the same rates after 6 time units::

    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, max_time=6.0)

If both conditions are given simultaneously, then tree growth will terminate when
*any* of the termination conditions (i.e., number of tips == `ntax`, or number
of tips == len(taxon_set) or maximum time = `max_time`) are met.

By default, a new |Taxon| object will be created and associated with each leaf (labeled "T1", "T2", etc.),  all belonging to a new |TaxonSet| object associated with the resulting tree.

You can pass in an explicit |TaxonSet| object using the "``taxon_set``" keyword.
For example, assuming "``ts``" is a pre-existing |TaxonSet|::

    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, taxon_set=ts)

In this case, the branching process underlying the tree generation will terminate when the number of leaves in the tree equals the number of taxa in the |TaxonSet| "``ts``", and the |Taxon| objects in "``ts``" will be randomly assigned to the leaves.

The "``taxon_set``" keyword can be combined with the "``ntax``" keyword::

    >>> t = treesim.birth_death(birth_rate=1.0, death_rate=0.5, ntax=8, taxon_set=ts)

Here, if the size of the |TaxonSet| object given by the ``taxon_set`` argument is greater than the specified target tree taxon number, then a random subset of |Taxon| object in the |TaxonSet| will be assigned to the leaves.
If the size of the |TaxonSet| object is less than the target taxon number, then new |Taxon| objects will be created as needed and added to the |TaxonSet| object as well as associated with the leaves.

With a non-zero death rate, it is possible for all lineages of a tree to go extinct before the termination conditions are reached.
In this case, by default a :class:`~dendropy.treesim.TreeSimTotalExtinctionException` will be raised.
If the keyword argument "``repeat_on_total_extinction``" is given, then instead of raising an exception, the process starts again and repeats until the termination condition is met.

Variable Birth and Death Rates
------------------------------

The same functions can also produce trees generated under variable birth and death rates.
The "``birth_rate_sd``" keyword argument specifies the standard deviation of the normally-distributed error of birth rates as they evolve from parent to child node, while the "``death_rate_sd``" keyword argument specifies the same of the the death rates.
For example, to get a 10-taxon tree generated under a birth- and death-rate that evolves with a standard deviation of 0.1::

    >>> t = treesim.birth_death(birth_rate=1.0,
                death_rate=0.5,
                birth_rate_sd=0.1,
                death_rate_sd=0.1,
                ntax=10)

Extending an Existing Tree
--------------------------

Both these functions also accept a Tree object (with valid branch lengths) as an argument passed using the keyword `tree`.
If given, then this tree will be used as the starting point; otherwise a new one will be created.



Star Tree
=========

The :func:`~dendropy.treesim.star_tree()` generates a simply polytomy tree, with a single node as the immediate ancestor to a set of leaves, with one leaf per |Taxon| in the |TaxonSet| object given by the `taxon_set` argument.
For example::

    >>> from dendropy import treesim




Population Genetic Tree
=======================

The :func:`~dendropy.treesim.pop_gen_tree()` function generates a tree with edges decorated with population sizes and leaf nodes decorated by the number of genes (samples or lineages) in each leaf.
This tree is useful for coalescent-simulations (see below).

Coalescent Tree
===============
The :func:`~dendropy.treesim.pure_kingman()` function simulates a tree under Kingman's n-coalescent (i.e., the pure, unconstrained coalescent process).


Censored Coalescent Tree
========================
The :func:`~dendropy.treesim.constrained_kingman()` function simulates a tree under the censored coalescent, i.e., the coalescent conditional or constrained by a containing species or population tree.

