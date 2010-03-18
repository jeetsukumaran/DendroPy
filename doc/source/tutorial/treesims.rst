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

The birth rate is specified by the `birth_rate` argument, and death rate specified by the `death_rate` argument.
A variable birth-rate or death-rate can be set by the `birth_rate_sd` or the `death_rate_sd`  arguments.
`birth_rate_sd` is the standard deviation of the normally-distributed mutation
added to the birth rate as it is inherited by daughter nodes; if 0, birth
rate does not evolve on the tree.
`death_rate_sd` is the standard deviation of the normally-distributed mutation
added to the death rate as it is inherited by daughter nodes; if 0, death
rate does not evolve on the tree.

Tree growth is controlled by one or more of the following arguments, of which
at least one must be specified:

    - If `ntax` is given as a keyword argument, tree is grown until the number of
      tips == ntax.
    - If `taxon_set` is given as a keyword argument, tree is grown until the
      number of tips == len(taxon_set), and the taxa are assigned randomly to the
      tips.
    - If 'max_time' is given as a keyword argument, tree is grown for
      a maximum of `max_time`.

If more than one of the above is given, then tree growth will terminate when
*any* of the termination conditions (i.e., number of tips == `ntax`, or number
of tips == len(taxon_set) or maximum time = `max_time`) are met.

These functions also accept a Tree object (with valid branch lengths) as an argument passed using the keyword `tree`: if given, then this tree will be used; otherwise
a new one will be created.

If `assign_taxa` is False, then taxa will *not* be assigned to the tips;
otherwise (default), taxa will be assigned. If `taxon_set` is given
(`tree.taxon_set`, if `tree` is given), and the final number of tips on the
tree after the termination condition is reached is less then the number of
taxa in `taxon_set` (as will be the case, for example, when
`ntax` < len(`taxon_set`)), then a random subset of taxa in `taxon_set` will
be assigned to the tips of tree. If the number of tips is more than the number
of taxa in the `taxon_set`, new Taxon objects will be created and added
to the `taxon_set` if the keyword argument `create_required_taxa` is not given as
False.

In addition, a Random() object or equivalent can be passed using the `rng` keyword;
otherwise a global random number generator will be used.

Star Tree
=========

The :func:`~dendropy.treesim.star_tree()` generates a simply polytomy tree, with a single node as the immediate ancestor to a set of leaves, with one leaf per |Taxon| in the |TaxonSet| object given by the `taxon_set` argument.

Population Genetic Tree
=======================

The :func:`pop_gen_tree()` function generates a tree with edges decorated with population sizes and leaf nodes decorated by the number of genes (samples or lineages) in each leaf.
This tree is useful for coalescent-simulations (see below).

Coalescent Tree
===============
The :func:`pure_kingman()` function simulates a tree under Kingman's n-coalescent (i.e., the pure, unconstrained coalescent process).


Censored Coalescent Tree
========================
The :func:`constrained_kingman()` function simulates a tree under the censored coalescent, i.e., the coalescent conditional or constrained by a containing species or population tree.

