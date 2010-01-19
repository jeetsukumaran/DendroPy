******************************
Tree Simulation and Generation
******************************

The :mod:`~dendropy.treesim` module provides functions for the simulation of trees under a variety of theoretical models.

Birth-Death Process Trees
=========================

There are three different birth-death process tree simulation routines in DendroPy:

.. SCRATCH:

    Returns a birth-death tree with birth rate specified by `birth_rate`, and
    death rate specified by `death_rate`, with edge lengths in continuous (real)
    units.

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

    Also accepts a Tree object (with valid branch lengths) as an argument passed
    using the keyword `tree`: if given, then this tree will be used; otherwise
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
    otherwise GLOBAL_RNG is used.
