#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

##############################################################################
## Parts of this code were adapted from:
##
##  -   TESS
##
##          https://github.com/hoehna/TESS
##
##          H{"o}hna S. 2013. Fast simulation of reconstructed phylogenies under
##          global time-dependent birth--death processes. Bioinformatics, 29(11)
##          1367-1374.
##
##          Copyright (c) 2012- Sebastian Hoehna
##
##          TESS is free software; you can redistribute it and/or modify
##          it under the terms of the GNU Lesser General Public License as
##          published by the Free Software Foundation; either version 2
##          of the License, or (at your option) any later version.
##
##############################################################################

"""
Models, modeling and model-fitting of birth-death processes.

-   Nee, S.  2001.  Inferring speciation rates from phylogenies.
    Evolution 55:661-668.
-   Yule, G. U. 1924. A mathematical theory of evolution based on the
    conclusions of Dr.  J. C. Willis. Phil. Trans. R. Soc. Lond. B
    213:21-87.
-   Hoehna, S. (2015). The time-dependent reconstructed evolutionary process
    with a key-role for mass-extinction events. Journal of theoretical biology,
    380, 321-331.

"""

import math
from dendropy.calculate import combinatorics
from dendropy.calculate import probability
from dendropy.utility import GLOBAL_RNG
from dendropy.utility.error import TreeSimTotalExtinctionException
from dendropy.utility import constants
from dendropy.utility import deprecate

import dendropy

def birth_death_tree(birth_rate, death_rate, birth_rate_sd=0.0, death_rate_sd=0.0, **kwargs):
    """
    Returns a birth-death tree with birth rate specified by ``birth_rate``, and
    death rate specified by ``death_rate``, with edge lengths in continuous (real)
    units.

    Tree growth is controlled by one or more of the following arguments, of which
    at least one must be specified:

        - If ``num_extant_tips`` is given as a keyword argument, tree is grown until the
          number of EXTANT tips equals this number.
        - If ``num_extinct_tips`` is given as a keyword argument, tree is grown until the
          number of EXTINCT tips equals this number.
        - If ``num_total_tips`` is given as a keyword argument, tree is grown until the
          number of EXTANT plus EXTINCT tips equals this number.
        - If 'max_time' is given as a keyword argument, tree is grown for
          a maximum of ``max_time``.
        - If ``gsa_ntax`` is given then the tree will be simulated up to this number of
          EXTANT tips (or 0 tips), then a tree will be randomly selected from the
          intervals which corresond to times at which the tree had exactly ``num_extant_tips``
          leaves. This allows for simulations according to the "General
          Sampling Approach" of Hartmann et al. (2010). If this option is
          specified, then ``num_extant_tips`` MUST be specified and
          ``num_extinct_tips`` and ``num_total_tips`` CANNOT be specified.

    If more than one of the above is given, then tree growth will terminate when
    *any* one of the termination conditions are met.

    Parameters
    ----------

    birth_rate : float
        The birth rate.
    death_rate : float
        The death rate.
    birth_rate_sd : float
        The standard deviation of the normally-distributed mutation added to
        the birth rate as it is inherited by daughter nodes; if 0, birth rate
        does not evolve on the tree.
    death_rate_sd : float
        The standard deviation of the normally-distributed mutation added to
        the death rate as it is inherited by daughter nodes; if 0, death rate
        does not evolve on the tree.

    Keyword Arguments
    -----------------

    num_extant_tips: int
        If specified, branching process is terminated when number of EXTANT
        tips equals this number.
    num_extinct_tips: int
        If specified, branching process is terminated when number of EXTINCT
        tips equals this number.
    num_total_tips: int
        If specified, branching process is terminated when number of EXTINCT
        plus EXTANT tips equals this number.
    max_time: float
        If specified, branching process is terminated when time reaches or
        exceeds this value.
    gsa_ntax: int
        The General Sampling Approach threshold for number of taxa. See above
        for details.
    tree : Tree instance
        If given, then this tree will be used; otherwise a new one will be created.
    taxon_namespace : TaxonNamespace instance
        If given, then this will be assigned to the new tree, and, in addition,
        taxa assigned to tips will be sourced from or otherwise created with
        reference to this.
    is_assign_extant_taxa : bool [default: True]
        If False, then taxa will not be assigned to extant tips. If True
        (default), then taxa will be assigned to extant tips. Taxa will be
        assigned from the specified ``taxon_namespace`` or
        ``tree.taxon_namespace``. If the number of taxa required exceeds the
        number of taxa existing in the taxon namespace, new |Taxon| objects
        will be created as needed and added to the taxon namespace.
    is_assign_extinct_taxa : bool [default: True]
        If False, then taxa will not be assigned to extant tips. If True
        (default), then taxa will be assigned to extant tips. Taxa will be
        assigned from the specified ``taxon_namespace`` or
        ``tree.taxon_namespace``. If the number of taxa required exceeds the
        number of taxa existing in the taxon namespace, new |Taxon| objects
        will be created as needed and added to the taxon namespace. Note that
        this option only makes sense if extinct tips are retained (specified via
        'is_retain_extinct_tips' option), and will otherwise be ignored.
    is_add_extinct_attr: bool [default: True]
        If True (default), add an boolean attribute indicating whether or not a
        node is an extinct tip or not. False will skip this. Name of attribute
        is set by 'extinct_attr_name' argument, defaulting to 'is_extinct'.
        Note that this option only makes sense if extinct tips are retained
        (specified via 'is_retain_extinct_tips' option), and will otherwise be
        ignored.
    extinct_attr_name: str [default: 'is_extinct']
        Name of attribute to add to nodes indicating whether or not tip is extinct.
        Note that this option only makes sense if extinct tips are retained
        (specified via 'is_retain_extinct_tips' option), and will otherwise be
        ignored.
    is_retain_extinct_tips : bool [default: False]
        If True, extinct tips will be retained on tree. Defaults to False:
        extinct lineages removed from tree.
    repeat_until_success: bool [default: True]
        Under some conditions, it is possible for all lineages on a tree to go
        extinct. In this case, if this argument is given as |True| (the
        default), then a new branching process is initiated. If |False|
        (default), then a TreeSimTotalExtinctionException is raised.
    rng: random.Random() or equivalent instance
        A Random() object or equivalent can be passed using the ``rng`` keyword;
        otherwise GLOBAL_RNG is used.

    References
    ----------

    Hartmann, Wong, and Stadler "Sampling Trees from Evolutionary Models" Systematic Biology. 2010. 59(4). 465-476

    """
    if "assign_taxa" in kwargs:
        deprecate.dendropy_deprecation_warning(
                message="Deprecated: 'assign_taxa' will no longer be supported as an argument to this function. Use 'is_assign_extant_taxa' and/or 'is_assign_extinct_taxa' instead",
                stacklevel=3)
        a = kwargs.pop("assign_taxa")
        kwargs["is_assign_extant_taxa"] = a
        kwargs["is_assign_extant_taxa"] = a
    if "ntax" in kwargs:
        deprecate.dendropy_deprecation_warning(
                message="Deprecated: 'ntax' is no longer supported as an argument to this function. Use one or more of the following instead: 'num_extant_tips', 'num_extinct_tips', 'num_total_tips', or 'max_time'",
                stacklevel=3)
        kwargs["num_extant_tips"] = kwargs.pop("ntax")
    if (("num_extant_tips" not in kwargs)
            and ("num_extinct_tips" not in kwargs)
            and ("num_total_tips" not in kwargs)
            and ("max_time" not in kwargs) ):
        # if "taxon_namespace" in kwargs:
        #     ### cannot support legacy approach, b/c ``taxon_namespace`` may grow during function, leading to unpredictable behavior
        #     # deprecate.dendropy_deprecation_warning(
        #     #         preamble="Deprecated: The 'taxon_namespace' argument can no longer be used to specify a termination condition as a side-effect. Use one or more of the following instead with the length of the taxon namespace instance as a value: 'num_extant_tips', 'num_extinct_tips', or 'num_total_tips'",
        #     #         old_construct="tree = birth_death_tree(\n    ...\n    taxon_namespace=taxon_namespace,\n    ...\n)",
        #     #         new_construct="tree = birth_death_tree(\n    ...\n    taxon_namespace=taxon_namespace,\n    num_extant_tips=len(taxon_namespace),\n    ...\n)")
        #     # kwargs["num_extant_tips"] = len(kwargs["taxon_namespace"])
        #     raise ValueError("The 'taxon_namespace' argument can no longer be used to specify a termination condition as a side-effect."
        #                      "Use one or more of the following instead with the length of the taxon namespace instance as a value: "
        #                      "'num_extant_tips', 'num_extinct_tips', or 'num_total_tips'.\n"
        #                      "That is, instead of:\n\n"
        #                      "    tree = birth_death_tree(\n        ...\n        taxon_namespace=taxon_namespace,\n        ...\n    )\n\n"
        #                      "Use:\n\n"
        #                      "    ntax = len(taxon_namespace)\n    tree = birth_death_tree(\n        ...\n        taxon_namespace=taxon_namespace,\n        num_extant_tips=ntax,\n        ...\n    )\n"
        #                      "\nOr (recommended):\n\n"
        #                      "    tree = birth_death_tree(\n        ...\n        taxon_namespace=taxon_namespace,\n        num_extant_tips=100,\n        ...\n    )\n"
        #                      "\nNote that the taxon namespace instance size may grow during any particular call of the function depending on taxon assignment/creation settings, so"
        #                      " for stable and predictable behavor it is important to take a snapshot of the desired taxon namespace size before any call of the function, or, better yet"
        #                      " simply pass in a constant value."
        #                      )
        #  else:
        raise ValueError("One or more of the following must be specified: 'num_extant_tips', 'num_extinct_tips', or 'max_time'")
    target_num_extant_tips = kwargs.pop("num_extant_tips", None)
    target_num_extinct_tips = kwargs.pop("num_extinct_tips", None)
    target_num_total_tips = kwargs.pop("num_total_tips", None)
    max_time = kwargs.pop('max_time', None)
    gsa_ntax = kwargs.pop('gsa_ntax', None)
    is_add_extinct_attr = kwargs.pop('is_add_extinct_attr', True)
    extinct_attr_name = kwargs.pop('extinct_attr_name', 'is_extinct')
    is_retain_extinct_tips = kwargs.pop('is_retain_extinct_tips', False)
    is_assign_extant_taxa = kwargs.pop('is_assign_extant_taxa', True)
    is_assign_extinct_taxa = kwargs.pop('is_assign_extinct_taxa', True)
    repeat_until_success = kwargs.pop('repeat_until_success', True)

    tree = kwargs.pop("tree", None)
    taxon_namespace = kwargs.pop("taxon_namespace", None)

    rng = kwargs.pop('rng', GLOBAL_RNG)

    ignore_unrecognized_keyword_arguments = kwargs.pop('ignore_unrecognized_keyword_arguments', False)
    if kwargs and not ignore_unrecognized_keyword_arguments:
        raise ValueError("Unsupported keyword arguments: {}".format(kwargs.keys()))

    terminate_at_full_tree = False

    if gsa_ntax is None:
        terminate_at_full_tree = True
        # gsa_ntax = 1 + target_num_taxa
    elif target_num_extant_tips is None:
        raise ValueError("If 'gsa_ntax' is specified, 'num_extant_tips' must be specified")
    elif target_num_extinct_tips is not None:
        raise ValueError("If 'gsa_ntax' is specified, 'num_extinct_tips' cannot be specified")
    elif target_num_total_tips is not None:
        raise ValueError("If 'gsa_ntax' is specified, 'num_total_tips' cannot be specified")
    elif gsa_ntax < target_num_extant_tips:
        raise ValueError("'gsa_ntax' ({}) must be greater than 'num_extant_tips' ({})".format(gsa_ntax, target_num_extant_tips))

    # initialize tree
    if tree is not None:
        if taxon_namespace is not None:
            assert tree.taxon_namespace is taxon_namespace
        else:
            taxon_namespace = tree.taxon_namespace
        extant_tips = []
        extinct_tips = []
        for nd in tree:
            if not nd._child_nodes:
                if not getattr(nd, extinct_attr_name, False):
                    extant_tips.append(nd)
                    if is_add_extinct_attr:
                        setattr(nd, extinct_attr_name, False)
                else:
                    extinct_tips.append(nd)
                    if is_add_extinct_attr:
                        setattr(nd, extinct_attr_name, True)
            elif is_add_extinct_attr:
                setattr(nd, extinct_attr_name, None)
    else:
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        tree = dendropy.Tree(taxon_namespace=taxon_namespace)
        tree.is_rooted = True
        tree.seed_node.edge.length = 0.0
        tree.seed_node.birth_rate = birth_rate
        tree.seed_node.death_rate = death_rate
        if is_add_extinct_attr:
            setattr(tree.seed_node, extinct_attr_name, False)
        extant_tips = [tree.seed_node]
        extinct_tips = []
    initial_extant_tip_set = list(extant_tips)
    initial_extinct_tip_set = list(extinct_tips)

    total_time = 0

    # for the GSA simulations targetted_time_slices is a list of tuple
    #   the first element in the tuple is the duration of the amount
    #   that the simulation spent at the (targetted) number of taxa
    #   and a list of edge information. The list of edge information includes
    #   a list of terminal edges in the tree and the length for that edge
    #   that marks the beginning of the time slice that corresponds to the
    #   targetted number of taxa.
    targetted_time_slices = []

    while True:
        if gsa_ntax is None:
            if target_num_extant_tips is not None and len(extant_tips) >= target_num_extant_tips:
                break
            if target_num_extinct_tips is not None and len(extinct_tips) >= target_num_extinct_tips:
                break
            if target_num_total_tips is not None and (len(extant_tips) + len(extinct_tips)) >= target_num_total_tips:
                break
            if max_time is not None and total_time >= max_time:
                break
        elif len(extant_tips) >= gsa_ntax:
            break

        # get vector of birth/death probabilities, and
        # associate with nodes/events
        event_rates = []
        event_nodes = []
        for nd in extant_tips:
            if not hasattr(nd, 'birth_rate'):
                nd.birth_rate = birth_rate
            if not hasattr(nd, 'death_rate'):
                nd.death_rate = death_rate
            event_rates.append(nd.birth_rate)
            event_nodes.append((nd, True)) # birth event = True
            event_rates.append(nd.death_rate)
            event_nodes.append((nd, False)) # birth event = False; i.e. death

        # get total probability of any birth/death
        rate_of_any_event = sum(event_rates)

        # waiting time based on above probability
        waiting_time = rng.expovariate(rate_of_any_event)

        if ( (gsa_ntax is not None)
                and (len(extant_tips) == target_num_extant_tips)
                ):
            edge_and_start_length = []
            for nd in extant_tips:
                e = nd.edge
                edge_and_start_length.append((e, e.length))
            targetted_time_slices.append((waiting_time, edge_and_start_length))
            if terminate_at_full_tree:
                break

        # add waiting time to nodes
        for nd in extant_tips:
            try:
                nd.edge.length += waiting_time
            except TypeError:
                nd.edge.length = waiting_time
        total_time += waiting_time

        # if event occurs within time constraints
        if max_time is None or total_time <= max_time:
            # normalize probability
            for i in range(len(event_rates)):
                event_rates[i] = event_rates[i]/rate_of_any_event
            # select node/event and process
            nd, birth_event = probability.weighted_choice(event_nodes, event_rates, rng=rng)
            extant_tips.remove(nd)
            if birth_event:
                if is_add_extinct_attr:
                    setattr(nd, extinct_attr_name, None)
                c1 = nd.new_child()
                c2 = nd.new_child()
                c1.edge.length = 0
                c2.edge.length = 0
                c1.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c1.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                c2.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c2.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                extant_tips.append(c1)
                extant_tips.append(c2)
            else:
                if len(extant_tips) > 0:
                    extinct_tips.append(nd)
                    if is_add_extinct_attr:
                        setattr(nd, extinct_attr_name, None)
                else:
                    # total extinction
                    if (gsa_ntax is not None):
                        if (len(targetted_time_slices) > 0):
                            break
                    if not repeat_until_success:
                        raise TreeSimTotalExtinctionException()
                    # We are going to basically restart the simulation because
                    # the tree has gone extinct (without reaching the specified
                    # ntax)
                    extant_tips = list(initial_extant_tip_set)
                    extinct_tips = list(initial_extinct_tip_set)
                    for nd in extant_tips:
                        if is_add_extinct_attr:
                            setattr(nd, extinct_attr_name, False)
                        nd.clear_child_nodes()
                    total_time = 0

    if gsa_ntax is not None:
        total_duration_at_target_n_tax = 0.0
        for i in targetted_time_slices:
            total_duration_at_target_n_tax += i[0]
        r = rng.random()*total_duration_at_target_n_tax
        selected_slice = None
        for n, i in enumerate(targetted_time_slices):
            r -= i[0]
            if r < 0.0:
                selected_slice = i
        assert(selected_slice is not None)
        edges_at_slice = selected_slice[1]
        last_waiting_time = selected_slice[0]

        for e, prev_length in edges_at_slice:
            daughter_nd = e.head_node
            for nd in daughter_nd.child_nodes():
                nd._parent_node = None
                try:
                    extinct_tips.remove(nd)
                except ValueError:
                    pass
                try:
                    extant_tips.remove(nd)
                except ValueError:
                    pass
                for desc in nd.preorder_iter():
                    try:
                        extant_tips.remove(desc)
                    except ValueError:
                        pass
            daughter_nd.clear_child_nodes()
            try:
                extinct_tips.remove(daughter_nd)
            except ValueError:
                pass
            extant_tips.append(daughter_nd)
            if is_add_extinct_attr:
                setattr(daughter_nd, extinct_attr_name, False)
            e.length = prev_length + last_waiting_time

    if not is_retain_extinct_tips:
        processed_nodes = set()
        for nd in list(extinct_tips):
            if nd in processed_nodes:
                continue
            processed_nodes.add(nd)
            try:
                extinct_tips.remove(nd)
            except ValueError:
                pass
            assert not nd._child_nodes
            while (nd.parent_node is not None) and (len(nd.parent_node._child_nodes) == 1):
                nd = nd.parent_node
                processed_nodes.add(nd)
            tree.prune_subtree(nd, suppress_unifurcations=False)
    tree.suppress_unifurcations()

    if is_assign_extant_taxa or is_assign_extinct_taxa:
        taxon_pool = [t for t in taxon_namespace]
        rng.shuffle(taxon_pool)
        taxon_pool_labels = set([t.label for t in taxon_pool])

        ### ONLY works if in GSA sub-section we remove ALL extant and
        ### extinct nodes beyond time slice: expensive
        ### Furthermore, main reason to use this approach is to have different
        ### label prefixes for extinct vs. extant lineages, but the second time
        ### this function is called with the same taxon namespace or any time
        ### this function is called with a populated taxon namespace, that
        ### aesthetic is lost.
        # node_pool_labels = ("T", "X")
        # for node_pool_idx, node_pool in enumerate((extant_tips, extinct_tips)):
        #     for node_idx, nd in enumerate(node_pool):
        #         if taxon_pool:
        #             taxon = taxon_pool.pop()
        #         else:
        #             taxon = taxon_namespace.require_taxon("{}{}".format(node_pool_labels[node_pool_idx], node_idx+1))
        #         nd.taxon = taxon
        #         assert not nd._child_nodes

        tlabel_counter = 0
        leaf_nodes = tree.leaf_nodes()
        rng.shuffle(leaf_nodes)
        for nd_idx, nd in enumerate(leaf_nodes):
            if not is_assign_extant_taxa and nd in extant_tips:
                continue
            if not is_assign_extant_taxa and nd in extinct_tips:
                continue
            if taxon_pool:
                taxon = taxon_pool.pop()
            else:
                while True:
                    tlabel_counter += 1
                    label = "{}{}".format("T", tlabel_counter)
                    if label not in taxon_pool_labels:
                        break
                taxon = taxon_namespace.require_taxon(label=label)
                taxon_pool_labels.add(label)
            nd.taxon = taxon
    return tree


#NICOLA: new version of birth-death tree simulation with linear time cost
#in number of tips instead of quadratic. To obtain this accelleration, however,
#I have removed variation in birth and death rates.
def fast_birth_death_tree(birth_rate, death_rate, **kwargs):
    #NICOLA: we are not allowing variation in birth and death rate so to have
    #increased efficiency.
    birth_rate_sd=0.0
    death_rate_sd=0.0

    """
    Returns a birth-death tree with birth rate specified by ``birth_rate``, and
    death rate specified by ``death_rate``, with edge lengths in continuous (real)
    units.

    Tree growth is controlled by one or more of the following arguments, of which
    at least one must be specified:

        - If ``num_extant_tips`` is given as a keyword argument, tree is grown until the
          number of EXTANT tips equals this number.
        - If ``num_extinct_tips`` is given as a keyword argument, tree is grown until the
          number of EXTINCT tips equals this number.
        - If ``num_total_tips`` is given as a keyword argument, tree is grown until the
          number of EXTANT plus EXTINCT tips equals this number.
        - If 'max_time' is given as a keyword argument, tree is grown for
          a maximum of ``max_time``.
        - If ``gsa_ntax`` is given then the tree will be simulated up to this number of
          EXTANT tips (or 0 tips), then a tree will be randomly selected from the
          intervals which corresond to times at which the tree had exactly ``num_extant_tips``
          leaves. This allows for simulations according to the "General
          Sampling Approach" of Hartmann et al. (2010). If this option is
          specified, then ``num_extant_tips`` MUST be specified and
          ``num_extinct_tips`` and ``num_total_tips`` CANNOT be specified.

    If more than one of the above is given, then tree growth will terminate when
    *any* one of the termination conditions are met.

    Parameters
    ----------

    birth_rate : float
        The birth rate.
    death_rate : float
        The death rate.
    birth_rate_sd : float
        The standard deviation of the normally-distributed mutation added to
        the birth rate as it is inherited by daughter nodes; if 0, birth rate
        does not evolve on the tree.
    death_rate_sd : float
        The standard deviation of the normally-distributed mutation added to
        the death rate as it is inherited by daughter nodes; if 0, death rate
        does not evolve on the tree.

    Keyword Arguments
    -----------------

    num_extant_tips: int
        If specified, branching process is terminated when number of EXTANT
        tips equals this number.
    num_extinct_tips: int
        If specified, branching process is terminated when number of EXTINCT
        tips equals this number.
    num_total_tips: int
        If specified, branching process is terminated when number of EXTINCT
        plus EXTANT tips equals this number.
    max_time: float
        If specified, branching process is terminated when time reaches or
        exceeds this value.
    gsa_ntax: int
        The General Sampling Approach threshold for number of taxa. See above
        for details.
    tree : Tree instance
        If given, then this tree will be used; otherwise a new one will be created.
    taxon_namespace : TaxonNamespace instance
        If given, then this will be assigned to the new tree, and, in addition,
        taxa assigned to tips will be sourced from or otherwise created with
        reference to this.
    is_assign_extant_taxa : bool [default: True]
        If False, then taxa will not be assigned to extant tips. If True
        (default), then taxa will be assigned to extant tips. Taxa will be
        assigned from the specified ``taxon_namespace`` or
        ``tree.taxon_namespace``. If the number of taxa required exceeds the
        number of taxa existing in the taxon namespace, new |Taxon| objects
        will be created as needed and added to the taxon namespace.
    is_assign_extinct_taxa : bool [default: True]
        If False, then taxa will not be assigned to extant tips. If True
        (default), then taxa will be assigned to extant tips. Taxa will be
        assigned from the specified ``taxon_namespace`` or
        ``tree.taxon_namespace``. If the number of taxa required exceeds the
        number of taxa existing in the taxon namespace, new |Taxon| objects
        will be created as needed and added to the taxon namespace. Note that
        this option only makes sense if extinct tips are retained (specified via
        'is_retain_extinct_tips' option), and will otherwise be ignored.
    is_add_extinct_attr: bool [default: True]
        If True (default), add an boolean attribute indicating whether or not a
        node is an extinct tip or not. False will skip this. Name of attribute
        is set by 'extinct_attr_name' argument, defaulting to 'is_extinct'.
        Note that this option only makes sense if extinct tips are retained
        (specified via 'is_retain_extinct_tips' option), and will otherwise be
        ignored.
    extinct_attr_name: str [default: 'is_extinct']
        Name of attribute to add to nodes indicating whether or not tip is extinct.
        Note that this option only makes sense if extinct tips are retained
        (specified via 'is_retain_extinct_tips' option), and will otherwise be
        ignored.
    is_retain_extinct_tips : bool [default: False]
        If True, extinct tips will be retained on tree. Defaults to False:
        extinct lineages removed from tree.
    repeat_until_success: bool [default: True]
        Under some conditions, it is possible for all lineages on a tree to go
        extinct. In this case, if this argument is given as |True| (the
        default), then a new branching process is initiated. If |False|
        (default), then a TreeSimTotalExtinctionException is raised.
    rng: random.Random() or equivalent instance
        A Random() object or equivalent can be passed using the ``rng`` keyword;
        otherwise GLOBAL_RNG is used.

    References
    ----------

    Hartmann, Wong, and Stadler "Sampling Trees from Evolutionary Models" Systematic Biology. 2010. 59(4). 465-476

    """
    if "assign_taxa" in kwargs:
        deprecate.dendropy_deprecation_warning(
                message="Deprecated: 'assign_taxa' will no longer be supported as an argument to this function. Use 'is_assign_extant_taxa' and/or 'is_assign_extinct_taxa' instead",
                stacklevel=3)
        a = kwargs.pop("assign_taxa")
        kwargs["is_assign_extant_taxa"] = a
        kwargs["is_assign_extant_taxa"] = a
    if "ntax" in kwargs:
        deprecate.dendropy_deprecation_warning(
                message="Deprecated: 'ntax' is no longer supported as an argument to this function. Use one or more of the following instead: 'num_extant_tips', 'num_extinct_tips', 'num_total_tips', or 'max_time'",
                stacklevel=3)
        kwargs["num_extant_tips"] = kwargs.pop("ntax")
    if (("num_extant_tips" not in kwargs)
            and ("num_extinct_tips" not in kwargs)
            and ("num_total_tips" not in kwargs)
            and ("max_time" not in kwargs) ):
        if "taxon_namespace" in kwargs:
            ### cannot support legacy approach, b/c ``taxon_namespace`` may grow during function, leading to unpredictable behavior
            # deprecate.dendropy_deprecation_warning(
            #         preamble="Deprecated: The 'taxon_namespace' argument can no longer be used to specify a termination condition as a side-effect. Use one or more of the following instead with the length of the taxon namespace instance as a value: 'num_extant_tips', 'num_extinct_tips', or 'num_total_tips'",
            #         old_construct="tree = birth_death_tree(\n    ...\n    taxon_namespace=taxon_namespace,\n    ...\n)",
            #         new_construct="tree = birth_death_tree(\n    ...\n    taxon_namespace=taxon_namespace,\n    num_extant_tips=len(taxon_namespace),\n    ...\n)")
            # kwargs["num_extant_tips"] = len(kwargs["taxon_namespace"])
            raise ValueError("The 'taxon_namespace' argument can no longer be used to specify a termination condition as a side-effect."
                             "Use one or more of the following instead with the length of the taxon namespace instance as a value: "
                             "'num_extant_tips', 'num_extinct_tips', or 'num_total_tips'.\n"
                             "That is, instead of:\n\n"
                             "    tree = birth_death_tree(\n        ...\n        taxon_namespace=taxon_namespace,\n        ...\n    )\n\n"
                             "Use:\n\n"
                             "    ntax = len(taxon_namespace)\n    tree = birth_death_tree(\n        ...\n        taxon_namespace=taxon_namespace,\n        num_extant_tips=ntax,\n        ...\n    )\n"
                             "\nOr (recommended):\n\n"
                             "    tree = birth_death_tree(\n        ...\n        taxon_namespace=taxon_namespace,\n        num_extant_tips=100,\n        ...\n    )\n"
                             "\nNote that the taxon namespace instance size may grow during any particular call of the function depending on taxon assignment/creation settings, so"
                             " for stable and predictable behavor it is important to take a snapshot of the desired taxon namespace size before any call of the function, or, better yet"
                             " simply pass in a constant value."
                             )
        else:
            raise ValueError("One or more of the following must be specified: 'num_extant_tips', 'num_extinct_tips', or 'max_time'")
    target_num_extant_tips = kwargs.pop("num_extant_tips", None)
    target_num_extinct_tips = kwargs.pop("num_extinct_tips", None)
    target_num_total_tips = kwargs.pop("num_total_tips", None)
    max_time = kwargs.pop('max_time', None)
    gsa_ntax = kwargs.pop('gsa_ntax', None)
    is_add_extinct_attr = kwargs.pop('is_add_extinct_attr', True)
    extinct_attr_name = kwargs.pop('extinct_attr_name', 'is_extinct')
    is_retain_extinct_tips = kwargs.pop('is_retain_extinct_tips', False)
    is_assign_extant_taxa = kwargs.pop('is_assign_extant_taxa', True)
    is_assign_extinct_taxa = kwargs.pop('is_assign_extinct_taxa', True)
    repeat_until_success = kwargs.pop('repeat_until_success', True)

    tree = kwargs.pop("tree", None)
    taxon_namespace = kwargs.pop("taxon_namespace", None)

    rng = kwargs.pop('rng', GLOBAL_RNG)

    ignore_unrecognized_keyword_arguments = kwargs.pop('ignore_unrecognized_keyword_arguments', False)
    if kwargs and not ignore_unrecognized_keyword_arguments:
        raise ValueError("Unsupported keyword arguments: {}".format(kwargs.keys()))

    terminate_at_full_tree = False

    if gsa_ntax is None:
        terminate_at_full_tree = True
        # gsa_ntax = 1 + target_num_taxa
    elif target_num_extant_tips is None:
        raise ValueError("If 'gsa_ntax' is specified, 'num_extant_tips' must be specified")
    elif target_num_extinct_tips is not None:
        raise ValueError("If 'gsa_ntax' is specified, 'num_extinct_tups' cannot be specified")
    elif target_num_total_tips is not None:
        raise ValueError("If 'gsa_ntax' is specified, 'num_total_tips' cannot be specified")
    elif gsa_ntax < target_num_extant_tips:
        raise ValueError("'gsa_ntax' ({}) must be greater than 'num_extant_tips' ({})".format(gsa_ntax, target_num_extant_tips))

    # initialize tree
    if tree is not None:
        if taxon_namespace is not None:
            assert tree.taxon_namespace is taxon_namespace
        else:
            taxon_namespace = tree.taxon_namespace
        extant_tips = []
        extinct_tips = []
        for nd in tree:
            if nd.edge.length is None:
                nd.edge.length = 0.0;
            if not nd._child_nodes:
                if not getattr(nd, extinct_attr_name, False):
                    extant_tips.append(nd)
                    #NICOLA: the reason for this will become clear later
                    nd.edge.length=-nd.edge.length
                    if is_add_extinct_attr:
                        setattr(nd, extinct_attr_name, False)
                else:
                    extinct_tips.append(nd)
                    if is_add_extinct_attr:
                        setattr(nd, extinct_attr_name, True)
            elif is_add_extinct_attr:
                setattr(nd, extinct_attr_name, None)
    else:
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        tree = dendropy.Tree(taxon_namespace=taxon_namespace)
        tree.is_rooted = True
        tree.seed_node.edge.length = 0.0
        tree.seed_node.birth_rate = birth_rate
        tree.seed_node.death_rate = death_rate
        if is_add_extinct_attr:
            setattr(tree.seed_node, extinct_attr_name, False)
        extant_tips = [tree.seed_node]
        extinct_tips = []
    initial_extant_tip_set = list(extant_tips)
    initial_extinct_tip_set = list(extinct_tips)
    #NICOLA: in case we need a restart, let's record the initial branch lengths.
    initial_lengths=[]
    for nd in initial_extant_tip_set:
        initial_lengths.append(nd.edge.length)

    total_time = 0

    # for the GSA simulations targetted_time_slices is a list of tuple
    #   the first element in the tuple is the duration of the amount
    #   that the simulation spent at the (targetted) number of taxa
    #   and a list of edge information. The list of edge information includes
    #   a list of terminal edges in the tree and the length for that edge
    #   that marks the beginning of the time slice that corresponds to the
    #   targetted number of taxa.
    targetted_time_slices = []

    while True:
        if gsa_ntax is None:
            if target_num_extant_tips is not None and len(extant_tips) >= target_num_extant_tips:
                #NICOLA: here and in the other cases, update the branch length so to
                #represent the true branch length (explained better below).
                for nd in extant_tips:
                        nd.edge.length=total_time-nd.edge.length
                break
            if target_num_extinct_tips is not None and len(extinct_tips) >= target_num_extinct_tips:
                for nd in extant_tips:
                        nd.edge.length=total_time-nd.edge.length
                break
            if target_num_total_tips is not None and (len(extant_tips) + len(extinct_tips)) >= target_num_total_tips:
                for nd in extant_tips:
                        nd.edge.length=total_time-nd.edge.length
                break
            if max_time is not None and total_time >= max_time:
                for nd in extant_tips:
                        nd.edge.length=total_time-nd.edge.length
                break
        elif len(extant_tips) >= gsa_ntax:
                for nd in extant_tips:
                        nd.edge.length=total_time-nd.edge.length
                break


        # get vector of birth/death probabilities, and
        # associate with nodes/events
        #NICOLA: replaced this
   #      event_rates = []
#         event_nodes = []
#         for nd in extant_tips:
#             if not hasattr(nd, 'birth_rate'):
#                 nd.birth_rate = birth_rate
#             if not hasattr(nd, 'death_rate'):
#                 nd.death_rate = death_rate
#             event_rates.append(nd.birth_rate)
#             event_nodes.append((nd, True)) # birth event = True
#             event_rates.append(nd.death_rate)
#             event_nodes.append((nd, False)) # birth event = False; i.e. death


        # get total probability of any birth/death
        #rate_of_any_event = sum(event_rates)
        #NICOLA: instead, assume that each node has the same birth and death rate, and sample
        #from the total rate and then sample one node at random from extant ones.
        #this has constant cost instead of linear in the number of extant taxa.
        rate_of_any_event=len(extant_tips)*(birth_rate+death_rate)

        # waiting time based on above probability
        waiting_time = rng.expovariate(rate_of_any_event)

        if ( (gsa_ntax is not None)
                and (len(extant_tips) == target_num_extant_tips)
                ):
            edge_and_start_length = []
            for nd in extant_tips:
                e = nd.edge
                #NICOLA: modified accordingly
                edge_and_start_length.append((e, total_time-e.length))
            targetted_time_slices.append((waiting_time, edge_and_start_length))
            if terminate_at_full_tree:
                #NICOLA: terminal update as before
                for nd in extant_tips:
                        nd.edge.length=total_time-nd.edge.length
                break

        # add waiting time to nodes
        #Nicola: this is the part that makes it quadratic in ntax.
        #instead, we use a trick of initializing nd.edge.length to the time of node creation,
        # and then we update it to the correct branch length when the node is closed or simulations are terminated.
        #for nd in extant_tips:
        #    try:
        #        nd.edge.length += waiting_time
        #    except TypeError:
        #        nd.edge.length = waiting_time
        total_time += waiting_time

        # if event occurs within time constraints
        if max_time is None or total_time <= max_time:
            # normalize probability
            #for i in range(len(event_rates)):
            #    event_rates[i] = event_rates[i]/rate_of_any_event
            # select node/event and process
            #nd, birth_event = probability.weighted_choice(event_nodes, event_rates, rng=rng)
            #NICOLA: instead of sampling from a categorical distribution (requiring linear time in ntax)
            #sample a random integer to represent the chosen taxon, and sample if birth or death event.
            taxI=rng.randint(0, len(extant_tips)-1)
            if rng.random()< birth_rate/(birth_rate+death_rate):
                birth_event=True
            else:
                birth_event=False
            nd=extant_tips[taxI]
            #extant_tips.remove(nd)
            if birth_event:
                if is_add_extinct_attr:
                    setattr(nd, extinct_attr_name, None)
                c1 = nd.new_child()
                c2 = nd.new_child()
                extant_tips[taxI]=c1
                #NICOLA: we initialize branch lengths to total time;
                #this may not make sense at first, but when we are done with a node, we
                #update the branch length to the new current time minus the time of node creation
                #which gives us the branch length. The advantage is that this takes linear time
                #instead of quadratic.
                c1.edge.length = total_time
                c2.edge.length = total_time
                nd.edge.length = total_time-nd.edge.length
                #c1.edge.length = 0
                #c2.edge.length = 0
                #NICOLA: here, for speed, we don't allow variation in birth rates.
                #It could be possible to allow it, but it would take ntax*log(ntax) time,
                #which wouldn't be bad, but it would also require a more complicated approach.
                #c1.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                #c1.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                #c2.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                #c2.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                c1.birth_rate = birth_rate
                c1.death_rate = death_rate
                c2.birth_rate = birth_rate
                c2.death_rate = death_rate
                #extant_tips.append(c1)
                extant_tips.append(c2)
            else:
                del extant_tips[taxI]
                if len(extant_tips) > 0:
                    extinct_tips.append(nd)
                    if is_add_extinct_attr:
                        setattr(nd, extinct_attr_name, None)
                else:
                    # total extinction
                    if (gsa_ntax is not None):
                        if (len(targetted_time_slices) > 0):
                            break
                    if not repeat_until_success:
                        raise TreeSimTotalExtinctionException()
                    # We are going to basically restart the simulation because
                    # the tree has gone extinct (without reaching the specified
                    # ntax)
                    extant_tips = list(initial_extant_tip_set)
                    extinct_tips = list(initial_extinct_tip_set)
                    ndIndex=0
                    for nd in extant_tips:
                        if is_add_extinct_attr:
                            setattr(nd, extinct_attr_name, False)
                        nd.clear_child_nodes()
                        #NICOLA: restore branch length to original value
                        nd.edge.length=initial_lengths[ndIndex]
                        ndIndex+=1
                    total_time = 0


    if gsa_ntax is not None:
        total_duration_at_target_n_tax = 0.0
        for i in targetted_time_slices:
            total_duration_at_target_n_tax += i[0]
        r = rng.random()*total_duration_at_target_n_tax
        selected_slice = None
        for n, i in enumerate(targetted_time_slices):
            r -= i[0]
            if r < 0.0:
                selected_slice = i
        assert(selected_slice is not None)
        edges_at_slice = selected_slice[1]
        last_waiting_time = selected_slice[0]

        for e, prev_length in edges_at_slice:
            daughter_nd = e.head_node
            for nd in daughter_nd.child_nodes():
                nd._parent_node = None
                try:
                    extinct_tips.remove(nd)
                except ValueError:
                    pass
                try:
                    extant_tips.remove(nd)
                except ValueError:
                    pass
                for desc in nd.preorder_iter():
                    try:
                        extant_tips.remove(desc)
                    except ValueError:
                        pass
            daughter_nd.clear_child_nodes()
            try:
                extinct_tips.remove(daughter_nd)
            except ValueError:
                pass
            extant_tips.append(daughter_nd)
            if is_add_extinct_attr:
                setattr(daughter_nd, extinct_attr_name, False)
            e.length = prev_length + last_waiting_time

    if not is_retain_extinct_tips:
        processed_nodes = set()
        for nd in list(extinct_tips):
            if nd in processed_nodes:
                continue
            processed_nodes.add(nd)
            try:
                extinct_tips.remove(nd)
            except ValueError:
                pass
            assert not nd._child_nodes
            while (nd.parent_node is not None) and (len(nd.parent_node._child_nodes) == 1):
                nd = nd.parent_node
                processed_nodes.add(nd)
            tree.prune_subtree(nd, suppress_unifurcations=False)
    tree.suppress_unifurcations()

    if is_assign_extant_taxa or is_assign_extinct_taxa:
        taxon_pool = [t for t in taxon_namespace]
        rng.shuffle(taxon_pool)
        taxon_pool_labels = set([t.label for t in taxon_pool])

        ### ONLY works if in GSA sub-section we remove ALL extant and
        ### extinct nodes beyond time slice: expensive
        ### Furthermore, main reason to use this approach is to have different
        ### label prefixes for extinct vs. extant lineages, but the second time
        ### this function is called with the same taxon namespace or any time
        ### this function is called with a populated taxon namespace, that
        ### aesthetic is lost.
        # node_pool_labels = ("T", "X")
        # for node_pool_idx, node_pool in enumerate((extant_tips, extinct_tips)):
        #     for node_idx, nd in enumerate(node_pool):
        #         if taxon_pool:
        #             taxon = taxon_pool.pop()
        #         else:
        #             taxon = taxon_namespace.require_taxon("{}{}".format(node_pool_labels[node_pool_idx], node_idx+1))
        #         nd.taxon = taxon
        #         assert not nd._child_nodes

        tlabel_counter = 0
        leaf_nodes = tree.leaf_nodes()
        rng.shuffle(leaf_nodes)
        for nd_idx, nd in enumerate(leaf_nodes):
            if not is_assign_extant_taxa and nd in extant_tips:
                continue
            if not is_assign_extant_taxa and nd in extinct_tips:
                continue
            if taxon_pool:
                taxon = taxon_pool.pop()
            else:
                while True:
                    tlabel_counter += 1
                    label = "{}{}".format("T", tlabel_counter)
                    if label not in taxon_pool_labels:
                        break
                taxon = taxon_namespace.require_taxon(label=label)
                taxon_pool_labels.add(label)
            nd.taxon = taxon
    return tree



def discrete_birth_death_tree(birth_rate, death_rate, birth_rate_sd=0.0, death_rate_sd=0.0, **kwargs):
    """
    Returns a birth-death tree with birth rate specified by ``birth_rate``, and
    death rate specified by ``death_rate``, with edge lengths in discrete (integer)
    units.

    ``birth_rate_sd`` is the standard deviation of the normally-distributed mutation
    added to the birth rate as it is inherited by daughter nodes; if 0, birth
    rate does not evolve on the tree.

    ``death_rate_sd`` is the standard deviation of the normally-distributed mutation
    added to the death rate as it is inherited by daughter nodes; if 0, death
    rate does not evolve on the tree.

    Tree growth is controlled by one or more of the following arguments, of which
    at least one must be specified:

        - If ``ntax`` is given as a keyword argument, tree is grown until the number of
          tips == ntax.
        - If ``taxon_namespace`` is given as a keyword argument, tree is grown until the
          number of tips == len(taxon_namespace), and the taxa are assigned randomly to the
          tips.
        - If 'max_time' is given as a keyword argument, tree is grown for ``max_time``
          number of generations.

    If more than one of the above is given, then tree growth will terminate when
    *any* of the termination conditions (i.e., number of tips == ``ntax``, or number
    of tips == len(taxon_namespace) or number of generations = ``max_time``) are met.

    Also accepts a Tree object (with valid branch lengths) as an argument passed
    using the keyword ``tree``: if given, then this tree will be used; otherwise
    a new one will be created.

    If ``assign_taxa`` is False, then taxa will *not* be assigned to the tips;
    otherwise (default), taxa will be assigned. If ``taxon_namespace`` is given
    (``tree.taxon_namespace``, if ``tree`` is given), and the final number of tips on the
    tree after the termination condition is reached is less then the number of
    taxa in ``taxon_namespace`` (as will be the case, for example, when
    ``ntax`` < len(``taxon_namespace``)), then a random subset of taxa in ``taxon_namespace`` will
    be assigned to the tips of tree. If the number of tips is more than the number
    of taxa in the ``taxon_namespace``, new Taxon objects will be created and added
    to the ``taxon_namespace`` if the keyword argument ``create_required_taxa`` is not given as
    False.

    Under some conditions, it is possible for all lineages on a tree to go extinct.
    In this case, if the keyword argument ``repeat_until_success`` is |True|, then a new
    branching process is initiated.
    If |False| (default), then a TreeSimTotalExtinctionException is raised.

    A Random() object or equivalent can be passed using the ``rng`` keyword;
    otherwise GLOBAL_RNG is used.
    """
    if 'ntax' not in kwargs \
        and 'taxon_namespace' not in kwargs \
        and 'max_time' not in kwargs:
            raise ValueError("At least one of the following must be specified: 'ntax', 'taxon_namespace', or 'max_time'")
    target_num_taxa = None
    taxon_namespace = None
    target_num_gens = kwargs.get('max_time', None)
    if 'taxon_namespace' in kwargs:
        taxon_namespace = kwargs.get('taxon_namespace')
        target_num_taxa = kwargs.get('ntax', len(taxon_namespace))
    elif 'ntax' in kwargs:
        target_num_taxa = kwargs['ntax']
    if taxon_namespace is None:
        taxon_namespace = dendropy.TaxonNamespace()
    repeat_until_success = kwargs.get('repeat_until_success', False)
    rng = kwargs.get('rng', GLOBAL_RNG)

    # grow tree
    if "tree" in kwargs:
        tree = kwargs['tree']
        if "taxon_namespace" in kwargs and kwargs['taxon_namespace'] is not tree.taxon_namespace:
            raise ValueError("Cannot specify both ``tree`` and ``taxon_namespace``")
    else:
        tree = dendropy.Tree(taxon_namespace=taxon_namespace)
        tree.is_rooted = True
        tree.seed_node.edge.length = 0
        tree.seed_node.birth_rate = birth_rate
        tree.seed_node.death_rate = death_rate
    leaf_nodes = tree.leaf_nodes()
    num_gens = 0
    while (target_num_taxa is None or len(leaf_nodes) < target_num_taxa) \
            and (target_num_gens is None or num_gens < target_num_gens):
        for nd in leaf_nodes:
            if not hasattr(nd, 'birth_rate'):
                nd.birth_rate = birth_rate
            if not hasattr(nd, 'death_rate'):
                nd.death_rate = death_rate
            try:
                nd.edge.length += 1
            except TypeError:
                nd.edge.length = 1
            u = rng.uniform(0, 1)
            if u < nd.birth_rate:
                c1 = nd.new_child()
                c2 = nd.new_child()
                c1.edge.length = 0
                c2.edge.length = 0
                c1.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c1.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                c2.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c2.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
            elif u > nd.birth_rate and u < (nd.birth_rate + nd.death_rate):
                if nd is not tree.seed_node:
                    tree.prune_subtree(nd)
                elif not repeat_until_success:
                    # all lineages are extinct: raise exception
                    raise TreeSimTotalExtinctionException()
                else:
                    # all lineages are extinct: repeat
                    num_gens = 0

        num_gens += 1
        leaf_nodes = tree.leaf_nodes()

    # If termination condition specified by ntax or taxon_namespace, then the last
    # split will have a daughter edges of length == 0;
    # so we continue growing the edges until the next birth/death event *or*
    # the max number of generations condition is given and met
    gens_to_add = 0
    while (target_num_gens is None or num_gens < target_num_gens):
        u = rng.uniform(0, 1)
        if u < (birth_rate + death_rate):
            break
        gens_to_add += 1
    for nd in tree.leaf_nodes():
        nd.edge.length += gens_to_add

    if kwargs.get("assign_taxa", True):
        tree.randomly_assign_taxa(create_required_taxa=True, rng=rng)

    # return
    return tree

def uniform_pure_birth_tree(taxon_namespace, birth_rate=1.0, rng=None):
    "Generates a uniform-rate pure-birth process tree. "
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    tree.seed_node.edge.length = 0.0
    leaf_nodes = tree.leaf_nodes()
    while len(leaf_nodes) < len(taxon_namespace):
        waiting_time = rng.expovariate(len(leaf_nodes)/birth_rate)
        for nd in leaf_nodes:
            nd.edge.length += waiting_time
        parent_node = rng.choice(leaf_nodes)
        c1 = parent_node.new_child()
        c2 = parent_node.new_child()
        c1.edge.length = 0.0
        c2.edge.length = 0.0
        leaf_nodes = tree.leaf_nodes()
    leaf_nodes = tree.leaf_nodes()
    waiting_time = rng.expovariate(len(leaf_nodes)/birth_rate)
    for nd in leaf_nodes:
        nd.edge.length += waiting_time
    for idx, leaf in enumerate(leaf_nodes):
        leaf.taxon = taxon_namespace[idx]
    tree.is_rooted = True
    return tree

def fit_pure_birth_model(**kwargs):
    r"""
    Calculates the maximum-likelihood estimate of the birth rate of a set of *internal* node ages under a Yule (pure-birth) model.

    Requires either a |Tree| object or an interable of *internal* node ages to be passed in via keyword arguments ``tree`` or ``internal_node_ages``, respectively. The former is more convenient when doing one-off calculations, while the latter is more efficient if the list of internal node ages needs to be used in other places and you already have it calculated and want to avoid re-calculating it here.

    Parameters
    ----------
    \*\*kwargs : keyword arguments, mandatory

        Exactly *one* of the following *must* be specified:

            tree : a |Tree| object.
                A |Tree| object. The tree needs to be ultrametric for the internal node ages (time from each internal node to the tips) to make sense. The precision by which the ultrametricity is checked can be specified using the ``ultrametricity_precision`` keyword argument (see below). If ``tree`` is given, then ``internal_node_ages`` cannot be given, and vice versa. If ``tree`` is not given, then ``internal_node_ages`` must be given.
            internal_node_ages : iterable (of numerical values)
                Iterable of node ages of the internal nodes of a tree, i.e., the list of sum of the edge lengths between each internal node and the tips of the tree. If ``internal_node_ages`` is given, then ``tree`` cannot be given, and vice versa. If ``internal_node_ages`` is not given, then ``tree`` must be given.

        The following are optional, and are only used if internal node ages are
        specified (i.e., 'internal_node_ages' are passed in):

            is_node_ages_presorted : bool
                By default, the vector of node ages are sorted. If this
                argument is specified as ``True``, then this sorting will be
                skipped, in which case it is the client code's responsibility
                to make sure that the node ages are given in REVERSE order
                (i.e., oldest nodes -- nodes closer to the root -- given
                first).
    Returns
    -------
    m : dictionary

    A dictionary with keys being parameter names and values being
    estimates:

        "birth_rate"
            The birth rate.
        "log_likelihood"
            The log-likelihood of the model and given birth rate.

    Examples
    --------

    Birth rates can be estimated by passing in trees directly::

        for idx, tree in enumerate(trees):
            m = birthdeath.fit_pure_birth_model(tree=tree)
            print("Tree {}: birth rate = {} (logL = {})".format(
                idx+1, m["birth_rate"], m["log_likelihood"]))

    Or by pre-calculating and passing in a list of node ages::

        for idx, tree in enumerate(trees):
            m = birthdeath.fit_pure_birth_model(
                    internal_node_ages=tree.internal_node_ages())
            print("Tree {}: birth rate = {} (logL = {})".format(
                idx+1, m["birth_rate"], m["log_likelihood"]))


    Notes
    -----
    Adapted from the laser package for R:

        -   Dan Rabosky and Klaus Schliep (2013). laser: Likelihood Analysis of
            Speciation/Extinction Rates from Phylogenies. R package version
            2.4-1. http://CRAN.R-project.org/package=laser

    """
    tree = kwargs.get("tree", None)
    if tree is not None:
        internal_node_ages = tree.internal_node_ages(ultrametricity_precision=kwargs.get("ultrametricity_precision", 1e-6))
    else:
        try:
            internal_node_ages = kwargs["internal_node_ages"]
        except KeyError:
            raise TypeError("Need to specify 'tree' or 'internal_node_ages'")
    x = internal_node_ages
    if tree is not None or not kwargs.get("is_node_ages_presorted", False):
        x = sorted(internal_node_ages, reverse=True)
    st1 = x[0]
    st2 = 0
    nvec = range(2, len(x)+2)
    nv = [i for i in x if (i < st1) and (i >= st2)]
    lo = max(nvec[idx] for idx, i in enumerate(x) if i >= st1)
    up = max(nvec[idx] for idx, i in enumerate(x) if i >= st2)
    if st1 <= x[0]:
        nv.insert(0, st1)
        nv = [i - st2 for i in nv]
    else:
        nv = [i - st2 for i in nv]
    t1 = (up-lo)
    t2 = (lo*nv[0])
    t3 = sum( nv[1:(up-lo+1)] )
    smax = t1/(t2 + t3)

    try:
        s1 = sum(map(math.log, range(lo,up)))
        s2 = (up-lo) * math.log(smax)
        s3 = lo - up
        lh = s1 + s2 + s3
    except ValueError:
        if kwargs.get("ignore_likelihood_calculation_failure", False):
            raise ValueError("Likelihood estimation failure")
        else:
            lh = float("-inf")

    result = {
        "birth_rate" : smax,
        "log_likelihood" : lh,
    }
    return result

def fit_pure_birth_model_to_tree(tree, ultrametricity_precision=constants.DEFAULT_ULTRAMETRICITY_PRECISION):
    """
    Calculates the maximum-likelihood estimate of the birth rate a tree under a
    Yule (pure-birth) model.

    Parameters
    ----------
    tree : |Tree| object
        A tree to be fitted.

    Returns
    -------
    m : dictionary

    A dictionary with keys being parameter names and values being
    estimates:

        -   "birth_rate"
            The birth rate.
        -   "log_likelihood"
            The log-likelihood of the model and given birth rate.

    Examples
    --------

    ::

        import dendropy
        from dendropy.model import birthdeath
        trees = dendropy.TreeList.get_from_path(
                "pythonidae.nex", "nexus")
        for idx, tree in enumerate(trees):
            m = birthdeath.fit_pure_birth_model_to_tree(tree)
            print("Tree {}: birth rate = {} (logL = {})".format(
                idx+1, m["birth_rate"], m["log_likelihood"]))

    """
    return fit_pure_birth_model(tree=tree, ultrametricity_precision=ultrametricity_precision)


def birth_death_likelihood(**kwargs):
    r"""
    Calculates the log-likelihood of a tree (or a set of internal nodes) under
    a birth death model.

    Requires either a |Tree| object or an interable of *internal* node
    ages to be passed in via keyword arguments ``tree`` or ``internal_node_ages``,
    respectively. The former is more convenient when doing one-off
    calculations, while the latter is more efficient if the list of internal
    node ages needs to be used in other places and you already have it
    calculated and want to avoid re-calculating it here.

    Parameters
    ----------
    \*\*kwargs : keyword arguments, mandatory

        Exactly *one* of the following *must* be specified:

            tree : a |Tree| object.
                A |Tree| object. The tree needs to be ultrametric for the
                internal node ages (time from each internal node to the tips)
                to make sense. The precision by which the ultrametricity is
                checked can be specified using the ``ultrametricity_precision`` keyword
                argument (see below). If ``tree`` is given, then
                ``internal_node_ages`` cannot be given, and vice versa. If ``tree``
                is not given, then ``internal_node_ages`` must be given.
            internal_node_ages : iterable (of numerical values)
                Iterable of node ages of the internal nodes of a tree, i.e., the
                list of sum of the edge lengths between each internal node and
                the tips of the tree. If ``internal_node_ages`` is given, then
                ``tree`` cannot be given, and vice versa. If ``internal_node_ages``
                is not given, then ``tree`` must be given.

        The following keyword parameters are mandatory:

            birth_rate : float
                The birth rate.
            death_rate : float
                The death rate.

        The following keyword parameters are optional:

            sampling_probability
                The probability for a species to be included in the sample. Defaults to 1.0 (all species sampled).
            sampling_strategy
                The strategy how samples were obtained. Options are: uniform|diversified|age.
            is_mrca_included
                Does the process start with the most recent common ancestor?
            condition_on : string
                Do we condition the process on: "time", "survival", or "taxa"?

        The following are optional, and are only used if internal node ages
        need to be calculated (i.e., 'tree' is passed in).

            ultrametricity_precision : float
                When calculating the node ages, an error will be raised if the tree in
                o ultrametric. This error may be due to floating-point or numerical
                imprecision. You can set the precision of the ultrametricity validation
                by setting the ``ultrametricity_precision`` parameter. E.g., use
                ``ultrametricity_precision=0.01`` for a more relaxed precision,
                down to 2 decimal places. Use ``ultrametricity_precision=False``
                to disable checking of ultrametricity precision.

            ignore_likelihood_calculation_failure: bool (default: False)
                In some cases (typically, abnormal trees, e.g., 1-tip), the
                likelihood estimation will fail. In this case a ValueError will
                be raised. If ``ignore_likelihood_calculation_failure`` is
                |True|, then the function call will still succeed, with the
                likelihood set to -``inf``.

        The following are optional, and are only used if internal node ages are
        specified (i.e., 'internal_node_ages' are passed in):

            is_node_ages_presorted : bool
                By default, the vector of node ages are sorted. If this
                argument is specified as ``True``, then this sorting will be
                skipped, in which case it is the client code's responsibility
                to make sure that the node ages are given in REVERSE order
                (i.e., oldest nodes -- nodes closer to the root -- given
                first).

    Notes
    -----
    Lifted directly from the (fantastic!) TESS package for R:

        H{"o}hna S. 2013. Fast simulation of reconstructed phylogenies under
        global time-dependent birth--death processes. Bioinformatics, 29(11)
        1367-1374.

    Returns
    -------
    lnl : float

    The log-likehood of the tree under the birth-death model.

    """
    tree = kwargs.get("tree", None)
    if tree is not None:
        internal_node_ages = tree.internal_node_ages(ultrametricity_precision=kwargs.get("ultrametricity_precision", 1e-6))
    else:
        try:
            internal_node_ages = kwargs["internal_node_ages"]
        except KeyError:
            raise TypeError("Need to specify 'tree' or 'internal_node_ages'")
    if tree is not None or not kwargs.get("is_node_ages_presorted", False):
        internal_node_ages = sorted(internal_node_ages, reverse=True)

    # check for sensible parameter values
    # if ( lambda <= 0 || mu < 0 || sampling_probability <= 0 || sampling_probability > 1.0) {
    #     stop("Invalid parameter values for lambda and mu!")
    # }

    # massExtinctionTimes=c(),
    # massExtinctionSurvivalProbabilities=c(),
    # missingSpecies = c(),
    # timesMissingSpecies = c(),
    # tess.likelihood 11
    # samplingProbability=1.0,
    # samplingStrategy="uniform",
    # MRCA=TRUE,
    # CONDITION="survival",
    birth_rate = kwargs.get("birth_rate")
    death_rate = kwargs.get("death_rate")
    massExtinctionTimes = None
    massExtinctionSurvivalProbabilities = None
    sampling_probability = kwargs.get("sampling_probability", 1.0)
    sampling_strategy = kwargs.get("sampling_strategy", "uniform")
    is_mrca_included = kwargs.get("is_mrca_included", True)
    condition_on = kwargs.get("condition_on", "survival")

    ntax = len(internal_node_ages) + 1
    PRESENT = max(internal_node_ages)
    times = [(PRESENT - t) for t in internal_node_ages]

    # if we condition on the MRCA, then we need to remove the root speciation event
    if is_mrca_included:
        times = times[1:]

    # set the uniform taxon sampling probability
    if sampling_strategy == "uniform":
        rho = sampling_probability
    else:
        rho = 1.0

    # initialize the log likelihood
    lnl = 0.0

    # what do we condition on?
    # did we condition on survival?
    if condition_on == "survival":
        lnl = - _p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,PRESENT)

    # death_rateltiply the probability of a descendant of the initial species
    lnl = lnl + _p1_constant(
            birth_rate=birth_rate,
            death_rate=death_rate,
            massExtinctionTimes=massExtinctionTimes,
            massExtinctionSurvivalProbabilities=massExtinctionSurvivalProbabilities,
            sampling_probability=rho,
            t=0,
            T=PRESENT)

    # add the survival of a second species if we condition on the MRCA
    if is_mrca_included:
        lnl = 2 * lnl

    # did we condition on observing n species today
    if condition_on == "taxa":
        lnl = lnl - _p_N_constant(
                birth_rate=birth_rate,
                death_rate=death_rate,
                massExtinctionTimes=massExtinctionTimes,
                massExtinctionSurvivalProbabilities=massExtinctionSurvivalProbabilities,
                sampling_probability=sampling_probability,
                i=ntax,
                s=0,
                t=PRESENT,
                SURVIVAL=False,
                MRCA=is_mrca_included)

    # if we assume diversified sampling, we need to death_rateltiply with the
    # probability that all missing species happened after the last speciation
    # event
    if sampling_strategy == "diversified":
        # We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
        lastEvent = times[-1]
        p_0_T = 1.0 - math.exp(_p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,0,PRESENT,PRESENT)) * math.exp((death_rate-birth_rate)*PRESENT)
        p_0_t = 1.0 - math.exp(_p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,lastEvent,PRESENT,PRESENT)) * math.exp((death_rate-birth_rate)*(PRESENT-lastEvent))
        F_t = p_0_t / p_0_T
        # get an estimate of the actual number of taxa
        m = round(float(ntax) / sampling_probability)
        # remove the number of species that we started with
        if is_mrca_included:
            k = 2
        else:
            k = 1
        lnl = lnl + (m-ntax) * math.log(F_t) + math.log(combinatorics.choose(int(m-k),int(ntax-k)))

    # multiply the probability for the missing species
    # if len(missing_species) > 0:
    #     # We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
    #     # now iterate over the vector of missing species per interval
    #     lastEvent = timesMissingSpecies
    #     p_0_T = 1.0 - math.exp( _p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,0,PRESENT,PRESENT,log=TRUE) + ((death_rate-birth_rate)*PRESENT) )
    #     p_0_t = 1.0 - math.exp( _p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,lastEvent,PRESENT,PRESENT,log=TRUE) + ((death_rate-birth_rate)*(PRESENT-lastEvent)) )
    #     log_F_t = math.log(p_0_t) - math.log(p_0_T)
    #     # get an estimate of the actual number of taxa
    #     m = missingSpecies
    #     # remove the number of species that we started with
    #     lnl = lnl + sum( m * log_F_t ) #+ lchoose(m-k,ntax-k)


    # multiply the probability for each speciation time
    if len(times) > 0:
        # lnl = lnl + len(times)*math.log(birth_rate) + sum(tess.equations.p1.constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,times,PRESENT,log=TRUE))
        lnl += (len(times) * math.log(birth_rate))
        for time in times:
            lnl += _p1_constant(
                    birth_rate=birth_rate,
                    death_rate=death_rate,
                    massExtinctionTimes=massExtinctionTimes,
                    massExtinctionSurvivalProbabilities=massExtinctionSurvivalProbabilities,
                    sampling_probability=rho,
                    t=time,
                    T=PRESENT)

    return lnl



################################################################################
# From: TESS
#
#   https://github.com/hoehna/TESS
#
#   H{"o}hna S. 2013. Fast simulation of reconstructed phylogenies under
#   global time-dependent birth--death processes. Bioinformatics, 29(11)
#   1367-1374.
#
#   Copyright (c) 2012- Sebastian Hoehna
#
#   TESS is free software; you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as
#   published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 (or 2) lineage at time s.
#
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    i                                             scalar        number of lineages
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        log-probability
#
################################################################################
# tess.equations.pN.constant = function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,i,s,t,SURVIVAL=FALSE,MRCA=FALSE,log=FALSE) {
def _p_N_constant(
        birth_rate,
        death_rate,
        massExtinctionTimes,
        massExtinctionSurvivalProbabilities,
        sampling_probability,
        i,
        s,
        t,
        SURVIVAL=False,
        MRCA=False,):
    if i < 1: # we assume conditioning on survival
        p = 0
    elif i == 1:
        if MRCA: # we assume conditioning on survival of the two species
            p = 0
        else:
            if SURVIVAL:
                p =  _p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,sampling_probability,s,t,t) + (death_rate-birth_rate)*(t-s) - math.log(sampling_probability)
            else:
                p =  2*_p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,sampling_probability,s,t,t) + (death_rate-birth_rate)*(t-s) - math.log(sampling_probability)
    else:
        p_s = math.exp(_p_survival_constant(birth_rate,death_rate,massExtinctionTimes,massExtinctionSurvivalProbabilities,sampling_probability,s,t,t))
        r   = (death_rate-birth_rate)*(t-s) - math.log(sampling_probability)
        # for (j in seq_len(length(massExtinctionTimes)) ) {
        #     cond =  (s < massExtinctionTimes[j]) & (t >= massExtinctionTimes[j])
        #     r  = r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
        e = p_s * math.exp(r)
        if e > 1:
            e = 1.0

        if not MRCA:
            if SURVIVAL:
                p   = math.log(p_s) + r + math.log( 1 - e) * (i-1)
            else:
                p   = 2*math.log(p_s) + r + math.log( 1 - e) * (i-1)
        else:
            if SURVIVAL:
                p = math.log(i-1) + 2*math.log(p_s) + 2*r + math.log( 1 - e) * (i-2)
            else:
                p = math.log(i-1) + 4*math.log(p_s) + 2*r + math.log( 1 - e) * (i-2)

    return p

################################################################################
# From: TESS
#
#   https://github.com/hoehna/TESS
#
#   H{"o}hna S. 2013. Fast simulation of reconstructed phylogenies under
#   global time-dependent birth--death processes. Bioinformatics, 29(11)
#   1367-1374.
#
#   Copyright (c) 2012- Sebastian Hoehna
#
#   TESS is free software; you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as
#   published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
# @brief  Calculate the probability of survival in the interval [t_low,t_high].
#
# See Equation (2) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
#
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    t_low                                         scalar        starting time
# @param    t_high                                        scalar        end time
# @param    T                                             scalar        present time (time goes forward and the origin/MRCA might be 0)
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        probability of survival in [t,tau]
#
################################################################################
# tess.equations.pSurvival.constant = function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,t_high,T,log=FALSE) {
def _p_survival_constant(
        birth_rate,
        death_rate,
        massExtinctionTimes,
        massExtinctionSurvivalProbabilities,
        sampling_probability,
        t_low,
        t_high,
        T):

    # compute the rate
    rate = death_rate - birth_rate

    # do the integration of int_{t_low}^{t_high} ( death_rate(s) exp(rate(t,s)) ds )
    # where rate(t,s) = int_{t}^{s} ( death_rate(x)-birth_rate(x) dx ) - sum_{for all t < m_i < s in massExtinctionTimes }( log(massExtinctionSurvivalProbability[i]) )

    # we compute the integral stepwise for each epoch between mass-extinction events
    # add mass-extinction
    accudeath_ratelatedMassExtinction = 1.0
    prev_time = t_low
    den = 1.0
    # if ( length(massExtinctionTimes) > 0 ) {
    #     for (j in 1:length(massExtinctionTimes) ) {
    #         cond =  (t_low < massExtinctionTimes[j]) & (t_high >= massExtinctionTimes[j])
    #         # compute the integral for this time episode until the mass-extinction event
    #     #       den = den + ifelse(cond, exp(-rate*t_low) * death_rate / (rate * accudeath_ratelatedMassExtinction ) * ( exp(rate* massExtinctionTimes[j]) - exp(rate*prev_time)) , 0 )
    #         den = den + cond * exp(-rate*t_low) * death_rate / (rate * accudeath_ratelatedMassExtinction ) * ( exp(rate* massExtinctionTimes[j]) - exp(rate*prev_time))
    #         # store the current time so that we remember from which episode we need to integrate next
    #     #       prev_time = ifelse(cond, massExtinctionTimes[j], prev_time)
    #         prev_time[cond] = massExtinctionTimes[j]
    #         accudeath_ratelatedMassExtinction = accudeath_ratelatedMassExtinction * ifelse(cond, massExtinctionSurvivalProbabilities[j], 1.0)
    #         # integrate over the tiny time interval of the mass-extinction event itself and add it to the integral
    #     #       den = den - ifelse(cond, (massExtinctionSurvivalProbabilities[j]-1) / accudeath_ratelatedMassExtinction * exp( rate*(massExtinctionTimes[j] - t_low) ), 0.0 )
    #         den = den - cond * (massExtinctionSurvivalProbabilities[j]-1) / accudeath_ratelatedMassExtinction * exp( rate*(massExtinctionTimes[j] - t_low) )
    #     }
    # }

    # add the integral of the final epoch until the present time
    den = den + math.exp(-rate*t_low) * death_rate / (rate * accudeath_ratelatedMassExtinction ) * ( math.exp(rate*t_high) - math.exp(rate*prev_time))

    # add sampling
    if (sampling_probability < 1) and (t_low < T) and (t_high >= T):
        accudeath_ratelatedMassExtinction = accudeath_ratelatedMassExtinction *  sampling_probability
        cond = 1
    else:
        cond = 0

    #   den = den - ifelse(cond, (sampling_probability-1)*math.exp( rate*(T-t_low) ) / accudeath_ratelatedMassExtinction, 0.0)
    den = den - cond * (sampling_probability-1)*math.exp( rate*(T-t_low) ) / accudeath_ratelatedMassExtinction

    res = 1.0 / den
    res = math.log(res)

    return (res)


################################################################################
# From: TESS
#
#   https://github.com/hoehna/TESS
#
#   H{"o}hna S. 2013. Fast simulation of reconstructed phylogenies under
#   global time-dependent birth--death processes. Bioinformatics, 29(11)
#   1367-1374.
#
#   Copyright (c) 2012- Sebastian Hoehna
#
#   TESS is free software; you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as
#   published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
# @brief  Calculate the probability of exactly 1 lineage surviving until time T
#         if we started with 1 lineage at time t.
#
# @see    Equation (4) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    t                                             scalar        time
# @param    T                                             scalar        present time
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        ln-probability
#
################################################################################
# tess.equations.p1.constant = function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,log=FALSE) {
def _p1_constant(
        birth_rate,
        death_rate,
        massExtinctionTimes,
        massExtinctionSurvivalProbabilities,
        sampling_probability,
        t,
        T):
    # compute the survival probability
    a = _p_survival_constant(
            birth_rate=birth_rate,
            death_rate=death_rate,
            massExtinctionTimes=massExtinctionTimes,
            massExtinctionSurvivalProbabilities=massExtinctionSurvivalProbabilities,
            sampling_probability=sampling_probability,
            t_low=t,
            t_high=T,
            T=T)
    # compute the rate
    rate = (death_rate - birth_rate)*(T-t)
    # add mass-extinction
    # for (j in range(len(massExtinctionTimes)) ) :
    #     rate = rate - ifelse( t < massExtinctionTimes[j] & T >= massExtinctionTimes[j], log(massExtinctionSurvivalProbabilities[j]), 0 )
    # add sampling
    rate = rate - math.log(sampling_probability)
    p = 2*a + rate
    return p
