#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tree simulation and generation.
"""

import sys
import copy
import math

from dendropy.utility import GLOBAL_RNG
from dendropy.utility import probability
from dendropy import coalescent
from dendropy import dataobject
from dendropy import treemanip
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)


class TreeSimTotalExtinctionException(Exception):
    """Exception to be raised when branching process results in all lineages going extinct."""
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

def star_tree(taxon_set):
    "Builds and returns a star tree from the given taxa block."
    star_tree = dataobject.Tree(taxon_set=taxon_set)
    for taxon in taxon_set:
        star_tree.seed_node.new_child(taxon=taxon)
    return star_tree

def birth_death(birth_rate, death_rate, birth_rate_sd=0.0, death_rate_sd=0.0, **kwargs):
    """
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
        - If `gsa_ntax` is given then the tree will be simulated up to this number of
          tips (or 0 tips), then a tree will be randomly selected from the
          intervals which corresond to times at which the tree had exactly `ntax`
          leaves (or len(taxon_set) tips). This allows for simulations according to
          the "General Sampling Approach" of [citeHartmannWS2010]_


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

    Under some conditions, it is possible for all lineages on a tree to go extinct.
    In this case, if the keyword argument `repeat_until_success` is `True` (the
    default), then a new branching process is initiated.
    If `False` (default), then a TreeSimTotalExtinctionException is raised.

    A Random() object or equivalent can be passed using the `rng` keyword;
    otherwise GLOBAL_RNG is used.

    .. [citeHartmannWS2010] Hartmann, Wong, and Stadler "Sampling Trees from Evolutionary Models" Systematic Biology. 2010. 59(4). 465-476

    """
    target_num_taxa = kwargs.get('ntax')
    max_time = kwargs.get('max_time')
    taxon_set = kwargs.get('taxon_set')
    if (target_num_taxa is None) and (taxon_set is not None):
        target_num_taxa = len(taxon_set)
    elif taxon_set is None:
        taxon_set = dataobject.TaxonSet()
    gsa_ntax = kwargs.get('gsa_ntax')
    terminate_at_full_tree = False
    if target_num_taxa is None:
        if gsa_ntax is not None:
            raise ValueError("When 'gsa_ntax' is used, either 'ntax' or 'taxon_set' must be used")
        if max_time is None:
            raise ValueError("At least one of the following must be specified: 'ntax', 'taxon_set', or 'max_time'")
    else:
        if gsa_ntax is None:
            terminate_at_full_tree = True
            gsa_ntax = 1 + target_num_taxa
        elif gsa_ntax < target_num_taxa:
            raise ValueError("gsa_ntax must be greater than target_num_taxa")
    repeat_until_success = kwargs.get('repeat_until_success', True)
    rng = kwargs.get('rng', GLOBAL_RNG)

    # initialize tree
    if "tree" in kwargs:
        tree = kwargs['tree']
        if "taxon_set" in kwargs and kwargs['taxon_set'] is not tree.taxon_set:
            raise ValueError("Cannot specify both `tree` and `taxon_set`")
    else:
        tree = dataobject.Tree(taxon_set=taxon_set)
        tree.is_rooted = True
        tree.seed_node.edge.length = 0.0
        tree.seed_node.birth_rate = birth_rate
        tree.seed_node.death_rate = death_rate

    # grow tree
    leaf_nodes = tree.leaf_nodes()
    #_LOG.debug("Will generate a tree with no more than %s leaves to get a tree of %s leaves" % (str(gsa_ntax), str(target_num_taxa)))
    curr_num_leaves = len(leaf_nodes)
    total_time = 0
    # for the GSA simulations targetted_time_slices is a list of tuple
    #   the first element in the tuple is the duration of the amount
    #   that the simulation spent at the (targetted) number of taxa
    #   and a list of edge information. The list of edge information includes
    #   a list of terminal edges in the tree and the length for that edge
    #   that marks the beginning of the time slice that corresponds to the
    #   targetted number of taxa.

    targetted_time_slices = []
    extinct_tips = []
    while True:
        if gsa_ntax is None:
            assert (max_time is not None)
            if total_time >= max_time:
                break
        elif curr_num_leaves >= gsa_ntax:
            break

        # get vector of birth/death probabilities, and
        # associate with nodes/events
        event_rates = []
        event_nodes = []
        for nd in leaf_nodes:
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
        #_LOG.debug("rate_of_any_event = %f" % (rate_of_any_event))
        waiting_time = rng.expovariate(rate_of_any_event)
        #_LOG.debug("Drew waiting time of %f from hazard parameter of %f" % (waiting_time, rate_of_any_event))

        if (gsa_ntax is not None) and (curr_num_leaves == target_num_taxa):
            edge_and_start_length = []
            for nd in leaf_nodes:
                e = nd.edge
                edge_and_start_length.append((e, e.length))
            targetted_time_slices.append((waiting_time, edge_and_start_length))
            #_LOG.debug("Recording slice with %d edges" % len(edge_and_start_length))
            if terminate_at_full_tree:
                break

        # add waiting time to nodes
        for nd in leaf_nodes:
            try:
                nd.edge.length += waiting_time
            except TypeError:
                nd.edge.length = waiting_time
        #_LOG.debug("Next waiting_time = %f" % waiting_time)
        total_time += waiting_time

        # if event occurs within time constraints
        if max_time is None or total_time <= max_time:

            # normalize probability
            for i in xrange(len(event_rates)):
                event_rates[i] = event_rates[i]/rate_of_any_event

            # select node/event and process
            nd, birth_event = probability.weighted_choice(event_nodes, event_rates, rng=rng)
            leaf_nodes.remove(nd)
            curr_num_leaves -= 1
            if birth_event:
                #_LOG.debug("Speciation")
                c1 = nd.new_child()
                c2 = nd.new_child()
                c1.edge.length = 0
                c2.edge.length = 0
                c1.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c1.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                c2.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c2.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                leaf_nodes.append(c1)
                leaf_nodes.append(c2)
                curr_num_leaves += 2
            else:
                #_LOG.debug("Extinction")
                if curr_num_leaves > 0:
                    #_LOG.debug("Will delete " + str(id(nd)) + " with parent = " + str(id(nd.parent_node)))
                    extinct_tips.append(nd)
                else:
                    if (gsa_ntax is not None):
                        if (len(targetted_time_slices) > 0):
                            break
                    if not repeat_until_success:
                        raise TreeSimTotalExtinctionException()
                    # We are going to basically restart the simulation because the tree has gone extinct (without reaching the specified ntax)
                    leaf_nodes = [tree.seed_node]
                    curr_num_leaves = 1
                    for nd in tree.seed_node.child_nodes():
                        treemanip.prune_subtree(tree, nd, delete_outdegree_one=False)
                    extinct_tips = []
                    total_time = 0
            assert(curr_num_leaves == len(leaf_nodes))
            #_LOG.debug("Current tree \n%s" % (tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True)))
    #tree._debug_tree_is_valid()
    #_LOG.debug("Terminated with %d leaves (%d, %d  according to len(leaf_nodes))" % (curr_num_leaves, len(leaf_nodes), len(tree.leaf_nodes())))
    if gsa_ntax is not None:
        total_duration_at_target_n_tax = 0.0
        for i in targetted_time_slices:
            total_duration_at_target_n_tax += i[0]
        r = rng.random()*total_duration_at_target_n_tax
        #_LOG.debug("Selected rng = %f out of (0, %f)" % (r, total_duration_at_target_n_tax))
        selected_slice = None
        for n, i in enumerate(targetted_time_slices):
            r -= i[0]
            if r < 0.0:
                selected_slice = i
        assert(selected_slice is not None)
        #_LOG.debug("Selected time slice index %d" % n)
        edges_at_slice = selected_slice[1]
        last_waiting_time = selected_slice[0]
        for e, prev_length in edges_at_slice:
            daughter_nd = e.head_node
            for nd in daughter_nd.child_nodes():
                treemanip.prune_subtree(tree, nd, delete_outdegree_one=False)
                #_LOG.debug("After pruning %s:\n%s" % (str(id(nd)), tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True)))
                try:
                    extinct_tips.remove(nd)
                except:
                    pass
            try:
                extinct_tips.remove(daughter_nd)
            except:
                pass
            e.length = prev_length + last_waiting_time



#     tree._debug_tree_is_valid()
#     for nd in extinct_tips:
#         _LOG.debug("Will be deleting " + str(id(nd)))

    for nd in extinct_tips:
        bef = len(tree.leaf_nodes())
        while (nd.parent_node is not None) and (len(nd.parent_node.child_nodes()) == 1):
            _LOG.debug("Will be pruning %d rather than its only child (%d)" % (id(nd.parent_node), id(nd)))
            nd = nd.parent_node
#         _LOG.debug("Deleting " + str(nd.__dict__) + '\n' + str(nd.edge.__dict__))
#         for n, pnd in enumerate(tree.postorder_node_iter()):
#             _LOG.debug("%d %s" % (n, repr(pnd)))
#        _LOG.debug("Before prune of %s:\n%s" % (str(id(nd)), tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True)))
        if nd.parent_node:
            treemanip.prune_subtree(tree, nd, delete_outdegree_one=False)
        _LOG.debug("After prune (went from %d to %d leaves):\n%s" % (bef, len(tree.leaf_nodes()), tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True)))
#         _LOG.debug("Deleted " + str(nd.__dict__))
#         for n, pnd in enumerate(tree.postorder_node_iter()):
#             _LOG.debug("%d %s" % (n, repr(pnd)))
#         tree._debug_tree_is_valid()
    tree.delete_outdegree_one_nodes()
#    tree._debug_tree_is_valid()
#    _LOG.debug("After deg2suppression:\n%s" % (tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True)))

    if kwargs.get("assign_taxa", True):
        tree.randomly_assign_taxa(create_required_taxa=True, rng=rng)


    # return
    return tree

def discrete_birth_death(birth_rate, death_rate, birth_rate_sd=0.0, death_rate_sd=0.0, **kwargs):
    """
    Returns a birth-death tree with birth rate specified by `birth_rate`, and
    death rate specified by `death_rate`, with edge lengths in discrete (integer)
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
        - If 'max_time' is given as a keyword argument, tree is grown for `max_time`
          number of generations.

    If more than one of the above is given, then tree growth will terminate when
    *any* of the termination conditions (i.e., number of tips == `ntax`, or number
    of tips == len(taxon_set) or number of generations = `max_time`) are met.

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

    Under some conditions, it is possible for all lineages on a tree to go extinct.
    In this case, if the keyword argument `repeat_until_success` is `True`, then a new
    branching process is initiated.
    If `False` (default), then a TreeSimTotalExtinctionException is raised.

    A Random() object or equivalent can be passed using the `rng` keyword;
    otherwise GLOBAL_RNG is used.
    """
    if 'ntax' not in kwargs \
        and 'taxon_set' not in kwargs \
        and 'max_time' not in kwargs:
            raise ValueError("At least one of the following must be specified: 'ntax', 'taxon_set', or 'max_time'")
    target_num_taxa = None
    taxon_set = None
    target_num_gens = kwargs.get('max_time', None)
    if 'taxon_set' in kwargs:
        taxon_set = kwargs.get('taxon_set')
        target_num_taxa = kwargs.get('ntax', len(taxon_set))
    elif 'ntax' in kwargs:
        target_num_taxa = kwargs['ntax']
    if taxon_set is None:
        taxon_set = dataobject.TaxonSet()
    repeat_until_success = kwargs.get('repeat_until_success', False)
    rng = kwargs.get('rng', GLOBAL_RNG)

    # grow tree
    if "tree" in kwargs:
        tree = kwargs['tree']
        if "taxon_set" in kwargs and kwargs['taxon_set'] is not tree.taxon_set:
            raise ValueError("Cannot specify both `tree` and `taxon_set`")
    else:
        tree = dataobject.Tree(taxon_set=taxon_set)
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
                    treemanip.prune_subtree(tree, nd)
                elif not repeat_until_success:
                    # all lineages are extinct: raise exception
                    raise TreeSimTotalExtinctionException()
                else:
                    # all lineages are extinct: repeat
                    num_gens = 0

        num_gens += 1
        leaf_nodes = tree.leaf_nodes()

    # If termination condition specified by ntax or taxon_set, then the last
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

def uniform_pure_birth(taxon_set, birth_rate=1.0, rng=None):
    "Generates a uniform-rate pure-birth process tree. "
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default
    tree = dataobject.Tree(taxon_set=taxon_set)
    tree.seed_node.edge.length = 0.0
    leaf_nodes = tree.leaf_nodes()
    while len(leaf_nodes) < len(taxon_set):
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
        leaf.taxon = taxon_set[idx]
    tree.is_rooted = True
    return tree

def pop_gen_tree(tree=None,
                 taxon_set=None,
                 ages=None,
                 num_genes=None,
                 pop_sizes=None,
                 num_genes_attr = 'num_genes',
                 pop_size_attr = 'pop_size',
                 rng=None):
    """
    This will simulate and return a tree with edges decorated with
    population sizes and leaf nodes decorated by the number of genes
    (samples or lineages) in each leaf.

    If `tree` is given, then this is used as the tree to be decorated.
    Otherwise, a Yule tree is generated based on the given taxon_set.
    Either `tree` or `taxon_set` must be given.

    The timing of the divergences can be controlled by specifying a vector of
    ages, `ages`. This should be sequences of values specifying the ages of the
    first, second, third etc. divergence events, in terms of time from the
    present, specified either in generations (if the `pop_sizes` vector is
    given) or population units (if the pop_size vector is not given).
    If an ages vector is given and there are less than num_pops-1 of these,
    then an exception is raised.

    The number of gene lineages per population can be specified through
    the 'num_genes', which can either be an scalar integer or a list.
    If it is an integer, all the population get the same number of
    genes. If it is a list, it must be at least as long as num_pops.

    The population sizes of each edge can be specified using the `pop_sizes`
    vector, which should be a sequence of values specifying the population
    sizes of the edges in postorder. If the pop_size vector is given, then it
    must be at least as long as there are branches on a tree, i.e. 2 * num_pops
    + 1, otherwise it is an error.  The population size should be the effective
    *haploid* population size; i.e., number of gene copies in the population: 2
    * N in a diploid population of N individuals, or N in a haploid population
    * of N individuals.

    If `pop_size` is 1 or 0 or None, then edge lengths of the tree are in
    haploid population units; i.e. where 1 unit of time equals 2N generations
    for a diploid population of size N, or N generations for a haploid
    population of size N. Otherwise edge lengths of the tree are in
    generations.

    This function first generates a tree using a pure-birth model with a
    uniform birth rate of 1.0. If an ages vector is given, it then sweeps
    through the internal nodes, assigning branch lengths such that the
    divergence events correspond to the ages in the vector. If a population
    sizes vector is given, it then visits all the edges in postorder, assigning
    population sizes to the attribute with the name specified in
    'pop_size_attr' (which is persisted as an annotation). During this, if an
    ages vector was *not* given, then the edge lengths are multiplied by the
    population size of the edge so the branch length units will be in
    generations. If an ages vector was given, then it is assumed that the ages
    are already in the proper scale/units.
    """

    # get our random number generator
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default

    # get a yule tree
    if not tree:
        if taxon_set:
            tree = uniform_pure_birth(taxon_set=taxon_set,
                                      rng=rng)
        else:
            raise Exception("Either tree or taxa block must be given")

    num_pops = len(tree.leaf_nodes())

    # basic idiot-checking
    if ages is not None and len(ages) < (num_pops - 1):
        msg = "Too few ages specified."
        raise Exception(msg)
    if num_genes is not None:
        if isinstance(num_genes, list):
            if len(num_genes) < num_pops:
                msg = "Too few number of gene samples specified"
                raise Exception(msg)
            else:
                samples = num_genes
        else:
            samples = [num_genes for tax in range(num_pops)]
    else:
        samples = None
    if pop_sizes is not None and len(pop_sizes) < (2 * num_pops + 1):
        msg = "Too few population sizes specified."
        raise Exception(msg)

    # set the ages
    if ages is not None:

        # get the internal nodes on the tree in reverse branching
        # order, so that newest nodes are returned first
        nodes = tree.nodes(cmp_fn = lambda x, y : \
                           int((y.distance_from_root()-x.distance_from_root())*10e+6),
                           filter_fn = lambda x : not x.is_leaf())
        # assign the ages
        for index, node in enumerate(nodes):
            for child in node.child_nodes():
                child.edge.length = ages[index] - child.distance_from_tip()

    # set the gene samples
    if samples is not None:
        for index, leaf in enumerate(tree.leaf_iter()):
            setattr(leaf, num_genes_attr, samples[index])
            leaf.annotate(num_genes_attr)

    # set the population sizes
    if pop_sizes is not None:
        index = 0
        for edge in tree.postorder_edge_iter():
            setattr(edge, pop_size_attr, pop_sizes[index])
            edge.annotate(pop_size_attr)
            if ages is None:
                edge.length = edge.length * getattr(edge, pop_size_attr)
            index = index + 1

    return tree

def contained_coalescent(containing_tree,
        gene_to_containing_taxon_map,
        edge_pop_size_attr="pop_size",
        default_pop_size=1,
        rng=None):
    """
    Returns a gene tree simulated under the coalescent contained within a
    population or species tree.

        `containing_tree`
            The population or species tree. If `edge_pop_size_map` is not None,
            and population sizes given are non-trivial (i.e., >1), then edge
            lengths on this tree are in units of generations. Otherwise edge
            lengths are in population units; i.e. 2N generations for diploid
            populations of size N, or N generations for diploid populations of
            size N.

        `gene_to_containing_taxon_map`
            A TaxonSetMapping object mapping Taxon objects in the
            `containing_tree` TaxonSet to corresponding Taxon objects in the
            resulting gene tree.

        `edge_pop_size_attr`
            Name of attribute of edges that specify population size. By default
            this is "pop_size". If this attribute does not exist,
            `default_pop_size` will be used.  The value for this attribute
            should be the haploid population size or the number of genes;
            i.e.  2N for a diploid population of N individuals, or N for a
            haploid population of N individuals. This value determines how
            branch length units are interpreted in the input tree,
            `containing_tree`.  If a biologically-meaningful value, then branch
            lengths on the `containing_tree` are properly read as generations.
            If not (e.g. 1 or 0), then they are in population units, i.e. where
            1 unit of time equals 2N generations for a diploid population of
            size N, or N generations for a haploid population of size N.
            Otherwise time is in generations. If this argument is None, then
            population sizes default to `default_pop_size`.

        `default_pop_size`
            Population size to use if `edge_pop_size_attr` is None or
            if an edge does not have the attribute. Defaults to 1.

    The returned gene tree will have the following extra attributes:

        `pop_node_genes`
            A dictionary with nodes of `containing_tree` as keys and a list of gene
            tree nodes that are uncoalesced as values.

    Note that this function does very much the same thing as
    `constrained_kingman()`, but provides a very different API.
    """

    if rng is None:
        rng = GLOBAL_RNG

    gene_tree_taxon_set = gene_to_containing_taxon_map.domain_taxon_set
    if gene_tree_taxon_set is None:
        gene_tree_taxon_set = dendropy.TaxonSet()
        for gene_taxa in pop_gene_taxa_map:
            for taxon in gene_taxa:
                gene_tree_taxon_set.add(taxon)
    gene_tree = dataobject.Tree(taxon_set=gene_tree_taxon_set)
    gene_tree.is_rooted = True

    pop_node_genes = {}
    pop_gene_taxa = gene_to_containing_taxon_map.reverse
    for nd in containing_tree.postorder_node_iter():
        if nd.taxon and nd.taxon in pop_gene_taxa:
            pop_node_genes[nd] = []
            gene_taxa = pop_gene_taxa[nd.taxon]
            for gene_taxon in gene_taxa:
                gene_node = dataobject.Node()
                gene_node.taxon = gene_taxon
                pop_node_genes[nd].append(gene_node)
            #gene_nodes = [dataobject.Node() for i in range(len(gene_taxa))]
            #for gidx, gene_node in enumerate(gene_nodes):
            #    gene_node.taxon = gene_taxa[gidx]
            #    pop_node_genes[nd].append(gene_node)

    for edge in containing_tree.postorder_edge_iter():

        if edge_pop_size_attr and hasattr(edge, edge_pop_size_attr):
            pop_size = getattr(edge, edge_pop_size_attr)
        else:
            pop_size = default_pop_size
        if edge.head_node.parent_node is None:
            if len(pop_node_genes[edge.head_node]) > 1:
                final = coalescent.coalesce(nodes=pop_node_genes[edge.head_node],
                                            pop_size=default_pop_size,
                                            period=None,
                                            rng=rng)
            else:
                final = pop_node_genes[edge.head_node]
            gene_tree.seed_node = final[0]
        else:
            uncoal = coalescent.coalesce(nodes=pop_node_genes[edge.head_node],
                                         pop_size=pop_size,
                                         period=edge.length,
                                         rng=rng)
            if edge.tail_node not in pop_node_genes:
                pop_node_genes[edge.tail_node] = []
            pop_node_genes[edge.tail_node].extend(uncoal)

    gene_tree.pop_node_genes = pop_node_genes
    return gene_tree

def pure_kingman(taxon_set, pop_size=1, rng=None):
    """
    Generates a tree under the unconstrained Kingman's coalescent process.
    """

    # get our random number generator
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default

    nodes = [dataobject.Node(taxon=t) for t in taxon_set]
    seed_node = coalescent.coalesce(nodes=nodes,
                                    pop_size=pop_size,
                                    period=None,
                                    rng=rng,
                                    use_expected_tmrca=True)[0]
    tree = dataobject.Tree(taxon_set=taxon_set, seed_node=seed_node)
    return tree

def mean_kingman(taxon_set, pop_size=1):
    """
    Returns a tree with coalescent intervals given by the expected times under
    Kingman's neutral coalescent.
    """

    # get our random number generator
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default

    nodes = [dataobject.Node(taxon=t) for t in taxon_set]
    seed_node = coalescent.coalesce(nodes=nodes,
                                    pop_size=pop_size,
                                    period=None,
                                    rng=rng,
                                    use_expected_tmrca=True)[0]
    tree = dataobject.Tree(taxon_set=taxon_set, seed_node=seed_node)
    return tree

def constrained_kingman(pop_tree,
                        gene_tree_list=None,
                        rng=None,
                        gene_node_label_func=None,
                        num_genes_attr='num_genes',
                        pop_size_attr='pop_size',
                        decorate_original_tree=False):
    """
    Given a population tree, `pop_tree` this will return a *pair of
    trees*: a gene tree simulated on this population tree based on
    Kingman's n-coalescent, and population tree with the additional
    attribute 'gene_nodes' on each node, which is a list of
    uncoalesced nodes from the gene tree associated with the given
    node from the population tree.

    `pop_tree` should be a DendroPy Tree object or an object
    of a class derived from this with the following attribute
    `num_genes` -- the number of gene samples from each population in the
    present.  Each edge on the tree should also have the attribute

    `pop_size_attr` is the attribute name of the edges of `pop_tree` that
    specify the population size. By default it is `pop_size`. The should
    specify the effective *haploid* population size; i.e., number of gene
    in the population: 2 * N in a diploid population of N individuals,
    or N in a haploid population of N individuals.

    If `pop_size` is 1 or 0 or None, then the edge lengths of `pop_tree` is
    taken to be in haploid population units; i.e. where 1 unit equals 2N
    generations for a diploid population of size N, or N generations for a
    haploid population of size N. Otherwise the edge lengths of `pop_tree` is
    taken to be in generations.

    If `gene_tree_list` is given, then the gene tree is added to the
    tree block, and the tree block's taxa block will be used to manage
    the gene tree's `taxa`.

    `gene_node_label_func` is a function that takes two arguments (a string
    and an integer, respectively, where the string is the containing species
    taxon label and the integer is the gene index) and returns a label for
    the corresponding the gene node.

    if `decorate_original_tree` is True, then the list of uncoalesced nodes at
    each node of the population tree is added to the original (input) population
    tree instead of a copy.

    Note that this function does very much the same thing as `contained_coalescent()`,
    but provides a very different API.
    """

    # get our random number generator
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default

    if gene_tree_list is not None:
        gtaxa = gene_tree_list.taxon_set
    else:
        gtaxa = dataobject.TaxonSet()

    if gene_node_label_func is None:
        gene_node_label_func = lambda x, y: "%s_%02d" % (x, y)

    # we create a set of gene nodes for each leaf node on the population
    # tree, and associate those gene nodes to the leaf by assignment
    # of 'taxon'.
    for leaf_count, leaf in enumerate(pop_tree.leaf_iter()):
        gene_nodes = []
        for gene_count in range(getattr(leaf, num_genes_attr)):
            gene_node = dataobject.Node()
            gene_node.taxon = gtaxa.require_taxon(label=gene_node_label_func(leaf.taxon.label, gene_count+1))
            gene_nodes.append(gene_node)
        leaf.gene_nodes = gene_nodes

    # We iterate through the edges of the population tree in post-order,
    # i.e., visiting child edges before we visit parent edges. For
    # each edge visited, we take the genes found in the child nodes,
    # and run the coalescent simulation on them attacheded by the length
    # of the edge. Any genes that have not yet coalesced at the end of
    # this period are added to the genes of the tail (parent) node of
    # the edge.

    if decorate_original_tree:
        working_poptree = pop_tree
    else:
        # start with a new (deep) copy of the population tree so as to not
        # to change the original tree
        working_poptree = copy.deepcopy(pop_tree)

    # start with a new tree
    gene_tree = dataobject.Tree()
    gene_tree.taxon_set = gtaxa
    for edge in working_poptree.postorder_edge_iter():

        # if mrca root, run unconstrained coalescent
        if edge.head_node.parent_node is None:
            if len(edge.head_node.gene_nodes) > 1:
                final = coalescent.coalesce(nodes=edge.head_node.gene_nodes,
                                            pop_size=pop_size,
                                            period=None,
                                            rng=rng)
            else:
                final = edge.head_node.gene_nodes
            gene_tree.seed_node = final[0]
        else:

            if hasattr(edge, pop_size_attr):
                pop_size = getattr(edge, pop_size_attr)
            else:
                # this means all our time will be in population units
                pop_size = 1

            uncoal = coalescent.coalesce(nodes=edge.head_node.gene_nodes,
                                         pop_size=pop_size,
                                         period=edge.length,
                                         rng=rng)
            if not hasattr(edge.tail_node, 'gene_nodes'):
                edge.tail_node.gene_nodes = []
            edge.tail_node.gene_nodes.extend(uncoal)

    gene_tree.is_rooted = True
    if gene_tree_list is not None:
        gene_tree_list.append(gene_tree)
        return gene_tree, working_poptree
    else:
        return gene_tree, working_poptree

