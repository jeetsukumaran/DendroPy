#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

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

def star_tree(taxon_set):
    "Builds and returns a star tree from the given taxa block."
    star_tree = dataobject.Tree(taxon_set=taxon_set)
    for taxon in taxon_set:
        star_tree.seed_node.new_child(taxon=taxon)
    return star_tree

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
        - If 'max_gens' is given as a keyword argument, tree is grown for `max_gens`
          number of generations.

    If more than one of the above is given, then tree growth will terminate when
    *any* of the termination conditions (i.e., number of tips == `ntax`, or number
    of tips == len(taxon_set) or number of generations = `max_gens`) are met.

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
    """
    if 'ntax' not in kwargs \
        and 'taxon_set' not in kwargs \
        and 'max_gens' not in kwargs:
            raise ValueError("At least one of the following must be specified: 'ntax', 'taxon_set', or 'max_gens'")
    target_num_taxa = None
    taxon_set = None
    target_num_gens = kwargs.get('max_gens', None)
    if 'taxon_set' in kwargs:
        taxon_set = kwargs.get('taxon_set')
        target_num_taxa = kwargs.get('ntax', len(taxon_set))
    elif 'ntax' in kwargs:
        target_num_taxa = kwargs['ntax']
    if taxon_set is None:
        taxon_set = dataobject.TaxonSet()
    rng = kwargs.get('rng', GLOBAL_RNG)

    # grow tree
    if "tree" in kwargs:
        tree = kwargs['tree']
        if "taxon_set" in kwargs and kwargs['taxon_set'] is not tree.taxon_set:
            raise ValueError("Cannot specify both `tree` and `taxon_set`")
    else:
        tree = dataobject.Tree(taxon_set=taxon_set)
        tree.seed_node.edge.length = 0
    leaf_nodes = tree.leaf_nodes()
    num_gens = 0
    while (target_num_taxa is None or len(leaf_nodes) < target_num_taxa) \
            and (target_num_gens is None or num_gens < target_num_gens):
        for nd in leaf_nodes:
            if not hasattr(nd, 'birth_rate'):
                nd.birth_rate = birth_rate
            if not hasattr(nd, 'death_rate'):
                nd.death_rate = death_rate
            nd.edge.length += 1
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
                else:
                    return tree
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
    """
    if 'ntax' not in kwargs \
        and 'taxon_set' not in kwargs \
        and 'max_time' not in kwargs:
            raise ValueError("At least one of the following must be specified: 'ntax', 'taxon_set', or 'max_time'")
    target_num_taxa = None
    taxon_set = None
    max_time = kwargs.get('max_time', None)
    if 'taxon_set' in kwargs:
        taxon_set = kwargs.get('taxon_set')
        target_num_taxa = kwargs.get('ntax', len(taxon_set))
    elif 'ntax' in kwargs:
        target_num_taxa = kwargs['ntax']
    if taxon_set is None:
        taxon_set = dataobject.TaxonSet()
    rng = kwargs.get('rng', GLOBAL_RNG)

    # initialize tree
    if "tree" in kwargs:
        tree = kwargs['tree']
        if "taxon_set" in kwargs and kwargs['taxon_set'] is not tree.taxon_set:
            raise ValueError("Cannot specify both `tree` and `taxon_set`")
    else:
        tree = dataobject.Tree(taxon_set=taxon_set)
        tree.seed_node.edge.length = 0
        tree.seed_node.birth_rate = birth_rate
        tree.seed_node.death_rate = death_rate

    # grow tree
    leaf_nodes = tree.leaf_nodes()
    total_time = 0
    while (target_num_taxa is None or len(leaf_nodes) < target_num_taxa) \
            and (max_time is None or total_time < max_time):

        # get vector of birth/death probabilities, and
        # associate with nodes/events
        event_probs = []
        event_nodes = []
        for nd in tree.leaf_nodes():
            if not hasattr(nd, 'birth_rate'):
                nd.birth_rate = birth_rate
            if not hasattr(nd, 'death_rate'):
                nd.death_rate = death_rate
            event_probs.append(nd.birth_rate)
            event_nodes.append((nd, True)) # birth event = True
            event_probs.append(nd.death_rate)
            event_nodes.append((nd, False)) # birth event = False; i.e. death

            # get total probability of any birth/death
            prob_of_event = sum(event_probs)

            # waiting time based on above probability
            waiting_time = rng.expovariate(1.0/prob_of_event)

        # add waiting time to nodes
        for nd in tree.leaf_nodes():
            nd.edge.length += waiting_time
        total_time += waiting_time

        # if event occurs within time constraints
        if max_time is None or total_time <= max_time:

            # normalize probability
            for i in xrange(len(event_probs)):
                event_probs[i] = event_probs[i]/prob_of_event

            # select node/event and process
            nd, birth_event = probability.lengthed_choice(event_nodes, event_probs)
            if birth_event:
                c1 = nd.new_child()
                c2 = nd.new_child()
                c1.edge.length = 0
                c2.edge.length = 0
                c1.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c1.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
                c2.birth_rate = nd.birth_rate + rng.gauss(0, birth_rate_sd)
                c2.death_rate = nd.death_rate + rng.gauss(0, death_rate_sd)
            else:
                if nd is not tree.seed_node:
                    treemanip.prune_subtree(tree, nd)
                else:
                    return tree

        leaf_nodes = tree.leaf_nodes()

    # If termination condition specified by ntax or taxon_set, then the last
    # split will have a daughter edges of length == 0;
    # so we continue growing the edges until the next birth/death event *or*
    # the max time is reached
    if max_time is None or total_time < max_time:
        for nd in tree.leaf_nodes():
            if not hasattr(nd, 'birth_rate'):
                nd.birth_rate = birth_rate
            if not hasattr(nd, 'death_rate'):
                nd.death_rate = death_rate
            event_probs.append(nd.birth_rate)
            event_probs.append(nd.death_rate)
            waiting_time = rng.expovariate(1/sum(event_probs))
            if max_time is None or (waiting_time + total_time) < max_time:
                remaining_time = waiting_time
            else:
                remaining_time = total_time - max_time

        for nd in tree.leaf_nodes():
            nd.edge.length += remaining_time

    if kwargs.get("assign_taxa", True):
        tree.randomly_assign_taxa(create_required_taxa=True, rng=rng)

    # return
    return tree

def uniform_pure_birth(taxon_set,
                       birth_rate=1.0,
                       ultrametricize=True,
                       rng=None):
    "Generates a uniform-rate pure-birth process tree. "
    if rng is None:
        rng = GLOBAL_RNG # use the global rng by default
    tree = dataobject.Tree(taxon_set=taxon_set)

    leaf_nodes = tree.leaf_nodes()
    count = 0
    while len(leaf_nodes) < len(taxon_set):
        parent_node = rng.choice(leaf_nodes)
        edge_length = rng.expovariate(len(leaf_nodes)/birth_rate)
        child1 = dataobject.Node()
        child2 = dataobject.Node()
        child1.node_id = 'n' + str(count+1)
        child2.node_id = 'n' + str(count+2)
        child1.edge.length = edge_length
        child2.edge.length = edge_length
        parent_node.add_child(child1)
        parent_node.add_child(child2)
        count = count + 2
        leaf_nodes = tree.leaf_nodes()
    leaf_nodes = tree.leaf_nodes()
    for idx, leaf in enumerate(leaf_nodes):
        leaf.taxon = taxon_set[idx]
    if ultrametricize:
        max_distance_from_root = max([node.distance_from_root() for node in leaf_nodes])
        for node in leaf_nodes:
            node.edge.length = node.edge.length + (max_distance_from_root - node.distance_from_root())
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

    The timing of the divergences can be
    controlled by specifying a vector of ages, `ages`. This should be
    sequences of values specifying the ages of the first, second,
    third etc. divergence events, in terms of time from the present,
    specified either in generations (if the pop_sizes vector is given)
    or population units (if the pop_size vector is not given). If an
    ages vector is given and there are less than num_pops-1 of these, then
    an exception is raised.

    The number of gene lineages per population can be specified through
    the 'num_genes', which can either be an scalar integer or a list.
    If it is an integer, all the population get the same number of
    genes. If it is a list, it must be at least as long as num_pops.

    The population sizes of each edge can be specified using the
    `pop_sizes` vector, which should be a sequence of values
    specifying the population sizes of the edges in postorder. If the
    pop_size vector is given, then it must be at least as long as
    there are branches on a tree (=2 * num_pops + 1), otherwise it is an
    error. If it is not given, then the branch lengths of the population
    trees will be in population units.

    This function first generates a tree using a pure-birth model with
    a uniform birth rate of 1.0. If an ages vector is given, it then
    sweeps through the internal nodes, assigning branch lengths such
    that the divergence events correspond to the ages in the
    vector. If a population sizes vector is given, it then visits all
    the edges in postorder, assigning population sizes to the
    attribute 'pop_size' (which is persisted as an annotation). During
    this, if an ages vector was *not* given, then the edge lengths are
    multiplied by the population size of the edge so the branch length
    units will be in generations. If an ages vector was given, then it
    is assumed that the ages are already in the proper scale/units.
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
                                    rng=rng)[0]
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
    `pop_size` -- the effective size of the population at this time.

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
        edge.head_node.gene_nodes = edge.head_node.gene_nodes

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

