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

"""
Statistics, metrics, measurements, and values calculated on (single) trees.
"""

import math
from dendropy.calculate import phylogeneticdistance
from dendropy.utility import deprecate

EULERS_CONSTANT = 0.5772156649015328606065120900824024310421


class PatristicDistanceMatrix(phylogeneticdistance.PhylogeneticDistanceMatrix):

    def __init__(self, tree):
        deprecate.dendropy_deprecation_warning(
            message="PatristicDistanceMatrix is deprecated since Dendropy 5. "
            "Use entropy. phylogenetic distance.PhylogeneticDistanceMatrix instead.",
        )
        phylogeneticdistance.PhylogeneticDistanceMatrix.__init__(self)
        self.compile_from_tree(tree=tree)

def patristic_distance(tree, taxon1, taxon2, is_bipartitions_updated=False):
    """
    Given a tree with bipartitions encoded, and two taxa on that tree, returns the
    patristic distance between the two. Much more inefficient than constructing
    a PhylogeneticDistanceMatrix object.
    """
    mrca = tree.mrca(taxa=[taxon1, taxon2], is_bipartitions_updated=is_bipartitions_updated)
    dist = 0
    n = tree.find_node(lambda x: x.taxon == taxon1)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    n = tree.find_node(lambda x: x.taxon == taxon2)
    while n != mrca:
        if n.edge.length is not None:
            dist += n.edge.length
        n = n.parent_node
    return dist

###########################################################################
### Metrics -- Unary

def B1(tree):
    """
    Returns the B1 statistic: the reciprocal of the sum of the maximum
    number of nodes between each interior node and tip over all internal
    nodes excluding root.
    """
    b1 = 0.0
    nd_mi = {}
    for nd in tree.postorder_node_iter():
        if nd._parent_node is None:
            continue
        child_nodes = nd._child_nodes
        if len(child_nodes) == 0:
            nd_mi[nd] = 0.0
            continue
        mi = max(nd_mi[ch] for ch in child_nodes)
        mi += 1
        nd_mi[nd] = mi
        b1 += 1.0/mi
    return b1

def colless_tree_imbalance(tree, normalize="max"):
    """
    Returns Colless' tree imbalance or I statistic: the sum of differences
    of numbers of children in left and right subtrees over all internal
    nodes. ``normalize`` specifies the normalization:

        * "max" or True [DEFAULT]
            normalized to maximum value for tree of
            this size
        * "yule"
            normalized to the Yule model
        * "pda"
            normalized to the PDA (Proportional to Distinguishable
            Arrangements) model
        * None or False
            no normalization

    """
    colless = 0.0
    num_leaves = 0
    subtree_leaves = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            subtree_leaves[nd] = 1
            num_leaves += 1
        else:
            total_leaves = 0
            if len(nd._child_nodes) > 2:
                raise TypeError("Colless' tree imbalance statistic requires strictly bifurcating trees")
            left = subtree_leaves[nd._child_nodes[0]]
            right = subtree_leaves[nd._child_nodes[1]]
            colless += abs(right-left)
            subtree_leaves[nd] = right + left
    if normalize == "yule":
        colless = float(colless - (num_leaves * math.log(num_leaves)) - (num_leaves * (EULERS_CONSTANT - 1.0 - math.log(2))))/num_leaves
    elif normalize == "pda":
        colless = colless / pow(num_leaves, 3.0/2)
    elif normalize is True or normalize == "max":
        ## note that Mooers 1995 (Evolution 49(2):379-384)
        ## remarks that the correct normalization factor is
        ## 2/((num_leaves - 1) * (num_leaves -2))
        colless = colless * (2.0/(num_leaves * (num_leaves-3) + 2))
    elif normalize is not None and normalize is not False:
        raise TypeError("``normalization`` accepts only None, True, False, 'yule' or 'pda' as argument values")
    return colless

def pybus_harvey_gamma(tree, prec=0.00001):
    """Returns the gamma statistic of Pybus and Harvey (2000). This statistic
    is used to test for constancy of birth and death rates over the course of
    a phylogeny.  Under the pure-birth process, the statistic should follow
    a standard Normal distibution: a Normal(mean=0, variance=1).

    If the lengths of different paths to the node differ by more than ``prec``,
        then a ValueError exception will be raised indicating deviation from
        ultrametricty.
    Raises a Value Error if the tree is not ultrametric, is non-binary, or has
        only 2 leaves.

    As a side effect a ``age`` attribute is added to the nodes of the tree.

    Pybus and Harvey. 2000. "Testing macro-evolutionary models using incomplete
    molecular phylogenies." Proc. Royal Society Series B: Biological Sciences.
    (267). 2267-2272
    """
    # the equation is given by:
    #   T = \sum_{j=2}^n (jg_j)
    #   C = T \sqrt{\frac{1}{12(n-2)}}
    #   C gamma = \frac{1}{n-2}\sum_{i=2}^{n-1} (\sum_{k=2}^i kg_k) - \frac{T}{2}
    # where n is the number of taxa, and g_2 ... g_n is the vector of waiting
    #   times between consecutive (in time, not along a branch) speciation times.
    node = None
    speciation_ages = []
    n = 0
    if tree.seed_node.age is None:
        tree.calc_node_ages(ultrametricity_precision=prec)
    for node in tree.postorder_node_iter():
        if len(node.child_nodes()) == 2:
            speciation_ages.append(node.age)
        else:
            n += 1
    if node is None:
        raise ValueError("Empty tree encountered")
    speciation_ages.sort(reverse=True)
    g = []
    older = speciation_ages[0]
    for age in speciation_ages[1:]:
        g.append(older - age)
        older = age
    g.append(older)
    if not g:
        raise ValueError("No internal nodes found (other than the root)")
    assert(len(g) == (n - 1))
    T = 0.0
    accum = 0.0
    for i in range(2, n):
        list_index = i - 2
        T += i * float(g[list_index])
        accum += T
    list_index = n - 2
    T += (n) * g[list_index]
    nmt = n - 2.0
    numerator = accum/nmt - T/2.0
    C = T*pow(1/(12*nmt), 0.5)
    return numerator/C

def N_bar(tree):
    """
    Returns the $\bar{N}$ statistic: the average number of nodes above a
    terminal node.
    """
    leaf_count = 0
    nbar = 0
    for leaf_node in tree.leaf_node_iter():
        leaf_count += 1
        for parent in leaf_node.ancestor_iter(inclusive=False):
            nbar += 1
    return float(nbar) / leaf_count

def sackin_index(tree, normalize=True):
    """
    Returns the Sackin's index: the sum of the number of ancestors for each
    tip of the tree. The larger the Sackin's index, the less balanced the
    tree. ``normalize`` specifies the normalization:

        * True [DEFAULT]
            normalized to number of leaves; this results in a value
            equivalent to that given by Tree.N_bar()
        * "yule"
            normalized to the Yule model
        * "pda"
            normalized to the PDA (Proportional to Distinguishable
            Arrangements) model
        * None or False
            no normalization

    """
    leaf_count = 0
    num_anc = 0
    for leaf_node in tree.leaf_node_iter():
        leaf_count += 1
        for parent in leaf_node.ancestor_iter(inclusive=False):
            num_anc += 1
    if normalize == "yule":
        x = sum(1.0/j for j in range(2, leaf_count+1))
        s = float(num_anc - (2 * leaf_count * x))/leaf_count
    elif normalize == "pda":
        s = float(num_anc)/(pow(leaf_count, 3.0/2))
    elif normalize is True:
        s = float(num_anc)/leaf_count
    elif normalize is None or normalize is False:
        s = float(num_anc)
    elif normalize is not None and normalize is not False:
        raise TypeError("``normalization`` accepts only None, True, False, 'yule' or 'pda' as argument values")
    return s

def treeness(tree):
    """
    Returns the proportion of total tree length that is taken up by
    internal branches.
    """
    internal = 0.0
    external = 0.0
    for nd in tree.postorder_node_iter():
        if not nd._parent_node:
            continue
        if nd.is_leaf():
            external += nd.edge.length
        else:
            internal += nd.edge.length
    return internal/(external + internal)

def node_ages(tree, is_internal_only=False):
    """
    Returns vector of branching events indexed in backward time.
    """
    return sorted(age for (node, age) in tree.resolve_node_ages().items() if (not is_internal_only) or node.is_internal())

def node_depths(tree, is_internal_only=False):
    """
    Returns vector of branching events indexed in forward time.
    """
    return sorted(age for (node, age) in tree.resolve_node_depths().items() if (not is_internal_only) or node.is_internal())

def coalescence_ages(tree):
    """
    Returns vector of splits indexed in backward time.
    """
    return node_ages(tree, is_internal_only=True)

def divergence_times(tree):
    """
    Returns vector of splits indexed in forward time.
    """
    return divergence_times(tree, is_internal_only=True)
