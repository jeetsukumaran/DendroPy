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
Functions and classes for manipulating or changing tree
structural relations.
"""

import sys
import copy
import math

from dendropy.utility import GLOBAL_RNG
from dendropy import treesplit
from dendropy import dataobject

def collapse_edge(edge):
    """
    Inserts all children of the head_node of edge as children of the
    tail_node of edge in the same place in the child_node list that head_node
    had occupied. The edge length and head_node will no longer be
    part of the tree.
    """
    edge.collapse()

def collapse_clade(node):
    """Collapses all internal edges that are descendants of node."""
    node.collapse_clade()

def prune_subtree(tree, node, delete_outdegree_one=True):
    tree.prune_subtree(node=node,
            delete_outdegree_one=delete_outdegree_one)
    return tree

def prune_leaves_without_taxa(tree, delete_outdegree_one=True):
    """
    Removes all terminal nodes that have their ``taxon`` attribute set to
    ``None``.
    """
    tree.prune_leaves_without_taxa(delete_outdegree_one=delete_outdegree_one)
    return tree

def prune_taxa(tree, taxa, delete_outdegree_one=True):
    """
    Removes terminal nodes associated with Taxon objects given by the container
    `taxa` (which can be any iterable, including a TaxonSet object) from `tree`.
    """
    tree.prune_taxa(taxa=taxa,
            delete_outdegree_one=delete_outdegree_one)
    return tree

def retain_taxa(tree, taxa, delete_outdegree_one=True):
    """
    Removes terminal nodes that are not associated with any
    of the Taxon objects given by ``taxa`` (which can be any iterable, including a
    TaxonSet object) from the ``tree``.
    """
    tree.retain_taxa(taxa=taxa,
            delete_outdegree_one=delete_outdegree_one)

def randomly_reorient_tree(tree, rng=None, splits=False):
    """
    Randomly picks a new rooting position and rotates the branches around all
    internal nodes in the `tree`. If `splits` is True, the the `split_bitmask`
    and `split_edges` attributes kept valid.
    """
    tree.randomly_reorient_tree(rng=rng, update_splits=splits)

def randomly_rotate(tree, rng=None):
    "Randomly rotates the branches around all internal nodes in the `tree`"
    tree.randomly_rotate(rng=rng)

def collapse_conflicting(subtree_root, split, split_bitmask):
    """
    Takes a node that is the root of a subtree.  Collapses every edge in the
    subtree that conflicts with split.  This can include the edge subtending
    subtree_root.
    """
    # we flip splits so that both the split and each edges split  have the
    # lowest bit of the clade mask set to one
    lb = treesplit.lowest_bit_only(split_bitmask)

    if lb & split:
        cropped_split = split & split_bitmask
    else:
        cropped_split = (~split) & split_bitmask

    to_collapse_head_nodes = []
    for nd in subtree_root.postorder_iter():
        if not nd.is_leaf():
            ncm = nd.edge.split_bitmask
            if lb & ncm:
                nd_split = ncm & split_bitmask
            else:
                nd_split = (~ncm) & split_bitmask

            cm_union = nd_split | cropped_split
            if (cm_union != nd_split) and (cm_union != cropped_split):
                to_collapse_head_nodes.append(nd)

    for nd in to_collapse_head_nodes:
        e = nd.edge
        e.collapse()

def scale_edges(tree, edge_len_multiplier):
    tree.scale_edges(edge_len_multiplier=edge_len_multiplier)

