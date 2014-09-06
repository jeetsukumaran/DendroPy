#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
All of these functions are now native methods of the Tree class. They are
maintained here only for background compatibility.
"""

from dendropy.utility import error

def collapse_edge(edge):
    error.dendropy_construct_migration_warning("treemanip.collapse_edge(edge)", "edge.collapse()")
    edge.collapse()

def collapse_clade(node):
    error.dendropy_construct_migration_warning("treemanip.collapse_clade(node)", "node.collapse_clade()")
    node.collapse_clade()

def prune_subtree(tree, node, delete_outdegree_one=True):
    error.dendropy_construct_migration_warning("treemanip.prune_subtree(tree, node)", "tree.prune_subtree(node)")
    tree.prune_subtree(node=node,
            delete_outdegree_one=delete_outdegree_one)
    return tree

def prune_leaves_without_taxa(tree, delete_outdegree_one=True):
    error.dendropy_construct_migration_warning("treemanip.prune_leaves_without_taxa(tree)", "tree.prune_leaves_without_taxa()")
    tree.prune_leaves_without_taxa(delete_outdegree_one=delete_outdegree_one)
    return tree

def prune_taxa(tree, taxa, delete_outdegree_one=True):
    error.dendropy_construct_migration_warning("treemanip.prune_taxa(tree, taxa)", "tree.prune_taxa(taxa)")
    tree.prune_taxa(taxa=taxa,
            delete_outdegree_one=delete_outdegree_one)
    return tree

def retain_taxa(tree, taxa, delete_outdegree_one=True):
    error.dendropy_construct_migration_warning("treemanip.retain_taxa(tree, taxa)", "tree.retain_taxa(taxa)")
    tree.retain_taxa(taxa=taxa,
            delete_outdegree_one=delete_outdegree_one)

def randomly_reorient_tree(tree, rng=None, splits=False):
    error.dendropy_construct_migration_warning("randomly_reorient_tree(tree)", "tree.randomly_reorient()")
    tree.randomly_reorient(rng=rng, update_splits=splits)

def randomly_rotate(tree, rng=None):
    error.dendropy_construct_migration_warning("treemanip.randomly_rotate(tree)", "tree.randomly_rotate()")
    tree.randomly_rotate(rng=rng)

def collapse_conflicting(subtree_root, split, split_bitmask):
    error.dendropy_construct_migration_warning("treemanip.collapse_conflicting(node)", "node.collapse_conflicting()")
    subtree_root.collapse_conflicting(split, split_bitmask)

def scale_edges(tree, edge_len_multiplier):
    error.dendropy_construct_migration_warning("treemanip.scale_edges(tree)", "tree.scale_edges()")
    tree.scale_edges(edge_len_multiplier=edge_len_multiplier)

