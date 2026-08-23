#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
All of these functions are now native methods of the Tree class. They are
maintained here only for background compatibility.
"""

from dendropy.utility import deprecate

def _deprecate_tree_manip(old, new):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: Tree structure manipulation and editing functionality are now native methods of the 'dendropy.Tree' class.",
            old_construct=old,
            new_construct=new,
            stacklevel=4)

def collapse_edge(edge):
    _deprecate_tree_manip("treemanip.collapse_edge(edge)", "edge.collapse()")
    edge.collapse()

def collapse_clade(node):
    _deprecate_tree_manip("treemanip.collapse_clade(node)", "node.collapse_clade()")
    node.collapse_clade()

def prune_subtree(tree, node, suppress_unifurcations=True):
    _deprecate_tree_manip("treemanip.prune_subtree(tree, node)", "tree.prune_subtree(node)")
    tree.prune_subtree(node=node,
            suppress_unifurcations=suppress_unifurcations)
    return tree

def prune_leaves_without_taxa(tree, suppress_unifurcations=True):
    _deprecate_tree_manip("treemanip.prune_leaves_without_taxa(tree)", "tree.prune_leaves_without_taxa()")
    tree.prune_leaves_without_taxa(suppress_unifurcations=suppress_unifurcations)
    return tree

def prune_taxa(tree, taxa, suppress_unifurcations=True):
    _deprecate_tree_manip("treemanip.prune_taxa(tree, taxa)", "tree.prune_taxa(taxa)")
    tree.prune_taxa(taxa=taxa,
            suppress_unifurcations=suppress_unifurcations)
    return tree

def retain_taxa(tree, taxa, suppress_unifurcations=True):
    _deprecate_tree_manip("treemanip.retain_taxa(tree, taxa)", "tree.retain_taxa(taxa)")
    tree.retain_taxa(taxa=taxa,
            suppress_unifurcations=suppress_unifurcations)

def randomly_reorient_tree(tree, rng=None, splits=False):
    _deprecate_tree_manip("randomly_reorient_tree(tree)", "tree.randomly_reorient()")
    tree.randomly_reorient(rng=rng, update_splits=splits)

def randomly_rotate(tree, rng=None):
    _deprecate_tree_manip("treemanip.randomly_rotate(tree)", "tree.randomly_rotate()")
    tree.randomly_rotate(rng=rng)

def collapse_conflicting(subtree_root, split, split_bitmask):
    _deprecate_tree_manip("treemanip.collapse_conflicting(node)", "node.collapse_conflicting()")
    subtree_root.collapse_conflicting(split, split_bitmask)

def scale_edges(tree, edge_len_multiplier):
    _deprecate_tree_manip("treemanip.scale_edges(tree)", "tree.scale_edges()")
    tree.scale_edges(edge_len_multiplier=edge_len_multiplier)

