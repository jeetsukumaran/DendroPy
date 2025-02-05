# !/usr/bin/env python

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
Generates lists of trees and dendropy.TreeList for tests.
"""

import dendropy

DEFAULT_NUM_TREES = 5
_TREE_COUNTER = 0

def get_tree(taxon_namespace=None,
        label=None,
        suppress_internal_node_taxa=False,
        suppress_leaf_node_taxa=False):
    global _TREE_COUNTER
    _TREE_COUNTER += 1
    if taxon_namespace is None:
        taxon_namespace = dendropy.TaxonNamespace()
    if label is None:
        label = "Tree{}".format(_TREE_COUNTER)
    t1 = dendropy.Tree(label=label,
            taxon_namespace=taxon_namespace)
    t1.seed_node.label = "i0"
    c1 = t1.seed_node.new_child(label="i1")
    c2 = t1.seed_node.new_child(label="i2")
    c1.new_child(label="t1")
    c1.new_child(label="t2")
    c2.new_child(label="t3")
    c2.new_child(label="t4")
    tax_labels = set()
    for nd in t1:
        is_leaf = nd.is_leaf()
        if is_leaf and not suppress_leaf_node_taxa:
            tax1 = t1.taxon_namespace.require_taxon(nd.label)
            nd.taxon = tax1
            tax_labels.add(nd.label)
        elif (not is_leaf) and not suppress_internal_node_taxa:
            tax1 = t1.taxon_namespace.require_taxon(nd.label)
            nd.taxon = tax1
            tax_labels.add(nd.label)
    t1.tax_labels = tax_labels
    try:
        t1.taxon_namespace.tax_labels.update(tax_labels)
    except AttributeError:
        t1.taxon_namespace.tax_labels = set(tax_labels)
    return t1

def get_trees(
        num_trees,
        taxon_namespace=None,
        label=None,
        suppress_internal_node_taxa=False,
        suppress_leaf_node_taxa=False):
    trees = []
    for idx in range(num_trees):
        t1 = get_tree(
                taxon_namespace=taxon_namespace,
                label=label,
                suppress_internal_node_taxa=suppress_internal_node_taxa,
                suppress_leaf_node_taxa=suppress_leaf_node_taxa)
        trees.append(t1)
    return trees

def get_tree_list(
        num_trees,
        taxon_namespace=None,
        label=None,
        suppress_internal_node_taxa=False,
        suppress_leaf_node_taxa=False):
    if taxon_namespace is None:
        taxon_namespace = dendropy.TaxonNamespace()
    tlist1 = dendropy.TreeList(label="1",
            taxon_namespace=taxon_namespace)
    for idx in range(num_trees):
        t1 = get_tree(
                taxon_namespace=taxon_namespace,
                label=label,
                suppress_internal_node_taxa=suppress_internal_node_taxa,
                suppress_leaf_node_taxa=suppress_leaf_node_taxa)
        assert t1.taxon_namespace is tlist1.taxon_namespace
        tlist1.append(t1)
    return tlist1

def get_tree_list_and_list_of_trees(
        num_trees,
        tree_list_taxon_namespace=None,
        list_of_trees_taxon_namespace=None):
    tlist = get_tree_list(
            num_trees=0,
            taxon_namespace=tree_list_taxon_namespace,
            label=None,
            suppress_internal_node_taxa=False,
            suppress_leaf_node_taxa=False)
    trees = get_trees(
            num_trees=num_trees,
            taxon_namespace=list_of_trees_taxon_namespace,
            label=None,
            suppress_internal_node_taxa=False,
            suppress_leaf_node_taxa=False)
    return tlist, trees

