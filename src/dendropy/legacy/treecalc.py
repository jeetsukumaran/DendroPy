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
DEPRECATED IN DENDROPY 4: USE
`dendropy.calculate.treecompare`,`dendropy.calculate.treemeasure`,
or `dendropy.calculate.treescore`,
"""

import dendropy
from dendropy.calculate import treecompare
from dendropy.calculate import treemeasure
from dendropy.utility import deprecate

##############################################################################
## dendropy.calculate.treemeasure

class PhylogeneticDistanceMatrix(dendropy.PhylogeneticDistanceMatrix):
    def __init__(self, tree=None):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.PhylogeneticDistanceMatrix' class has moved to 'dendropy.PhylogeneticDistanceMatrix'.",
                old_construct="from dendropy import treecalc\npdm = treecalc.PhylogeneticDistanceMatrix(...)",
                new_construct="import dendropy\npdm = dendropy.PhylogeneticDistanceMatrix(...)")
        dendropy.PhylogeneticDistanceMatrix.__init__(self, tree=tree)

def patristic_distance(tree, taxon1, taxon2):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.patristic_distance()' function has moved to 'dendropy.calculate.treemeasure.patristic_distance()'.",
            old_construct="from dendropy import treecalc\npdm = treecalc.patristic_distance(...)",
            new_construct="from dendropy.calculate import treemeasure\npdm = treemeasure.patristic_distance(...)")
    return treemeasure.patristic_distance(
            tree=tree,
            taxon1=taxon1,
            taxon2=taxon2)

##############################################################################
## dendropy.calculate.treecompare

def robinson_foulds_distance(tree1,
        tree2,
        edge_length_attr="length"):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.robinson_foulds_distance()' function has moved to 'dendropy.calculate.treecompare.weighted_robinson_foulds_distance()'.",
            old_construct="from dendropy import treecalc\nd = treecalc.robinson_foulds_distance(...)",
            new_construct="from dendropy.calculate import treecompare\nd = treecompare.weighted_robinson_foulds_distance(...)")
    return treecompare.weighted_robinson_foulds_distance(
            tree1=tree1,
            tree2=tree2,
            edge_weight_attr=edge_length_attr)

def euclidean_distance(tree1,
        tree2,
        edge_length_attr="length",
        value_type=float):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.euclidean_distance()' function has moved to 'dendropy.calculate.treecompare.euclidean_distance()'.",
            old_construct="from dendropy import treecalc\nd = treecalc.euclidean_distance(...)",
            new_construct="from dendropy.calculate import treecompare\nd = treecompare.euclidean_distance(...)")
    return treecompare.euclidean_distance(
            tree1=tree1,
            tree2=tree2,
            edge_weight_attr=edge_length_attr,
            value_type=value_type)

def false_positives_and_negatives(reference_tree, test_tree):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.false_positives_and_negatives()' function has moved to 'dendropy.calculate.treecompare.false_positives_and_negatives()'.",
            old_construct="from dendropy import treecalc\nd = treecalc.false_positives_and_negatives(...)",
            new_construct="from dendropy.calculate import treecompare\nd = treecompare.false_positives_and_negatives(...)")
    return treecompare.false_positives_and_negatives(
            reference_tree=reference_tree,
            comparison_tree=test_tree)

def symmetric_difference(tree1, tree2):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.symmetric_difference()' function has moved to 'dendropy.calculate.treecompare.symmetric_difference()'.",
            old_construct="from dendropy import treecalc\nd = treecalc.symmetric_difference(...)",
            new_construct="from dendropy.calculate import treecompare\nd = treecompare.symmetric_difference(...)")
    return treecompare.symmetric_difference(
            tree1=tree1,
            tree2=tree2)

def find_missing_splits(reference_tree, test_tree):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.find_missing_splits()' function has moved to 'dendropy.calculate.treecompare.find_missing_splits()'.",
            old_construct="from dendropy import treecalc\nd = treecalc.find_missing_splits(...)",
            new_construct="from dendropy.calculate import treecompare\nd = treecompare.find_missing_splits(...)")
    return treecompare.find_missing_splits(
            reference_tree=reference_tree,
            comparison_tree=test_tree)

##############################################################################
## dendropy.calculate.treescore

# def fitch_up_pass(preorder_node_list, attr_name="state_sets", taxa_to_state_set_map=None):
#     deprecate.dendropy_deprecation_warning(
#             preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.euclidean_distance()' function has moved to 'dendropy.calculate.treecompare.euclidean_distance()'.",
#             old_construct="from dendropy import treecalc\nd = treecalc.euclidean_distance(...)",
#             new_construct="from dendropy.calculate import treecompare\nd = treecompare.euclidean_distance(...)")
#     return treecompare.euclidean_distance(
#             tree1=tree1,
#             tree2=tree2,
#             edge_weight_attr=edge_length_attr,
#             value_type=value_type)

# def fitch_down_pass(postorder_node_list, attr_name="state_sets", weight_list=None, taxa_to_state_set_map=None):
#     deprecate.dendropy_deprecation_warning(
#             preamble="Deprecated since DendroPy 4: The 'dendropy.treecalc.euclidean_distance()' function has moved to 'dendropy.calculate.treecompare.euclidean_distance()'.",
#             old_construct="from dendropy import treecalc\nd = treecalc.euclidean_distance(...)",
#             new_construct="from dendropy.calculate import treecompare\nd = treecompare.euclidean_distance(...)")
#     return treecompare.euclidean_distance(
#             tree1=tree1,
#             tree2=tree2,
#             edge_weight_attr=edge_length_attr,
#             value_type=value_type)

# def mason_gamer_kellogg_score(tree1, tree2):

