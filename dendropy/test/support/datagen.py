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
Test data generation and verification.
"""

import dendropy

###############################################################################
## A simple four taxon tree

def four_taxon_tree1():
    taxa = dendropy.TaxonNamespace(['A', 'B', 'C', 'D'])
    tree = dendropy.Tree(taxon_namespace=taxa)
    assert tree.taxon_set == taxa
    tree.seed_node.oid = 'root'
    tree.seed_node.label = 'root'
    tree.seed_node.edge.length = 2.0
    i1 = tree.seed_node.new_child(oid='i1', label='i1')
    i1.edge.length = 2.5
    a = i1.new_child(oid='a', taxon=taxa.require_taxon(label='A'))
    a.edge.length = 3.5
    b = i1.new_child(oid='b', taxon=taxa.require_taxon(label='B'))
    b.edge.length = 3.5
    i2 = tree.seed_node.new_child(oid='i2', label='i2')
    i2.edge.length = 4.0
    c = i2.new_child(oid='c', taxon=taxa.require_taxon(label='C'))
    c.edge.length = 2.0
    d = i2.new_child(oid='d', taxon=taxa.require_taxon(label='D'))
    d.edge.length = 2.0
    return tree


###############################################################################
## NodeRelationship

class NodeRelationship(object):

    def from_tree(tree):
        return [NodeRelationship(node) for node in tree]
    from_tree = staticmethod(from_tree)

    def from_node(node):
        ndr = NodeRelationship(None, None, None, None)
        if node.parent_node is not None:
            ndr.parent_label = node.parent_node.label
        else:
            ndr.parent_label = 'None'
        ndr.child_labels = [cnd.label for cnd in node.child_nodes()]
        ndr.edge_length = node.edge.length
        if node.taxon is not None:
            ndr.taxon_label = node.taxon.label
        else:
            ndr.taxon_label = 'None'
        return ndr
    from_node = staticmethod(from_node)

    def __init__(self, parent_label, child_labels, edge_length, taxon_label):
        self.parent_label = parent_label
        self.child_labels = child_labels
        self.edge_length = edge_length
        self.taxon_label = taxon_label

    def test_node(self, testcase, node):
        if self.parent_label is not None:
            testcase.assertTrue(node.parent_node is not None)
            testcase.assertEqual(self.parent_label, node.parent_node.label)
        else:
            testcase.assertTrue(node.parent_node is None or node.parent_node.label is None)
        testcase.assertEqual(self.child_labels, [cnd.label for cnd in node.child_nodes()])
        if self.edge_length is not None:
            testcase.assertTrue(node.edge.length is not None)
            testcase.assertAlmostEqual(self.edge_length, node.edge.length)
        else:
            testcase.assertTrue(node.edge.length is None)
        if self.taxon_label is not None:
            testcase.assertTrue(node.taxon is not None)
            testcase.assertEqual(self.taxon_label, node.taxon.label)
        else:
            testcase.assertTrue(node.taxon is None or node.taxon.label is None)
