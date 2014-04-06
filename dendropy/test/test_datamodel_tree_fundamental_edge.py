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
Tests basic Edge operations.
"""

import copy
import unittest
from dendropy.test.support import compare_and_validate
import dendropy

class EdgeIdentity(unittest.TestCase):

    def setUp(self):
        self.e1 = dendropy.Edge(label="a")
        self.e2 = dendropy.Edge(label="a")

    def test_equal(self):
        # two distinct :class:`Edge` objects are never equal, even if all
        # member values are the same.
        self.assertNotEqual(self.e1, self.e2)

    def test_hash_dict_membership(self):
        k = {}
        k[self.e1] = 1
        k[self.e2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.e1], 1)
        self.assertEqual(k[self.e2], 2)
        self.assertIn(self.e1, k)
        self.assertIn(self.e2, k)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.e1)
        k.add(self.e2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.e1, k)
        self.assertIn(self.e2, k)

class EdgeCloning(compare_and_validate.Comparator, unittest.TestCase):

    def setUp(self):
        self.taxa = [dendropy.Taxon(label=label) for label in ["a", "b", "c", "d"]]
        self.n0 = dendropy.Node(label="0", taxon=self.taxa[0])
        self.c1 = dendropy.Node(label="1", taxon=None)
        self.c2 = dendropy.Node(label=None, taxon=self.taxa[1])
        self.c3 = dendropy.Node(label=None, taxon=None)
        self.c3 = dendropy.Node(label=None, taxon=self.taxa[2])
        self.p1 = dendropy.Node(label="-1", taxon=self.taxa[3])
        self.n0.parent_node = self.p1
        self.n0.set_child_nodes([self.c1, self.c2])
        self.c2.set_child_nodes([self.c3])
        self.nodes = [self.n0, self.c1, self.c2, self.c3, self.p1]
        for idx, nd in enumerate(self.nodes):
            if idx % 2 == 0:
                nd.edge.label = "E{}".format(idx)
                nd.edge.length = idx
            an1 = nd.annotations.add_new("a{}".format(idx),
                    "{}{}{}".format(nd.label, nd.taxon, idx))
            an2 = nd.annotations.add_bound_attribute("label")
            an3 = an1.annotations.add_bound_attribute("name")
            ae1 = nd.edge.annotations.add_new("a{}".format(idx),
                    "{}{}".format(nd.edge.label, idx))
            ae2 = nd.edge.annotations.add_bound_attribute("label")
            ae3 = ae1.annotations.add_bound_attribute("name")
        self.e0 = self.n0._edge

    def test_unsupported_copy(self):
        with self.assertRaises(TypeError):
            self.e0.clone(0)
        with self.assertRaises(TypeError):
            copy.copy(self.e0)
        with self.assertRaises(TypeError):
            self.e0.clone(1)
        with self.assertRaises(TypeError):
            self.e0.taxon_namespace_scoped_copy()

    def test_deepcopy(self):
        for clone in (
                self.e0.clone(2),
                copy.deepcopy(self.e0),
                ):
            self.compare_distinct_nodes(
                    clone._head_node, self.n0,
                    distinct_taxon_objects=True,
                    compare_annotations=True)

if __name__ == "__main__":
    unittest.main()
