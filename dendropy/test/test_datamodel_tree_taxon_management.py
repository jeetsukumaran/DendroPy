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
Tests Tree taxon management
"""

import unittest
import dendropy
import copy
from dendropy.test.support import datagen_curated_test_tree
from dendropy.test.support import compare_and_validate

class TestTreeTaxonManagement(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        self.tree1, self.anodes1, self.lnodes1, self.inodes1 = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True)
        self.expected_labels = set([nd.label for nd in self.anodes1 if nd.label is not None])
        self.expected_taxa = set()
        for nd in self.tree1:
            if nd.label is not None:
                nd.taxon = dendropy.Taxon(label=nd.label)
                self.expected_taxa.add(nd.taxon)

    # def test_infer_taxa(self):
    #     tns1 = self.tree1.taxon_namespace
    #     self.assertEqual(len(tns1), 0)
    #     tns2 = self.tree1.infer_taxa()
    #     self.assertNotEqual(len(tns2), 0)
    #     self.assertIs(self.tree1.taxon_namespace, tns2)
    #     self.assertIsNot(self.tree1.taxon_namespace, tns1)
    #     self.assertEqual(set(tns2._taxa), self.expected_taxa)

    def test_update_taxon_namespace(self):
        tns1 = self.tree1.taxon_namespace
        self.assertEqual(len(tns1), 0)
        tns2 = self.tree1.update_taxon_namespace()
        self.assertIs(self.tree1.taxon_namespace, tns1)
        self.assertEqual(set(tns2._taxa), self.expected_taxa)

    def test_reindex_subcomponent_taxa(self):
        tns1 = self.tree1.taxon_namespace
        self.assertEqual(len(tns1), 0)
        self.tree1.reindex_subcomponent_taxa()
        tns2 = self.tree1.taxon_namespace
        # self.assertIsNot(tns1, tns2) tns1 *IS* tns2
        self.assertEqual(len(tns2), len(self.expected_taxa))
        self.assertEqual(set(t.label for t in self.tree1.taxon_namespace),
                self.expected_labels)

    def test_unassign_taxa(self):
        pass

    def test_randomly_assign_taxa(self):
        pass

if __name__ == "__main__":
    unittest.main()
