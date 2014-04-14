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
        self.tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True)

    def test_infer_taxa(self):
        pass

    def test_update_taxon_namespace(self):
        pass

    def test_reindex_subcomponent_taxa(self):
        pass

    def test_unassign_taxa(self):
        pass

    def test_randomly_assign_taxa(self):
        pass

if __name__ == "__main__":
    unittest.main()
