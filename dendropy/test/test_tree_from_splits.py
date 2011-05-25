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
NEXUS data read/write parse/format tests.
"""

import unittest
from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
from dendropy.utility.messaging import get_logger
from dendropy import treesplit
import dendropy

_LOG = get_logger(__name__)

class NexusTreeListReaderTest(datatest.DataObjectVerificationTestCase):

    def testReferenceTree(self):
        ref_tree_list = datagen.reference_tree_list()
        t_tree_list = dendropy.TreeList()
        for ref_tree in ref_tree_list:
            treesplit.encode_splits(ref_tree)
            splits = ref_tree.split_edges.keys()
            t_tree = treesplit.tree_from_splits(splits=splits,
                    taxon_set=ref_tree_list.taxon_set,
                    is_rooted=ref_tree.is_rooted)
            self.assertEqual(ref_tree.symmetric_difference(t_tree), 0)

    def testUltrametricTrees(self):
        tree_files = [
                "pythonidae.beast.summary.tre",
                "primates.beast.mcct.medianh.tre"
                ]

        for tree_file in tree_files:
            ref_tree = dendropy.Tree.get_from_path(pathmap.tree_source_path(tree_file),
                    "nexus",
                    as_rooted=True)
            treesplit.encode_splits(ref_tree)
            splits = ref_tree.split_edges.keys()
            t_tree = treesplit.tree_from_splits(splits=splits,
                    taxon_set=ref_tree.taxon_set,
                    is_rooted=ref_tree.is_rooted)
            treesplit.encode_splits(t_tree)
            self.assertEqual(ref_tree.symmetric_difference(t_tree), 0)

if __name__ == "__main__":
    unittest.main()
