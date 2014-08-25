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
from dendropy.utility.messaging import get_logger
from dendropy import treesplit
import dendropy

_LOG = get_logger(__name__)

class TreeFromSplitsTest(unittest.TestCase):

    def testTrees(self):
        tree_files = [
                ("pythonidae.beast.summary.tre", "force-rooted", True),
                ("primates.beast.mcct.medianh.tre", "force-rooted", True),
                ("dendropy-test-trees-n33-unrooted-x100a.nexus", "force-unrooted", False),
                ("dendropy-test-trees-multifurcating-unrooted.nexus", "force-unrooted", False),
                ]
        for tree_file, rooting, is_rooted in tree_files:
            ref_tree = dendropy.Tree.get_from_path(pathmap.tree_source_path(tree_file),
                    "nexus",
                    rooting=rooting)
            treesplit.encode_splits(ref_tree)
            splits = ref_tree.split_edges.keys()
            t_tree = treesplit.tree_from_splits(
                    splits=splits,
                    taxon_namespace=ref_tree.taxon_namespace,
                    is_rooted=ref_tree.is_rooted)
            treesplit.encode_splits(t_tree)
            self.assertEqual(ref_tree.symmetric_difference(t_tree), 0)

if __name__ == "__main__":
    unittest.main()
