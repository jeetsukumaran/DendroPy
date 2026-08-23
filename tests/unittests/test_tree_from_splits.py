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
NEXUS data read/write parse/format tests.
"""

import unittest
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import pathmap
from dendropy.utility.messaging import get_logger
from dendropy.calculate import treecompare
import dendropy

_LOG = get_logger(__name__)

class TreeFromBipartitionsTest(unittest.TestCase):
    def testTrees(self):
        tree_files = [
                ("dendropy-test-trees-n33-unrooted-x100a.nexus", "force-unrooted", False),
                ("dendropy-test-trees-multifurcating-unrooted.nexus", "force-unrooted", False),
                ("pythonidae.beast.summary.tre", "force-rooted", True),
                ("primates.beast.mcct.medianh.tre", "force-rooted", True),
                ]
        for tree_file, rooting, is_rooted in tree_files:
            ref_tree = dendropy.Tree.get_from_path(pathmap.tree_source_path(tree_file),
                    "nexus",
                    rooting=rooting)
            bipartition_encoding = ref_tree.encode_bipartitions()
            t_tree = dendropy.Tree.from_bipartition_encoding(
                    bipartition_encoding,
                    taxon_namespace=ref_tree.taxon_namespace,
                    is_rooted=ref_tree.is_rooted)
            # t_tree.encode_bipartitions()
            _LOG.debug("--\n       File: {} ({})".format(tree_file, ref_tree.is_rooted))
            _LOG.debug("     Original: {}".format(ref_tree.as_string("newick")))
            _LOG.debug("Reconstructed: {}".format(t_tree.as_string("newick")))
            self.assertEqual(treecompare.symmetric_difference(ref_tree, t_tree), 0)

if __name__ == "__main__":
    unittest.main()
