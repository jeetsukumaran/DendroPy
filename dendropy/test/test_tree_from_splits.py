#! /usr/bin/env python

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
NEXUS data read/write parse/format tests.
"""

import unittest
from dendropy.test.support import pathmap
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
