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
Tests for NEWICK writing.
"""

import collections
import unittest
import dendropy
import re
from dendropy.test.support import pathmap
from dendropy.test.support import compare_and_validate
from dendropy.test.support import dendropytest
from dendropy.test.support import curated_test_tree

class NexusTreeWriterTests(
        curated_test_tree.CuratedTestTree,
        compare_and_validate.ValidateWriteable,
        unittest.TestCase):

    def test_simple(self):
        for suppress_internal_node_taxa in (True, False):
            for suppress_leaf_node_taxa in (True, False):
                tree1, all_nodes, leaf_nodes, internal_nodes = self.get_tree(
                            suppress_internal_node_taxa=suppress_internal_node_taxa,
                            suppress_leaf_node_taxa=suppress_leaf_node_taxa
                        )
                kwargs = {}
                s = self.write_out_validate_equal_and_return(
                        tree1, "nexus", kwargs)
                print(s)


if __name__ == "__main__":
    unittest.main()
