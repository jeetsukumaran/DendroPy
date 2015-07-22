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
Tests of ascii tree plots.
"""

import unittest
from dendropy.test.support import curated_test_tree
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

class AsciiTreeTest(
        curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.tree = tree1

    def test_plot_by_depth(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='depth'))

    def test_plot_by_level(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='level'))

    def test_plot_by_age(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='age'))

    def test_plot_by_length(self):
        _LOG.debug(self.tree.as_ascii_plot(plot_metric='length'))

if __name__ == "__main__":
    unittest.main()
