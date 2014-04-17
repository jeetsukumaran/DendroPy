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
Tests for general NEWICK reading.
"""

import sys
import os
import unittest
import dendropy
from dendropy.test.support import datagen_standard_file_test_trees
from dendropy.test.support import datagen_curated_test_tree
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

class NewickTreeReaderBasic(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def test_default(self):
        s = self.get_newick_string()
        t = dendropy.Tree.get_from_string(s, "newick")
        self.verify_curated_tree(t,
                suppress_internal_taxa=True,
                suppress_external_taxa=False,
                suppress_edge_lengths=False)

    def test_basic_options(self):
        s = self.get_newick_string()
        for suppress_internal_taxa in (False, True):
            for suppress_external_taxa in (False, True):
                for suppress_edge_lengths in (False, True):
                    t = dendropy.Tree.get_from_string(s,
                            "newick",
                            suppress_internal_node_taxa=suppress_internal_taxa,
                            suppress_external_node_taxa=suppress_external_taxa,
                            suppress_edge_lengths=suppress_edge_lengths)
                    self.verify_curated_tree(t,
                            suppress_internal_taxa=suppress_internal_taxa,
                            suppress_external_taxa=suppress_external_taxa,
                            suppress_edge_lengths=suppress_edge_lengths)



if __name__ == "__main__":
    unittest.main()
