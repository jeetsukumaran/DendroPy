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
Tests for NEWICK reading.
"""

import sys
import os
import unittest
import dendropy
from dendropy.test.support import standard_test_tree_data
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

class NewickTreeReaderTest(unittest.TestCase):


    def test_treelist_standard_newick_reader(self):
        for tree_filename in standard_test_tree_data.newick_tree_filenames:
            tree_file_title = os.path.splitext(os.path.basename(tree_filename))[0]
            trees = dendropy.TreeList.get_from_path(
                    pathmap.tree_source_path(tree_filename),
                    "newick",
                    suppress_internal_node_taxa=True,
                    suppress_external_node_taxa=True)
            self.assertEqual(len(trees), standard_test_tree_data.expected_number_of_trees[tree_file_title])
            for tree_idx, (tree, check_tree) in enumerate(zip(trees, standard_test_tree_data.tree_directory[tree_file_title])):
                _LOG.debug("{}: {}".format(tree_file_title, tree_idx))
                self.assertEqual(tree.comments, check_tree.comments)
                for nd in tree:
                    _LOG.debug("{}: {}: {}".format(tree_file_title, tree_idx, nd.label))
                    check_node = standard_test_tree_data.node_directory[(tree_file_title, tree_idx, nd.label)]
                    assert check_node.label == nd.label
                    self.assertEqual(nd.comments, check_node.comments)


if __name__ == "__main__":
    unittest.main()
