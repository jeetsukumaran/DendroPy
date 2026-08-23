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
Tests for NeXML tree writing.
"""

import collections
import unittest
import dendropy
import re
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import pathmap
from support import compare_and_validate
from support import dendropytest
from support import curated_test_tree
from support import standard_file_test_trees

class NexmlTreeWriterTests(
        curated_test_tree.CuratedTestTree,
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    def test_simple(self):
        tree1, all_nodes, leaf_nodes, internal_nodes = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False
                )
        s = tree1.as_string("nexml")
        tree2 = dendropy.Tree.get_from_string(
                s, "nexml",
                )
        self.verify_curated_tree(tree2)

class NexmlTreeWriterDefaultTest(
        standard_file_test_trees.NexmlTestTreesChecker,
        compare_and_validate.ValidateWriteable,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NexmlTestTreesChecker.create_class_fixtures(cls)

    def test_roundtrip_full(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-annotated-x10a'
        tree_filepath = standard_file_test_trees._TREE_FILEPATHS["nexml"][tree_file_title]
        tree1 = dendropy.Tree.get_from_path(
                tree_filepath,
                "nexml",
        )
        s = tree1.as_string("nexml")
        tree2 = dendropy.Tree.get_from_string(
                s,
                "nexml",
        )
        self.compare_to_reference_by_title_and_index(
            tree=tree2,
            tree_file_title=tree_file_title,
            reference_tree_idx=0)

class NexmlStandardTreeListWriterTestCase(
        compare_and_validate.ValidateWriteable,
        standard_file_test_trees.NexmlTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NexmlTestTreesChecker.create_class_fixtures(cls)

    def test_annotated_tree_list_writing(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-annotated-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_list1 = dendropy.TreeList.get_from_path(
                tree_filepath,
                "nexml")
        s = self.write_out_validate_equal_and_return(
                tree_list1, "nexml", {})
        tree_list2 = dendropy.TreeList.get_from_string(s,
                "nexml")
        self.verify_standard_trees(
                tree_list=tree_list2,
                tree_file_title=tree_file_title,
                tree_offset=0)

if __name__ == "__main__":
    unittest.main()
