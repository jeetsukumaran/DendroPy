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
Tests for NEXUS tree list writing.
"""

import unittest
import dendropy
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import dendropytest
from support import standard_file_test_trees
from support import compare_and_validate

class NexusStandardTreeListWriterTestCase(
        compare_and_validate.ValidateWriteable,
        standard_file_test_trees.NexusTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NexusTestTreesChecker.create_class_fixtures(cls)

    def test_annotated_tree_list_writing(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-annotated-x10a'
        tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
        expected_non_metadata_comments = tree_reference["tree_list_comments"]
        expected_metadata_comments = tree_reference["tree_list_metadata"]
        expected_metadata = tree_reference["tree_list_metadata"]
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_list1 = dendropy.TreeList.get_from_path(
                tree_filepath,
                "nexus",
                extract_comment_metadata=True)
        s = self.write_out_validate_equal_and_return(
                tree_list1, "nexus", {})
        tree_list2 = dendropy.TreeList.get_from_string(s,
                "nexus",
                extract_comment_metadata=True)
        self.verify_standard_trees(
                tree_list=tree_list2,
                tree_file_title=tree_file_title,
                tree_offset=0)

if __name__ == "__main__":
    unittest.main()
