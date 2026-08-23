# !/usr/bin/env python

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
Tests for general NEXML tree list reading.
"""

import copy
import sys
import os
import unittest
import dendropy
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from support import dendropytest
from support import standard_file_test_trees
import tests.unittests.base_standard_trees_parsing_test_cases as base_standard_trees_parsing_test_cases
from support import curated_test_tree
from support import pathmap

class NexmlStandardTreeParsingTestCase(
        base_standard_trees_parsing_test_cases.StandardTreesParsingTestCase,
        standard_file_test_trees.NexmlTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NexmlTestTreesChecker.create_class_fixtures(cls)

    ## NOTE: many tests are in standard_file_test_trees.StandardTreeParsingTestCase !! ##

    def test_collection_comments_and_annotations(self):
        for tree_file_title in (
                "dendropy-test-trees-multifurcating-rooted-annotated",
                "dendropy-test-trees-n33-unrooted-annotated-x10a",
                ):
            tree_reference = copy.deepcopy(dict(standard_file_test_trees._TREE_REFERENCES[tree_file_title]))
            expected_non_metadata_comments = tree_reference["tree_list_comments"]
            expected_metadata = tree_reference["tree_list_metadata"]
            tree_filepath = self.schema_tree_filepaths[tree_file_title]
            tree_list = dendropy.TreeList.get_from_path(
                    tree_filepath,
                    "nexml")
            expected_comments = expected_non_metadata_comments
            self.compare_annotations_to_json_metadata_dict(
                    tree_list,
                    expected_metadata)
            if self.__class__.is_check_comments:
                self.assertEqual(len(tree_list.comments), len(expected_comments))
                self.assertEqual(set(tree_list.comments), set(expected_comments))
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    tree_offset=0)

if __name__ == "__main__":
    unittest.main()
