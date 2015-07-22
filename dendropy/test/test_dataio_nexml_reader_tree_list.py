# !/usr/bin/env python

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
Tests for general NEXML tree list reading.
"""

import sys
import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import standard_file_test_trees
from dendropy.test import base_standard_trees_parsing_test_cases
from dendropy.test.support import curated_test_tree
from dendropy.test.support import pathmap
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

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
            tree_reference = dict(standard_file_test_trees._TREE_REFERENCES[tree_file_title])
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
