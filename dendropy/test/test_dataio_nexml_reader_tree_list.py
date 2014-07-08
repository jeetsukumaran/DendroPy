# !/usr/bin/env python

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
Tests for general NEXML tree list reading.
"""

import sys
import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import standard_file_test_trees
from dendropy.test.support import curated_test_tree
from dendropy.test.support import pathmap
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

class NexmlStandardTreeListReaderTestCase(
        standard_file_test_trees.StandardTreeListReaderTestCase,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.create_nexml_checker_class_fixtures(
                cls,
                schema="nexml",
                is_distinct_nodes_and_edges_representation=True,
                is_distinct_taxa_and_labels_on_tree=True,
                is_taxa_managed_separately_from_tree=False,
                is_coerce_metadata_values_to_string=False,
                is_check_comments=False,
                )

    ## NOTE: many tests are in standard_file_test_trees.StandardTreeListReaderTestCase !! ##

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
                    expected_metadata,
                    is_coerce_metadata_values_to_string=self.__class__.is_coerce_metadata_values_to_string)
            if self.__class__.is_check_comments:
                self.assertEqual(len(tree_list.comments), len(expected_comments))
                self.assertEqual(set(tree_list.comments), set(expected_comments))
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    tree_offset=0,
                    suppress_internal_node_taxa=True,
                    suppress_leaf_node_taxa=False,
                    is_metadata_extracted=True,
                    is_coerce_metadata_values_to_string=True,
                    is_distinct_nodes_and_edges_representation=False,
                    is_taxa_managed_separately_from_tree=True,
                    is_check_comments=self.__class__.is_check_comments)

if __name__ == "__main__":
    unittest.main()
