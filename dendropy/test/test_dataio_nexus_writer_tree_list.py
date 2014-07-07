#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Tests for NEXUS tree list writing.
"""

import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import standard_file_test_trees
from dendropy.test.support import compare_and_validate

class NexusStandardTreeListWriterTestCase(
        compare_and_validate.ValidateWriteable,
        standard_file_test_trees.StandardTestTreeChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        cls.schema = "nexus"
        cls.schema_tree_filepaths = dict(standard_file_test_trees.tree_filepaths[cls.schema])
        cls.is_taxa_managed_separately_from_tree = True

    def test_annotated_tree_list_writing(self):
        tree_file_title = 'standard-test-trees-n33-annotated'
        tree_reference = standard_file_test_trees.tree_references[tree_file_title]
        expected_non_metadata_comments = tree_reference["tree_list_comments"]
        expected_metadata_comments = tree_reference["tree_list_metadata_comments"]
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
                tree_offset=0,
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=False,
                metadata_extracted=True,
                is_coerce_metadata_values_to_string=True,
                is_distinct_nodes_and_edges_representation=False,
                is_taxa_managed_separately_from_tree=True)

if __name__ == "__main__":
    unittest.main()
