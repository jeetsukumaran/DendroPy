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
Tests for general NEXUS tree list reading.
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

class NexusStandardTreeParsingTestCase(
        standard_file_test_trees.StandardTreeParsingTestCase,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.StandardTreeParsingTestCase.create_class_fixtures(
                cls,
                schema="nexus",
                is_taxa_managed_separately_from_tree=False)

    ## NOTE: many tests are in standard_file_test_trees.StandardTreeParsingTestCase !! ##

    def test_collection_comments_and_annotations(self):
        tree_file_title = 'standard-test-trees-n33-annotated'
        tree_reference = standard_file_test_trees.tree_references[tree_file_title]
        expected_non_metadata_comments = tree_reference["tree_list_comments"]
        expected_metadata_comments = tree_reference["tree_list_metadata_comments"]
        expected_metadata = tree_reference["tree_list_metadata"]
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        for extract_comment_metadata in (True, False):
            tree_list = dendropy.TreeList.get_from_path(
                    tree_filepath,
                    "nexus",
                    extract_comment_metadata=extract_comment_metadata)
            if extract_comment_metadata:
                expected_comments = expected_non_metadata_comments
                self.compare_annotations_to_json_metadata_dict(
                        tree_list,
                        expected_metadata,
                        is_coerce_metadata_values_to_string=True)
            else:
                expected_comments = ["&" + ",".join(s for s in expected_metadata_comments)] + expected_non_metadata_comments
                # expected_comments = [expected_metadata_comments] + expected_non_metadata_comments
                tree_list_metadata = tree_list.annotations.values_as_dict()
                self.assertEqual(tree_list_metadata, {})
            self.assertEqual(len(tree_list.comments), len(expected_comments))
            self.assertEqual(set(tree_list.comments), set(expected_comments))
            self.verify_standard_trees(
                    tree_list=tree_list,
                    tree_file_title=tree_file_title,
                    tree_offset=0,
                    suppress_internal_node_taxa=True,
                    suppress_leaf_node_taxa=False,
                    is_metadata_extracted=extract_comment_metadata,
                    is_coerce_metadata_values_to_string=True,
                    is_distinct_nodes_and_edges_representation=False,
                    is_taxa_managed_separately_from_tree=True)

class NexusMultiTreeListTestCase(dendropytest.ExtendedTestCase):

    def test_multiple_trees1(self):
        src_filename = "multitreeblocks.nex"
        src_path = pathmap.tree_source_path(src_filename)
        trees = dendropy.TreeList.get_from_path(src_path, "nexus")
        self.assertEqual(len(trees), 9)

    def test_multiple_trees2(self):
        src_filename = "multitreeblocks2.nex"
        src_path = pathmap.tree_source_path(src_filename)
        trees = dendropy.TreeList.get_from_path(src_path, "nexus")
        self.assertEqual(len(trees), 4)
        labels = ["x2.1","x2.2","x2.3","x2.4"]
        # self.assertEqual(len(trees.taxon_namespace), len(labels))
        self.assertEqual([t.label for t in trees.taxon_namespace], labels)
        for tree in trees:
            self.assertIs(tree.taxon_namespace, trees.taxon_namespace)
            seen_taxa = 0
            for nd in tree:
                if nd.taxon is not None:
                    seen_taxa += 1
                    self.assertIn(nd.taxon, tree.taxon_namespace)
            self.assertEqual(seen_taxa, len(tree.taxon_namespace))

if __name__ == "__main__":
    unittest.main()
