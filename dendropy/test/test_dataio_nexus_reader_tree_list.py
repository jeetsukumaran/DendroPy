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
Tests for general NEXUS tree list reading.
"""

import sys
import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import standard_file_test_trees
from dendropy.test.support import curated_test_tree
from dendropy.test.support import pathmap
from dendropy.test import base_standard_trees_parsing_test_cases
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

class NexusStandardTreeParsingTestCase(
        base_standard_trees_parsing_test_cases.StandardTreesParsingTestCase,
        standard_file_test_trees.NexusTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NexusTestTreesChecker.create_class_fixtures(cls)

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
                    "nexus")
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

class NexusStandardTreesWithTranslateBlockButNoTaxaBlockTestCase(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def test_with_translate_but_no_taxa_block(self):
        src_filename = "curated-with-translate-block-and-no-taxa-block-and-untranslated-internal-taxa.nex"
        src_path = pathmap.tree_source_path(src_filename)
        tree_list = dendropy.TreeList.get_from_path(src_path, "nexus")
        tree_labels = ("1", "2", "3")
        self.assertEqual(len(tree_list), len(tree_labels))
        for tree_idx, (tree, label) in enumerate(zip(tree_list, tree_labels)):
            self.assertEqual(tree.label, label)
            self.verify_curated_tree(tree=tree)

if __name__ == "__main__":
    unittest.main()
