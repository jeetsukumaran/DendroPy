#! /usr/bin/env python

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
Tests for general NEWICK tree iteration reading.
"""

import sys
import unittest
import dendropy
from dendropy.test.support import dendropytest
from dendropy.test.support import standard_file_test_trees
from dendropy.test.support import pathmap

if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open

class NexusTreeYielderDefaultTestCase(
        standard_file_test_trees.NexusTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NexusTestTreesChecker.create_class_fixtures(cls)

    def test_basic(self):
        tree_file_titles = [
            # "dendropy-test-trees-multifurcating-rooted-annotated",
            # "dendropy-test-trees-multifurcating-rooted",
            # "dendropy-test-trees-multifurcating-unrooted",
            # "dendropy-test-trees-n10-rooted-treeshapes",
            "dendropy-test-trees-n12-x2",
            "dendropy-test-trees-n33-unrooted-x10a",
            "dendropy-test-trees-n33-unrooted-x10b",
            "dendropy-test-trees-n33-unrooted-annotated-x10a",
            "dendropy-test-trees-n33-unrooted-annotated-x10a",
        ]
        expected_file_names = []
        expected_tree_references = []
        tree_files = []
        for file_idx, tree_file_title in enumerate(tree_file_titles):
            tree_filepath = self.schema_tree_filepaths[tree_file_title]
            if False and idx % 2 == 0:
                tree_files.append(open(tree_filepath, "r"))
            else:
                tree_files.append(tree_filepath)
            num_trees = self.tree_references[tree_file_title]["num_trees"]
            for tree_idx in range(num_trees):
                expected_file_names.append(tree_filepath)
                expected_tree_references.append(self.tree_references[tree_file_title][str(tree_idx)])
        collected_trees = []
        tns = dendropy.TaxonNamespace()
        # for f in tree_files:
        #     dendropy.TreeList.get_from_path(f, "nexus")
        tree_sources = dendropy.Tree.yield_from_files(
                files=tree_files,
                schema="nexus",
                taxon_namespace=tns)
        for tree_idx, tree in enumerate(tree_sources):
            self.assertEqual(tree_sources.current_file_name, expected_file_names[tree_idx])
            tree.current_file_name = tree_sources.current_file_name
            collected_trees.append(tree)
        self.assertEqual(len(collected_trees), len(expected_tree_references))
        for tree, ref_tree in zip(collected_trees, expected_tree_references):
            self.assertIs(tree.taxon_namespace, tns)
            self.compare_to_reference_tree(tree, ref_tree)

## TODO:
# - test multiple trees blocks
# - mix of newick/nexus

if __name__ == "__main__":
    unittest.main()

