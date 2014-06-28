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

class NexusStandardTreeListReaderTestCase(
        standard_file_test_trees.StandardTreeListReaderTestCase,
        dendropytest.ExtendedTestCase):

    @classmethod
    def build(cls):
        standard_file_test_trees.StandardTreeListReaderTestCase.build(schema="nexus",
                taxa_on_tree_equal_taxa_in_taxon_namespace=False)

    @classmethod
    def setUpClass(cls):
        cls.build()

    ## NOTE: many tests are in standard_file_test_trees.StandardTreeListReaderTestCase !! ##

    def test_collection_comments_and_annotations(self):
        tree_file_title = 'standard-test-trees-n33-annotated'
        expected_non_metadata_comments = ["esquireship2232",
                "puborectalis8922 nullity2204 belligerence4219 curb6499 amphanthium2186 montilla265",
                "grotesqueness1554 urd3366 stomachic6982 bibliology967 strippler64"]
        expected_metadata_comments = [
                '&!color="#0000ff",word1="saccharifier noninterpolation voiturette hypaton",word2="suspend physic tigerhearted",word3="reinstruct antitartaric dilli",word4="bartending discursativeness",a_95%_HPD={-79.7997141181,-44.405279378},a_mean=27.1368888529,b_95%_HPD={-27.033564206,-25.9638065766},b_mean=-58.0985938457'
                ]
        expected_metadata = {
                '!color': '"#0000ff"',
                'word1': '"saccharifier noninterpolation voiturette hypaton"',
                'word2': '"suspend physic tigerhearted"',
                'word3': '"reinstruct antitartaric dilli"',
                'word4': '"bartending discursativeness"',
                'a_95%_HPD': ["-79.7997141181","-44.405279378"],
                'a_mean': "27.1368888529",
                'b_95%_HPD': ["-27.033564206","-25.9638065766"],
                'b_mean': "-58.0985938457"
                }
        tree_reference = standard_file_test_trees.tree_references[tree_file_title]
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        for extract_comment_metadata in (True, False):
            tree_list = dendropy.TreeList.get_from_path(
                    tree_filepath,
                    "nexus",
                    extract_comment_metadata=extract_comment_metadata)
            if extract_comment_metadata:
                expected_comments = expected_non_metadata_comments
                tree_list_metadata = tree_list.annotations.values_as_dict()
                self.assertEqual(set(tree_list_metadata.keys()), set(expected_metadata.keys()))
                # for annote in tree_list.annotations:
                #     print("{}: {}".format(annote.name, annote.value))
                # for key in set(tree_list_metadata.keys()):
                #     if tree_list_metadata[key] != expected_metadata[key]:
                #         print("**** {}:\t\t{}\t\t{}".format(
                #             key,
                #             tree_list_metadata[key],
                #             expected_metadata[key]))
                self.assertEqual(tree_list_metadata, expected_metadata)
            else:
                expected_comments = expected_non_metadata_comments + expected_metadata_comments
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
                    metadata_extracted=extract_comment_metadata,
                    distinct_nodes_and_edges=False,
                    taxa_on_tree_equal_taxa_in_taxon_namespace=True)

class NexusMultiTreeListTestCase(dendropytest.ExtendedTestCase):

    def test_multiple_trees(self):
        src_filename = "multitreeblocks.nex"
        src_path = pathmap.tree_source_path(src_filename)
        trees = dendropy.TreeList.get_from_path(src_path, "nexus")
        self.assertEqual(len(trees), 9)

if __name__ == "__main__":
    unittest.main()
