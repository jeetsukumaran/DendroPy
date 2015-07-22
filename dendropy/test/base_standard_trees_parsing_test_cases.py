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

import random
import dendropy
from dendropy.test.support import standard_file_test_trees

class StandardTreesParsingTestCase(standard_file_test_trees.StandardTestTreesChecker):

    def test_default_get(self):
        for tree_file_title in [
            "dendropy-test-trees-multifurcating-rooted",
            "dendropy-test-trees-multifurcating-unrooted",
            "dendropy-test-trees-n10-rooted-treeshapes",
            "dendropy-test-trees-n14-unrooted-treeshapes",
                ]:
            tree_filepath = self.schema_tree_filepaths[tree_file_title]
            with open(tree_filepath, "r") as src:
                tree_string = src.read()
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        {"path": tree_filepath},
                        {"file": tree_stream},
                        {"data": tree_string},
                        )
                for approach_kwargs in approaches:
                    approach_kwargs["schema"] = self.__class__.schema
                    tree_list = dendropy.TreeList.get(**approach_kwargs)
                    self.verify_standard_trees(tree_list=tree_list,
                            tree_file_title=tree_file_title)

    def test_default_read(self):
        preloaded_tree_file_title = "dendropy-test-trees-n33-unrooted-x10a"
        preloaded_tree_reference = self.tree_references[preloaded_tree_file_title]
        tree_file_title = "dendropy-test-trees-n33-unrooted-x10a"
        tree_reference = self.tree_references[tree_file_title]
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    {"path": tree_filepath},
                    {"file": tree_stream},
                    {"data": tree_string},
                    )
            for approach_kwargs in approaches:
                # prepopulate
                tree_list = dendropy.TreeList.get(
                        path=self.schema_tree_filepaths[preloaded_tree_file_title],
                        schema=self.__class__.schema)
                # check to make sure trees were loaded
                old_len = len(tree_list)
                self.assertEqual(old_len, len(tree_list._trees))
                self.assertEqual(old_len, preloaded_tree_reference["num_trees"])
                self.verify_standard_trees(tree_list, preloaded_tree_file_title)

                # load
                old_id = id(tree_list)
                approach_kwargs["schema"] = self.__class__.schema
                trees_read = tree_list.read(**approach_kwargs)
                new_id = id(tree_list)
                self.assertEqual(old_id, new_id)

                # make sure new trees added
                new_len = len(tree_list)
                self.assertEqual(new_len, len(tree_list._trees))
                expected_number_of_trees = tree_reference["num_trees"]
                self.assertEqual(old_len + expected_number_of_trees, new_len)
                self.assertEqual(trees_read, expected_number_of_trees)

                # check new trees
                for tree_idx, tree in enumerate(tree_list[old_len:]):
                    self.compare_to_reference_by_title_and_index(
                            tree=tree,
                            tree_file_title=tree_file_title,
                            reference_tree_idx=tree_idx)

                # make sure old ones still intact
                for tree_idx, tree in enumerate(tree_list[:old_len]):
                    self.compare_to_reference_by_title_and_index(
                            tree=tree,
                            tree_file_title=preloaded_tree_file_title,
                            reference_tree_idx=tree_idx)

    def test_tree_offset_get(self):
        tree_file_title = "dendropy-test-trees-n33-unrooted-x100a"
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        tree_offsets = set([0, expected_number_of_trees-1, -1, -expected_number_of_trees])
        while len(tree_offsets) < 8:
            tree_offsets.add(random.randint(1, expected_number_of_trees-2))
        while len(tree_offsets) < 12:
            tree_offsets.add(random.randint(-expected_number_of_trees-2, -2))
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for tree_offset in tree_offsets:
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.TreeList.get_from_path, tree_filepath),
                        (dendropy.TreeList.get_from_stream, tree_stream),
                        (dendropy.TreeList.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree_list = method(
                            src,
                            self.__class__.schema,
                            collection_offset=0,
                            tree_offset=tree_offset)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title,
                            tree_offset=tree_offset)

    def test_tree_offset_read(self):
        tree_file_title = "dendropy-test-trees-n33-unrooted-x100a"
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        tree_offsets = set([0, expected_number_of_trees-1, -1, -expected_number_of_trees])
        while len(tree_offsets) < 8:
            tree_offsets.add(random.randint(1, expected_number_of_trees-2))
        while len(tree_offsets) < 12:
            tree_offsets.add(random.randint(-expected_number_of_trees-2, -2))
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for tree_offset in tree_offsets:
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        ("read_from_path", tree_filepath),
                        ("read_from_stream", tree_stream),
                        ("read_from_string", tree_string),
                        )
                for method, src in approaches:
                    tree_list = dendropy.TreeList()
                    f = getattr(tree_list, method)
                    trees_read = f(src,
                            self.__class__.schema,
                            # collection_offset=0,
                            tree_offset=tree_offset)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title,
                            tree_offset=tree_offset)

    def test_out_of_range_tree_offset_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, self.__class__.schema, collection_offset=0, tree_offset=expected_number_of_trees)

    def test_out_of_range_tree_offset_read(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(IndexError):
                    f(src, self.__class__.schema, collection_offset=0, tree_offset=expected_number_of_trees)

    def test_out_of_range_collection_offset_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, self.__class__.schema, collection_offset=1, tree_offset=0)

    def test_out_of_range_collection_offset_read(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(IndexError):
                    f(src, self.__class__.schema, collection_offset=1, tree_offset=0)

    def test_unsupported_keyword_arguments_get(self):
        tree_file_title = 'dendropy-test-trees-n12-x2'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(TypeError):
                    method(src,
                            self.__class__.schema,
                            suppress_internal_taxa=True,  # should be suppress_internal_node_taxa
                            gobbledegook=False,
                            )

    def test_unsupported_keyword_arguments_read(self):
        tree_file_title = 'dendropy-test-trees-n12-x2'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(TypeError):
                    f(src,
                      self.__class__.schema,
                      suppress_internal_taxa=True,  # should be suppress_internal_node_taxa
                      gobbledegook=False,
                    )

