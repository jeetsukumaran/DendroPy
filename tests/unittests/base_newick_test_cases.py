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

import dendropy
from support import standard_file_test_trees

class NewickTreeListReaderTaxaManagementBaseTestCase(standard_file_test_trees.NewickTestTreesChecker):

    def test_get(self):
        tree_file_title = "dendropy-test-trees-n12-x2"
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        kwargs = {
                "suppress_internal_node_taxa": self.__class__.suppress_internal_node_taxa,
                "suppress_leaf_node_taxa": self.__class__.suppress_leaf_node_taxa,
        }
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.TreeList.get_from_path, tree_filepath),
                        (dendropy.TreeList.get_from_stream, tree_stream),
                        (dendropy.TreeList.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree_list = method(src,
                            self.__class__.schema,
                            **kwargs)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title)

    def test_selective_taxa_read(self):
        tree_file_title = "dendropy-test-trees-n12-x2"
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        kwargs = {
            "suppress_internal_node_taxa": self.__class__.suppress_internal_node_taxa,
            "suppress_leaf_node_taxa": self.__class__.suppress_leaf_node_taxa,
        }
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
                    old_id = id(tree_list)
                    f = getattr(tree_list, method)
                    f(src, self.__class__.schema, **kwargs)
                    new_id = id(tree_list)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title)

