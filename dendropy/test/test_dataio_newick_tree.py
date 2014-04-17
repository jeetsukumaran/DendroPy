#! /usr/bin/env python

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
Tests for general NEWICK reading.
"""

import sys
import os
import unittest
import dendropy
from dendropy.test.support import datagen_standard_file_test_trees
from dendropy.test.support import datagen_curated_test_tree
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)
if sys.hexversion < 0x03040000:
    from dendropy.utility.filesys import pre_py34_open as open

class NewickTreeReaderBasic(
        datagen_curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def test_default(self):
        s = self.get_newick_string()
        t = dendropy.Tree.get_from_string(s, "newick")
        self.verify_curated_tree(t,
                suppress_internal_taxa=True,
                suppress_external_taxa=False,
                suppress_edge_lengths=False)

    def test_basic_options(self):
        s = self.get_newick_string()
        for suppress_internal_taxa in (False, True):
            for suppress_external_taxa in (False, True):
                for suppress_edge_lengths in (False, True):
                    t = dendropy.Tree.get_from_string(s,
                            "newick",
                            suppress_internal_node_taxa=suppress_internal_taxa,
                            suppress_external_node_taxa=suppress_external_taxa,
                            suppress_edge_lengths=suppress_edge_lengths)
                    self.verify_curated_tree(t,
                            suppress_internal_taxa=suppress_internal_taxa,
                            suppress_external_taxa=suppress_external_taxa,
                            suppress_edge_lengths=suppress_edge_lengths)

    def test_unsupported_keyword_arguments(self):
        s = self.get_newick_string()
        with self.assertRaises(TypeError):
            t = dendropy.Tree.get_from_string(s,
                    "newick",
                    suppress_internal_taxa=True,
                    suppress_external_taxa=False)

    def test_rooting_interpretation(self):
        rooting_tokens = ("", "[&R]", "[&U]", "[&r]", "[&u]", "[&0]", "[&invalid]", "[R]", "[U]", "[&]")
        rooting_interpretations = ("force-rooted", "force-unrooted", "default-rooted", "default-unrooted", None)
        for rooting_token in rooting_tokens:
            for rooting_interpretation in rooting_interpretations:
                if rooting_interpretation == "force-rooted":
                    expected_is_rooted = True
                    expected_is_unrooted = False
                    expected_is_rootedness_undefined = False
                elif rooting_interpretation == "force-unrooted":
                    expected_is_rooted = False
                    expected_is_unrooted = True
                    expected_is_rootedness_undefined = False
                elif rooting_interpretation == "default-rooted":
                    if rooting_token.upper() == "[&R]":
                        expected_is_rooted = True
                        expected_is_unrooted = False
                        expected_is_rootedness_undefined = False
                    elif rooting_token.upper() == "[&U]":
                        expected_is_rooted = False
                        expected_is_unrooted = True
                        expected_is_rootedness_undefined = False
                    else:
                        expected_is_rooted = True
                        expected_is_unrooted = False
                        expected_is_rootedness_undefined = False
                elif rooting_interpretation == "default-unrooted":
                    if rooting_token.upper() == "[&R]":
                        expected_is_rooted = True
                        expected_is_unrooted = False
                        expected_is_rootedness_undefined = False
                    elif rooting_token.upper() == "[&U]":
                        expected_is_rooted = False
                        expected_is_unrooted = True
                        expected_is_rootedness_undefined = False
                    else:
                        expected_is_rooted = False
                        expected_is_unrooted = True
                        expected_is_rootedness_undefined = False
                else:
                    if rooting_token.upper() == "[&R]":
                        expected_is_rooted = True
                        expected_is_unrooted = False
                        expected_is_rootedness_undefined = False
                    elif rooting_token.upper() == "[&U]":
                        expected_is_rooted = False
                        expected_is_unrooted = True
                        expected_is_rootedness_undefined = False
                    else:
                        expected_is_rooted = None
                        expected_is_unrooted = None
                        expected_is_rootedness_undefined = True
                _LOG.info("Rooting token = '{}', Rooting interpretation = '{}'".format(rooting_token, rooting_interpretation))
                s = self.get_newick_string(rooting_token=rooting_token)
                t = dendropy.Tree.get_from_string(s, "newick",
                        rooting=rooting_interpretation)
                self.assertIs(t.is_rooted, expected_is_rooted)
                self.assertIs(t.is_unrooted, expected_is_unrooted)
                self.assertIs(t.is_rootedness_undefined, expected_is_rootedness_undefined)

        with self.assertRaises(TypeError):
            t = dendropy.Tree.get_from_string(s,
                    "newick",
                    suppress_internal_taxa=True,
                    suppress_external_taxa=False)


if __name__ == "__main__":
    unittest.main()
