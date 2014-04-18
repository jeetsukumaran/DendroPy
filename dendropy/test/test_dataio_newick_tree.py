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
from dendropy.utility import error
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
                elif rooting_interpretation is None:
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
                else:
                    raise Exception("Unexpecged rooting interpretation: '{}'".format(rooting_interpretation))
                _LOG.info("Rooting token = '{}', Rooting interpretation = '{}'".format(rooting_token, rooting_interpretation))
                s = self.get_newick_string(rooting_token=rooting_token)
                _LOG.debug(s)
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

class NewickTreeInvalidStatements(unittest.TestCase):

    def test_invalid_trees(self):
        invalid_tree_statements = (
            "(a,(b,c))a",
            "(a,(b,c)) (b,(a,c))",
            "(a,(b,c)) (d,(e,f))",
            "(a,(b,c)),",
            "(a,(b,c)z1)z2,",
            "(a,(b,c)))",
            "(a,(b,c)):",
            "(a,(b,c))(",
            )
        for s in invalid_tree_statements:
            with self.assertRaises(Exception):
                t = dendropy.Tree.get_from_string(s, "newick")

class NewickTreeQuotedLabels(unittest.TestCase):

    def test_edge_lengths1(self):
        tree = dendropy.Tree.get_from_string(
                """
                ((('T1 = 1.242e-10':1.242e-10,
                'T2 is 213.31e-4':213.31e-4)i1:3.44e-3,
                ('T3 is a (nice) taxon':3.3e7,
                T4:4.4e+8)'this is an internal node called "i2"':4.0e+1)i3:4.0E-4,
                (T5:6.7E+2,
                'and this so-called ''node'' is ("T6" with a length of ''7.2E-9'')':7.2E-9)i4:4.0E8)'this is the ''root\'\'\':7.0;
                """,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        expected_edge_lens = {
            'T1 = 1.242e-10': 1.242e-10,
            'T2 is 213.31e-4': 213.31e-4,
            'i1': 3.44e-3,
            'T3 is a (nice) taxon': 3.3e7,
            'T4': 4.4e+8,
            'this is an internal node called "i2"': 4.0e+1,
            'i3': 4.0e-4,
            'T5': 6.7e+2,
            "and this so-called 'node' is (\"T6\" with a length of '7.2E-9')": 7.2e-9,
            'i4': 4.0e8,
            "this is the 'root'": 7.0,
        }
        for nd in tree.postorder_node_iter():
            self.assertAlmostEqual(nd.edge.length, expected_edge_lens[nd.label])


class CommentReadingTests(unittest.TestCase):

    def testSimplePostNodeComments(self):
        s = "((A[A]:1,B[B]:1)AB[AB]:1,(C[C]:1,D[D]:1)CD[CD]:1)Root[Root]:1;"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 1)
            self.assertEqual(nd.comments[0], nd.label)
            self.assertEqual(nd.edge.length, 1)

    def testSimplePostEdgeLengthComments(self):
        s = "((A:1[A],B:1[B])AB:1[AB],(C:1[C],D:1[D])CD:1[CD])Root:1[Root];"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 1)
            self.assertEqual(nd.comments[0], nd.label)

    def testPostNodeAndEdgeLengthComments(self):
        s = "((A[A]:1[A],B[B]:1[B])AB[AB]:1[AB],(C[C]:1[C],D[D]:1[D])CD[CD]:1[CD])Root[Root]:1[Root];"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 2)
            self.assertEqual(nd.comments[0], nd.label)
            self.assertEqual(nd.comments[1], nd.label)

    def testMultiPositionComments(self):
        s = """(([xxx]A[A][A]:[A][A]1[A][A],
                 [xxx]B[B][B]:[B][B]1[B][B])
                 [xxx]AB[AB][AB]:[AB][AB]1[AB][AB],
                ([xxx]C[C][C]:[C][C]1[C][C],
                 [xxx]D[D][D]:[D][D]1[D][D])
                 [xxx]CD[CD][CD]:[CD][CD]1[CD][CD])
                 [xxx]Root[Root][Root]:[Root][Root]1[Root][Root];"""
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        for nd in tree.postorder_node_iter():
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 7)
            self.assertEqual(nd.comments[0], 'xxx')
            for i in range(1,7):
                self.assertEqual(nd.comments[i], nd.label)

    def testIncompleteMetadata(self):
        s = """[&color=blue](A[&region=Asia,id=00012][cryptic],(B[&region=Africa],C[&region=Madagascar,id=19391][two of three]));"""
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                extract_comment_metadata=True,
                )
        self.assertEqual(tree.annotations.values_as_dict(), {'color': 'blue'})
        expected = [ {'region': 'Asia', 'id': '00012'},
                {'region': 'Africa'},
                {'region': 'Madagascar', 'id': '19391'},
                {},
                {},]
        for idx, nd in enumerate(tree.postorder_node_iter()):
            self.assertEqual(nd.annotations.values_as_dict(), expected[idx])

if __name__ == "__main__":
    unittest.main()
