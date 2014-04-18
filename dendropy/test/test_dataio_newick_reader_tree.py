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

    def test_basic_parsing(self):
        tree_string = self.get_newick_string()
        reader_kwargs = {}
        with pathmap.SandboxedFile() as tempf:
            tempf.write(tree_string)
            tempf.flush()
            tree_filepath = tempf.name
            for suppress_internal_node_taxa in (None, False, True):
                if suppress_internal_node_taxa is None:
                    expected_suppress_internal_node_taxa = True
                    reader_kwargs.pop("suppress_internal_node_taxa", None)
                else:
                    expected_suppress_internal_node_taxa = suppress_internal_node_taxa
                    reader_kwargs["suppress_internal_node_taxa"] = suppress_internal_node_taxa
                for suppress_external_node_taxa in (None, False, True):
                    if suppress_external_node_taxa is None:
                        expected_suppress_external_node_taxa = False
                        reader_kwargs.pop("suppress_external_node_taxa", None)
                    else:
                        expected_suppress_external_node_taxa = suppress_external_node_taxa
                        reader_kwargs["suppress_external_node_taxa"] = suppress_external_node_taxa
                    for suppress_edge_lengths in (None, False, True):
                        if suppress_edge_lengths is None:
                            expected_suppress_edge_lengths = False
                            reader_kwargs.pop("suppress_edge_lengths", None)
                        else:
                            expected_suppress_edge_lengths = suppress_edge_lengths
                            reader_kwargs["suppress_edge_lengths"] = suppress_edge_lengths
                        # print("--")
                        # print(suppress_internal_node_taxa,
                        #         suppress_external_node_taxa,
                        #         suppress_edge_lengths)
                        # print(reader_kwargs)
                        # print(expected_suppress_internal_node_taxa,
                        #         expected_suppress_external_node_taxa,
                        #         expected_suppress_edge_lengths)
                        with open(tree_filepath, "r") as tree_stream:
                            approaches = (
                                    (dendropy.Tree.get_from_path, tree_filepath),
                                    (dendropy.Tree.get_from_stream, tree_stream),
                                    (dendropy.Tree.get_from_string, tree_string),
                                    )
                            for method, src in approaches:
                                t = method(src,
                                        "newick",
                                        **reader_kwargs
                                        )
                                self.verify_curated_tree(t,
                                        suppress_internal_node_taxa=expected_suppress_internal_node_taxa,
                                        suppress_external_node_taxa=expected_suppress_external_node_taxa,
                                        suppress_edge_lengths=expected_suppress_edge_lengths)
                        with open(tree_filepath, "r") as tree_stream:
                            approaches = (
                                    ("read_from_path", tree_filepath),
                                    ("read_from_stream", tree_stream),
                                    ("read_from_string", tree_string),
                                    )
                            for method, src in approaches:
                                t = dendropy.Tree()
                                f = getattr(t, method)
                                f(src, "newick", **reader_kwargs)
                                self.verify_curated_tree(t,
                                        suppress_internal_node_taxa=expected_suppress_internal_node_taxa,
                                        suppress_external_node_taxa=expected_suppress_external_node_taxa,
                                        suppress_edge_lengths=expected_suppress_edge_lengths)

    def test_unsupported_keyword_arguments(self):
        s = self.get_newick_string()
        with self.assertRaises(TypeError):
            t = dendropy.Tree.get_from_string(s,
                    "newick",
                    suppress_internal_taxa=True,  # should be suppress_internal_node_taxa
                    gobbledegook=False,
                    )

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
                    raise Exception("Unexpected rooting interpretation: '{}'".format(rooting_interpretation))
                _LOG.info("Rooting token = '{}', Rooting interpretation = '{}'".format(rooting_token, rooting_interpretation))
                s = self.get_newick_string(rooting_token=rooting_token)
                _LOG.debug(s)
                t = dendropy.Tree.get_from_string(s, "newick",
                        rooting=rooting_interpretation)
                self.assertIs(t.is_rooted, expected_is_rooted)
                self.assertIs(t.is_unrooted, expected_is_unrooted)
                self.assertIs(t.is_rootedness_undefined, expected_is_rootedness_undefined)

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

    def test_simple_post_node_comments(self):
        s = "((A[A]:1,B[B]:1)AB[AB]:1,(C[C]:1,D[D]:1)CD[CD]:1)Root[Root]:1;"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        for nd in tree:
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 1)
            self.assertEqual(nd.comments[0], nd.label)
            self.assertEqual(nd.edge.length, 1)

    def test_simple_post_edge_length_comments(self):
        s = "((A:1[A],B:1[B])AB:1[AB],(C:1[C],D:1[D])CD:1[CD])Root:1[Root];"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        for nd in tree:
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 1)
            self.assertEqual(nd.comments[0], nd.label)

    def test_post_node_and_edge_comments(self):
        s = "((A[A]:1[A],B[B]:1[B])AB[AB]:1[AB],(C[C]:1[C],D[D]:1[D])CD[CD]:1[CD])Root[Root]:1[Root];"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                )
        for nd in tree:
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 2)
            self.assertEqual(nd.comments[0], nd.label)
            self.assertEqual(nd.comments[1], nd.label)

    def test_multi_position_comments(self):
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
        for nd in tree:
            _LOG.info("%s: %s" % (nd.label, nd.comments))
            self.assertEqual(len(nd.comments), 7)
            self.assertEqual(nd.comments[0], 'xxx')
            for i in range(1,7):
                self.assertEqual(nd.comments[i], nd.label)

    def test_comment_association(self):
        tree_strings = [
                "([a1][a2]a[a3]:[a4]1[a5],[h1][h2][h3]([b1]b[b2]:[b3][b4]2[b5],[g1][g2]([c1]c[c2]:[c3]3[c4][c5],[f1]([d1]d[d2][d3]:[d4]4[d5],[e1]e[e2]:5[e3][e4][e5])[f2]f[f3]:[f4]6[f5])[g3][g4]g:7[g5])[h4]h[h5]:8)[i1][i2]i[i3]:[i4]9[i5];",
            ]
        for tree_string in tree_strings:
            tree = dendropy.Tree.get_from_string(
                    tree_string,
                    "newick",
                    suppress_external_node_taxa=True)
            for nd in tree:
                exp_brlen = ord(nd.label[0]) - ord('a') + 1
                self.assertEqual(nd.edge.length, exp_brlen)
                self.assertEqual(len(nd.comments), 5)
                for comment in nd.comments:
                    self.assertEqual(comment[0], nd.label)

    def test_anonymous_node_comment_association(self):
        tree_string1 = "[x1]([x2],[x3]([x4],[x5]([x6],[x7],[x8],[x9])[x10])[x11])[x12];"
        tree_string2 = "[x1](a[x2],[x3]([x4]b,[x5]([x6]c,[x7]d,[x8]e,[x9]f)g[x10])h[x11])i[x12];"
        tree1 = dendropy.Tree.get_from_string(tree_string1, "newick", suppress_external_node_taxa=True)
        tree2 = dendropy.Tree.get_from_string(tree_string2, "newick", suppress_external_node_taxa=True)
        expected_comments = {
                "a": ["x2",],
                "b": ["x4",],
                "c": ["x6",],
                "d": ["x7",],
                "e": ["x8",],
                "f": ["x9",],
                "g": ["x5", "x10"],
                "h": ["x3", "x11"],
                "i": ["x12"],
                }
        nodes1 = [nd for nd in tree1]
        nodes2 = [nd for nd in tree2]
        for nd1, nd2 in zip(nodes1, nodes2):
            self.assertEqual(nd2.comments, expected_comments[nd2.label])
            self.assertEqual(nd1.comments, expected_comments[nd2.label])


class CommentMetaDataTests(unittest.TestCase):
    figtree_metadata_str = """[&Tree1=1,Tree2=2, Tree3={1,2,3}](([xxx]A[&A1=1,A2=2,A3={1,2,3},  ,][A]:[A][A]1[A][A],
                 [xxx]B[&B1=1,B2=2,B3={1,2,3}][B]:[B][B]1[B][B])
                 [xxx]AB[&AB1=1,AB2=2,AB3={1,2,3}][AB]:[AB][AB]1[AB][AB],
                ([xxx]C[&C1=1,C2=2,C3={1,2,3}][C]:[C][C]1[C][C],
                 [xxx]D[&D1=1,D2=2,D3={1,2,3}][D]:[D][D]1[D][D])
                 [xxx]CD[&CD1=1, CD2=2, CD3={1,2,3}][CD]:[CD][CD]1[CD][CD])
                 [xxx]Root[&Root1=1, Root2=2, Root3={1,2,3}][Root]:[Root][Root]1[Root][Root];"""

    nhx_metadata_str = """[&Tree1=1,Tree2=2, Tree3={1,2,3}](([xxx]A[&&A1=1:A2=2:A3={1,2,3}][A]:[A][A]1[A][A],
                 [xxx]B[&&B1=1:B2=2:B3={1,2,3}][B]:[B][B]1[B][B])
                 [xxx]AB[&&AB1=1:AB2=2:AB3={1,2,3}][AB]:[AB][AB]1[AB][AB],
                ([xxx]C[&&C1=1:C2=2:C3={1,2,3}][C]:[C][C]1[C][C],
                 [xxx]D[&&D1=1:D2=2:D3={1,2,3}][D]:[D][D]1[D][D])
                 [xxx]CD[&&CD1=1: CD2=2: CD3={1,2,3}][CD]:[CD][CD]1[CD][CD])
                 [xxx]Root[&&Root1=1: Root2=2: Root3={1,2,3}][Root]:[Root][Root]1[Root][Root];"""

    def check_results(self, tree):
        metadata = tree.annotations.values_as_dict()
        self.assertEqual(metadata, {'Tree1': '1', 'Tree2': '2', 'Tree3':['1','2','3']})
        for nd in tree.postorder_node_iter():
            metadata = nd.annotations.values_as_dict()
            #print("%s: %s => %s" % (nd.label, nd.comments, metadata))
            self.assertEqual(len(metadata), 3)
            values = ["1", "2", ["1","2","3"]]
            for i in range(3):
                key = "{}{}".format(nd.label, i+1)
                self.assertTrue(key in metadata)
                self.assertEqual(metadata[key], values[i])

    def testFigtreeStyleBasic(self):
        s = self.figtree_metadata_str
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                extract_comment_metadata=True)
        self.check_results(tree)

    def testNHXBasic(self):
        s = self.nhx_metadata_str
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_external_node_taxa=True,
                extract_comment_metadata=True)
        self.check_results(tree)

    def test_incomplete_metadata(self):
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

class NewickTreeTaxonHandlingTest(unittest.TestCase):

    def test_basic_taxa(self):
        s = "(a1:3.14e-2,(b2:1.2,(c3:0.5,d4:0.7)e5:111)f6:222)g7:333;"
        tree = dendropy.Tree.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=False)
        expected = {
                "a1": 3.14e-2,
                "b2": 1.2,
                "c3": 0.5,
                "d4": 0.7,
                "e5": 111,
                "f6": 222,
                "g7": 333,
                }
        tns = tree.taxon_namespace
        self.assertEqual(len(tns), len(expected))
        labels = set([t.label for t in tns])
        self.assertEqual(labels, set(expected.keys()))
        for nd in tree:
            self.assertEqual(nd.edge.length, expected[nd.taxon.label])

    def test_quoted_underscores(self):
        s = "('a_1':3.14e-2,('b_2':1.2,('c_3':0.5,'d_4':0.7)'e_5':111)'f_6':222)'g_7':333;"
        tree = dendropy.Tree.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=False)
        expected = {
                "a_1": 3.14e-2,
                "b_2": 1.2,
                "c_3": 0.5,
                "d_4": 0.7,
                "e_5": 111,
                "f_6": 222,
                "g_7": 333,
                }
        tns = tree.taxon_namespace
        self.assertEqual(len(tns), len(expected))
        labels = set([t.label for t in tns])
        self.assertEqual(labels, set(expected.keys()))
        for nd in tree:
            self.assertEqual(nd.edge.length, expected[nd.taxon.label])

    def test_unquoted_underscores(self):
        s = "(a_1:3.14e-2,(b_2:1.2,(c_3:0.5,d_4:0.7)e_5:111)f_6:222)g_7:333;"
        tree = dendropy.Tree.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=False)
        expected = {
                "a 1": 3.14e-2,
                "b 2": 1.2,
                "c 3": 0.5,
                "d 4": 0.7,
                "e 5": 111,
                "f 6": 222,
                "g 7": 333,
                }
        tns = tree.taxon_namespace
        self.assertEqual(len(tns), len(expected))
        labels = set([t.label for t in tns])
        self.assertEqual(labels, set(expected.keys()))
        for nd in tree:
            self.assertEqual(nd.edge.length, expected[nd.taxon.label])

    def test_unquoted_underscores_preserved(self):
        s = "(a_1:3.14e-2,(b_2:1.2,(c_3:0.5,d_4:0.7)e_5:111)f_6:222)g_7:333;"
        tree = dendropy.Tree.get_from_string(s,
                "newick",
                suppress_internal_node_taxa=False,
                preserve_underscores=True)
        expected = {
                "a_1": 3.14e-2,
                "b_2": 1.2,
                "c_3": 0.5,
                "d_4": 0.7,
                "e_5": 111,
                "f_6": 222,
                "g_7": 333,
                }
        tns = tree.taxon_namespace
        self.assertEqual(len(tns), len(expected))
        labels = set([t.label for t in tns])
        self.assertEqual(labels, set(expected.keys()))
        for nd in tree:
            self.assertEqual(nd.edge.length, expected[nd.taxon.label])

if __name__ == "__main__":
    unittest.main()
