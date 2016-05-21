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
Tests for general NEWICK reading.
"""

import sys
import os
import unittest
import itertools
import collections
import random
import dendropy
from dendropy.utility import error
from dendropy.dataio import newickreader
from dendropy.test.support import dendropytest
from dendropy.test.support import compare_and_validate
from dendropy.test.support import standard_file_test_trees
from dendropy.test.support import curated_test_tree
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)
if sys.hexversion < 0x03040000:
    from dendropy.utility.filesys import pre_py34_open as open

class NewickTreeReaderBasic(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

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
                for suppress_leaf_node_taxa in (None, False, True):
                    if suppress_leaf_node_taxa is None:
                        expected_suppress_leaf_node_taxa = False
                        reader_kwargs.pop("suppress_leaf_node_taxa", None)
                    else:
                        expected_suppress_leaf_node_taxa = suppress_leaf_node_taxa
                        reader_kwargs["suppress_leaf_node_taxa"] = suppress_leaf_node_taxa
                    for suppress_edge_lengths in (None, False, True):
                        if suppress_edge_lengths is None:
                            expected_suppress_edge_lengths = False
                            reader_kwargs.pop("suppress_edge_lengths", None)
                        else:
                            expected_suppress_edge_lengths = suppress_edge_lengths
                            reader_kwargs["suppress_edge_lengths"] = suppress_edge_lengths
                        with open(tree_filepath, "r") as tree_stream:
                            approaches = (
                                    {"path": tree_filepath},
                                    {"file": tree_stream},
                                    {"data": tree_string},
                                    )
                            for approach_kwargs in approaches:
                                approach_kwargs.update(reader_kwargs)
                                approach_kwargs["schema"] = "newick"
                                t = dendropy.Tree.get(**approach_kwargs)
                                self.verify_curated_tree(t,
                                        suppress_internal_node_taxa=expected_suppress_internal_node_taxa,
                                        suppress_leaf_node_taxa=expected_suppress_leaf_node_taxa,
                                        suppress_edge_lengths=expected_suppress_edge_lengths)
                            # approaches = (
                            #         (dendropy.Tree.get_from_path, tree_filepath),
                            #         (dendropy.Tree.get_from_stream, tree_stream),
                            #         (dendropy.Tree.get_from_string, tree_string),
                            #         )
                            # for method, src in approaches:
                            #     t = method(src,
                            #             "newick",
                            #             **reader_kwargs
                            #             )
                            #     self.verify_curated_tree(t,
                            #             suppress_internal_node_taxa=expected_suppress_internal_node_taxa,
                            #             suppress_leaf_node_taxa=expected_suppress_leaf_node_taxa,
                            #             suppress_edge_lengths=expected_suppress_edge_lengths)
                        # with open(tree_filepath, "r") as tree_stream:
                        #     approaches = (
                        #             ("read_from_path", tree_filepath),
                        #             ("read_from_stream", tree_stream),
                        #             ("read_from_string", tree_string),
                        #             )
                        #     for method, src in approaches:
                        #         t = dendropy.Tree.get_from_string("(zzz1,(zzz2,(zzz3,zzz4)));",
                        #                 "newick",
                        #                 suppress_internal_node_taxa=False,
                        #                 suppress_leaf_node_taxa=False,)
                        #         tns0 = t.taxon_namespace
                        #         self.assertIs(t.taxon_namespace, tns0)
                        #         f = getattr(t, method)
                        #         self.assertIs(t.taxon_namespace, tns0)
                        #         f(src, "newick", **reader_kwargs)
                        #         self.verify_curated_tree(t,
                        #                 suppress_internal_node_taxa=expected_suppress_internal_node_taxa,
                        #                 suppress_leaf_node_taxa=expected_suppress_leaf_node_taxa,
                        #                 suppress_edge_lengths=expected_suppress_edge_lengths)

    def test_rooting_weighting_and_tree_metadata_handling(self):
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
                weighting_tokens = ("", "[&w 0.25]", "[ &W 0.25]", "[&w 1/4]", "[&W 1/4]")
                for weighting_token in weighting_tokens:
                    for store_tree_weights in (False, True):
                        for extract_comment_metadata in (False, True):
                            # pre_tree_token_candidates = ["[&color=blue]", weighting_token, rooting_token]
                            # for token_combination in itertools.permutations(pre_tree_token_candidates):
                            token_combination = ["[&color=blue]", "[&Why=42]", weighting_token, "[&rate=0.15]", rooting_token, "[&lnL=-3.14]"]
                            token_str = "".join(token_combination)
                            _LOG.debug("Token = '{}', Rooting interpretation = '{}'".format(token_str, rooting_interpretation))
                            s = self.get_newick_string(tree_preamble_tokens=token_str)
                            _LOG.debug(s)
                            t = dendropy.Tree.get(
                                    data=s,
                                    schema="newick",
                                    rooting=rooting_interpretation,
                                    store_tree_weights=store_tree_weights,
                                    extract_comment_metadata=extract_comment_metadata)
                            self.assertIs(t.is_rooted, expected_is_rooted)
                            self.assertIs(t.is_unrooted, expected_is_unrooted)
                            self.assertIs(t.is_rootedness_undefined, expected_is_rootedness_undefined)
                            if store_tree_weights:
                                if weighting_token:
                                    self.assertEqual(t.weight, 0.25)
                                else:
                                    self.assertEqual(t.weight, 1.0)
                            else:
                                self.assertIs(t.weight, None)
                            if extract_comment_metadata:
                                self.assertEqual(t.annotations.get_value("color", None), "blue")
                                self.assertEqual(t.annotations.get_value("Why", None), "42")

class NewickTreeMultifurcatingtree(dendropytest.ExtendedTestCase):

    def test_multifurcating(self):
        s = """\
                ([p]([a]a:1[a],[b]b:2[b],[c]c:3[c],[s]([d]d:4[d],
         [e]e:5[e],[f]f:6[f],[g]g:7[g])[s]s:19[s])[p]p:16[p],[w]
         ([t]t:20[t],[u]u:21[u],[v]v:22[v])[w]w:23[w],[q]
         ([h]h:8[h],[i]i:9[i],[j]j:10[j],[k]k:11[k],[o]([l]l:12[l],[m]m:13[m],[n]n:14[n])[o]o:15[o])[q]q:17[q])[r]r:18[r][r];
        """
        tree = dendropy.Tree.get(data=s,
                schema="newick",
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True)
        expected_children = {
            'a': [],
            'b': [],
            'c': [],
            'd': [],
            'e': [],
            'f': [],
            'g': [],
            'h': [],
            'i': [],
            'j': [],
            'k': [],
            'l': [],
            'm': [],
            'n': [],
            'o': ['l','m','n'],
            'p': ['s', 'a', 'b', 'c'],
            'q': ['o', 'h', 'i', 'j','k'],
            'r': ['q','p', 'w'],
            's': ['d', 'e', 'f', 'g'],
            't': [],
            'u': [],
            'v': [],
            'w': ['t', 'u', 'v'],
        }
        expected_parent = {
            'a': 'p',
            'b': 'p',
            'c': 'p',
            'd': 's',
            'e': 's',
            'f': 's',
            'g': 's',
            'h': 'q',
            'i': 'q',
            'j': 'q',
            'k': 'q',
            'l': 'o',
            'm': 'o',
            'n': 'o',
            'o': 'q',
            'p': 'r',
            'q': 'r',
            'r': None,
            's': 'p',
            't': 'w',
            'u': 'w',
            'v': 'w',
            'w': 'r',
        }
        for nd in tree:
            children = [ch.label for ch in nd.child_node_iter()]
            self.assertCountEqual(children, expected_children[nd.label])
            if nd.parent_node is not None:
                self.assertEqual(nd.parent_node.label, expected_parent[nd.label])
            else:
                self.assertIs(expected_parent[nd.label], None)
            if nd.is_leaf():
                self.assertEqual(len(nd.comments), 2)
            else:
                self.assertEqual(len(nd.comments), 3)
            for comment in nd.comments:
                self.assertEqual(comment, nd.label)
            self.assertEqual(nd.edge.length, ord(nd.label) - ord('a') + 1)

class NewickTreeInvalidStatements(dendropytest.ExtendedTestCase):

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
            "(e,(c,(d,e)a)b;(b,(a,e)c)d;",
            )
        for s in invalid_tree_statements:
            # t = dendropy.Tree.get_from_string(s, "newick")
            with self.assertRaises(error.DataParseError):
                t = dendropy.Tree.get(data=s, schema="newick")

class NewickTreeDuplicateTaxa(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def test_duplicate_taxa1(self):
        tree_statements = (
            "((a,b)c,(b,c)a)d;",
            "((_,_)_,(_,_)_)_;",
        )
        expected_labels = (
            ("a","b","c","b","c","a","d"),
            (" "," "," "," "," "," "," "),
        )
        for sidx, s in enumerate(tree_statements):
            with self.assertRaises(newickreader.NewickReader.NewickReaderDuplicateTaxonError):
                tree = dendropy.Tree.get(data=s, schema="newick")
            tree = dendropy.Tree.get(data=s,
                    schema="newick",
                    suppress_internal_node_taxa=True,
                    suppress_leaf_node_taxa=True)
            labels = [nd.label for nd in tree]
            self.assertCountEqual(labels, expected_labels[sidx])

class NewickTreeAnonymousTaxa(dendropytest.ExtendedTestCase):

    def test_anonymous_taxa_no_error(self):
        s = "((,),(,(,(,))));"
        tree = dendropy.Tree.get(data=s,
                schema="newick")

    def test_anonymous_taxa(self):
        s = "((:1[a],:2[b])[c]:3,(:4[d],([e]:5,([f]:6,:7[g]):8[h])[i]:9)[j]:10):11[k];"
        tree = dendropy.Tree.get(
                data=s,
                schema="newick")
        self.assertEqual(len(tree.taxon_namespace), 0)
        anodes = [nd for nd in tree]
        leaves = [nd for nd in tree.leaf_node_iter()]
        internal = [nd for nd in tree.postorder_internal_node_iter()]
        self.assertEqual(len(anodes), 11)
        self.assertEqual(len(leaves), 6)
        self.assertEqual(len(internal), 5)
        leaf_labels = [nd.comments[0] for nd in leaves]
        internal_labels = [nd.comments[0] for nd in internal]
        self.assertCountEqual(leaf_labels, ('a','b','d','e','f','g'))
        self.assertCountEqual(internal_labels, ('c','h','i','j','k'))
        for nd in tree:
            x = nd.comments[0]
            k = ord(x) - ord('a') + 1
            self.assertEqual(nd.edge.length, k)

class NewickTreeUnsupportedKeywordArguments(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def test_unsupported_keyword_arguments(self):
        tree_filepath = pathmap.tree_source_path('dendropy-test-trees-n12-x2.newick')
        tree_string = self.get_newick_string()
        reader_kwargs = {
                "suppress_internal_taxa": True,  # should be suppress_internal_node_taxa
                "gobbledegook": False,
        }
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.Tree.get_from_path, tree_filepath),
                    (dendropy.Tree.get_from_stream, tree_stream),
                    (dendropy.Tree.get_from_string, tree_string),
            )
            for method, src in approaches:
                with self.assertRaises(TypeError):
                    t = method(src, "newick", **reader_kwargs)
        # with open(tree_filepath, "r") as tree_stream:
        #     approaches = (
        #             ("read_from_path", tree_filepath),
        #             ("read_from_stream", tree_stream),
        #             ("read_from_string", tree_string),
        #     )
        #     for method, src in approaches:
        #         t = dendropy.Tree()
        #         tns0 = t.taxon_namespace
        #         self.assertIs(t.taxon_namespace, tns0)
        #         f = getattr(t, method)
        #         self.assertIs(t.taxon_namespace, tns0)
        #         with self.assertRaises(TypeError):
        #             f(src, "newick", **reader_kwargs)


class NewickTreeQuotedLabels(dendropytest.ExtendedTestCase):

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
                suppress_leaf_node_taxa=True,
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


class CommentReadingTests(dendropytest.ExtendedTestCase):

    def test_simple_post_node_comments(self):
        s = "((A[A]:1,B[B]:1)AB[AB]:1,(C[C]:1,D[D]:1)CD[CD]:1)Root[Root]:1;"
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True,
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
                suppress_leaf_node_taxa=True,
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
                suppress_leaf_node_taxa=True,
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
                suppress_leaf_node_taxa=True,
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
                    suppress_leaf_node_taxa=True)
            for nd in tree:
                exp_brlen = ord(nd.label[0]) - ord('a') + 1
                self.assertEqual(nd.edge.length, exp_brlen)
                self.assertEqual(len(nd.comments), 5)
                for comment in nd.comments:
                    self.assertEqual(comment[0], nd.label)

    def test_anonymous_node_comment_association(self):
        tree_string1 = "[x1]([x2],[x3]([x4],[x5]([x6],[x7],[x8],[x9])[x10])[x11])[x12];"
        tree_string2 = "[x1](a[x2],[x3]([x4]b,[x5]([x6]c,[x7]d,[x8]e,[x9]f)g[x10])h[x11])i[x12];"
        tree1 = dendropy.Tree.get_from_string(tree_string1, "newick", suppress_leaf_node_taxa=True)
        tree2 = dendropy.Tree.get_from_string(tree_string2, "newick", suppress_leaf_node_taxa=True)
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


class CommentMetaDataTests(dendropytest.ExtendedTestCase):
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
                suppress_leaf_node_taxa=True,
                extract_comment_metadata=True)
        self.check_results(tree)

    def testNHXBasic(self):
        s = self.nhx_metadata_str
        _LOG.info("Tree = %s" % s)
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True,
                extract_comment_metadata=True)
        self.check_results(tree)

    def test_incomplete_metadata(self):
        s = """[&color=blue](A[&region=Asia,id=00012][cryptic],(B[&region=Africa],C[&region=Madagascar,id=19391][two of three]));"""
        tree = dendropy.Tree.get_from_string(
                s,
                "newick",
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True,
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

# class NewickTreeTaxonNamespaceTest(dendropytest.ExtendedTestCase):

#     def test_namespace_passing(self):
#         tns1 = dendropy.TaxonNamespace()
#         s1 = "(a,(b,c));"
#         tree = dendropy.Tree.get_from_string(
#                 s1, "newick", taxon_namespace=tns1)
#         self.assertIs(tree.taxon_namespace, tns1)
#         self.assertEqual(len(tns1), 3)
#         s2 = "((e,f),(g,h));"
#         tree.read_from_string(
#                 s2, "newick")
#         self.assertIs(tree.taxon_namespace, tns1)
#         self.assertEqual(len(tns1), 7)
#         tns2 = dendropy.TaxonNamespace()
#         s3 = "((j,k),(l,m));"
#         with self.assertRaises(TypeError):
#             tree.read_from_string(
#                     s3, "newick",
#                     taxon_namespace=tns2)

class NewickTreeLabelParsingTest(dendropytest.ExtendedTestCase):

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

class NewickTreeReaderOffsetTreeTest(
        standard_file_test_trees.NewickTestTreesChecker,
        dendropytest.ExtendedTestCase):

    @classmethod
    def setUpClass(cls):
        standard_file_test_trees.NewickTestTreesChecker.create_class_fixtures(cls)

    def test_tree_offset_newick_get(self):
        tree_file_title = "dendropy-test-trees-n33-unrooted-x100a"
        tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
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
            tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
            expected_number_of_trees = tree_reference["num_trees"]
            if tree_offset < 0:
                if abs(tree_offset) > expected_number_of_trees:
                    tree_offset = 0
                else:
                    tree_offset = expected_number_of_trees + tree_offset
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.Tree.get_from_path, tree_filepath),
                        (dendropy.Tree.get_from_stream, tree_stream),
                        (dendropy.Tree.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree = method(
                            src,
                            "newick",
                            collection_offset=0,
                            tree_offset=tree_offset,
                            suppress_internal_node_taxa=True,
                            suppress_leaf_node_taxa=False,
                            rooting="default-unrooted")
                    reference_tree_idx = tree_offset
                    self.compare_to_reference_by_title_and_index(
                            tree=tree,
                            tree_file_title=tree_file_title,
                            reference_tree_idx=tree_offset)

    def test_offset_get_with_redundant_semicolons(self):
        s = """\
            ;;;;(a,(b,c)d)e;;;;(e,(c,a)d)b;;;;(b,(a,e)c)d;;;;
            """
        expected_roots = {
            0 : 'e',
            1 : 'b',
            2 : 'd',
        }
        expected_leaves = {
            0 : ['a', 'b', 'c'],
            1 : ['e', 'c', 'a'],
            2 : ['b', 'a', 'e'],
        }
        for idx in range(3):
            tree = dendropy.Tree.get_from_string(
                    s, "newick",
                    collection_offset=0,
                    tree_offset=idx,
                    suppress_internal_node_taxa=True,
                    suppress_leaf_node_taxa=True)
            self.assertEqual(tree.seed_node.label, expected_roots[idx])
            leaves = [nd.label for nd in tree.leaf_node_iter()]
            self.assertCountEqual(leaves, expected_leaves[idx])

    # def test_tree_offset_newick_read(self):
    #     tree_file_title = "dendropy-test-trees-n33-unrooted-x100a"
    #     tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
    #     expected_number_of_trees = tree_reference["num_trees"]
    #     tree_offsets = set([0, expected_number_of_trees-1, -1, -expected_number_of_trees])
    #     while len(tree_offsets) < 8:
    #         tree_offsets.add(random.randint(1, expected_number_of_trees-2))
    #     while len(tree_offsets) < 12:
    #         tree_offsets.add(random.randint(-expected_number_of_trees-2, -2))
    #     tree_filepath = self.schema_tree_filepaths[tree_file_title]
    #     with open(tree_filepath, "r") as src:
    #         tree_string = src.read()
    #     for tree_offset in tree_offsets:
    #         tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
    #         expected_number_of_trees = tree_reference["num_trees"]
    #         if tree_offset < 0:
    #             if abs(tree_offset) > expected_number_of_trees:
    #                 tree_offset = 0
    #             else:
    #                 tree_offset = expected_number_of_trees + tree_offset
    #         with open(tree_filepath, "r") as tree_stream:
    #             approaches = (
    #                     ("read_from_path", tree_filepath),
    #                     ("read_from_stream", tree_stream),
    #                     ("read_from_string", tree_string),
    #                     )
    #             for method, src in approaches:
    #                 tree = dendropy.Tree()
    #                 tns0 = tree.taxon_namespace
    #                 f = getattr(tree, method)
    #                 trees_read = f(src,
    #                         "newick",
    #                         collection_offset=0,
    #                         tree_offset=tree_offset,
    #                         suppress_internal_node_taxa=True,
    #                         suppress_leaf_node_taxa=False,
    #                         rooting="default-unrooted")
    #                 self.assertIs(tree.taxon_namespace, tns0)
    #                 reference_tree_idx = tree_offset
    #                 self.compare_to_reference_by_title_and_index(
    #                         tree=tree,
    #                         tree_file_title=tree_file_title,
    #                         reference_tree_idx=tree_offset)

    def test_tree_offset_without_collection_offset_newick_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.Tree.get_from_path, tree_filepath),
                    (dendropy.Tree.get_from_stream, tree_stream),
                    (dendropy.Tree.get_from_string, tree_string),
                    )
            for approach in approaches:
                tree_offset = 2
                tree = approach[0](approach[1], "newick", tree_offset=tree_offset)
                reference_tree_idx = tree_offset
                self.compare_to_reference_by_title_and_index(
                        tree=tree,
                        tree_file_title=tree_file_title,
                        reference_tree_idx=tree_offset)


    # def test_tree_offset_without_collection_offset_newick_read(self):
    #     tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
    #     tree_filepath = self.schema_tree_filepaths[tree_file_title]
    #     approaches = (
    #             "read_from_path",
    #             "read_from_stream",
    #             "read_from_string",
    #             )
    #     for approach in approaches:
    #         tree = dendropy.Tree()
    #         f = getattr(tree, approach)
    #         with self.assertRaises(TypeError):
    #             f(tree_filepath, "newick", collection_offset=None, tree_offset=0)

    def test_out_of_range_tree_offset_newick_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.Tree.get_from_path, tree_filepath),
                    (dendropy.Tree.get_from_stream, tree_stream),
                    (dendropy.Tree.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, "newick", collection_offset=0, tree_offset=expected_number_of_trees)

    # def test_out_of_range_tree_offset_newick_read(self):
    #     tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
    #     tree_filepath = self.schema_tree_filepaths[tree_file_title]
    #     tree_reference = standard_file_test_trees._TREE_REFERENCES[tree_file_title]
    #     expected_number_of_trees = tree_reference["num_trees"]
    #     with open(tree_filepath, "r") as src:
    #         tree_string = src.read()
    #     with open(tree_filepath, "r") as tree_stream:
    #         approaches = (
    #                 ("read_from_path", tree_filepath),
    #                 ("read_from_stream", tree_stream),
    #                 ("read_from_string", tree_string),
    #                 )
    #         for method, src in approaches:
    #             tree = dendropy.Tree()
    #             f = getattr(tree, method)
    #             with self.assertRaises(IndexError):
    #                 f(src, "newick", collection_offset=0, tree_offset=expected_number_of_trees)

    def test_out_of_range_collection_offset_newick_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.Tree.get_from_path, tree_filepath),
                    (dendropy.Tree.get_from_stream, tree_stream),
                    (dendropy.Tree.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, "newick", collection_offset=1, tree_offset=0)

    # def test_out_of_range_collection_offset_newick_read(self):
    #     tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
    #     tree_filepath = self.schema_tree_filepaths[tree_file_title]
    #     with open(tree_filepath, "r") as src:
    #         tree_string = src.read()
    #     with open(tree_filepath, "r") as tree_stream:
    #         approaches = (
    #                 ("read_from_path", tree_filepath),
    #                 ("read_from_stream", tree_stream),
    #                 ("read_from_string", tree_string),
    #                 )
    #         for method, src in approaches:
    #             tree = dendropy.Tree()
    #             f = getattr(tree, method)
    #             with self.assertRaises(IndexError):
    #                 f(src, "newick", collection_offset=1, tree_offset=0)

class JplaceParsingTest(dendropytest.ExtendedTestCase):

    def test_basic(self):
        s = "((A:.01[e]{0}, B:.02{1})D:.3{3}[g], C:.04{4}[h]) {5};"
        tree = dendropy.Tree.get_from_string(s,
                                             "newick",
                                             is_parse_jplace_tokens=True)
        self.assertEqual([0.01, 0.02, 0.3, 0.04, None],
                         [e.length for e in tree.edge_index])
        self.assertEqual('((A:0.01,B:0.02)D:0.3,C:0.04)',
                         str(tree))
        self.assertEqual([0,1,3,4,5],
                         [e.edge_number for e in tree.edge_index])

    def test_do_not_parse_by_default(self):
        s = "((A:.01[e]{0}, B:.02{1})D:.3{3}[g], C:.04{4}[h]) {5};"
        with self.assertRaises(newickreader.NewickReader.NewickReaderMalformedStatementError):
            try:
                tree = dendropy.Tree.get_from_string(s, "newick")
            except newickreader.NewickReader.NewickReaderMalformedStatementError as e:
                self.assertEqual(
                    "Expecting ':', ')', ',' or ';' after reading label but found '{'",
                    e.message)
                raise e

    def test_malformed_error_message_parsing_jplace(self):
        s = "((A}:.01[e]{0}, B:.02{1})D:.3{3}[g], C:.04{4}[h]) {5};"
        with self.assertRaises(newickreader.NewickReader.NewickReaderMalformedStatementError):
            try:
                tree = dendropy.Tree.get_from_string(s, "newick", is_parse_jplace_tokens=True)
            except newickreader.NewickReader.NewickReaderMalformedStatementError as e:
                self.assertEqual(
                    # extra 2nd element of list
                    "Expecting ':', '{', ')', ',' or ';' after reading label but found '}'",
                    e.message)
                raise e

class NewickInternalLabelAssociationTest(unittest.TestCase):

    def test_assign_to_edges(self):
        s = "((C:1.3,D:4.0)34:0.034,(A:1.1,(B:1.2,X:1.6)26:0.026)12:0.0126,E:1.5)seed;"
        expected_labels = {
                0 :"seed",
                60 :"34",
                28 :"12",
                24 :"26",
                }
        for labels_to_edges in (True, False):
            for rooting in ("force-rooted", "force-unrooted"):
                tree = dendropy.Tree.get(
                        data=s,
                        schema="newick",
                        is_assign_internal_labels_to_edges=labels_to_edges)
                tree.encode_bipartitions()
                for nd in tree:
                    if nd.is_leaf():
                        continue
                    expected_label = expected_labels[int(nd.bipartition)]
                    if labels_to_edges:
                        self.assertIs(nd.label, None)
                        self.assertEquals(nd.edge.label, expected_label)
                    else:
                        self.assertEquals(nd.label, expected_label)
                        self.assertIs(nd.edge.label, None)

if __name__ == "__main__":
    unittest.main()
