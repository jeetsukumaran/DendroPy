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
Tests for dendropy.TreeList.
"""

import copy
import sys
import unittest
import collections
import dendropy
import random
from dendropy.test.support import dendropytest
from dendropy.test.support import curated_test_tree
from dendropy.test.support import curated_test_tree_list
from dendropy.test.support import compare_and_validate

class TestTreeListBasicOperations(dendropytest.ExtendedTestCase):

    def test_insert_simple_list_foreign_namespace(self):
        for idx in range(6):
            tlist = curated_test_tree_list.get_tree_list(5)
            self.assertEqual(len(tlist), 5)
            self.assertEqual(len(tlist._trees), 5)
            original_tns = tlist.taxon_namespace
            tree = curated_test_tree_list.get_tree()
            tlist.insert(idx, tree)
            self.assertEqual(len(tlist), 6)
            self.assertEqual(len(tlist._trees), 6)
            self.assertIs(tlist.taxon_namespace, original_tns)
            self.assertIn(tree, tlist)
            self.assertIs(tree.taxon_namespace, tlist.taxon_namespace)
            self.assertEqual(len(tlist.taxon_namespace), 7)
            for t1 in tlist:
                self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
                for nd in t1:
                    self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_insert_simple_list_native_namespace(self):
        for idx in range(6):
            tns = dendropy.TaxonNamespace()
            tlist = curated_test_tree_list.get_tree_list(5, taxon_namespace=tns)
            self.assertEqual(len(tlist), 5)
            self.assertEqual(len(tlist._trees), 5)
            original_tns = tlist.taxon_namespace
            tree = curated_test_tree_list.get_tree(taxon_namespace=tns)
            tlist.insert(idx, tree)
            self.assertEqual(len(tlist), 6)
            self.assertEqual(len(tlist._trees), 6)
            self.assertIs(tlist.taxon_namespace, original_tns)
            self.assertIn(tree, tlist)
            self.assertIs(tree.taxon_namespace, tlist.taxon_namespace)
            self.assertEqual(len(tlist.taxon_namespace), 7)
            for t1 in tlist:
                self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
                for nd in t1:
                    self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_append_simple_list_foreign_namespace(self):
        tlist, trees = curated_test_tree_list.get_tree_list_and_list_of_trees(num_trees=curated_test_tree_list.DEFAULT_NUM_TREES)
        original_tns = tlist.taxon_namespace
        for t in trees:
            tlist.append(t)
        self.assertEqual(len(tlist), curated_test_tree_list.DEFAULT_NUM_TREES)
        self.assertIs(tlist.taxon_namespace, original_tns)
        # self.assertEqual(len(tlist.taxon_namespace), len(tlist[0].tax_labels))
        self.assertEqual(len(tlist.taxon_namespace), 7)
        for t1, t2 in zip(tlist, trees):
            self.assertIs(t1, t2)
            self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
            for nd in t1:
                self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_append_simple_list_same_namespace(self):
        tns = dendropy.TaxonNamespace()
        tlist, trees = curated_test_tree_list.get_tree_list_and_list_of_trees(
                num_trees=curated_test_tree_list.DEFAULT_NUM_TREES,
                tree_list_taxon_namespace=tns,
                list_of_trees_taxon_namespace=tns)
        original_tns = tlist.taxon_namespace
        for t in trees:
            tlist.append(t)
        self.assertEqual(len(tlist), curated_test_tree_list.DEFAULT_NUM_TREES)
        self.assertIs(tlist.taxon_namespace, original_tns)
        # self.assertEqual(len(tlist.taxon_namespace), len(tlist[0].tax_labels))
        self.assertEqual(len(tlist.taxon_namespace), 7)
        for t1, t2 in zip(tlist, trees):
            self.assertIs(t1, t2)
            self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
            for nd in t1:
                self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_iadd_from_another_tree_list_different_namespace(self):
        tlist = curated_test_tree_list.get_tree_list(num_trees=3)
        original_tns = tlist.taxon_namespace
        original_tlist_len = len(tlist)
        original_tree_labels = [t.label for t in tlist]
        self.assertEqual(len(original_tree_labels), len(tlist))
        self.assertEqual(original_tlist_len, 3)

        tlist_source = curated_test_tree_list.get_tree_list(num_trees=5)
        self.assertEqual(len(tlist_source), 5)
        source_tree_labels = [t.label for t in tlist_source]
        self.assertEqual(len(source_tree_labels), len(tlist_source))

        tlist += tlist_source

        self.assertEqual(len(tlist), original_tlist_len + len(tlist_source))
        self.assertIs(tlist.taxon_namespace, original_tns)
        # self.assertEqual(len(tlist.taxon_namespace), len(tlist[0].tax_labels))
        self.assertEqual(len(tlist.taxon_namespace), 7)
        expected_tree_labels = original_tree_labels + source_tree_labels
        self.assertEqual(len(tlist), len(expected_tree_labels))
        for t1, tlabel in zip(tlist, expected_tree_labels):
            self.assertIn(t1, tlist)
            self.assertNotIn(t1, tlist_source)
            self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
            self.assertEqual(t1.label, tlabel)
            for nd in t1:
                self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_iadd_from_list_of_trees_different_namespace(self):
        tlist = curated_test_tree_list.get_tree_list(num_trees=3)
        original_tns = tlist.taxon_namespace
        original_tlist_len = len(tlist)
        original_tree_labels = [t.label for t in tlist]
        self.assertEqual(len(original_tree_labels), len(tlist))
        self.assertEqual(original_tlist_len, 3)
        source_trees = curated_test_tree_list.get_trees(
                num_trees=5,
                taxon_namespace=None,
                label=None,
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.assertEqual(len(source_trees), 5)
        source_tree_labels = [t.label for t in source_trees]
        self.assertEqual(len(source_tree_labels), len(source_trees))

        tlist += source_trees

        self.assertEqual(len(tlist), original_tlist_len + len(source_trees))
        self.assertIs(tlist.taxon_namespace, original_tns)
        # self.assertEqual(len(tlist.taxon_namespace), len(tlist[0].tax_labels))
        self.assertEqual(len(tlist.taxon_namespace), 7)
        expected_tree_labels = original_tree_labels + source_tree_labels
        self.assertEqual(len(tlist), len(expected_tree_labels))
        for t1, tlabel in zip(tlist, expected_tree_labels):
            self.assertIn(t1, tlist)
            if tlabel in source_tree_labels:
                self.assertIn(t1, source_trees)
            else:
                self.assertNotIn(t1, source_trees)
            self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
            self.assertEqual(t1.label, tlabel)
            for nd in t1:
                self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_add_from_another_tree_list_different_namespace(self):
        tlist_source1 = curated_test_tree_list.get_tree_list(num_trees=3)
        original_tns = tlist_source1.taxon_namespace
        source1_tree_labels = [t.label for t in tlist_source1]
        self.assertEqual(len(source1_tree_labels), len(tlist_source1))
        self.assertEqual(len(tlist_source1), 3)

        tlist_source2 = curated_test_tree_list.get_trees(num_trees=5)
        self.assertEqual(len(tlist_source2), 5)
        source2_tree_labels = [t.label for t in tlist_source2]
        self.assertEqual(len(source2_tree_labels), len(tlist_source2))

        tlist = tlist_source1 + tlist_source2

        self.assertEqual(len(tlist_source1), 3)
        self.assertEqual(len(tlist_source2), 5)

        self.assertEqual(len(tlist), len(tlist_source1) + len(tlist_source2))
        self.assertIs(tlist.taxon_namespace, original_tns)
        self.assertEqual(len(tlist.taxon_namespace), 7)
        expected_tree_labels = source1_tree_labels + source2_tree_labels
        self.assertEqual(len(tlist), len(expected_tree_labels))
        for t1, tlabel in zip(tlist, expected_tree_labels):
            self.assertIn(t1, tlist)
            self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
            self.assertEqual(t1.label, tlabel)
            if t1.label in source1_tree_labels:
                self.assertNotIn(t1, tlist_source1)
                self.assertNotIn(t1, tlist_source2)
            else:
                self.assertNotIn(t1, tlist_source1)
                self.assertIn(t1, tlist_source2)
            for nd in t1:
                self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_contains(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        self.assertEqual(len(tlist._trees), len(tlist))
        self.assertEqual(len(tlist), 5)
        trees = curated_test_tree_list.get_trees(5)
        self.assertEqual(len(trees), 5)
        for t in tlist:
            self.assertTrue(t in tlist._trees)
            self.assertTrue(t in tlist)
        for t in trees:
            self.assertFalse(t in tlist._trees)
            self.assertFalse(t in tlist)
        tlist += trees
        for t in trees:
            self.assertTrue(t in tlist._trees)
            self.assertTrue(t in tlist)

    def test_delitem(self):
        tsize = 5
        for del_idx in range(-tsize, tsize):
            tlist = curated_test_tree_list.get_tree_list(tsize)
            original_trees = list(tlist._trees)
            self.assertIn(original_trees[del_idx], tlist._trees)
            del tlist[del_idx]
            self.assertNotIn(original_trees[del_idx], tlist._trees)
            self.assertEqual(len(tlist), tsize - 1)
            del original_trees[del_idx]
            self.assertEqual(tlist._trees, original_trees)

    def test_iter(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        self.assertEqual(len(tlist), 5)
        self.assertEqual(len(tlist._trees), len(tlist))
        for t1, t2 in zip(tlist, tlist._trees):
            self.assertIs(t1, t2)

    def test_reversed(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        self.assertEqual(len(tlist), 5)
        self.assertEqual(len(tlist._trees), len(tlist))
        for t1, t2 in zip(reversed(tlist), reversed(tlist._trees)):
            self.assertIs(t1, t2)

    def test_getitem_simple(self):
        tsize = 5
        tlist = curated_test_tree_list.get_tree_list(tsize)
        self.assertEqual(len(tlist), tsize)
        self.assertEqual(len(tlist._trees), len(tlist))
        for idx in range(-tsize, tsize):
            self.assertIs(tlist[idx], tlist._trees[idx])
            self.assertTrue(isinstance(tlist[idx], dendropy.Tree))

    def test_getitem_slice(self):
        tsize = 5
        tlist = curated_test_tree_list.get_tree_list(tsize)
        self.assertEqual(len(tlist), tsize)
        self.assertEqual(len(tlist._trees), len(tlist))
        for a in range(-tsize, tsize):
            for b in range(-tsize, tsize):
                for step in range(-tsize, tsize):
                    if step == 0:
                        continue
                    tt = tlist[a:b:step]
                    k = tlist._trees[a:b:step]
                    self.assertEqual(len(tt), len(k))
                    for t1, t2 in zip(tt, k):
                        self.assertIn(t1, tlist)
                        self.assertIn(t1, tlist._trees)
                        self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)

    def test_setitem_simple(self):
        tsize = 5
        for idx in range(-tsize, tsize):
            tlist = curated_test_tree_list.get_tree_list(tsize)
            self.assertEqual(len(tlist), tsize)
            self.assertEqual(len(tlist._trees), len(tlist))
            old_tree = tlist[idx]
            new_tree = curated_test_tree_list.get_tree()
            tlist[idx] = new_tree
            self.assertIs(tlist[idx], new_tree)
            self.assertIsNot(tlist[idx], old_tree)
            self.assertIn(new_tree, tlist)
            self.assertNotIn(old_tree, tlist)
            self.assertIs(new_tree.taxon_namespace,
                    tlist.taxon_namespace)
            self.assertEqual(len(tlist.taxon_namespace), 7)
            for tree in tlist:
                for nd in tree:
                    self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_setitem_slice_from_list(self):
        tsize = 5
        for a in range(-tsize, tsize):
            for b in range(-tsize, tsize):
                for step in range(-tsize, tsize):
                    if step == 0:
                        continue
                    slice_obj = slice(a, b, step)
                    slice_len = len(range(*slice_obj.indices(tsize)))
                    if slice_len <= 0:
                        continue
                    tlist = curated_test_tree_list.get_tree_list(tsize)
                    self.assertEqual(len(tlist), tsize)
                    self.assertEqual(len(tlist._trees), len(tlist))
                    copy_list = list(tlist._trees)
                    source = curated_test_tree_list.get_trees(slice_len)
                    tlist[a:b:step] = source
                    copy_list[a:b:step] = source
                    expected_tree_labels = [t.label for t in copy_list]
                    self.assertEqual(len(tlist), len(copy_list))
                    self.assertEqual(len(tlist), len(tlist._trees))
                    self.assertEqual(len(tlist.taxon_namespace), 7)
                    for t1, t2, tlabel in zip(tlist, copy_list, expected_tree_labels):
                        self.assertIs(t1, t2)
                        self.assertIn(t1, tlist)
                        self.assertIn(t1, tlist._trees)
                        self.assertEqual(t1.label, tlabel)
                        self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
                        for nd in t1:
                            self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_setitem_slice_from_tree_list(self):
        tsize = 5
        for a in range(-tsize, tsize):
            for b in range(-tsize, tsize):
                for step in range(-tsize, tsize):
                    if step == 0:
                        continue
                    slice_obj = slice(a, b, step)
                    slice_indexes = list(range(*slice_obj.indices(tsize)))
                    slice_len = len(slice_indexes)
                    if slice_len <= 0:
                        continue
                    tlist = curated_test_tree_list.get_tree_list(tsize)
                    self.assertEqual(len(tlist), tsize)
                    self.assertEqual(len(tlist._trees), len(tlist))
                    copy_list = list(tlist._trees)
                    source = curated_test_tree_list.get_tree_list(slice_len)
                    copy_list[a:b:step] = source._trees
                    tlist[a:b:step] = source
                    expected_tree_labels = [t.label for t in copy_list]
                    self.assertEqual(len(tlist), len(copy_list))
                    self.assertEqual(len(tlist), len(tlist._trees))
                    self.assertEqual(len(tlist.taxon_namespace), 7)
                    for idx, (t1, t2, tlabel) in enumerate(zip(tlist, copy_list, expected_tree_labels)):
                        self.assertIn(t1, tlist)
                        self.assertIn(t1, tlist._trees)
                        self.assertEqual(t1.label, tlabel)
                        self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
                        if idx in slice_indexes:
                            self.assertIsNot(t1, t2)
                            self.assertIn(t1, tlist)
                            self.assertIn(t1, tlist._trees)
                            self.assertNotIn(t2, tlist)
                            self.assertNotIn(t2, tlist._trees)
                            self.assertNotIn(t1, source)
                            self.assertNotIn(t1, source._trees)
                        else:
                            self.assertIs(t1, t2)
                        for nd in t1:
                            self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_clear(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        self.assertEqual(len(tlist._trees), 5)
        tlist.clear()
        self.assertEqual(len(tlist), 0)
        self.assertEqual(len(tlist._trees), 0)

    def test_extend_from_another_tree_list_different_namespace(self):
        tlist = curated_test_tree_list.get_tree_list(num_trees=3)
        original_tns = tlist.taxon_namespace
        original_tlist_len = len(tlist)
        original_tree_labels = [t.label for t in tlist]
        self.assertEqual(len(original_tree_labels), len(tlist))
        self.assertEqual(original_tlist_len, 3)

        tlist_source = curated_test_tree_list.get_tree_list(num_trees=5)
        self.assertEqual(len(tlist_source), 5)
        source_tree_labels = [t.label for t in tlist_source]
        self.assertEqual(len(source_tree_labels), len(tlist_source))

        tlist.extend(tlist_source)

        self.assertEqual(len(tlist), original_tlist_len + len(tlist_source))
        self.assertIs(tlist.taxon_namespace, original_tns)
        # self.assertEqual(len(tlist.taxon_namespace), len(tlist[0].tax_labels))
        self.assertEqual(len(tlist.taxon_namespace), 7)
        expected_tree_labels = original_tree_labels + source_tree_labels
        self.assertEqual(len(tlist), len(expected_tree_labels))
        for t1, tlabel in zip(tlist, expected_tree_labels):
            self.assertIn(t1, tlist)
            self.assertNotIn(t1, tlist_source)
            self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
            self.assertEqual(t1.label, tlabel)
            for nd in t1:
                self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_extend_from_list_of_trees_different_namespace(self):
        tlist = curated_test_tree_list.get_tree_list(num_trees=3)
        original_tns = tlist.taxon_namespace
        original_tlist_len = len(tlist)
        original_tree_labels = [t.label for t in tlist]
        self.assertEqual(len(original_tree_labels), len(tlist))
        self.assertEqual(original_tlist_len, 3)
        source_trees = curated_test_tree_list.get_trees(
                num_trees=5,
                taxon_namespace=None,
                label=None,
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.assertEqual(len(source_trees), 5)
        source_tree_labels = [t.label for t in source_trees]
        self.assertEqual(len(source_tree_labels), len(source_trees))

        tlist.extend(source_trees)

        self.assertEqual(len(tlist), original_tlist_len + len(source_trees))
        self.assertIs(tlist.taxon_namespace, original_tns)
        # self.assertEqual(len(tlist.taxon_namespace), len(tlist[0].tax_labels))
        self.assertEqual(len(tlist.taxon_namespace), 7)
        expected_tree_labels = original_tree_labels + source_tree_labels
        self.assertEqual(len(tlist), len(expected_tree_labels))
        for t1, tlabel in zip(tlist, expected_tree_labels):
            self.assertIn(t1, tlist)
            if tlabel in source_tree_labels:
                self.assertIn(t1, source_trees)
            else:
                self.assertNotIn(t1, source_trees)
            self.assertIs(t1.taxon_namespace, tlist.taxon_namespace)
            self.assertEqual(t1.label, tlabel)
            for nd in t1:
                self.assertIn(nd.taxon, tlist.taxon_namespace)

    def test_index(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        for idx, t in enumerate(tlist):
            self.assertIs(t, tlist[idx])
            self.assertEqual(tlist.index(t), idx)

    def test_pop1(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        k = tlist[-1]
        t = tlist.pop()
        self.assertIs(t, k)
        self.assertEqual(len(tlist), 4)
        self.assertNotIn(t, tlist)

    def test_pop2(self):
        for idx in range(5):
            tlist = curated_test_tree_list.get_tree_list(5)
            k = tlist[idx]
            t = tlist.pop(idx)
            self.assertIs(t, k)
            self.assertEqual(len(tlist), 4)
            self.assertNotIn(t, tlist)

    def test_remove(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        t = tlist[0]
        tlist.remove(t)
        self.assertEqual(len(tlist), 4)
        self.assertNotIn(t, tlist)

    def test_remove(self):
        tlist = curated_test_tree_list.get_tree_list(5)
        clist = list(tlist._trees)
        tlist.reverse()
        clist.reverse()
        for t1, t2 in zip(tlist, clist):
            self.assertIs(t1, t2)

    def test_sort(self):
        for r in (True, False):
            tlist = curated_test_tree_list.get_tree_list(5)
            clist = list(tlist._trees)
            tlist.sort(key=lambda x: x.label, reverse=r)
            clist.sort(key=lambda x: x.label, reverse=r)
            for t1, t2 in zip(tlist, clist):
                self.assertIs(t1, t2)

class TreeListCreatingAndCloning(
        curated_test_tree.CuratedTestTree,
        compare_and_validate.Comparator,
        unittest.TestCase):

    def add_tree_annotations(self, tree):
        for idx, nd in enumerate(tree):
            if idx % 2 == 0:
                nd.edge.label = "E{}".format(idx)
                nd.edge.length = idx
            an1 = nd.annotations.add_new("a{}".format(idx),
                    "{}{}{}".format(nd.label, nd.taxon, idx))
            an2 = nd.annotations.add_bound_attribute("label")
            an3 = an1.annotations.add_bound_attribute("name")
            ae1 = nd.edge.annotations.add_new("a{}".format(idx),
                    "{}{}".format(nd.edge.label, idx))
            ae2 = nd.edge.annotations.add_bound_attribute("label")
            ae3 = ae1.annotations.add_bound_attribute("name")
        tree.annotations.add_new("a", 0)
        tree.label = "hello"
        b = tree.annotations.add_bound_attribute("label")
        b.annotations.add_new("c", 3)

    def add_tree_list_annotations(self, tree_list):
        tree_list.annotations.add_new("a", 0)
        tree_list.label = "hello"
        b = tree_list.annotations.add_bound_attribute("label")
        b.annotations.add_new("c", 3)

    def add_taxon_namespace_annotations(self, tns):
        for idx, taxon in enumerate(tns):
            a = taxon.annotations.add_new("!color", str(idx))
            a.annotations.add_new("setbytest", "a")

    def setUp(self):
        self.num_trees = 5
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False)
        self.original_taxon_labels = [t.label for t in tree1.taxon_namespace]
        assert len(self.original_taxon_labels) == len(anodes1)

    def get_tree_list(self):
        tlist1 = dendropy.TreeList()
        self.num_trees = 5
        for idx in range(self.num_trees):
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False,
                    taxon_namespace=tlist1.taxon_namespace)
            self.add_tree_annotations(tree1)
            tlist1.append(tree1)
        self.add_tree_list_annotations(tlist1)
        self.add_taxon_namespace_annotations(tlist1.taxon_namespace)
        return tlist1

    def test_shallow_copy_with_initializer_list(self):
        tlist1 = self.get_tree_list()
        trees = tlist1._trees
        tlist2 = dendropy.TreeList(trees)
        self.assertEqual(len(tlist2), self.num_trees)
        for tcopy, toriginal in zip(tlist2, trees):
            self.assertIs(tcopy, toriginal)
            self.assertIs(tcopy.taxon_namespace, tlist2.taxon_namespace)

    def test_clone0(self):
        tlist1 = self.get_tree_list()
        for tlist2 in (
                tlist1.clone(0),
                ):
            self.assertIs(tlist2.taxon_namespace, tlist1.taxon_namespace)
            self.assertEqual(len(tlist2), self.num_trees)
            for tcopy, toriginal in zip(tlist2, tlist1):
                self.assertIs(tcopy, toriginal)
                self.assertIs(tcopy.taxon_namespace, tlist2.taxon_namespace)

    def test_taxon_namespace_scoped_copy(self):
        tlist1 = self.get_tree_list()
        for tlist2 in (
                tlist1.clone(1),
                dendropy.TreeList(tlist1),
                tlist1.taxon_namespace_scoped_copy(),):
            self.compare_distinct_tree_list(tlist2, tlist1,
                    taxon_namespace_scoped=True,
                    compare_tree_annotations=True,
                    compare_taxon_annotations=True)

    def test_deepcopy_including_namespace(self):
        tlist1 = self.get_tree_list()
        for idx, tlist2 in enumerate((
                tlist1.clone(2),
                copy.deepcopy(tlist1),
                )):
            self.compare_distinct_tree_list(tlist2, tlist1,
                    taxon_namespace_scoped=False,
                    compare_tree_annotations=True,
                    compare_taxon_annotations=True)

    def test_deepcopy_excluding_namespace(self):
        tlist1 = self.get_tree_list()
        tlist2 = dendropy.TreeList(tlist1,
                taxon_namespace=dendropy.TaxonNamespace())
        self.compare_distinct_tree_list(tlist2, tlist1,
                taxon_namespace_scoped=False,
                compare_tree_annotations=True,
                compare_taxon_annotations=False)

class TreeListIdentity(unittest.TestCase):

    def setUp(self):
        self.tns = dendropy.TaxonNamespace()
        self.t1 = dendropy.TreeList(label="a", taxon_namespace=self.tns)
        self.t2 = dendropy.TreeList(label="a", taxon_namespace=self.tns)
        self.t3 = dendropy.TreeList(label="a")

    def test_equal(self):
        # two distinct |TreeList| objects are equal
        # if they have the same namespace and trees
        trees = [dendropy.Tree() for i in range(5)]
        for tree in trees:
            self.t1._trees.append(tree)
            self.t2._trees.append(tree)
        self.assertEqual(self.t1, self.t2)

    def test_unequal1(self):
        # two distinct |TreeList| objects are equal
        # if they have the same namespace and trees
        trees1 = [dendropy.Tree() for i in range(5)]
        for tree in trees1:
            self.t1._trees.append(tree)
        trees2 = [dendropy.Tree() for i in range(5)]
        for tree in trees2:
            self.t2._trees.append(tree)
        self.assertNotEqual(self.t1, self.t2)

    def test_unequal2(self):
        # two distinct |TreeList| objects are equal
        # if they have the same namespace and trees
        trees1 = [dendropy.Tree() for i in range(5)]
        for tree in trees1:
            self.t1._trees.append(tree)
            self.t3._trees.append(tree)
        self.assertNotEqual(self.t1, self.t3)

    def test_hash_dict_membership(self):
        k = {}
        k[self.t1] = 1
        k[self.t2] = 2
        self.assertEqual(len(k), 2)
        self.assertEqual(k[self.t1], 1)
        self.assertEqual(k[self.t2], 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        del k[self.t1]
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

    def test_hash_set_membership(self):
        k = set()
        k.add(self.t1)
        k.add(self.t2)
        self.assertEqual(len(k), 2)
        self.assertIn(self.t1, k)
        self.assertIn(self.t2, k)
        k.discard(self.t1)
        self.assertNotIn(self.t1, k)
        self.assertIn(self.t2, k)
        self.assertEqual(len(k), 1)
        k1 = {self.t1: 1}
        k2 = {self.t2: 1}
        self.assertIn(self.t1, k1)
        self.assertIn(self.t2, k2)
        self.assertNotIn(self.t2, k1)
        self.assertNotIn(self.t1, k2)

class TestTreeListUpdateTaxonNamespace(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        trees = []
        for idx in range(5):
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    suppress_internal_node_taxa=True,
                    suppress_leaf_node_taxa=True)
            trees.append(tree1)
        self.expected_labels = set()
        self.expected_taxa = set()
        node_label_to_taxon_label_map = {
            "a" : "z01",
            "b" : "<NONE>",
            "c" : "z03",
            "e" : "z04",
            "f" : "z05",
            "g" : "z06",
            "h" : None,
            "i" : None,
            "j" : "z09",
            "k" : "z10",
            "l" : "z11",
            "m" : "<NONE>",
            "n" : None,
            "o" : "z14",
            "p" : "z15",
                }
        registry = {}
        for tree_idx, tree in enumerate(trees):
            for nd in tree:
                if nd.label is not None:
                    if tree_idx > 3:
                        nd.label = node_label_to_taxon_label_map[nd.label]
                    if nd.label == "<NONE>":
                        try:
                            t = registry[None]
                        except KeyError:
                            t = dendropy.Taxon(label=None)
                            registry[None] = t
                        self.expected_labels.add(None)
                    else:
                        try:
                            t = registry[nd.label]
                        except KeyError:
                            t = dendropy.Taxon(label=nd.label)
                            registry[nd.label] = t
                        self.expected_labels.add(nd.label)
                    nd.taxon = t
                    self.expected_taxa.add(nd.taxon)
        self.tree_list = dendropy.TreeList()
        self.tree_list._trees = trees

    def test_noop_update_with_no_taxa(self):
        trees = []
        tns = dendropy.TaxonNamespace()
        for idx in range(5):
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    suppress_internal_node_taxa=True,
                    suppress_leaf_node_taxa=True,
                    taxon_namespace=tns)
            trees.append(tree1)
        tlst = dendropy.TreeList(taxon_namespace=tns)
        tlst._trees = trees
        original_tns = tlst.taxon_namespace
        self.assertEqual(len(original_tns), 0)
        tlst.update_taxon_namespace()
        self.assertIs(tlst.taxon_namespace, original_tns)
        for tree in tlst:
            self.assertIs(tree.taxon_namespace, tlst.taxon_namespace)
        self.assertEqual(len(original_tns), 0)

    def test_update(self):
        original_tns = self.tree_list.taxon_namespace
        self.assertEqual(len(original_tns), 0)
        self.tree_list.update_taxon_namespace()
        self.tree_list.update_taxon_namespace()
        self.tree_list.update_taxon_namespace()
        for tree in self.tree_list:
            self.assertIs(tree.taxon_namespace, self.tree_list.taxon_namespace)
        self.assertIs(self.tree_list.taxon_namespace, original_tns)
        new_taxa = [t for t in original_tns]
        new_labels = [t.label for t in original_tns]
        self.assertCountEqual(new_taxa, self.expected_taxa)
        self.assertCountEqual(new_labels, self.expected_labels)

class TestTreeListMigrateAndReconstructTaxonNamespace(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        tns = dendropy.TaxonNamespace()
        trees = []
        for idx in range(8):
            tree, anodes, lnodes, inodes = self.get_tree(
                    suppress_internal_node_taxa=True,
                    suppress_leaf_node_taxa=True,
                    taxon_namespace=tns)
            trees.append(tree)
        self.node_label_to_taxon_label_map = {
            "a" : "a",
            "b" : "a",
            "c" : "2",
            "e" : "2",
            "f" : "b",
            "g" : "B",
            "h" : "B",
            "i" : "h",
            "j" : "H",
            "k" : "h",
            "l" : None,
            "m" : None,
            "n" : "H",
            "o" : "J",
            "p" : "j",
                }
        self.original_taxa = []
        registry = {}
        for tree in trees:
            for idx, nd in enumerate(tree):
                try:
                    t = registry[nd.label]
                except KeyError:
                    taxon_label = self.node_label_to_taxon_label_map[nd.label]
                    t = dendropy.Taxon(label=taxon_label)
                    registry[nd.label] = t
                    self.original_taxa.append(t)
                tree.taxon_namespace.add_taxon(t)
                nd.taxon = t
                nd.original_taxon = t
        assert len(tree.taxon_namespace) == len(self.node_label_to_taxon_label_map)
        assert len(tree.taxon_namespace) == len(self.original_taxa)
        self.tree_list = dendropy.TreeList(taxon_namespace=tns)
        self.tree_list._trees = trees

    def verify_taxon_namespace_reconstruction(self,
            unify_taxa_by_label=False,
            case_sensitive_label_mapping=True,
            original_tns=None,
            redundant_taxa=False):
        if unify_taxa_by_label:
            if not case_sensitive_label_mapping:
                expected_labels = []
                for label in self.node_label_to_taxon_label_map.values():
                    if label is None:
                        expected_labels.append(label)
                    else:
                        label = label.upper()
                        if label not in expected_labels:
                            expected_labels.append(label)
            else:
                expected_labels = list(set(label for label in self.node_label_to_taxon_label_map.values()))
        else:
            expected_labels = [label for label in self.node_label_to_taxon_label_map.values()]
        for tree in self.tree_list:
            seen_taxa = []
            self.assertIs(tree.taxon_namespace, self.tree_list.taxon_namespace)
            for nd in tree:
                self.assertIsNot(nd.taxon, nd.original_taxon)
                if not case_sensitive_label_mapping and nd.taxon.label is not None:
                    self.assertEqual(nd.taxon.label.upper(), nd.original_taxon.label.upper())
                    self.assertEqual(self.node_label_to_taxon_label_map[nd.label].upper(), nd.taxon.label.upper())
                else:
                    self.assertEqual(nd.taxon.label, nd.original_taxon.label)
                    self.assertEqual(self.node_label_to_taxon_label_map[nd.label], nd.taxon.label)
                self.assertNotIn(nd.original_taxon, tree.taxon_namespace)
                self.assertIn(nd.original_taxon, self.original_taxa)
                self.assertIn(nd.taxon, tree.taxon_namespace)
                self.assertNotIn(nd.taxon, self.original_taxa)
                if original_tns is not None:
                    self.assertNotIn(nd.taxon, original_tns)
                if nd.taxon not in seen_taxa:
                    seen_taxa.append(nd.taxon)
                else:
                    self.assertTrue(unify_taxa_by_label or redundant_taxa)
                    if not case_sensitive_label_mapping:
                        self.assertIn(nd.taxon.label, [t.label for t in seen_taxa])
                    else:
                        if nd.taxon.label is None:
                            self.assertIs(nd.original_taxon.label, None)
                            self.assertEqual([t.label for t in seen_taxa].count(None), 1)
                        else:
                            x1 = [t.label.upper() for t in seen_taxa if t.label is not None]
                            self.assertIn(nd.taxon.label.upper(), x1)
            self.assertEqual(len(seen_taxa), len(tree.taxon_namespace))
            if not case_sensitive_label_mapping:
                seen_labels = [(t.label.upper() if t.label is not None else None) for t in seen_taxa]
            else:
                seen_labels = [t.label for t in seen_taxa]
            c1 = collections.Counter(expected_labels)
            c2 = collections.Counter(seen_labels)
            self.assertEqual(c1, c2)

        self.assertEqual(len(tree.taxon_namespace), len(expected_labels))
        if not unify_taxa_by_label and not redundant_taxa:
            self.assertEqual(len(tree.taxon_namespace), len(self.node_label_to_taxon_label_map))

    def test_basic_reconstruction(self):
        tns = dendropy.TaxonNamespace()
        trees = []
        for idx in range(5):
            tree, anodes, lnodes, inodes = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False,
                    taxon_namespace=tns)
            trees.append(tree)
        tree_list = dendropy.TreeList(taxon_namespace=tns)
        tree_list._trees = trees
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        tree_list.taxon_namespace = new_tns
        tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=False)
        self.assertIsNot(tree_list.taxon_namespace, tns)
        self.assertIs(tree_list.taxon_namespace, new_tns)
        self.assertEqual(len(tree_list.taxon_namespace), len(tns))
        original_labels = [t.label for t in tns]
        new_labels = [t.label for t in new_tns]
        self.assertCountEqual(new_labels, original_labels)
        for tree in tree_list:
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            for nd in tree:
                if nd.taxon is not None:
                    self.assertIn(nd.taxon, tree.taxon_namespace)
                    self.assertNotIn(nd.taxon, tns)

    def test_reconstruct_taxon_namespace_non_unifying(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree_list._taxon_namespace = new_tns
        self.assertEqual(len(self.tree_list.taxon_namespace), 0)
        self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=False)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=False,
                case_sensitive_label_mapping=True)

    def test_reconstruct_taxon_namespace_unifying_case_sensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree_list._taxon_namespace = new_tns
        self.assertEqual(len(self.tree_list.taxon_namespace), 0)
        self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=True)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=True,
                original_tns=original_tns)

    def test_reconstruct_taxon_namespace_unifying_case_insensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = False
        self.tree_list._taxon_namespace = new_tns
        self.assertEqual(len(self.tree_list.taxon_namespace), 0)
        self.tree_list.reconstruct_taxon_namespace(unify_taxa_by_label=True)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=False,
                original_tns=original_tns)

    def test_basic_migration(self):
        tns = dendropy.TaxonNamespace()
        trees = []
        for idx in range(5):
            tree, anodes, lnodes, inodes = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False,
                    taxon_namespace=tns)
            trees.append(tree)
        tree_list = dendropy.TreeList(taxon_namespace=tns)
        tree_list._trees = trees
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        tree_list.taxon_namespace = new_tns
        tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=False)
        self.assertIsNot(tree_list.taxon_namespace, tns)
        self.assertIs(tree_list.taxon_namespace, new_tns)
        self.assertEqual(len(tree_list.taxon_namespace), len(tns))
        original_labels = [t.label for t in tns]
        new_labels = [t.label for t in new_tns]
        self.assertCountEqual(new_labels, original_labels)
        for tree in tree_list:
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            for nd in tree:
                if nd.taxon is not None:
                    self.assertIn(nd.taxon, tree.taxon_namespace)
                    self.assertNotIn(nd.taxon, tns)

    def test_migrate_taxon_namespace_non_unifying(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=False)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=False,
                case_sensitive_label_mapping=True,
                original_tns=original_tns)

    def test_migrate_taxon_namespace_unifying_case_sensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = True
        self.tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=True)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=True,
                original_tns=original_tns)

    def test_migrate_taxon_namespace_unifying_case_insensitive(self):
        original_tns = self.tree_list.taxon_namespace
        new_tns = dendropy.TaxonNamespace()
        new_tns.is_case_sensitive = False
        self.tree_list.migrate_taxon_namespace(
                new_tns,
                unify_taxa_by_label=True)
        self.assertIsNot(self.tree_list.taxon_namespace, original_tns)
        self.assertIs(self.tree_list.taxon_namespace, new_tns)
        self.verify_taxon_namespace_reconstruction(
                unify_taxa_by_label=True,
                case_sensitive_label_mapping=False,
                original_tns=original_tns)

class TestTreeListAppend(
        curated_test_tree.CuratedTestTree,
        unittest.TestCase):

    def setUp(self):
        self.native_tns = dendropy.TaxonNamespace()
        self.tree_list = dendropy.TreeList(taxon_namespace=self.native_tns)
        self.foreign_tns = dendropy.TaxonNamespace()
        self.foreign_tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False,
                taxon_namespace=self.foreign_tns)
        for nd in self.foreign_tree:
            nd.original_taxon = nd.taxon
        self.check_tns = dendropy.TaxonNamespace()
        self.check_tree, anodes, lnodes, inodes = self.get_tree(
                suppress_internal_node_taxa=False,
                suppress_leaf_node_taxa=False,
                taxon_namespace=self.check_tns)

    def test_append_default(self):
        self.assertIsNot(self.tree_list.taxon_namespace, self.foreign_tree.taxon_namespace)
        self.tree_list.append(self.foreign_tree)
        self.assertEqual(len(self.tree_list), 1)
        self.assertIn(self.foreign_tree, self.tree_list)
        self.assertIs(self.foreign_tree, self.tree_list[0])
        self.assertIs(self.tree_list.taxon_namespace, self.native_tns)
        self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
        self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
        for nd in self.foreign_tree:
            if nd.taxon:
                self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
                self.assertIsNot(nd.taxon, nd.original_taxon)
                self.assertIn(nd.original_taxon, self.foreign_tns)
                self.assertNotIn(nd.original_taxon, self.tree_list.taxon_namespace)
                self.assertEqual(nd.taxon.label, nd.original_taxon.label)

    def test_append_migrate_matching_labels(self):
        kwargs_groups = [
                {"taxon_import_strategy": "migrate", "unify_taxa_by_label": True},
                {"taxon_import_strategy": "migrate", "unify_taxa_by_label": False},
                {"taxon_import_strategy": "add", },
        ]
        for kwargs in kwargs_groups:
            self.setUp()
            self.assertEqual(len(self.tree_list.taxon_namespace), 0)
            native_tree, anodes, lnodes, inodes = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False,
                    taxon_namespace=self.native_tns)
            self.assertEqual(len(self.tree_list.taxon_namespace), len(self.postorder_sequence))
            self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
            original_tns_len = len(self.tree_list.taxon_namespace)
            self.tree_list.append(self.foreign_tree, **kwargs)
            self.assertEqual(len(self.tree_list), 1)
            self.assertIn(self.foreign_tree, self.tree_list)
            self.assertIs(self.foreign_tree, self.tree_list[0])
            self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
            if kwargs["taxon_import_strategy"] == "add":
                self.assertEqual(len(self.tree_list.taxon_namespace),
                        original_tns_len + len(self.foreign_tns))
                for nd in self.foreign_tree:
                    self.assertIn(nd.taxon, self.foreign_tns)
                    self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
            else:
                if "unify_taxa_by_label" not in kwargs or not kwargs["unify_taxa_by_label"]:
                    self.assertEqual(len(self.tree_list.taxon_namespace),
                            original_tns_len + len(self.foreign_tns))
                else:
                    self.assertEqual(len(self.tree_list.taxon_namespace), original_tns_len)
                for nd in self.foreign_tree:
                    self.assertNotIn(nd.taxon, self.foreign_tns)
                    self.assertIn(nd.taxon, self.tree_list.taxon_namespace)

    def test_append_add(self):
        self.assertIsNot(self.tree_list.taxon_namespace, self.foreign_tree.taxon_namespace)
        self.tree_list.append(self.foreign_tree,
                taxon_import_strategy="add")
        self.assertEqual(len(self.tree_list), 1)
        self.assertIn(self.foreign_tree, self.tree_list)
        self.assertIs(self.foreign_tree, self.tree_list[0])
        self.assertIs(self.tree_list.taxon_namespace, self.native_tns)
        self.assertIs(self.foreign_tree.taxon_namespace, self.tree_list.taxon_namespace)
        self.assertEqual(len(self.tree_list.taxon_namespace), len(self.foreign_tns))
        for nd in self.foreign_tree:
            if nd.taxon:
                self.assertIn(nd.taxon, self.tree_list.taxon_namespace)
                self.assertIs(nd.taxon, nd.original_taxon)
                self.assertIn(nd.original_taxon, self.foreign_tns)
                self.assertIn(nd.original_taxon, self.tree_list.taxon_namespace)

class TestTreeListTaxa(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.tree_list = dendropy.TreeList()
        self.expected_taxa = None
        for i in range(10):
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    taxon_namespace=self.tree_list.taxon_namespace,
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False)
            self.tree_list.append(tree1)
            if self.expected_taxa is None:
                self.expected_taxa = set([nd.taxon for nd in anodes1 if nd.taxon is not None])

    def test_basic_taxa(self):
        self.assertEqual(self.tree_list.poll_taxa(), self.expected_taxa)

class TestTreeListPurgeTaxonNamespace(
        curated_test_tree.CuratedTestTree,
        dendropytest.ExtendedTestCase):

    def setUp(self):
        self.tree_list = dendropy.TreeList()
        self.expected_taxa = None
        for i in range(10):
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    taxon_namespace=self.tree_list.taxon_namespace,
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False)
            self.tree_list.append(tree1)
            if self.expected_taxa is None:
                self.expected_taxa = set([nd.taxon for nd in anodes1 if nd.taxon is not None])

    def test_noop_purge(self):
        self.assertEqual(set(self.tree_list.taxon_namespace), self.expected_taxa)
        self.tree_list.purge_taxon_namespace()
        self.assertEqual(set(self.tree_list.taxon_namespace), self.expected_taxa)

    def test_basic_purge(self):
        self.assertEqual(set(self.tree_list.taxon_namespace), self.expected_taxa)
        added_taxa = set(self.expected_taxa)
        for label in ("z1", "z2", "z3", "z4"):
            t = self.tree_list.taxon_namespace.new_taxon(label=label)
            added_taxa.add(t)
        self.assertEqual(set(self.tree_list.taxon_namespace), added_taxa)
        self.tree_list.purge_taxon_namespace()
        self.assertEqual(set(self.tree_list.taxon_namespace), self.expected_taxa)

class TreeListCreation(unittest.TestCase):

    def test_create_with_taxon_namespace(self):
        tns = dendropy.TaxonNamespace()
        tt = dendropy.TreeList(label="a", taxon_namespace=tns)
        self.assertEqual(tt.label, "a")
        self.assertIs(tt.taxon_namespace, tns)

class TestSpecialTreeListConstruction(
        unittest.TestCase):

    def test_construction_from_another_tree_different_label(self):
        tlist1 = dendropy.TreeList()
        tlist1.label = "tlist1"
        self.assertEqual(tlist1.label, "tlist1")
        tlist2 = dendropy.TreeList(tlist1, label="tlist2")
        self.assertEqual(tlist2.label, "tlist2")
        self.assertNotEqual(tlist1.label, "tlist2")
        self.assertNotEqual(tlist1.label, tlist2.label)

if __name__ == "__main__":
    unittest.main()
