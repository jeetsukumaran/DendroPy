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
Tests native tree structuring routines.
"""

import sys
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
import unittest
from dendropy.test.support import curated_test_tree
from dendropy.test.support import compare_and_validate
from dendropy.test.support import pathmap
from dendropy.test.support import dendropytest
from dendropy.utility import messaging
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.test.support.mockrandom import MockRandom
import dendropy
from dendropy.calculate import treecompare
from dendropy.utility.textprocessing import StringIO
import re

_LOG = messaging.get_logger(__name__)

class ScaleTest(unittest.TestCase):

    def testScaleEdgesNoLens(self):
        newick_list = ['(5,((4,3),2),1);',
            '(5,(4,3,2),1);',
            '(5,((4,3),2),1);',
            '(5,(4,3),2,1);',
            '(5,((4,3),2),1);',
            '(5,4,3,2,1);']
        tree_list = dendropy.TreeList.get_from_stream(
                        StringIO("""%s""" % "\n".join(newick_list)),
                        schema="newick")
        for n, tree in enumerate(tree_list):
            tree.scale_edges(2.0)
            self.assertEqual(newick_list[n], "%s;" % tree._as_newick_string())
    def testScaleEdgesRealTest(self):
        newick_list = ['(5:3,((4:1,3:1):1.5,2:3),1:0);',
            '(5:7.5,4:1,3:-2,2:4,1:.1);']
        doubled = ['(5:6.0,((4:2.0,3:2.0):3.0,2:6.0),1:0.0);',
                    '(5:15.0,4:2.0,3:-4.0,2:8.0,1:0.2);']
        as_f = ['(5:3.0,((4:1.0,3:1.0):1.5,2:3.0),1:0.0);',
            '(5:7.5,4:1.0,3:-2.0,2:4.0,1:0.1);']
        tree_list = dendropy.TreeList.get_from_stream(
                        StringIO("""%s""" % "\n".join(newick_list)),
                        schema="newick")
        for n, tree in enumerate(tree_list):
            tree.scale_edges(2)
            self.assertEqual(doubled[n], "%s;" % tree._as_newick_string())
        for n, tree in enumerate(tree_list):
            tree.scale_edges(.5)
            self.assertEqual(as_f[n], "%s;" % tree._as_newick_string())

class RandomlyRotateTest(unittest.TestCase):

    def runTest(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        trees = dendropy.TreeList.get_from_stream(StringIO(n+n), schema="newick")
        ref = trees[0]
        changing = trees[1]
        rng = MockRandom()
        ref.encode_bipartitions()
        changing.encode_bipartitions()
        orig_root = changing.seed_node
        for i in range(50):
            changing.randomly_rotate(rng=rng)
            self.assertNotEqual(str(changing), n)
            self.assertEqual(orig_root, changing.seed_node)
            changing._debug_check_tree(logger_obj=_LOG, check_bipartitions=True)
            if treecompare.symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

class RandomlyReorientTest(unittest.TestCase):

    def runTest(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        k = dendropy.TreeList.get_from_stream(StringIO(n), schema="newick")[0]
        trees = dendropy.TreeList.get_from_stream(StringIO(n+n), schema="newick", taxon_namespace=k.taxon_namespace)
        ref = trees[0]
        changing = trees[1]
        rng = MockRandom()
        for i in range(50):
            changing.randomly_reorient(rng=rng, update_bipartitions=True)
            self.assertNotEqual(str(changing), n)
            changing._debug_check_tree(logger_obj=_LOG, check_bipartitions=True)
            d = treecompare.symmetric_difference(ref, changing, is_bipartitions_updated=False)
            if d != 0:
                self.fail("\n{}\n!=\n{}\nRF={}".format(str(ref), str(changing), d))

class CollapseConflictingTest(unittest.TestCase):

    def runTest(self):

        taxon_namespace = dendropy.TaxonNamespace([str(i+1) for i in range(5)])
        tree_list = dendropy.TreeList.get_from_stream(
            StringIO("""
            (5,((4,3),2),1);
            (5,(4,3,2),1);
            (5,((4,3),2),1);
            (5,(4,3),2,1);
            (5,((4,3),2),1);
            (5,4,3,2,1);
            """),
            schema="newick",
            taxon_namespace=taxon_namespace)

        tree = tree_list[0]
        expected_tree = tree_list[1]
        tree.encode_bipartitions()
        tree_leafset_bitmask = tree.seed_node.edge.bipartition._leafset_bitmask
        bipartition_to_target = dendropy.Bipartition(
                bitmask=0xA,
                tree_leafset_bitmask=tree_leafset_bitmask,
                compule_bitmasks=True)
        assert bipartition_to_target._lowest_relevant_bit is not None
        tree.seed_node.collapse_conflicting(bipartition_to_target)
        tree.encode_bipartitions()
        expected_tree.encode_bipartitions()
        self.assertEqual(treecompare.symmetric_difference(tree, expected_tree), 0)

        tree = tree_list[2]
        expected_tree = tree_list[3]
        tree.encode_bipartitions()
        tree_leafset_bitmask = tree.seed_node.edge.bipartition._leafset_bitmask
        bipartition_to_target = dendropy.Bipartition(bitmask=0x3,
                tree_leafset_bitmask=tree_leafset_bitmask,
                compile_bipartition=True)
        tree.seed_node.collapse_conflicting(bipartition_to_target)
        tree.encode_bipartitions()
        expected_tree.encode_bipartitions()
        self.assertEqual(treecompare.symmetric_difference(tree, expected_tree), 0)

        tree = tree_list[4]
        expected_tree = tree_list[5]
        tree.encode_bipartitions()
        tree_leafset_bitmask = tree.seed_node.edge.bipartition._leafset_bitmask
        bipartition_to_target = dendropy.Bipartition(bitmask=0x5,
                tree_leafset_bitmask=tree_leafset_bitmask,
                compile_bipartition=True)
        tree.seed_node.collapse_conflicting(bipartition_to_target)
        tree.encode_bipartitions()
        expected_tree.encode_bipartitions()
        self.assertEqual(treecompare.symmetric_difference(tree, expected_tree), 0)

class PruneTest(unittest.TestCase):

    def check(self,
            title,
            src_prefix,
            to_retain=False):
        input_ds = dendropy.DataSet.get_from_path(
                src=pathmap.tree_source_path(src_prefix + ".pre-pruned.nex"),
                schema='nexus')
        tns1 = dendropy.TaxonNamespace()
        input_ds.attach_taxon_namespace(tns1)
        input_taxa = input_ds.taxon_namespaces[0]
        output_ds = dendropy.DataSet.get_from_path(
                src=pathmap.tree_source_path(src_prefix + ".paup-pruned.nex"),
                schema='nexus',
                taxon_namespace=input_taxa)
        tns2 = dendropy.TaxonNamespace()
        output_ds.attach_taxon_namespace(tns2)
        if to_retain:
            taxf = open(pathmap.tree_source_path(src_prefix + ".retained_taxa.txt"), "r")
        else:
            taxf = open(pathmap.tree_source_path(src_prefix + ".pruned_taxa.txt"), "r")
        rows = taxf.readlines()
        taxon_idxs_list = [ [int(i) for i in row.split()] for row in rows ]
        for set_idx, src_trees in enumerate(input_ds.tree_lists):
            src_trees = input_ds.tree_lists[set_idx]
            ref_trees = output_ds.tree_lists[set_idx]
            taxon_idxs = taxon_idxs_list[set_idx]
            sub_taxa = [src_trees.taxon_namespace[i] for i in taxon_idxs]
            for tree_idx, src_tree in enumerate(src_trees):
                _LOG.debug("%s Set %d/%d, Tree %d/%d" % (title, set_idx+1, len(input_ds.tree_lists), tree_idx+1, len(src_trees)))
                ref_tree = ref_trees[tree_idx]
                if to_retain:
                    src_tree.retain_taxa(sub_taxa)
                else:
                    src_tree.prune_taxa(sub_taxa)
                # tree_dist = paup.symmetric_difference(src_tree, ref_tree)
                self.assertEqual(treecompare.symmetric_difference(src_tree, ref_tree), 0)
        taxf.close()

    def testPruneTaxaUnrooted(self):
        self.check("Unrooted", "prune_unrooted", False)

    def testPruneTaxaRooted(self):
        self.check("Rooted", "prune_rooted", False)

    def testRetainTaxaUnrooted(self):
        self.check("Unrooted", "prune_unrooted", True)

    def testRetainTaxaRooted(self):
        self.check("Rooted", "prune_rooted", True)

class TruncateTree(unittest.TestCase):

    def setUp(self):
        self.trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path("pythonidae.reference-trees.nexus"), "nexus")

    def check_ultrametric_tree(self, tree, dist):
        self.assertTrue(tree._debug_tree_is_valid())
        tree.calc_node_root_distances()
        for nd in tree.leaf_node_iter():
            self.assertAlmostEqual(nd.root_distance, nd.distance_from_root())
            self.assertAlmostEqual(dist, nd.root_distance)

    def test_truncate_ultrametric(self):
        for tree in self.trees:
            dists = tree.calc_node_root_distances()
            min_dist, max_dist = tree.minmax_leaf_distance_from_root()
            trunc_dists = [(max_dist * f) for f in (0.25, 0.5, 0.75, 0.90)]
            for td in trunc_dists:
                working = dendropy.Tree(tree)
                working.truncate_from_root(td)
                for idx, leaf in enumerate(working.leaf_node_iter()):
                    if leaf.label is None and leaf.taxon is None:
                        leaf.taxon = dendropy.Taxon(label="t%s" % (idx+1))
                self.check_ultrametric_tree(working, td)

class TestTreeLadderization(unittest.TestCase):

    def setUp(self):
        self.tree_str = "[&R] ((A, (B, (C, (D, E)))),(F, (G, H)));"

    def clean_newick_str(self, s):
        """
        Strips out everything but the core tree statement characters from a
        NEWICK string.
        """
        return re.sub(r'[^()A-H,]', '', s)

    def testLadderizeLeft(self):
        tree = dendropy.Tree.get_from_string(self.tree_str, "newick")
        tree.ladderize(ascending=True)
        self.assertEqual(self.clean_newick_str(tree.as_string("newick")),
                self.clean_newick_str("[&R] ((F,(G,H)),(A,(B,(C,(D,E)))));"))

    def testLadderizeRight(self):
        tree = dendropy.Tree.get_from_string(self.tree_str, "newick")
        tree.ladderize(ascending=False)
        self.assertEqual(self.clean_newick_str(tree.as_string("newick")),
               self.clean_newick_str("[&R] (((((D,E),C),B),A),((G,H),F));"))

class TreeMidpointRootingTest(ExtendedTestCase):

    def testMidpointRooting(self):
        taxa = dendropy.TaxonNamespace()
        test_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path('pythonidae.random.bd0301.randomly-rooted.tre'),
                "nexus",
                taxon_namespace=taxa,
                rooting="force-rooted")
        expected_trees = dendropy.TreeList.get_from_path(pathmap.tree_source_path('pythonidae.random.bd0301.midpoint-rooted.tre'),
                "nexus",
                taxon_namespace=taxa,
                rooting="force-rooted")
        for idx, test_tree in enumerate(test_trees):
            expected_tree = expected_trees[idx]
            test_tree.reroot_at_midpoint(update_bipartitions=True)
            self.assertEqual(treecompare.symmetric_difference(test_tree, expected_tree), 0)
            for bipartition in test_tree.bipartition_encoding:
                if test_tree.bipartition_edge_map[bipartition].head_node is test_tree.seed_node:
                    continue
                # self.assertAlmostEqual(bipartition.edge.length, expected_tree.split_bitmask_edge_map[bipartition.split_bitmask].length, 3)
                self.assertAlmostEqual(test_tree.bipartition_edge_map[bipartition].length,
                        expected_tree.bipartition_edge_map[bipartition].length,
                        3)

class TreeRerootingTests(dendropytest.ExtendedTestCase):
    #                  a
    #                 / \
    #                /   \
    #               /     \
    #              /       \
    #             /         \
    #            /           \
    #           /             c
    #          b             / \
    #         / \           /   \
    #        /   e         /     f
    #       /   / \       /     / \
    #      /   /   \     g     /   h
    #     /   /     \   / \   /   / \
    #    i    j     k  l  m  n   o   p

    # change in orientation == exchange in branch length
    # removing of outdegreee one
    # infomation is as follows
    # 'reseed' : { 'x1' : ('z1', 'z2' ..), 'x2': ('z1',) }
    #    - reseed: the new root
    #    - x1: the edge which has its length changed
    #    - z1, z2, z3 ... : the edges for which lengths, when summed, give the new length for x1
    reseeded_with_unifurcations_suppressed_edge_length_updates = {
            'a': {},
            'b': {'c': ('c', 'b'), },
            'c': {'b': ('c', 'b'), },
            'e': {'b': ('e'), 'c': ('b', 'c'), },
            'f': {'c': ('f'), 'b': ('c', 'b'), },
            'g': {'c': ('g'), 'b': ('c', 'b'), },
            'h': {'f': ('h'), 'c': ('f'), 'b': ('c', 'b'), },
            'i': {'b': ('i'), 'c': ('b', 'c'), },
            'j': {'e': ('j'), 'b': ('e'), 'c':('b','c'), },
            'k': {'e': ('k'), 'b': ('e'), 'c':('b','c'), },
            'l': {'g': ('l'), 'c': ('g'), 'b':('c','b'), },
            'm': {'g': ('m'), 'c': ('g'), 'b':('c','b'), },
            'n': {'f': ('n'), 'c': ('f'), 'b':('c','b'), },
            'o': {'h': ('o'), 'f': ('h'), 'c': ('f'), 'b':('c','b'), },
            'p': {'h': ('p'), 'f': ('h'), 'c': ('f'), 'b':('c','b'), },
    }
    reseeded_with_unifurcations_preserved_edge_length_updates = {
            'a': {},
            'b': {'a':('b'), },
            'c': {'a':('c'), },
            'e': {'b': ('e'), 'a': ('b'), 'c': ('c'),},
            'f': {'c': ('f'), 'a': ('c'), },
            'g': {'c': ('g'), 'a': ('c'), },
            'h': {'f': ('h'), 'a': ('c'), 'c':('f'), },
            'i': {'b': ('i'), 'a': ('b'), },
            'j': {'e': ('j'), 'a': ('b'), 'b': ('e'), },
            'k': {'e': ('k'), 'a': ('b'), 'b': ('e'), },
            'l': {'g': ('l'), 'c': ('g'), 'a':('c'), },
            'm': {'g': ('m'), 'c': ('g'), 'a':('c'), },
            'n': {'f': ('n'), 'c': ('f'), 'a':('c'), },
            'o': {'h': ('o'), 'f': ('h'), 'c': ('f'), 'a':('c'), },
            'p': {'h': ('p'), 'f': ('h'), 'c': ('f'), 'a':('c'), },
    }

    reseeded_with_unifurcations_suppressed_preorder_labels = {
            'a': ['a', 'b', 'i', 'e', 'j', 'k', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'b': ['b', 'i', 'e', 'j', 'k', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'c': ['c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p', 'b', 'i', 'e', 'j', 'k'],
            'e': ['e', 'j', 'k', 'b', 'i', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'f': ['f', 'n', 'h', 'o', 'p', 'c', 'g', 'l', 'm', 'b', 'i', 'e', 'j', 'k'],
            'g': ['g', 'l', 'm', 'c', 'f', 'n', 'h', 'o', 'p', 'b', 'i', 'e', 'j', 'k'],
            'h': ['h', 'o', 'p', 'f', 'n', 'c', 'g', 'l', 'm', 'b', 'i', 'e', 'j', 'k'],
            'i': ['i', 'e', 'j', 'k', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'j': ['j', 'k', 'b', 'i', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'k': ['k', 'j', 'b', 'i', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'l': ['l', 'm', 'c', 'f', 'n', 'h', 'o', 'p', 'b', 'i', 'e', 'j', 'k'],
            'm': ['m', 'l', 'c', 'f', 'n', 'h', 'o', 'p', 'b', 'i', 'e', 'j', 'k'],
            'n': ['n', 'h', 'o', 'p', 'c', 'g', 'l', 'm', 'b', 'i', 'e', 'j', 'k'],
            'o': ['o', 'p', 'f', 'n', 'c', 'g', 'l', 'm', 'b', 'i', 'e', 'j', 'k'],
            'p': ['p', 'o', 'f', 'n', 'c', 'g', 'l', 'm', 'b', 'i', 'e', 'j', 'k']
    }
    reseeded_with_unifurcations_preserved_preorder_labels = {
            'a': ['a', 'b', 'i', 'e', 'j', 'k', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'b': ['b', 'i', 'e', 'j', 'k', 'a', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'c': ['c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p', 'a', 'b', 'i', 'e', 'j', 'k'],
            'e': ['e', 'j', 'k', 'b', 'i', 'a', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'f': ['f', 'n', 'h', 'o', 'p', 'c', 'g', 'l', 'm', 'a', 'b', 'i', 'e', 'j', 'k'],
            'g': ['g', 'l', 'm', 'c', 'f', 'n', 'h', 'o', 'p', 'a', 'b', 'i', 'e', 'j', 'k'],
            'h': ['h', 'o', 'p', 'f', 'n', 'c', 'g', 'l', 'm', 'a', 'b', 'i', 'e', 'j', 'k'],
            'i': ['i', 'b', 'e', 'j', 'k', 'a', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'j': ['j', 'e', 'k', 'b', 'i', 'a', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'k': ['k', 'e', 'j', 'b', 'i', 'a', 'c', 'g', 'l', 'm', 'f', 'n', 'h', 'o', 'p'],
            'l': ['l', 'g', 'm', 'c', 'f', 'n', 'h', 'o', 'p', 'a', 'b', 'i', 'e', 'j', 'k'],
            'm': ['m', 'g', 'l', 'c', 'f', 'n', 'h', 'o', 'p', 'a', 'b', 'i', 'e', 'j', 'k'],
            'n': ['n', 'f', 'h', 'o', 'p', 'c', 'g', 'l', 'm', 'a', 'b', 'i', 'e', 'j', 'k'],
            'o': ['o', 'h', 'p', 'f', 'n', 'c', 'g', 'l', 'm', 'a', 'b', 'i', 'e', 'j', 'k'],
            'p': ['p', 'h', 'o', 'f', 'n', 'c', 'g', 'l', 'm', 'a', 'b', 'i', 'e', 'j', 'k']
    }
    reseeded_with_unifurcations_suppressed_relations = {
        'a': {'a': ('b', 'c'),
            'b': ('i', 'e'),
            'c': ('g', 'f'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'b': {'b': ('i', 'e', 'c'),
            'c': ('g', 'f'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'c': {'b': ('i', 'e'),
            'c': ('g', 'f', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'e': {'b': ('i', 'c'),
            'c': ('g', 'f'),
            'e': ('j', 'k', 'b'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'f': {'b': ('i', 'e'),
            'c': ('g', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'h', 'c'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'g': {'b': ('i', 'e'),
            'c': ('f', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm', 'c'),
            'h': ('o', 'p')},
        'h': {'b': ('i', 'e'),
            'c': ('g', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'c'),
            'g': ('l', 'm'),
            'h': ('o', 'p', 'f')},
        'i': {'c': ('g', 'f'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'i': ('e', 'c')},
        'j': {'b': ('i', 'c'),
            'c': ('g', 'f'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'j': ('k', 'b')},
        'k': {'b': ('i', 'c'),
            'c': ('g', 'f'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'k': ('j', 'b')},
        'l': {'b': ('i', 'e'),
            'c': ('f', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'h': ('o', 'p'),
            'l': ('m', 'c')},
        'm': {'b': ('i', 'e'),
            'c': ('f', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'h': ('o', 'p'),
            'm': ('l', 'c')},
        'n': {'b': ('i', 'e'),
            'c': ('g', 'b'),
            'e': ('j', 'k'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'n': ('h', 'c')},
        'o': {'b': ('i', 'e'),
            'c': ('g', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'c'),
            'g': ('l', 'm'),
            'o': ('p', 'f')},
        'p': {'b': ('i', 'e'),
            'c': ('g', 'b'),
            'e': ('j', 'k'),
            'f': ('n', 'c'),
            'g': ('l', 'm'),
            'p': ('o', 'f')}}
    reseeded_with_unifurcations_preserved_relations = {
        'a': {'a': ('b', 'c'),
            'b': ('i', 'e'),
            'c': ('g', 'f'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'b': {'a': ('c',),
            'b': ('i', 'e', 'a'),
            'c': ('g', 'f'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'c': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('g', 'f', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'e': {'a': ('c',),
            'b': ('i', 'a'),
            'c': ('g', 'f'),
            'e': ('j', 'k', 'b'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'f': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('g', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'h', 'c'),
            'g': ('l', 'm'),
            'h': ('o', 'p')},
        'g': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('f', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm', 'c'),
            'h': ('o', 'p')},
        'h': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('g', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'c'),
            'g': ('l', 'm'),
            'h': ('o', 'p', 'f')},
        'i': {'a': ('c',),
            'b': ('e', 'a'),
            'c': ('g', 'f'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'i': ('b',)},
        'j': {'a': ('c',),
            'b': ('i', 'a'),
            'c': ('g', 'f'),
            'e': ('k', 'b'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'j': ('e',)},
        'k': {'a': ('c',),
            'b': ('i', 'a'),
            'c': ('g', 'f'),
            'e': ('j', 'b'),
            'f': ('n', 'h'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'k': ('e',)},
        'l': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('f', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('m', 'c'),
            'h': ('o', 'p'),
            'l': ('g',)},
        'm': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('f', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'h'),
            'g': ('l', 'c'),
            'h': ('o', 'p'),
            'm': ('g',)},
        'n': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('g', 'a'),
            'e': ('j', 'k'),
            'f': ('h', 'c'),
            'g': ('l', 'm'),
            'h': ('o', 'p'),
            'n': ('f',)},
        'o': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('g', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'c'),
            'g': ('l', 'm'),
            'h': ('p', 'f'),
            'o': ('h',)},
        'p': {'a': ('b',),
            'b': ('i', 'e'),
            'c': ('g', 'a'),
            'e': ('j', 'k'),
            'f': ('n', 'c'),
            'g': ('l', 'm'),
            'h': ('o', 'f'),
            'p': ('h',)}}

    def test_reseeding(self):
        # import pprint
        curated_tree_gen = curated_test_tree.CuratedTestTree()
        # x1 = {}
        # x2 = {}
        # for reseed_at_label in curated_tree_gen.internal_labels:
        # to do: test (1) update bipartitions, (2) edge length
        ref_tree, all_nodes, leaf_nodes, internal_nodes = curated_tree_gen.get_tree(
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=True)
        ref_edge_lengths = {}
        for nd in ref_tree:
            ref_edge_lengths[nd.label] = nd.edge.length
        for is_rooted in (True, False):
            for suppress_unifurcations, expected_preorder_labels, expected_relations, edge_length_updates in (
                    (True, self.reseeded_with_unifurcations_suppressed_preorder_labels, self.reseeded_with_unifurcations_suppressed_relations, self.reseeded_with_unifurcations_suppressed_edge_length_updates),
                    (False, self.reseeded_with_unifurcations_preserved_preorder_labels, self.reseeded_with_unifurcations_preserved_relations, self.reseeded_with_unifurcations_preserved_edge_length_updates),
                    ):
                for reseed_at_label in curated_tree_gen.all_labels:

                    if reseed_at_label == "a":
                        continue

                    # identify
                    # print("\n### is_rooted = {}, suppress_unifurcations = {}, reseed_at = '{}'".format(is_rooted, suppress_unifurcations, reseed_at_label))

                    # get tree
                    tree, all_nodes, leaf_nodes, internal_nodes = curated_tree_gen.get_tree(
                            suppress_internal_node_taxa=False,
                            suppress_leaf_node_taxa=False)

                    # label nodes
                    for nd in tree:
                        nd.label = nd.taxon.label
                    # set rooting
                    tree.is_rooted = is_rooted

                    # calc bipartitions
                    tree.encode_bipartitions(suppress_unifurcations=False, collapse_unrooted_basal_bifurcation=False)

                    # save old root
                    old_root = tree.seed_node

                    # find new root
                    new_root = tree.find_node_with_label(reseed_at_label)

                    # reroot it
                    tree.reseed_at(
                            new_root,
                            suppress_unifurcations=suppress_unifurcations,
                            collapse_unrooted_basal_bifurcation=False,
                            update_bipartitions=True)

                    # check that new root is as expected
                    self.assertEqual(tree.seed_node.label, reseed_at_label)

                    # check old root integrity
                    if old_root is not new_root:
                        old_root_still_in_tree = False
                        for nd in tree:
                            if nd is old_root:
                                old_root_still_in_tree = True
                                break
                        if old_root_still_in_tree:
                            self.assertIsNot(old_root._parent_node, None)
                            self.assertIs(old_root._edge._head_node, old_root)
                            self.assertIs(old_root._edge.tail_node, old_root._parent_node)
                            found_parent = None
                            for node in all_nodes:
                                if old_root in node._child_nodes:
                                    self.assertIs(found_parent, None)
                                    found_parent = node
                            self.assertIsNot(found_parent, None)
                            self.assertIs(found_parent, old_root._parent_node)
                        else:
                            self.assertEqual(len(old_root._child_nodes), 1)

                    # check new root integrity
                    self.assertIs(new_root._parent_node, None)
                    self.assertIs(new_root._edge._head_node, new_root)
                    self.assertIs(new_root._edge.tail_node, None)
                    for node in all_nodes:
                        self.assertNotIn(new_root, node._child_nodes)

                    # check that rooting state is as expected
                    self.assertEqual(tree.is_rooted, is_rooted)

                    # check for structural integrity, including bipartitions
                    if is_rooted and suppress_unifurcations:
                        tree._debug_check_tree(
                                logger_obj=_LOG,
                                check_bipartitions=True,
                                unique_bipartition_edge_mapping=True)
                    else:
                        tree._debug_check_tree(
                                logger_obj=_LOG,
                                check_bipartitions=True,
                                unique_bipartition_edge_mapping=False)

                    # check that traversal is as expected
                    preorder_labels = [nd.label for nd in tree]
                    self.assertEqual(preorder_labels, expected_preorder_labels[reseed_at_label])

                    # check that parent-child relationships are as expected
                    relations = {}
                    for nd in tree:
                        child_labels = tuple(ch.label for ch in nd._child_nodes)
                        if child_labels:
                            relations[nd.label] = child_labels
                    self.assertEqual(relations, expected_relations[reseed_at_label])

                    # check that edge lengths are correct
                    updated_edge_lengths = edge_length_updates[reseed_at_label]
                    for nd in tree:
                        if nd.label == reseed_at_label:
                            self.assertEqual(nd.edge.length, ref_edge_lengths['a'])
                        elif nd.label not in updated_edge_lengths:
                            self.assertEqual(nd.edge.length, ref_edge_lengths[nd.label], "New seed: {}, Current node: {}".format(reseed_at_label, nd.label))
                        else:
                            new_length = sum([ref_edge_lengths[clabel] for clabel in updated_edge_lengths[nd.label]])
                            self.assertEqual(nd.edge.length, new_length, "New seed: {}, Current node: {}".format(reseed_at_label, nd.label))

class ResolvePolytomiesTestCase(dendropytest.ExtendedTestCase):

    def verify_resolve_polytomies(self, tree_string, rng):
        tree = dendropy.Tree.get_from_string(tree_string, "newick")
        if "&U" in tree_string:
            assert not tree.is_rooted
        else:
            assert tree.is_rooted
        for nd in tree:
            nd.edge.length = 100
        tree.resolve_polytomies(rng=rng)
        tree.encode_bipartitions()
        tree._debug_check_tree(
                check_bipartitions=True,
                unique_bipartition_edge_mapping=True)
        for nd in tree:
            if nd is tree.seed_node and not tree.is_rooted:
                self.assertEqual(len(nd._child_nodes), 3)
            elif len(nd._child_nodes) > 0:
                self.assertEqual(len(nd._child_nodes), 2)
        tree2 = dendropy.Tree.get_from_string(tree_string, "newick", taxon_namespace=tree.taxon_namespace)
        self.assertNotEqual(treecompare.symmetric_difference(tree, tree2), 0)
        tree.collapse_unweighted_edges()
        self.assertEqual(treecompare.symmetric_difference(tree, tree2), 0)

    def test_resolve_polytomies_at_root(self):
        for tree_string in (
                "(a,b,c,d)e;",
                ):
            for rooting in ("[&R]", "[&U]"):
                tree_string2 = rooting + " " +  tree_string
                # cycle through rng period
                self.verify_resolve_polytomies(tree_string2, None)
                for x in range(1001):
                    rng = MockRandom()
                    for i in range(x):
                        rng.uniform(0, 1)
                    self.verify_resolve_polytomies(tree_string2, rng)

    def test_resolve_polytomies(self):
        for tree_string in (
                "((((Homo:0.21,Bogus1:0.23,Pongo:0.21)N1:0.28,Bogus2:0.49,Macaca:0.49)N2:0.13,Bogus3:0.62,Ateles:0.62)N3:0.38,Galago:1.00)N4:0.0;",
                "(((t52,t62,(t2,t58,(t32,(t55,t28,t39,t17,t4,t44,t25)internal6,t26,t9,t48,(t41,t45)internal7)internal5,t56)internal4,t54,((t18,t14,t34)internal9,(t49,t22,t50,t27,t16,t40,t6,t19)internal10,t13,(t51,t35,t61,t53,t43)internal11)internal8,((t42,t5,t7,t33,t30,t21,t47,t38)internal13,t23,t1,t11,t46)internal12,(t63,t3,t37,t59)internal14)internal3,t57,t64,t31,(t12,(t60,t24,t10,t15,t20)internal16,t36)internal15,t8)internal2,t29)internal1;",
                "((t13,t37,t19,((t21,t44)internal4,(t61,t46,t4,t8,t63)internal5,t23,t28)internal3,t52,t64,(t39,t40,t24)internal6,(t54,t62,t15,t55,t51)internal7)internal2,(t3,(t53,t33,(t47,t9,t25)internal10,t45,t18,t27,t17)internal9,(t10,t22)internal11,(t59,t20,t12,t57,t56,t38,t7,t11)internal12,((t1,t31,t43,t36,t34,t14,t49)internal14,t2,t41,t50,((t30,t32,t58,t6)internal16,t60,t16)internal15,t26,t48,t35)internal13,t42,t29,t5)internal8)internal1;",
                "((t7,t13,(t42,t51,t20,t26,(t21,t18,t16)internal4)internal3,(t6,t48,(t23,t33,t34,t15,t2,t25)internal6,(t64,t45,t49,t3,t55,t31,(t19,t47,t38,t35)internal8,t14)internal7,t9,(t5,(t62,t50)internal10,t54,t32,(t40,t8,t58,t60,t10,t30)internal11)internal9)internal5)internal2,t28,(t4,t57,t52,t43,t46,t39)internal12,t63,t11,(t44,(t12,t22,t36,t29,t24,t1,t17,t56)internal14,(t41,t59,t53,t61)internal15,t37)internal13,t27)internal1;",
                "(t10,t42,(t54,t12,t40,(t30,(t55,t3,t22,t56)internal4,t39)internal3,t51,((t41,(t28,t52,t24,t14,t49,t38,t36,(t35,t34,t13,t9,t59,t58)internal8)internal7,t5,t45,t17,(t23,(t11,t53,t57,t19,t26)internal10)internal9,t63,(t62,t29,t18,t20,t27,t43)internal11)internal6,t4,(t25,t7,(((t47,t61,t21,t64)internal15,t44,t46,t15,t6,t37,t48)internal14,t50,t60)internal13)internal12)internal5,t31)internal2,t1,t32,t33,(t2,t8)internal16,t16)internal1;",
                "(t60,(((t29,t35,t5,t59,t4,t9,(t6,t25,t37,t44,t54)internal5)internal4,(t53,t23,t28,t48,((t16,t46,t26,t10,t52)internal8,t30,t17,t51,t40)internal7,t31,t3)internal6,t14,t7,t2,t22,(t45,t56,t20,t36,t43,t47)internal9,t8)internal3,t18,(t39,(t61,t27,t21,t58,t24)internal11,t50,(t64,t38,(t32,t11,t49,t63)internal13)internal12,t13,((t57,t34,t15,t55,t19)internal15,t1,t42)internal14,t33)internal10,(t62,t12,t41)internal16)internal2)internal1;",
                "((t12,t16,(t51,t24,t27,(t38,t44,t52,t6,t9,t53,t20,t18)internal4,(t31,t47,t56,t60)internal5)internal3,t34,t42)internal2,t5,t62,((t61,t64,t59,t15,t48,t2,t35)internal7,t57,(t14,t28,t40,t22,t58,t7,(t25,t4)internal9,((t8,t29,t21)internal11,t17,t1,t46,(t33,t50,t11)internal12,t13,t45,t37)internal10)internal8,(t26,t10,(t3,t30,t32,(t19,t36,t39,t43)internal15,t23,t54,t55)internal14)internal13,t41,t49,t63)internal6)internal1;",
                "(((t31,t19,t12,(t55,t25)internal4,t11,t60,(t46,t8,t56)internal5)internal3,(t52,t18)internal6,(((t24,t16,t7)internal9,t57,t40,(t2,t9,t50,t37,t43,t20,t15,t22)internal10,t42,(t47,t28,t58,t10)internal11)internal8,t21,(t14,(t51,t53,t26,t35,t49)internal13,t33,t62,t34,t54,t17)internal12,t5,t61,t4)internal7,t3,t36,t45,t63)internal2,(t27,t32,t30,t48)internal14,(t13,t64,(t23,t59,t41,t38,t44,t1)internal16,t29,t39,t6)internal15)internal1;",
                "(((t8,t40,t25,t36,t37,t11)internal3,t29,t31,t51,t3,((t61,t43,t63,t50,t23,t52,t24,t30)internal5,t39,t44,t58,((t15,t64,t9,t28,(t5,t34,t38,t22,(t33,t12,t35,t42,t10,t2,t27,t45)internal9)internal8,t7,t21,t60)internal7,t41,(t48,t17,t14)internal10)internal6,t13,t18)internal4)internal2,((t47,t16,t49,(t53,t26)internal13,(t19,t20,(t55,t54)internal15,t32,t4)internal14,t59,t57,t6)internal12,t46,t56,t1)internal11,t62)internal1;",
                "(t33,t39,t18,t35,((t19,t62,t55,(t41,t14,(t1,t36,t16,t38,t25,t59,t34)internal5)internal4,t61,t50)internal3,(((t2,t48,t22)internal8,t28,t37,t47,(t60,t30)internal9,t27)internal7,t12,(t31,t21,(t3,t5,t45,t53)internal11,(t23,t54,t20,t4,t64,t56,(t58,t13,t26,t11,t57,t44,t42,t46)internal13)internal12,t40)internal10)internal6,t7,(t52,(t43,t24)internal15,(t49,t29,t63,t32)internal16,t9,t6,t8,t15)internal14,t51)internal2,t10,t17)internal1;",
                "(t35,((t60,t47,t58,t26,t9)internal3,t3,((t13,t25)internal5,t1,((t6,t32,(t53,t7,t64,t34,t18,t23,t30,t33)internal8,t55,(t48,t20,t12,t4,t38,t28,t36)internal9)internal7,t8,((t57,t40,t52,t31,t43,t39)internal11,t59,(t37,t16,t27,(t44,t41,t10)internal13,t50)internal12,t24,t63)internal10,(t5,t56,t61,t29)internal14,t21,t49)internal6)internal4,((t17,t42,t62,t15)internal16,(t19,t2,t51,(t22,t14,(t45,t54,t46)internal19)internal18)internal17,t11)internal15)internal2)internal1;",
                ):
            for rooting in ("[&R]", "[&U]"):
                tree_string2 = rooting + " " +  tree_string
                for rng in (MockRandom(), None):
                    self.verify_resolve_polytomies(tree_string2, rng)

class TestStructureExtraction(
        curated_test_tree.CuratedTestTree,
        compare_and_validate.Comparator,
        unittest.TestCase):

    def setUp(self):
        self.taxon_label_inclusion_sets = (
            ("i", "k", "l", "n", "p"),
            ("i", "j", "k", "l", "m", "n", "p"),
            ("j", "k"),
            ("i", "k"),
            ("i", "p"),
            ("i", "j", "n"),
        )

    def test_full_extract(self):
        tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                suppress_internal_node_taxa=True,
                suppress_leaf_node_taxa=False)
        tree2 = tree1.extract_tree()
        self.compare_distinct_trees(tree1, tree2,
                taxon_namespace_scoped=True,
                compare_tree_annotations=False,
                compare_taxon_annotations=False)

    def test_filtered_leaf_extract(self):
        for suppress_internal_node_taxa in (True, ):#False):
            for inclusion_set in self.taxon_label_inclusion_sets:
                tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                        suppress_internal_node_taxa=suppress_internal_node_taxa,
                        suppress_leaf_node_taxa=False)
                to_include = set([tree1.taxon_namespace.get_taxon(label) for label in inclusion_set])
                assert None not in to_include
                node_filter_fn = lambda nd: nd.is_internal() or nd.taxon in to_include
                tree2 = tree1.extract_tree(node_filter_fn=node_filter_fn)
                tree1.retain_taxa(to_include)
                self.compare_distinct_trees(tree1, tree2,
                        taxon_namespace_scoped=True,
                        compare_tree_annotations=False,
                        compare_taxon_annotations=False)

    def test_extract_with_taxa(self):
        for taxon_label_inclusion_set in self.taxon_label_inclusion_sets:
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False)
            to_include_taxa = set([tree1.taxon_namespace.get_taxon(label) for label in taxon_label_inclusion_set])
            assert None not in to_include_taxa
            tree1.retain_taxa(to_include_taxa)
            tree2 = tree1.extract_tree_with_taxa_labels(taxon_label_inclusion_set)
            self.compare_distinct_trees(tree1, tree2,
                    taxon_namespace_scoped=True,
                    compare_tree_annotations=False,
                    compare_taxon_annotations=False)
            tree3 = tree1.extract_tree_with_taxa(to_include_taxa)
            self.compare_distinct_trees(tree1, tree3,
                    taxon_namespace_scoped=True,
                    compare_tree_annotations=False,
                    compare_taxon_annotations=False)

    def test_extract_without_taxa(self):
        for taxon_label_inclusion_set in self.taxon_label_inclusion_sets:
            tree1, anodes1, lnodes1, inodes1 = self.get_tree(
                    suppress_internal_node_taxa=False,
                    suppress_leaf_node_taxa=False)
            to_exclude_taxa = set([t for t in tree1.taxon_namespace if t.label not in taxon_label_inclusion_set])
            assert None not in to_exclude_taxa
            tree1.prune_taxa(to_exclude_taxa)
            to_exclude_taxa_labels = set([t.label for t in to_exclude_taxa])
            tree2 = tree1.extract_tree_without_taxa(to_exclude_taxa)
            self.compare_distinct_trees(tree1, tree2,
                    taxon_namespace_scoped=True,
                    compare_tree_annotations=False,
                    compare_taxon_annotations=False)
            tree3 = tree1.extract_tree_without_taxa_labels(to_exclude_taxa_labels)
            self.compare_distinct_trees(tree1, tree3,
                    taxon_namespace_scoped=True,
                    compare_tree_annotations=False,
                    compare_taxon_annotations=False)

    def test_special_case1(self):
        original_tree_str = """\
        [&R] ((((e1:4.25978504749,a0:4.25978504749):9.75100657322,(e5:11.2557415909,c9:11.2557415909):2.75505002977):5.25672273638,(c5:17.0225375511,e6:17.0225375511):2.24497680601):20.9755404109,(((c7:0.0433876754663,e4:0.0433876754663):16.2031718648,(b1:14.1628944123,d7:14.1628944123):2.08366512802):14.3825543479,((((d1:13.4235384066,(d4:7.64533761739,c3:7.64533761739):5.77820078917):2.00948796838,((d8:3.10025757397,b5:3.10025757397):5.07496414931,a4:8.17522172328):7.25780465166):4.52823355379,((((((a7:8.94718577977,(((a1:2.04048640276,c2:2.04048640276):1.45629935083,(e0:0.408302025932,b6:0.408302025932):3.08848372766):3.77714533326,(((c6:2.1238494561,(e8:2.03255428077,d6:2.03255428077):0.0912951753249):2.91822700988,a5:5.04207646598):1.92173681425,((a2:3.43218264885,(b8:0.515232535857,a9:0.515232535857):2.91695011299):1.6832785054,b4:5.11546115425):1.84835212598):0.310117806629):1.67325469292):0.613875266884,(d9:8.93428444448,(c1:5.91732320427,c8:5.91732320427):3.0169612402):0.626776602178):3.65721021136,((c0:3.99662328128,d2:3.99662328128):1.90572648225,(e9:1.84550535315,(b9:0.803660457957,e3:0.803660457957):1.0418448952):4.05684441038):7.31592149449):1.63573163655,d5:14.8540028946):2.25255068893,d3:17.1065535835):2.27795405888,(((a3:4.85356559967,(c4:3.08209866724,d0:3.08209866724):1.77146693244):2.74425153816,e7:7.59781713783):4.22596432824,(b2:2.86856170856,e2:2.86856170856):8.9552197575):7.56072617631):0.576752286338):1.05878660087,((b0:1.6464852541,b7:1.6464852541):10.0630186678,(a8:7.31781944487,(a6:7.13495568605,b3:7.13495568605):0.182863758824):4.39168447703):9.31054260769):9.60906735859):9.61394087983):6.65318140005;
        """
        expected_tree_strs = """\
        [&R] (a0:40.243054768,((a4:19.9612599287,((a7:8.94718577977,(a1:7.27393108685,(a5:6.96381328023,(a2:3.43218264885,a9:3.43218264885):3.53163063138):0.310117806629):1.67325469292):10.4373218626,a3:19.3845076424):0.576752286338):1.05878660087,(a8:7.31781944487,a6:7.31781944487):13.7022270847):19.2230082384):6.65318140005;
        [&R] (b1:30.6291138882,((b5:19.9612599287,(((b6:7.27393108685,(b8:5.11546115425,b4:5.11546115425):2.15846993261):5.94434017116,b9:13.218271258):6.16623638436,b2:19.3845076424):0.576752286338):1.05878660087,((b0:1.6464852541,b7:1.6464852541):10.0630186678,b3:11.7095039219):9.31054260769):9.60906735859):16.2671222799;
        [&R] ((c9:19.2675143571,c5:19.2675143571):20.9755404109,(c7:30.6291138882,(c3:19.9612599287,((((c2:7.27393108685,c6:7.27393108685):2.2871299598,(c1:5.91732320427,c8:5.91732320427):3.64373784238):3.65721021136,c0:13.218271258):6.16623638436,c4:19.3845076424):0.576752286338):10.6678539595):9.61394087983):6.65318140005;
        [&R] (d7:30.6291138882,(((d1:13.4235384066,d4:13.4235384066):2.00948796838,d8:15.4330263749):4.52823355379,(((((d6:9.56106104665,d9:9.56106104665):3.65721021136,d2:13.218271258):1.63573163655,d5:14.8540028946):2.25255068893,d3:17.1065535835):2.27795405888,d0:19.3845076424):0.576752286338):10.6678539595):16.2671222799;
        [&R] (((e1:14.0107916207,e5:14.0107916207):5.25672273638,e6:19.2675143571):20.9755404109,(e4:30.6291138882,(((e0:7.27393108685,e8:7.27393108685):5.94434017116,(e9:1.84550535315,e3:1.84550535315):11.3727659049):6.16623638436,(e7:11.8237814661,e2:11.8237814661):7.56072617631):11.2446062458):9.61394087983):6.65318140005;
        """

        tns = dendropy.TaxonNamespace()
        source_tree1 = dendropy.Tree.get(
                data=original_tree_str,
                schema="newick",
                taxon_namespace=tns)
        source_tree2 = dendropy.Tree.get(
                data=original_tree_str,
                schema="newick",
                taxon_namespace=tns)
        self.assertEqual(treecompare.weighted_robinson_foulds_distance(source_tree1, source_tree2), 0.0)
        group_ids = ("a", "b", "c", "d", "e")
        expected_induced_trees = dendropy.TreeList.get(
                data=expected_tree_strs,
                schema="newick",
                taxon_namespace=tns)
        assert len(expected_induced_trees) == len(group_ids)
        for group_id, expected_induced_tree in zip(group_ids, expected_induced_trees):
            extracted_tree = source_tree1.extract_tree(
                    node_filter_fn=lambda node: node.taxon.label.startswith(group_id),
                    is_apply_filter_to_leaf_nodes=True,
                    is_apply_filter_to_internal_nodes=False)
            for leaf_nd in extracted_tree.leaf_node_iter():
                self.assertTrue(leaf_nd.taxon.label.startswith(group_id))
            for leaf_nd in expected_induced_tree.leaf_node_iter():
                assert leaf_nd.taxon.label.startswith(group_id)

            # self.assertEqual(treecompare.weighted_robinson_foulds_distance(source_tree1, source_tree2), 0.0)
            self.assertEqual(treecompare.unweighted_robinson_foulds_distance(extracted_tree, expected_induced_tree), 0)
            self.assertAlmostEqual(treecompare.weighted_robinson_foulds_distance(extracted_tree, expected_induced_tree), 0.0)

    def test_special_case2(self):
        original_tree_str = """\
        [&R] ((a1,(a2,(a3,a4)a0)),(b1,(b2,(b3,b4))));
        """
        source_tree1 = dendropy.Tree.get(
                data=original_tree_str,
                schema="newick",
                )
        expected_tree_str = """\
        [&R] ((a3,a4)a0);
        """
        extracted_tree = source_tree1.extract_tree(
                # node_filter_fn=lambda node: node.taxon is not None and node.taxon.label.startswith("a"),
                node_filter_fn=lambda node: node.taxon.label in set(["a3", "a4"]),
                is_apply_filter_to_leaf_nodes=True,
                is_apply_filter_to_internal_nodes=False,
                )
        expected_tree = dendropy.Tree.get(
                data=expected_tree_str,
                schema="newick",
                taxon_namespace=source_tree1.taxon_namespace)
        self.assertEqual(treecompare.unweighted_robinson_foulds_distance(extracted_tree, expected_tree), 0.0)

class TreeRestructuring(dendropytest.ExtendedTestCase):

    def test_collapse_basal_bifurcation(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_reseed_at(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_to_outgroup_position(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_reroot_at_node(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_reroot_at_edge(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_suppress_unifurcations(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_collapse_unweighted_edges(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_resolve_polytomies(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_prune_subtree(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_prune_leaves_without_taxa(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_prune_taxa(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_prune_nodes(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_prune_taxa_with_labels(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_retain_taxa(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_retain_taxa_with_labels(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_randomly_reorient_tree(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_randomly_rotate(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_ladderize(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_truncate_from_root(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_scale_edges(self):
        self.assertFalse(self.fail_incomplete_tests())

    def test_set_edge_lengths_from_node_ages(self):
        self.assertFalse(self.fail_incomplete_tests())

if __name__ == "__main__":
    unittest.main()

