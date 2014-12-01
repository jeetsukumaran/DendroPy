
#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Tests native tree structuring routines.
"""

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
import sys
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
import unittest
from dendropy.test.support import curated_test_tree
from dendropy.test.support import pathmap
from dendropy.test.support import dendropytest
from dendropy.utility import messaging
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.test.support.mockrandom import MockRandom
import dendropy
from dendropy.calculate import treecompare
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
                tree_leafset_bitmask=tree_leafset_bitmask)
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
                tree_leafset_bitmask=tree_leafset_bitmask)
        tree.seed_node.collapse_conflicting(bipartition_to_target)
        tree.encode_bipartitions()
        expected_tree.encode_bipartitions()
        self.assertEqual(treecompare.symmetric_difference(tree, expected_tree), 0)

        tree = tree_list[4]
        expected_tree = tree_list[5]
        tree.encode_bipartitions()
        tree_leafset_bitmask = tree.seed_node.edge.bipartition._leafset_bitmask
        bipartition_to_target = dendropy.Bipartition(bitmask=0x5,
                tree_leafset_bitmask=tree_leafset_bitmask)
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
        for nd in tree:
            nd.edge.length = 100
        tree.resolve_polytomies(rng=rng)
        tree._debug_check_tree()
        for nd in tree:
            if len(nd._child_nodes) > 0:
                self.assertEqual(len(nd._child_nodes), 2)
        tree2 = dendropy.Tree.get_from_string(tree_string, "newick", taxon_namespace=tree.taxon_namespace)
        self.assertNotEqual(treecompare.symmetric_difference(tree, tree2), 0)
        tree.collapse_unweighted_edges()
        self.assertEqual(treecompare.symmetric_difference(tree, tree2), 0)

    def test_resolve_polytomies_at_root(self):
        for tree_string in (
                "(a,b,c,d)e;",
                ):
            # cycle through rng period
            self.verify_resolve_polytomies(tree_string, None)
            for x in range(1001):
                rng = MockRandom()
                for i in range(x):
                    rng.uniform(0, 1)
                self.verify_resolve_polytomies(tree_string, rng)

    def test_resolve_polytomies(self):
        for tree_string in (
                "(((t52:9.00000000E-04,t62:3.00000000E-04,(t2:1.98838500E+02,t58:5.50000000E-03,(t32:1.00000000E-04,(t55:1.60000000E-03,t28:3.40000000E-03,t39:9.67173000E+01,t17:6.10000000E-03,t4:7.30000000E-03,t44:9.80000000E-03,t25:5.60000000E-03)internal6:7.00000000E-03,t26:6.20000000E-03,t9:5.80000000E-03,t48:5.00000000E-04,(t41:2.32791700E+02,t45:2.77974400E+02)internal7:2.89388700E+02)internal5:5.89266100E+02,t56:3.70000000E-03)internal4:5.05355800E+02,t54:2.10000000E-03,((t18:7.29127100E+02,t14:3.79456100E+02,t34:6.50000000E-03)internal9:6.10000000E-03,(t49:7.80000000E-03,t22:3.21013000E+01,t50:2.70000000E-03,t27:9.42909800E+02,t16:3.20000000E-03,t40:4.00000000E-03,t6:4.60000000E-03,t19:2.46628300E+02)internal10:6.90000000E-03,t13:2.70000000E-03,(t51:2.20000000E-03,t35:5.10000000E-03,t61:4.71174000E+01,t53:6.27446500E+02,t43:4.30000000E-03)internal11:9.00000000E-03)internal8:9.36654700E+02,((t42:8.70000000E-03,t5:7.20722100E+02,t7:5.40000000E-03,t33:6.40962200E+02,t30:4.34765900E+02,t21:9.60000000E-03,t47:2.70000000E-03,t38:1.80000000E-03)internal13:5.30000000E-03,t23:8.80000000E-03,t1:6.38949900E+02,t11:1.60000000E-03,t46:5.40000000E-03)internal12:3.81620000E+02,(t63:3.24156800E+02,t3:9.29098700E+02,t37:8.40000000E-03,t59:6.00000000E-04)internal14:5.40000000E-03)internal3:6.80000000E-03,t57:9.50000000E-03,t64:4.85991000E+02,t31:7.60602500E+02,(t12:5.50000000E-03,(t60:2.20000000E-03,t24:7.30000000E-03,t10:3.11717000E+02,t15:6.50000000E-03,t20:5.20000000E-03)internal16:8.72433200E+02,t36:2.24698200E+02)internal15:1.30000000E-03,t8:5.90000000E-03)internal2:2.23211600E+02,t29:2.20218200E+02)internal1:3.00000000E-04;",
                "((t13:8.00000000E-04,t37:6.68978200E+02,t19:1.32312800E+02,((t21:2.00000000E-03,t44:4.31051800E+02)internal4:7.90000000E-03,(t61:6.80000000E-03,t46:1.00000000E-03,t4:3.40000000E-03,t8:2.50000000E-03,t63:4.48614100E+02)internal5:4.70000000E-03,t23:2.80000000E-03,t28:9.23265700E+02)internal3:5.80000000E-03,t52:8.70000000E-03,t64:6.00000000E-04,(t39:9.30000000E-03,t40:1.70000000E-03,t24:2.20000000E-03)internal6:8.40000000E-03,(t54:9.90000000E-03,t62:7.84070900E+02,t15:4.30000000E-03,t55:1.00000000E-02,t51:7.20000000E-03)internal7:6.00000000E-04)internal2:9.05420200E+02,(t3:5.80000000E-03,(t53:5.84178100E+02,t33:8.60000000E-03,(t47:1.85825800E+02,t9:6.80000000E-03,t25:1.19887500E+02)internal10:9.60779000E+02,t45:2.50000000E-03,t18:6.20000000E-03,t27:5.90000000E-03,t17:9.40000000E-03)internal9:7.50000000E-03,(t10:3.95786500E+02,t22:3.00000000E-03)internal11:7.16192100E+02,(t59:4.58286100E+02,t20:1.00000000E-02,t12:2.13155100E+02,t57:9.33259500E+02,t56:8.80000000E-03,t38:1.60000000E-03,t7:7.10000000E-03,t11:9.90000000E-03)internal12:7.60000000E-03,((t1:9.40000000E-03,t31:1.15429600E+02,t43:5.53224100E+02,t36:6.04830300E+02,t34:2.10000000E-03,t14:2.70000000E-03,t49:9.10000000E-03)internal14:3.00000000E-03,t2:1.00000000E-03,t41:2.80000000E-03,t50:7.71119500E+02,((t30:1.00000000E-04,t32:8.83106600E+02,t58:5.50000000E-03,t6:5.90000000E-03)internal16:5.60000000E-03,t60:1.27446400E+02,t16:9.00000000E-03)internal15:2.70000000E-03,t26:8.70000000E-03,t48:2.20000000E-03,t35:1.02794600E+02)internal13:1.00000000E-04,t42:8.90000000E-03,t29:6.30000000E-03,t5:9.29881100E+02)internal8:2.96708600E+02)internal1:6.31103400E+02;",
                "((t7:7.36564800E+02,t13:9.40000000E-03,(t42:2.70000000E-03,t51:1.10000000E-03,t20:8.60000000E-03,t26:8.16586800E+02,(t21:2.30000000E-03,t18:1.93130600E+02,t16:8.70000000E-03)internal4:3.10000000E-03)internal3:8.70000000E-03,(t6:4.00000000E-03,t48:5.40000000E-03,(t23:1.80000000E-03,t33:2.70000000E-03,t34:4.34564400E+02,t15:3.20000000E-03,t2:5.20000000E-03,t25:5.80000000E-03)internal6:1.20000000E-03,(t64:9.25575300E+02,t45:7.20000000E-03,t49:6.80000000E-03,t3:7.00000000E-04,t55:3.40000000E-03,t31:8.50000000E-03,(t19:4.10000000E-03,t47:3.00000000E-03,t38:4.20447000E+02,t35:6.80000000E-03)internal8:3.10000000E-03,t14:6.20000000E-03)internal7:7.08785600E+02,t9:5.50000000E-03,(t5:5.80000000E-03,(t62:2.20000000E-03,t50:9.10000000E-03)internal10:4.70000000E-03,t54:1.70000000E-03,t32:5.15452500E+02,(t40:7.60000000E-03,t8:2.30000000E-03,t58:2.44264000E+01,t60:4.75136900E+02,t10:8.00000000E-04,t30:6.30000000E-03)internal11:3.40000000E-03)internal9:2.86914500E+02)internal5:2.80000000E-03)internal2:8.90000000E-03,t28:6.96354600E+02,(t4:1.00000000E-04,t57:7.80000000E-03,t52:4.25146800E+02,t43:9.57966100E+02,t46:6.00000000E-04,t39:8.48336500E+02)internal12:2.50000000E-03,t63:8.10000000E-03,t11:9.90000000E-03,(t44:6.20000000E-03,(t12:5.50000000E-03,t22:7.50000000E-03,t36:2.60000000E-03,t29:6.37661800E+02,t24:5.30000000E-03,t1:2.80000000E-03,t17:2.85728900E+02,t56:3.19710300E+02)internal14:5.10000000E-03,(t41:6.93950200E+02,t59:7.00000000E-04,t53:5.50000000E-03,t61:2.10000000E-03)internal15:1.40000000E-03,t37:9.10000000E-03)internal13:9.60000000E-03,t27:7.00000000E-03)internal1:9.80000000E-03;",
                "(t10:1.00000000E-04,t42:7.60000000E-03,(t54:7.50000000E-03,t12:6.10000000E-03,t40:2.19422500E+02,(t30:6.80000000E-03,(t55:1.30000000E-03,t3:3.00000000E-04,t22:5.70000000E-03,t56:6.42750100E+02)internal4:1.70000000E-03,t39:4.61699000E+02)internal3:3.00000000E-04,t51:3.79104500E+02,((t41:8.40000000E-03,(t28:5.39419600E+02,t52:3.50000000E-03,t24:7.90000000E-03,t14:7.60000000E-03,t49:2.00000000E-03,t38:1.52383200E+02,t36:6.69459200E+02,(t35:7.70000000E-03,t34:6.07247900E+02,t13:1.20000000E-03,t9:9.70000000E-03,t59:2.56794000E+01,t58:6.80000000E-03)internal8:2.20000000E-03)internal7:8.19072000E+01,t5:4.00000000E-03,t45:8.00000000E-04,t17:6.30000000E-03,(t23:6.10000000E-03,(t11:3.50000000E-03,t53:3.80000000E-03,t57:3.50000000E-03,t19:8.30000000E-03,t26:9.60787700E+02)internal10:9.83844400E+02)internal9:7.72481200E+02,t63:8.30000000E-03,(t62:9.70000000E-03,t29:8.08199500E+02,t18:4.90000000E-03,t20:7.40000000E-03,t27:8.51713400E+02,t43:9.00000000E-04)internal11:4.40000000E-03)internal6:3.80000000E-03,t4:2.50000000E-03,(t25:3.00000000E-04,t7:1.90000000E-03,(((t47:9.48062700E+02,t61:7.20075000E+02,t21:7.60000000E-03,t64:6.50000000E-03)internal15:2.80000000E-03,t44:8.00000000E-03,t46:2.21597200E+02,t15:3.10000000E-03,t6:4.80000000E-03,t37:4.30000000E-03,t48:3.40000000E-03)internal14:8.90000000E-03,t50:2.80000000E-03,t60:1.20656800E+02)internal13:7.97832600E+02)internal12:6.20000000E-03)internal5:3.26846500E+02,t31:1.19555700E+02)internal2:9.60000000E-03,t1:7.70000000E-03,t32:2.16385500E+02,t33:7.30000000E-03,(t2:1.00527100E+02,t8:2.27484000E+02)internal16:5.30000000E-03,t16:1.41494000E+01)internal1:7.70000000E-03;",
                "(t60:5.60000000E-03,(((t29:4.50000000E-03,t35:5.13777000E+01,t5:4.80000000E-03,t59:4.10000000E-03,t4:6.29446100E+02,t9:1.49198500E+02,(t6:1.20000000E-03,t25:6.10000000E-03,t37:2.30000000E-03,t44:4.60000000E-03,t54:8.70000000E-03)internal5:3.10000000E-03)internal4:5.80000000E-03,(t53:6.10000000E-03,t23:9.50000000E-03,t28:2.11026000E+02,t48:1.60000000E-03,((t16:3.00000000E-03,t46:7.06956200E+02,t26:4.60000000E-03,t10:9.30000000E-03,t52:6.30000000E-03)internal8:4.50504300E+02,t30:9.40000000E-03,t17:5.50000000E-03,t51:9.10000000E-03,t40:8.00000000E-04)internal7:7.50658000E+01,t31:3.07612600E+02,t3:5.70000000E-03)internal6:3.10000000E-03,t14:1.24354600E+02,t7:7.00000000E-03,t2:5.10000000E-03,t22:9.00000000E-04,(t45:2.60000000E-03,t56:9.61911200E+02,t20:5.80000000E-03,t36:1.00000000E-02,t43:2.70000000E-03,t47:7.56269100E+02)internal9:4.32029300E+02,t8:6.60000000E-03)internal3:5.00000000E-03,t18:1.90000000E-03,(t39:1.00000000E-03,(t61:8.60000000E-03,t27:7.00000000E-03,t21:9.45193600E+02,t58:5.60000000E-03,t24:3.20000000E-03)internal11:3.40000000E-03,t50:9.80000000E-03,(t64:6.60000000E-03,t38:4.20000000E-03,(t32:6.30000000E-03,t11:2.10000000E-03,t49:9.30000000E-03,t63:7.00000000E-03)internal13:3.61780000E+02)internal12:5.20000000E-03,t13:9.73146900E+02,((t57:1.00000000E-03,t34:9.20000000E-03,t15:1.20000000E-03,t55:5.00000000E-03,t19:5.10000000E-03)internal15:5.50839500E+02,t1:8.34724200E+02,t42:2.50000000E-03)internal14:2.40000000E-03,t33:3.90000000E-03)internal10:6.40000000E-03,(t62:8.70000000E-03,t12:7.90006400E+02,t41:9.40000000E-03)internal16:5.10000000E-03)internal2:1.30000000E-03)internal1:2.80000000E-03;",
                "((t12:1.00000000E-02,t16:1.50000000E-03,(t51:5.50000000E-03,t24:3.30000000E-03,t27:6.97498800E+02,(t38:5.00000000E-04,t44:7.10000000E-03,t52:3.90000000E-03,t6:8.30000000E-03,t9:5.00000000E-03,t53:5.02291700E+02,t20:8.70000000E-03,t18:4.50000000E-03)internal4:2.40000000E-03,(t31:4.10000000E-03,t47:1.60000000E-03,t56:9.70000000E-03,t60:7.00000000E-03)internal5:4.60000000E-03)internal3:3.50000000E-03,t34:8.60000000E-03,t42:8.60000000E-03)internal2:2.98556600E+02,t5:3.20000000E-03,t62:7.60000000E-03,((t61:1.00000000E-02,t64:4.40000000E-03,t59:6.33748200E+02,t15:4.50000000E-03,t48:9.10000000E-03,t2:7.96143700E+02,t35:3.74841800E+02)internal7:6.31161400E+02,t57:5.31166800E+02,(t14:7.89693000E+01,t28:6.20000000E-03,t40:9.12829200E+02,t22:4.61150500E+02,t58:2.55327500E+02,t7:8.04633300E+02,(t25:4.41730400E+02,t4:5.90000000E-03)internal9:6.80000000E-03,((t8:3.84693900E+02,t29:2.40000000E-03,t21:5.80000000E-03)internal11:8.45304100E+02,t17:3.00000000E-03,t1:6.60000000E-03,t46:5.65929500E+02,(t33:6.06406900E+02,t50:9.10000000E-03,t11:7.98860400E+02)internal12:7.50000000E-03,t13:3.60000000E-03,t45:6.30000000E-03,t37:7.30000000E-03)internal10:4.30000000E-03)internal8:8.00000000E-03,(t26:9.10000000E-03,t10:3.10000000E-03,(t3:1.00000000E-03,t30:7.47481200E+02,t32:1.32160700E+02,(t19:9.20000000E-03,t36:2.60000000E-03,t39:4.90000000E-03,t43:7.50000000E-03)internal15:9.80000000E-03,t23:1.20000000E-03,t54:1.50000000E-03,t55:8.70000000E-03)internal14:5.80000000E-03)internal13:8.40000000E-03,t41:9.80000000E-03,t49:8.10000000E-03,t63:8.00000000E-03)internal6:4.00000000E-04)internal1:8.30000000E-03;",
                "(((t31:9.50000000E-03,t19:2.33822700E+02,t12:8.84235000E+02,(t55:6.32180000E+01,t25:9.10000000E-03)internal4:9.20000000E-03,t11:4.10000000E-03,t60:9.46261600E+02,(t46:8.87259500E+02,t8:4.53644400E+02,t56:7.50000000E-03)internal5:5.00000000E-03)internal3:5.90000000E-03,(t52:1.60000000E-03,t18:1.00000000E-03)internal6:4.20000000E-03,(((t24:8.60000000E-03,t16:6.85680600E+02,t7:9.80000000E-03)internal9:4.50000000E-03,t57:4.00000000E-03,t40:1.22160600E+02,(t2:6.50000000E-03,t9:3.00000000E-04,t50:2.50000000E-03,t37:5.64578500E+02,t43:7.65157700E+02,t20:2.15952200E+02,t15:3.30000000E-03,t22:9.00531200E+02)internal10:4.60000000E-03,t42:8.58406300E+02,(t47:1.74498000E+02,t28:3.00000000E-04,t58:7.89984900E+02,t10:3.23772200E+02)internal11:1.29993100E+02)internal8:4.00000000E-04,t21:5.24000000E+01,(t14:4.80000000E-03,(t51:5.03839200E+02,t53:5.00000000E-04,t26:8.70000000E-03,t35:4.60000000E-03,t49:7.00000000E-04)internal13:5.20000000E-03,t33:4.10000000E-03,t62:5.00000000E-03,t34:8.00000000E-04,t54:6.08297900E+02,t17:2.75016800E+02)internal12:5.30000000E-03,t5:5.50000000E-03,t61:1.00000000E-02,t4:4.60000000E-03)internal7:4.10000000E-03,t3:1.00000000E-03,t36:8.60000000E-03,t45:7.70000000E-03,t63:2.20000000E-03)internal2:6.80000000E-03,(t27:4.60000000E-03,t32:1.00000000E-03,t30:6.70000000E-03,t48:1.60000000E-03)internal14:2.30000000E-03,(t13:1.00000000E-03,t64:7.38796300E+02,(t23:5.70000000E-03,t59:1.30000000E-03,t41:6.70000000E-03,t38:8.70000000E-03,t44:9.70000000E-03,t1:3.60000000E-03)internal16:5.60000000E-03,t29:2.20000000E-03,t39:2.30000000E-03,t6:8.45373600E+02)internal15:7.00000000E-04)internal1:5.36572800E+02;",
                "(((t8:6.20000000E-03,t40:3.00000000E-04,t25:8.30000000E-03,t36:1.80416900E+02,t37:3.10000000E-03,t11:1.00000000E-04)internal3:9.80000000E-03,t29:5.70000000E-03,t31:1.50000000E-03,t51:4.00000000E-04,t3:2.20000000E-03,((t61:5.70000000E-03,t43:5.26098000E+01,t63:6.17838600E+02,t50:2.80000000E-03,t23:4.90000000E-03,t52:2.80000000E-03,t24:1.20000000E-03,t30:2.90000000E-03)internal5:7.30000000E-03,t39:4.90000000E-03,t44:5.00000000E-04,t58:6.00000000E-03,((t15:3.28123100E+02,t64:7.47014400E+02,t9:3.80208100E+02,t28:5.10000000E-03,(t5:8.61764600E+02,t34:1.87426000E+01,t38:8.70000000E-03,t22:5.80000000E-03,(t33:2.08576000E+01,t12:8.10000000E-03,t35:8.40000000E-03,t42:9.00000000E-04,t10:2.50000000E-03,t2:5.30000000E-03,t27:3.20000000E-03,t45:3.40000000E-03)internal9:4.20000000E-03)internal8:8.10000000E-03,t7:5.10483800E+02,t21:5.09952600E+02,t60:3.50000000E-03)internal7:1.37235700E+02,t41:8.20000000E-03,(t48:1.97272700E+02,t17:7.70000000E-03,t14:1.80000000E-03)internal10:2.40000000E-03)internal6:3.01420100E+02,t13:5.00000000E-03,t18:7.60000000E-03)internal4:3.50000000E-03)internal2:8.10000000E-03,((t47:8.70000000E-03,t16:1.60000000E-03,t49:4.42208700E+02,(t53:5.07337400E+02,t26:4.60000000E-03)internal13:2.80000000E-03,(t19:8.64250200E+02,t20:6.20000000E-03,(t55:8.00000000E-03,t54:4.50000000E-03)internal15:4.70000000E-03,t32:6.60000000E-03,t4:4.80000000E-03)internal14:8.31696800E+02,t59:6.18747000E+01,t57:8.99127000E+02,t6:7.20000000E-03)internal12:6.40000000E-03,t46:1.80000000E-03,t56:4.37758900E+02,t1:5.30000000E-03)internal11:9.30000000E-03,t62:3.72852700E+02)internal1:7.70000000E-03;",
                "(t33:6.00000000E-03,t39:6.50000000E-03,t18:9.45069000E+01,t35:2.90000000E-03,((t19:5.90000000E-03,t62:3.49821500E+02,t55:7.00000000E-03,(t41:3.10000000E-03,t14:8.00000000E-03,(t1:7.40000000E-03,t36:5.30000000E-03,t16:4.30000000E-03,t38:3.70000000E-03,t25:2.14193500E+02,t59:4.90000000E-03,t34:1.37566300E+02)internal5:6.80000000E-03)internal4:6.00000000E-03,t61:8.44428600E+02,t50:7.70875000E+02)internal3:8.80000000E-03,(((t2:1.70000000E-03,t48:7.80000000E-03,t22:1.70000000E-03)internal8:1.40000000E-03,t28:4.20000000E-03,t37:2.40000000E-03,t47:2.90000000E-03,(t60:3.10000000E-03,t30:4.70000000E-03)internal9:4.50000000E-03,t27:3.70000000E-03)internal7:3.36765100E+02,t12:9.90000000E-03,(t31:5.30000000E-03,t21:8.20000000E-03,(t3:6.41179400E+02,t5:5.60000000E-03,t45:8.45892500E+02,t53:2.90000000E-03)internal11:1.40000000E-03,(t23:1.50000000E-03,t54:9.50000000E-03,t20:4.51128600E+02,t4:3.00000000E-04,t64:5.02007700E+02,t56:9.94525400E+02,(t58:8.40000000E-03,t13:8.00000000E-03,t26:2.86880200E+02,t11:7.00000000E-03,t57:7.05536500E+02,t44:1.00000000E-04,t42:2.55924700E+02,t46:5.50000000E-03)internal13:3.00000000E-04)internal12:2.80000000E-03,t40:5.30000000E-03)internal10:2.90000000E-03)internal6:8.90000000E-03,t7:2.88102300E+02,(t52:8.00000000E-03,(t43:6.00000000E-03,t24:5.15376300E+02)internal15:1.20000000E-03,(t49:6.20000000E-03,t29:4.16568200E+02,t63:9.00000000E-04,t32:7.00000000E-04)internal16:1.43577700E+02,t9:5.10000000E-03,t6:5.60000000E-03,t8:4.80000000E-03,t15:9.90000000E-03)internal14:5.00000000E-04,t51:1.62433000E+01)internal2:6.60000000E-03,t10:3.50000000E-03,t17:6.30000000E-03)internal1:4.08770000E+02;",
                "(t35:5.00000000E-04,((t60:7.00000000E-03,t47:8.20000000E-03,t58:7.76190700E+02,t26:3.00000000E-03,t9:1.99470700E+02)internal3:7.79142900E+02,t3:9.33371100E+02,((t13:2.00000000E-03,t25:7.90000000E-03)internal5:2.60000000E-03,t1:3.40000000E-03,((t6:3.90000000E-03,t32:5.00000000E-03,(t53:8.10000000E-03,t7:8.00000000E-04,t64:8.30000000E-03,t34:3.30000000E-03,t18:6.60919600E+02,t23:6.10000000E-03,t30:6.90000000E-03,t33:5.00000000E-04)internal8:9.30000000E-03,t55:3.60000000E-03,(t48:9.50000000E-03,t20:6.37513000E+02,t12:8.10000000E-03,t4:9.60000000E-03,t38:6.07573300E+02,t28:4.00000000E-04,t36:7.78515600E+02)internal9:2.80000000E-03)internal7:7.50000000E-03,t8:3.90000000E-03,((t57:7.70000000E-03,t40:5.30000000E-03,t52:6.50000000E-03,t31:1.13449500E+02,t43:5.00000000E-03,t39:5.70000000E-03)internal11:9.50000000E-03,t59:2.60290600E+02,(t37:9.16479000E+02,t16:1.30000000E-03,t27:4.60000000E-03,(t44:6.81929600E+02,t41:3.20000000E-03,t10:7.00000000E-04)internal13:4.50000000E-03,t50:9.00000000E-03)internal12:1.00000000E-02,t24:7.20000000E-03,t63:6.40000000E-03)internal10:6.10000000E-03,(t5:2.00981200E+02,t56:7.90373300E+02,t61:5.54021900E+02,t29:8.10000000E-03)internal14:5.10000000E-03,t21:6.20000000E-03,t49:3.09291500E+02)internal6:3.68019900E+02)internal4:9.90000000E-03,((t17:8.00000000E-04,t42:6.50000000E-03,t62:6.90000000E-03,t15:8.80000000E-03)internal16:7.70000000E-03,(t19:9.66264700E+02,t2:9.20000000E-03,t51:1.24162000E+02,(t22:7.10000000E-03,t14:9.50000000E-03,(t45:6.35057500E+02,t54:1.30000000E-03,t46:5.40000000E-03)internal19:8.50000000E-03)internal18:1.64658000E+02)internal17:6.40000000E-03,t11:1.70000000E-03)internal15:7.20000000E-03)internal2:8.20000000E-03)internal1:7.40000000E-03;",
                ):
            for rng in (MockRandom(), None):
                self.verify_resolve_polytomies(tree_string, rng)

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

    def test_reroot_at_midpoint(self):
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

