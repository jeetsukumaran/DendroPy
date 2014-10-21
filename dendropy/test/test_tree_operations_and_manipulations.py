
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
                if bipartition.edge.head_node is test_tree.seed_node:
                    continue
                # self.assertAlmostEqual(bipartition.edge.length, expected_tree.split_bitmask_edge_map[bipartition.split_bitmask].length, 3)
                self.assertAlmostEqual(bipartition.edge.length, expected_tree.bipartition_edge_map[bipartition].length, 3)

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

                    # identify
                    # print("\n### is_rooted = {}, suppress_unifurcations = {}, reseed at = {}".format(is_rooted, suppress_unifurcations, reseed_at_label))

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
                    tree._debug_check_tree(logger_obj=_LOG, check_bipartitions=True)

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

