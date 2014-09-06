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
from dendropy.test.support import pathmap
from dendropy.utility import messaging
from dendropy.test.support.dendropytest import ExtendedTestCase
from dendropy.test.support.mockrandom import MockRandom
import dendropy
from dendropy.calculate import treesplit
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
        treesplit.encode_splits(ref)
        treesplit.encode_splits(changing)
        orig_root = changing.seed_node
        for i in range(50):
            changing.randomly_rotate(rng=rng)
            self.assertNotEqual(str(changing), n)
            self.assertEqual(orig_root, changing.seed_node)
            changing._debug_check_tree(logger_obj=_LOG, splits=True)
            if treecompare.symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

class RandomlyReorientTest(unittest.TestCase):

    def runTest(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        k = dendropy.TreeList.get_from_stream(StringIO(n), schema="newick")[0]
        trees = dendropy.TreeList.get_from_stream(StringIO(n+n), schema="newick", encode_splits=True, taxon_namespace=k.taxon_namespace)
        ref = trees[0]
        changing = trees[1]
        rng = MockRandom()
        for i in range(50):
            changing.randomly_reorient(rng=rng, update_splits=True)
            self.assertNotEqual(str(changing), n)
            changing._debug_check_tree(logger_obj=_LOG, splits=True)
            if treecompare.symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

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
        treesplit.encode_splits(tree)
        all_cm = tree.seed_node.edge.split_bitmask
        split_to_target = 0xA
        tree.seed_node.collapse_conflicting(split_to_target, all_cm)
        treesplit.encode_splits(tree)
        treesplit.encode_splits(expected_tree)
        self.assertEqual(treecompare.symmetric_difference(tree, expected_tree), 0)

        tree = tree_list[2]
        expected_tree = tree_list[3]
        treesplit.encode_splits(tree)
        all_cm = tree.seed_node.edge.split_bitmask
        split_to_target = 0x3
        tree.seed_node.collapse_conflicting(split_to_target, all_cm)
        treesplit.encode_splits(tree)
        treesplit.encode_splits(expected_tree)
        self.assertEqual(treecompare.symmetric_difference(tree, expected_tree), 0)

        tree = tree_list[4]
        expected_tree = tree_list[5]
        treesplit.encode_splits(tree)
        all_cm = tree.seed_node.edge.split_bitmask
        split_to_target = 0x5
        tree.seed_node.collapse_conflicting(split_to_target, all_cm)
        treesplit.encode_splits(tree)
        treesplit.encode_splits(expected_tree)
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
            test_tree.reroot_at_midpoint(update_splits=True)
            self.assertEqual(treecompare.symmetric_difference(test_tree, expected_tree), 0)
            for split in test_tree.split_edges:
                if test_tree.split_edges[split].head_node is test_tree.seed_node:
                    continue
                self.assertAlmostEqual(test_tree.split_edges[split].length, expected_tree.split_edges[split].length, 3)

if __name__ == "__main__":
    unittest.main()

