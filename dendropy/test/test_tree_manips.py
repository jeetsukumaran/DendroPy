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
Tests of tree structural manipulations.
"""

from cStringIO import StringIO
import unittest

import random
import dendropy
from dendropy import treecalc
from dendropy import treemanip
from dendropy import treesplit
from dendropy.test.support.datagen import RepeatedRandom
from dendropy.test.support import pathmap
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

class ScaleTest(unittest.TestCase):

    def testScaleEdgesNoLens(self):
        newick_list = ['(5,((4,3),2),1);',
            '(5,(4,3,2),1);',
            '(5,((4,3),2),1);',
            '(5,(4,3),2,1);',
            '(5,((4,3),2),1);',
            '(5,4,3,2,1);']
        tree_list = dendropy.TreeList(
                        stream=StringIO("""%s""" % "\n".join(newick_list)),
                        schema="newick")
        for n, tree in enumerate(tree_list):
            treemanip.scale_edges(tree, 2.0)
            self.assertEqual(newick_list[n], "%s;" % tree.as_newick_string())
    def testScaleEdgesRealTest(self):
        newick_list = ['(5:3,((4:1,3:1):1.5,2:3),1:0);',
            '(5:7.5,4:1,3:-2,2:4,1:.1);']
        doubled = ['(5:6.0,((4:2.0,3:2.0):3.0,2:6.0),1:0.0);',
                    '(5:15.0,4:2.0,3:-4.0,2:8.0,1:0.2);']
        as_f = ['(5:3.0,((4:1.0,3:1.0):1.5,2:3.0),1:0.0);',
            '(5:7.5,4:1.0,3:-2.0,2:4.0,1:0.1);']
        tree_list = dendropy.TreeList(
                        stream=StringIO("""%s""" % "\n".join(newick_list)),
                        schema="newick")
        for n, tree in enumerate(tree_list):
            treemanip.scale_edges(tree, 2)
            self.assertEqual(doubled[n], "%s;" % tree.as_newick_string())
        for n, tree in enumerate(tree_list):
            treemanip.scale_edges(tree, .5)
            self.assertEqual(as_f[n], "%s;" % tree.as_newick_string())


class RandomlyRotateTest(unittest.TestCase):

    def runTest(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        trees = dendropy.TreeList(stream=StringIO(n+n), schema="newick")
        ref = trees[0]
        changing = trees[1]
        rng = RepeatedRandom()
        treesplit.encode_splits(ref)
        treesplit.encode_splits(changing)
        orig_root = changing.seed_node
        for i in xrange(50):
            treemanip.randomly_rotate(changing, rng=rng)
            self.assertNotEqual(str(changing), n)
            self.assertEqual(orig_root, changing.seed_node)
            changing.debug_check_tree(logger_obj=_LOG, splits=True)
            if treecalc.symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

class RandomlyReorientTest(unittest.TestCase):

    def runTest(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        k = dendropy.TreeList(stream=StringIO(n), schema="newick")[0]
        trees = dendropy.TreeList(stream=StringIO(n+n), schema="newick", encode_splits=True, taxon_set=k.taxon_set)
        ref = trees[0]
        changing = trees[1]
        rng = RepeatedRandom()
        for i in xrange(50):
            treemanip.randomly_reorient_tree(changing, rng=rng, splits=True)
            self.assertNotEqual(str(changing), n)
            changing.debug_check_tree(logger_obj=_LOG, splits=True)
            if treecalc.symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

class CollapseConflictingTest(unittest.TestCase):

    def runTest(self):

        taxon_set = dendropy.TaxonSet([str(i+1) for i in range(5)])
        tree_list = dendropy.TreeList(
            stream=StringIO("""
            (5,((4,3),2),1);
            (5,(4,3,2),1);
            (5,((4,3),2),1);
            (5,(4,3),2,1);
            (5,((4,3),2),1);
            (5,4,3,2,1);
            """),
            schema="newick",
            taxon_set=taxon_set)
        tree = tree_list[0]
        expected_tree = tree_list[1]
        treesplit.encode_splits(tree)
        all_cm = tree.seed_node.edge.split_bitmask
        split_to_target = 0xA
        treemanip.collapse_conflicting(tree.seed_node, split_to_target, all_cm)
        treesplit.encode_splits(tree)
        treesplit.encode_splits(expected_tree)
        self.assertEqual(treecalc.symmetric_difference(tree, expected_tree), 0)

        tree = tree_list[2]
        expected_tree = tree_list[3]
        treesplit.encode_splits(tree)
        all_cm = tree.seed_node.edge.split_bitmask
        split_to_target = 0x3
        treemanip.collapse_conflicting(tree.seed_node, split_to_target, all_cm)
        treesplit.encode_splits(tree)
        treesplit.encode_splits(expected_tree)
        self.assertEqual(treecalc.symmetric_difference(tree, expected_tree), 0)

        tree = tree_list[4]
        expected_tree = tree_list[5]
        treesplit.encode_splits(tree)
        all_cm = tree.seed_node.edge.split_bitmask
        split_to_target = 0x5
        treemanip.collapse_conflicting(tree.seed_node, split_to_target, all_cm)
        treesplit.encode_splits(tree)
        treesplit.encode_splits(expected_tree)
        self.assertEqual(treecalc.symmetric_difference(tree, expected_tree), 0)

class PruneTest(unittest.TestCase):

    def check(self,
            title,
            src_prefix,
            to_retain=False):
        input_ds = dendropy.DataSet.get_from_path(
                src=pathmap.tree_source_path(src_prefix + ".pre-pruned.nex"),
                schema='nexus',
                attach_taxon_set=True)
        input_taxa = input_ds.taxon_sets[0]
        output_ds = dendropy.DataSet.get_from_path(
                src=pathmap.tree_source_path(src_prefix + ".paup-pruned.nex"),
                schema='nexus',
                attach_taxon_set=True,
                taxon_set=input_taxa)
        if to_retain:
            taxf = open(pathmap.tree_source_path(src_prefix + ".retained_taxa.txt"), "rU")
        else:
            taxf = open(pathmap.tree_source_path(src_prefix + ".pruned_taxa.txt"), "rU")
        rows = taxf.readlines()
        taxon_idxs_list = [ [int(i) for i in row.split()] for row in rows ]
        for set_idx, src_trees in enumerate(input_ds.tree_lists):
            src_trees = input_ds.tree_lists[set_idx]
            ref_trees = output_ds.tree_lists[set_idx]
            taxon_idxs = taxon_idxs_list[set_idx]
            sub_taxa = [src_trees.taxon_set[i] for i in taxon_idxs]
            for tree_idx, src_tree in enumerate(src_trees):
                _LOG.debug("%s Set %d/%d, Tree %d/%d" % (title, set_idx+1, len(input_ds.tree_lists), tree_idx+1, len(src_trees)))
                ref_tree = ref_trees[tree_idx]
                if to_retain:
                    src_tree.retain_taxa(sub_taxa)
                else:
                    src_tree.prune_taxa(sub_taxa)
                # tree_dist = paup.symmetric_difference(src_tree, ref_tree)
                self.assertEqual(src_tree.symmetric_difference(ref_tree), 0)

    def testPruneTaxaUnrooted(self):
        self.check("Unrooted", "prune_unrooted", False)

    def testPruneTaxaRooted(self):
        self.check("Rooted", "prune_rooted", False)

    def testRetainTaxaUnrooted(self):
        self.check("Unrooted", "prune_unrooted", True)

    def testRetainTaxaRooted(self):
        self.check("Rooted", "prune_rooted", True)

if __name__ == "__main__":
    unittest.main()
