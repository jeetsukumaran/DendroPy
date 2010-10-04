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

import dendropy
from dendropy import treecalc
from dendropy import treemanip
from dendropy import treesplit
from dendropy.test.support.datagen import RepeatedRandom
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

    def testPruneNodes(self):
        """NOT IMPLEMENTED YET: PRIORITY TODO!!!"""
        pass

    def testPruneTaxa(self):
        """NOT IMPLEMENTED YET: PRIORITY TODO!!!"""
        pass


if __name__ == "__main__":
    unittest.main()
