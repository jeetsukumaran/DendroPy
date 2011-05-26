#!/usr/bin/env python
############################################################################
##  test_tree_mods.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Tests of the strict consensus merger
"""
import os
import sys
import unittest
import copy
import logging
from cStringIO import StringIO

from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

from dendropy import TaxonSet, TreeList, DataSet
from dendropy.dataobject.tree import format_split
from dendropy.treesplit import encode_splits, count_bits, lowest_bit_only
from dendropy.treemanip import collapse_clade, collapse_edge, randomly_reorient_tree
from dendropy.treecalc import symmetric_difference
from dendropy.test.support.datagen import RepeatedRandom
from dendropy.test.support.runlevel import is_test_enabled
from dendropy.test.support import runlevel

verbose = False
IS_DEBUG_LOGGING = _LOG.isEnabledFor(logging.DEBUG)

# from dendropy.splits import encode_splits
# from dendropy.dataio import trees_from_newick



### MODULE THAT WE ARE TESTING ###
from dendropy.scm import inplace_strict_consensus_merge
### MODULE THAT WE ARE TESTING ###

def trees_from_newick_str_list(newick_list):
    all_tree_str = " ".join(newick_list)
    return TreeList(stream=StringIO(all_tree_str), taxon_set=TaxonSet(), schema="NEWICK")
_counter = 0

class SCMTest(unittest.TestCase):
    def kernelOfTest(self, trees):
        expected = trees[-1]
        input = trees[:-1]
        _LOG.debug('input = %s' % str(input))
        output = inplace_strict_consensus_merge(input)
        encode_splits(output)
        encode_splits(expected)
        if symmetric_difference(expected, output) != 0:
            self.fail("\n%s\n!=\n%s" % (str(output), str(expected)))
    def testThree(self):
        o = [
            '(Athrotaxi,(Liriodchi,Nelumbo),Sagittari);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,((((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend),Nelumbo),Liriodchi),Sagittari);',
            '(Athrotaxi,((((Liriodchi,Annona),Gyrocarpu),Illicium),Nelumbo),((((Ravenala,Calathea),Tacca),Calochort),Sagittari));',
            ]
        expected = '((Athrotaxi,(Callitris,(Juniperusc,Libocedrus))),(((((((Basichlsac,(Mougeotisp,Lamprothma)),Thuidium),(Petalaphy,Haplomitr2)),((Botrychbit,(Vittarifle,((Dicksonant,((Polypodapp,Oleandrapi),Dennstasam)),Azollacaro))),Angiopteri)),Isoetesmel),((Sagittari,(Calochort,(Tacca,(Calathea,Ravenala)))),((Nelumbo,((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend)),(((Liriodchi,Annona),Gyrocarpu),Illicium)))),(Pseudotsu,(Agathisova,Agathismac))));'
        n = o + [expected]
        trees = trees_from_newick_str_list(n)
        self.kernelOfTest(trees)
        return
        o.reverse()
        n = o + [expected]
        trees = trees_from_newick_str_list(n)
        self.kernelOfTest(trees)


    def testPolytomy(self):
        trees = trees_from_newick_str_list([
            '(Athrotaxi,(Liriodchi,Nelumbo2),Sagittari2);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,Liriodchi,Nelumbo2,Sagittari2,Basichlsac,Lamprothma,Mougeotisp,Haplomitr2,Petalaphy,Angiopteri,Azollacaro,Dennstasam,Oleandrapi,Polypodapp,Dicksonant,Vittarifle,Botrychbit,Isoetesmel,Agathismac,Agathisova,Pseudotsu,Libocedrus,Juniperusc,Callitris,Nelumbo,Sagittari,Thuidium);',
            ])
        self.kernelOfTest(trees)

    def dofour_five_compat(self, four_taxon_newick, five_taxon_newick):
        #sys.stdout.write("\n4 taxon:%s\n" % four_taxon_newick)
        #sys.stdout.write("5 taxon:%s\n" % five_taxon_newick)
        trees = trees_from_newick_str_list([
            five_taxon_newick,
            four_taxon_newick,
            five_taxon_newick
            ])
        self.kernelOfTest(trees)
        # make sure that the behavior is not order dependent
        trees = trees_from_newick_str_list([
            four_taxon_newick,
            five_taxon_newick,
            five_taxon_newick
            ])
        self.kernelOfTest(trees)

    def testSimple(self):
        if not is_test_enabled(runlevel.SLOW, _LOG, module_name=__name__, message="skipping all rotation scm tests"):
            return
        clades = ['A', 'D', None, None]
        for m in [0, 1]:
            four_sisters = ['B', 'C']
            clades[2] = '(%s,%s)' % (four_sisters[m], four_sisters[1-m])
            for i in [0, 1, 2]:
                for j in [0, 1, 2]:
                    if i == j:
                        continue
                    k = [0, 1, 2]
                    k.remove(i)
                    k.remove(j)
                    k = k[0]
                    base = [clades[i], clades[j], clades[k]]
                    four_taxon_newick = '(%s);' % ','.join(base)
                    for ii in [0, 1, 2]:
                        for jj in [0, 1, 2]:
                            if ii == jj:
                                continue
                            kk = [0, 1, 2]
                            kk.remove(ii)
                            kk.remove(jj)
                            kk = kk[0]
                            base = [clades[ii], clades[jj], clades[kk]]
                            for n in range(4):
                                if n == 0:
                                    for p in range(4):
                                        c = copy.copy(base)
                                        c.insert(p, 'E')
                                        five_taxon_newick = '(%s);' % ','.join(c)
                                        self.dofour_five_compat(four_taxon_newick, five_taxon_newick)
                                elif n == 1:
                                    for p in range(3):
                                        c = copy.copy(base)
                                        sisters = [base[p], 'E']
                                        for q in [0,1]:
                                            c[p] = '(%s,%s)' % (sisters[q], sisters[1-q])
                                            five_taxon_newick = '(%s);' % ','.join(c)
                                            self.dofour_five_compat(four_taxon_newick, five_taxon_newick)
                                elif n == 2:
                                    sc = copy.copy(clades)
                                    for p in range(3):
                                        upc = copy.copy(four_sisters)
                                        upc.insert(p, 'E')
                                        sc[2] = '(%s)' % ','.join(upc)
                                        nsc = [sc[ii], sc[jj], sc[kk]]
                                        five_taxon_newick = '(%s);' % ','.join(nsc)
                                        self.dofour_five_compat(four_taxon_newick, five_taxon_newick)
                                else:
                                    for p in range(2):
                                        for q in range(2):
                                            upc = copy.copy(four_sisters)
                                            s = upc[p]
                                            ns = [s, 'E']
                                            upc[p] = '(%s,%s)' % (ns[q], ns[1-q])
                                            for r in range(2):
                                                third = '(%s,%s)' % (upc[r], upc[1-r])
                                                nc = [clades[0], clades[1], third]
                                                five_taxon_newick = '(%s);' % ','.join(nc)
                                                self.dofour_five_compat(four_taxon_newick, five_taxon_newick)

    def testConflict(self):
        taxa = TaxonSet([str(i) for i in xrange(1,7)])
        o = ['(1,5,(2,((3,6),4)));', '(2,1,(3,(6,4)));', ]
        m = [o[0], o[1], '(1,5,(2,(3,6,4)));']
        trees = trees_from_newick_str_list(m)
        self.kernelOfTest(trees)
        rng = RepeatedRandom()
        for i in xrange(50):
            trees = trees_from_newick_str_list(m)
            for t in trees:
                randomly_reorient_tree(t, rng=rng)
            self.kernelOfTest(trees)


    def testInsertPath(self):
        trees = trees_from_newick_str_list([
            '(((1,2),3),4,5);',
            '(1,2,(3,(7,(8,(9,(4,5))))));',
            '(1,2,(3,(7,(8,(9,(4,5))))));',
            ])
        self.kernelOfTest(trees)



    def testOrderDependent(self):
        o = ['(1,5,(2,(3,4)));', '(2,4,(3,(6,7)));', '(3,4,(6,(7,8)));']
        n = [o[0], o[2], o[1], '(1,2,3,4,5,6,7,8);']

        trees = trees_from_newick_str_list(n)
        self.kernelOfTest(trees)

        expected = '(1,5,(2,((3,(6,(7,8))),4)));'
        trees = trees_from_newick_str_list(o + [expected])
        self.kernelOfTest(trees)

        o.reverse()
        trees = trees_from_newick_str_list(o + [expected])
        self.kernelOfTest(trees)

        o = ['(1,5,(3,((2,6),4)));', '(2,1,(3,(6,4)));', ]
        n = [o[0], o[1], '((1,5),2,3,6,4);']
        trees = trees_from_newick_str_list(n)
        self.kernelOfTest(trees)

    def testMultiEdgeCollision(self):
        trees = trees_from_newick_str_list([
            '(1,2,(3,(4,(5,6))));',
            '(1,2,(3,(7,(8,6))));',
            '(1,2,(3,(4,5,6,7,8)));',
            ])
        self.kernelOfTest(trees)
if __name__ == "__main__":
    unittest.main()

