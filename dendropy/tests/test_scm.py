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
from dendropy import get_logger
from dendropy.taxa import TaxaBlock
from dendropy.splits import encode_splits
from dendropy.dataio import trees_from_newick
from dendropy.tests import is_test_enabled, TestLevel
from dendropy.treedists import symmetric_difference
from dendropy.tests.debugging_random import DebuggingRandom
from dendropy.treegen import randomly_reorient_tree
_LOG = get_logger("StrictConsensusMerger")

### MODULE THAT WE ARE TESTING ###
from dendropy.scripts.strict_consensus_merge import strict_consensus_merge
### MODULE THAT WE ARE TESTING ###  

_counter = 0
class SCMTest(unittest.TestCase):
    def kernelOfTest(self, trees):
        expected = trees[-1]
        input = trees[:-1]
        output = strict_consensus_merge(input)
        encode_splits(output)
        encode_splits(expected)
        if symmetric_difference(expected, output) != 0:
            self.fail("\n%s\n!=\n%s" % (str(output), str(expected)))

    def testOrderDependent(self):
        o = ['(1,5,(2,(3,4)))', '(2,4,(3,(6,7)))', '(3,4,(6,(7,8)))']
        n = [o[0], o[2], o[1], '(1,2,3,4,5,6,7,8)']
        dataset = trees_from_newick(n)
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)
    
        expected = '(1,5,(2,((3,(6,(7,8))),4)))'
        dataset = trees_from_newick(o + [expected])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)
        
        o.reverse()
        dataset = trees_from_newick(o + [expected])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)

    def testConflict(self):
        o = ['(1,5,(2,((3,6),4)))', '(2,1,(3,(6,4)))', ]
        m = [o[0], o[1], '(1,5,(2,(3,6,4)))']
        n = list(m)
        dataset = trees_from_newick(n, taxa_block=TaxaBlock([str(i) for i in xrange(1,7)]))
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)
        
        rng = DebuggingRandom()
        for i in xrange(50):
            n = list(m)
            dataset = trees_from_newick(n, taxa_block=TaxaBlock([str(i) for i in xrange(1,7)]))
            trees = [i[0] for i in dataset.trees_blocks]
            for t in trees:
                randomly_reorient_tree(t, rng=rng)
            self.kernelOfTest(trees)

        
        
        o = ['(1,5,(3,((2,6),4)))', '(2,1,(3,(6,4)))', ]
        n = [o[0], o[1], '((1,5),2,3,6,4)']
        dataset = trees_from_newick(n, taxa_block=TaxaBlock([str(i) for i in xrange(1,7)]))
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)


    def testMultiEdgeCollision(self):
        dataset = trees_from_newick([
            '(1,2,(3,(4,(5,6))));',
            '(1,2,(3,(7,(8,6))));',
            '(1,2,(3,(4,5,6,7,8)));',
            ])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)

    def testInsertPath(self):
        dataset = trees_from_newick([
            '(((1,2),3),4,5);',
            '(1,2,(3,(7,(8,(9,(4,5))))));',
            '(1,2,(3,(7,(8,(9,(4,5))))));',
            ])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)

    def testPolytomy(self):
        dataset = trees_from_newick([
            '(Athrotaxi,(Liriodchi,Nelumbo2),Sagittari2);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,Liriodchi,Nelumbo2,Sagittari2,Basichlsac,Lamprothma,Mougeotisp,Haplomitr2,Petalaphy,Angiopteri,Azollacaro,Dennstasam,Oleandrapi,Polypodapp,Dicksonant,Vittarifle,Botrychbit,Isoetesmel,Agathismac,Agathisova,Pseudotsu,Libocedrus,Juniperusc,Callitris,Nelumbo,Sagittari,Thuidium);',
            ])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)

    def dofour_five_compat(self, four_taxon_newick, five_taxon_newick):
        #sys.stdout.write("\n4 taxon:%s\n" % four_taxon_newick)
        #sys.stdout.write("5 taxon:%s\n" % five_taxon_newick)
        dataset = trees_from_newick([
            five_taxon_newick,
            four_taxon_newick,
            five_taxon_newick
            ])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)
        # make sure that the behavior is not order dependent
        dataset = trees_from_newick([
            four_taxon_newick,
            five_taxon_newick,
            five_taxon_newick
            ])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)
        
    def testSimple(self):
        if not is_test_enabled(TestLevel.SLOW, _LOG, module_name=__name__, message="skipping all rotation scm tests"):
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
                    four_taxon_newick = '(%s)' % ','.join(base)
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
                                        five_taxon_newick = '(%s)' % ','.join(c)
                                        self.dofour_five_compat(four_taxon_newick, five_taxon_newick)
                                elif n == 1:
                                    for p in range(3):
                                        c = copy.copy(base)
                                        sisters = [base[p], 'E']
                                        for q in [0,1]:
                                            c[p] = '(%s,%s)' % (sisters[q], sisters[1-q])
                                            five_taxon_newick = '(%s)' % ','.join(c)
                                            self.dofour_five_compat(four_taxon_newick, five_taxon_newick)
                                elif n == 2:
                                    sc = copy.copy(clades)
                                    for p in range(3):
                                        upc = copy.copy(four_sisters)
                                        upc.insert(p, 'E')
                                        sc[2] = '(%s)' % ','.join(upc)
                                        nsc = [sc[ii], sc[jj], sc[kk]]
                                        five_taxon_newick = '(%s)' % ','.join(nsc)
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
                                                five_taxon_newick = '(%s)' % ','.join(nc)
                                                self.dofour_five_compat(four_taxon_newick, five_taxon_newick)
                                

    def testThree(self):
        o = [
            '(Athrotaxi,(Liriodchi,Nelumbo),Sagittari);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,((((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend),Nelumbo),Liriodchi),Sagittari);',
            '(Athrotaxi,((((Liriodchi,Annona),Gyrocarpu),Illicium),Nelumbo),((((Ravenala,Calathea),Tacca),Calochort),Sagittari));',
            ]
        expected = '((Athrotaxi,(Callitris,(Juniperusc,Libocedrus))),(((((((Basichlsac,(Mougeotisp,Lamprothma)),Thuidium),(Petalaphy,Haplomitr2)),((Botrychbit,(Vittarifle,((Dicksonant,((Polypodapp,Oleandrapi),Dennstasam)),Azollacaro))),Angiopteri)),Isoetesmel),((Sagittari,(Calochort,(Tacca,(Calathea,Ravenala)))),((Nelumbo,((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend)),(((Liriodchi,Annona),Gyrocarpu),Illicium)))),(Pseudotsu,(Agathisova,Agathismac))));'
        n = o + [expected]
        dataset = trees_from_newick(n)
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)
        o.reverse()
        n = o + [expected]
        dataset = trees_from_newick(n)
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees)
    
if __name__ == "__main__":
    unittest.main()

