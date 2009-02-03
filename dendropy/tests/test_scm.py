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
import unittest
from dendropy import get_logger
from dendropy.splits import encode_splits
from dendropy.dataio import trees_from_newick
import dendropy.tests
from dendropy.treedists import symmetric_difference
_LOG = get_logger("StrictConsensusMerger")

### MODULE THAT WE ARE TESTING ###
from dendropy.scripts.strict_consensus_merge import strict_consensus_merge
### MODULE THAT WE ARE TESTING ###  


class SCMTest(unittest.TestCase):
    def testOne(self):
        dataset = trees_from_newick([
            '(Athrotaxi,(Liriodchi,Nelumbo2),Sagittari2);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,Liriodchi,Nelumbo2,Sagittari2,Basichlsac,Lamprothma,Mougeotisp,Haplomitr2,Petalaphy,Angiopteri,Azollacaro,Dennstasam,Oleandrapi,Polypodapp,Dicksonant,Vittarifle,Botrychbit,Isoetesmel,Agathismac,Agathisova,Pseudotsu,Libocedrus,Juniperusc,Callitris,Nelumbo,Sagittari,Thuidium);',
            ])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees, dataset.taxa_blocks[0])

    def kernelOfTest(self, trees, taxa_block):
        if not "TESTING_DENDROPY_SCM" in os.environ:
            return
        expected = trees[-1]
        input = trees[:-1]
        output = strict_consensus_merge(input, taxa_block=taxa_block)
        encode_splits(expected)
        self.assertEqual(symmetric_difference(expected, output), 0)

    def testTwo(self):
        return
        dataset = trees_from_newick([
            '(Athrotaxi,(Liriodchi,Nelumbo),Sagittari);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,((((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend),Nelumbo),Liriodchi),Sagittari);',
            '(Athrotaxi,((((Liriodchi,Annona),Gyrocarpu),Illicium),Nelumbo),((((Ravenala,Calathea),Tacca),Calochort),Sagittari));',
            '((Athrotaxi,(Callitris,(Juniperusc,Libocedrus))),(((((((Basichlsac,(Mougeotisp,Lamprothma)),Thuidium),(Petalaphy,Haplomitr2)),((Botrychbit,(Vittarifle,((Dicksonant,((Polypodapp,Oleandrapi),Dennstasam)),Azollacaro))),Angiopteri)),Isoetesmel),((Sagittari,(Calochort,(Tacca,(Calathea,Ravenala)))),((Nelumbo,((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend)),(((Liriodchi,Annona),Gyrocarpu),Illicium)))),(Pseudotsu,(Agathisova,Agathismac))));'
            ])
        trees = [i[0] for i in dataset.trees_blocks]
        self.kernelOfTest(trees, dataset.taxa_blocks[0])
    
if __name__ == "__main__":
    unittest.main()

