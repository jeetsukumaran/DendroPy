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

import unittest
from dendropy import get_logger
from dendropy.dataio import trees_from_newick
import dendropy.tests
_LOG = get_logger("StrictConsensusMerger")

### MODULE THAT WE ARE TESTING ###
from dendropy.scripts import strict_consensus_merge
### MODULE THAT WE ARE TESTING ###  


class CalcTreeTest(unittest.TestCase):
    def testOne(self):
        return
        trees = trees_from_newick([
            '(Athrotaxi,(Liriodchi,Nelumbo),Sagittari);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,((((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend),Nelumbo),Liriodchi),Sagittari);',
            '(Athrotaxi,((((Liriodchi,Annona),Gyrocarpu),Illicium),Nelumbo),((((Ravenala,Calathea),Tacca),Calochort),Sagittari));',
            '((Athrotaxi,(Callitris,(Juniperusc,Libocedrus))),(((((((Basichlsac,(Mougeotisp,Lamprothma)),Thuidium),(Petalaphy,Haplomitr2)),((Botrychbit,(Vittarifle,((Dicksonant,((Polypodapp,Oleandrapi),Dennstasam)),Azollacaro))),Angiopteri)),Isoetesmel),((Sagittari,(Calochort,(Tacca,(Calathea,Ravenala)))),((Nelumbo,((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend)),(((Liriodchi,Annona),Gyrocarpu),Illicium)))),(Pseudotsu,(Agathisova,Agathismac))));'
            ])
        self.kernelOfTest(trees)

    def kernelOfTest(self, trees):
        expected = trees[-1]
        input = trees[:-1]
        output = strict_consensus_merge(input)
        self.assertEqual(symmetric_difference(expected, output), 0)

    def testTwo(self):
        return
        trees = trees_from_newick([
            '(Athrotaxi,(Liriodchi,Nelumbo),Sagittari);',
            '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));',
            '(Athrotaxi,((((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend),Nelumbo),Liriodchi),Sagittari);',
            '(Athrotaxi,((((Liriodchi,Annona),Gyrocarpu),Illicium),Nelumbo),((((Ravenala,Calathea),Tacca),Calochort),Sagittari));',
            '((Athrotaxi,(Callitris,(Juniperusc,Libocedrus))),(((((((Basichlsac,(Mougeotisp,Lamprothma)),Thuidium),(Petalaphy,Haplomitr2)),((Botrychbit,(Vittarifle,((Dicksonant,((Polypodapp,Oleandrapi),Dennstasam)),Azollacaro))),Angiopteri)),Isoetesmel),((Sagittari,(Calochort,(Tacca,(Calathea,Ravenala)))),((Nelumbo,((((((Verbena,((Thunbergi,Acanthus),(Proboscid,Harpogoph))),Asclepias),Menyanthe),(Phyllonom,(Chamaedap,Pyrola))),((((Mirabilus,Pisum),Circaea),((Rheinward,Octomeles),Greyia)),Dudleya)),Phoradend)),(((Liriodchi,Annona),Gyrocarpu),Illicium)))),(Pseudotsu,(Agathisova,Agathismac))));'
            ])
        self.kernelOfTest(trees)
    
if __name__ == "__main__":
    unittest.main()

