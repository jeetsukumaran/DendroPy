#! /usr/bin/env python

############################################################################
##  test_paup.py
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
Tests paup wrapper
"""

from dendropy import taxa

import unittest
from dendropy import get_logger
import dendropy.tests
_LOG = get_logger("TreeGenerationAndSimulation")

### MODULE THAT WE ARE TESTING ###
from dendropy import paup
### MODULE THAT WE ARE TESTING ###


class TreeDistTest(unittest.TestCase):

    def check_taxa_block(self, filename, taxlabels):
        p = paup.Paup()
        commands = []             
        commands.extend(p.compose_load_trees([dendropy.tests.data_source_path(filename)]))
        commands.extend(p.compose_list_taxa())
        print commands
        results = p.run(commands) 
        taxa_block = p.parse_taxa_block(results)
        assert len(taxa_block) == len(taxlabels)
        for i, t in enumerate(taxa_block):
            assert t.label == taxlabels[i]
                       
    def testTaxaBlock(self):
        test_cases = (
            ("feb032009.tre",("T01", "T02", "T03", "T04", "T05", "T06",
            "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14",
            "T15", "T16", "T17", "T18", "T19", "T20", "T21", "T22",
            "T23", "T24", "T25", "T26", "T27", "T28", "T29", "T30",
            "T31", "T32", "T33", "T34", "T35", "T36", "T37", "T38",
            "T39", "T40", "T41", "T42", "T43", "T44", "T45", "T46",
            "T47", "T48", "T49", "T50", "T51", "T52", "T53", "T54",
            "T55", "T56", "T57", "T58", "T59")),
            ("primates.chars.nexus", ("Lemur_catta", "Homo_sapiens",
            "Pan", "Gorilla", "Pongo", "Hylobates", "Macaca_fuscata",
            "Macaca_mulatta", "Macaca_fascicularis", "Macaca_sylvanus",
            "Saimiri_sciureus", "Tarsius_syrichta", ))
        )
        
        for i in test_cases:
            self.check_taxa_block(i[0], i[1])

if __name__ == "__main__":
    unittest.main()
    