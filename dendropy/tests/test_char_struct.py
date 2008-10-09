#! /usr/bin/env python

############################################################################
##  test_char_struct.py
##
##  Part of the DendroPy phylogenetic computation library.
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
Tests character block manipulations.
"""

import unittest
from dendropy import get_logger
from dendropy import get_logging_level
import dendropy.tests
from dendropy import datasets
_LOG = get_logger("Character Block Structure")

class CharStructTest(unittest.TestCase):
    
    def testCharBlockMerger(self):
        ds1 = datasets.Dataset()
        tb1 = ds1.add_taxa_block(label="Dataset 1, Taxa Block 1")
        for i in range(1,11):
            tb1.add_taxon(label="T%02d" % i)
        ds1 = datasets.Dataset()
        tb1 = ds1.add_taxa_block(label="Dataset 2, Taxa Block 1")
        for i in range(1,21):
            tb1.add_taxon(label="T%02d" % i)           
                   
if __name__ == "__main__":
    unittest.main()
