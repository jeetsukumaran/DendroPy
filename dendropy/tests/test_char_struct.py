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
from copy import deepcopy
from dendropy import get_logger
from dendropy import get_logging_level
import dendropy.tests
from dendropy import datasets
_LOG = get_logger("Character Block Structure")

from dendropy import characters

class CharStructTest(unittest.TestCase):
    
    def testCharBlockMerge(self):
        ds1 = datasets.Dataset()
        tb1 = ds1.add_taxa_block(label="Dataset 1, Taxa Block 1")
        for i in range(1,11):
            tb1.add_taxon(label="T%02d" % i)
            
        cb1 = ds1.add_char_block(char_block=characters.DnaCharactersBlock(label="Dataset 2, Taxa Block 1"))
        for t in tb1:
            cb1.append_taxon_sequence(t, state_symbols="AAAAAAAAAA")
            
        ds2 = datasets.Dataset()
        tb2 = ds2.add_taxa_block(label="Dataset 2, Taxa Block 1")
        for i in range(1,21):
            tb2.add_taxon(label="T%02d" % i)  
            
        cb2 = ds2.add_char_block(char_block=characters.DnaCharactersBlock(label="Dataset 2, Taxa Block 1"))
        for t in tb2:
            cb2.append_taxon_sequence(t, state_symbols="CCCCC")            
            
        ds1b = deepcopy(ds1)
        cb = ds1b.char_blocks[0]
        ntax_pre = len(cb)
        nchars_pre = len(cb.values()[0])
        cb.extend_characters(ds2.char_blocks[0])
        self.failIf(len(cb) != ntax_pre, 
                    "Number of taxa have changed after from %d to %d" % (ntax_pre, len(cb)))
        for t in cb:
            _LOG.info("\n%s: %s" \
                % (str(t), cb[t].values_as_string()))
            print cb[t], len(cb[t])
            self.failIf(len(cb[t]) != 15,
                "Data vector is incorrect length (%d):\n%s: %s" \
                % (len(cb[t]), str(t), cb[t].values_as_string()))
            self.failIf(cb[t].values_as_string() != "AAAAAAAAAACCCCC",
                "Incorrect sequence:\n%s: %s" % (str(t), cb[t].values_as_string()))
if __name__ == "__main__":
    unittest.main()
