#! /usr/bin/env python

############################################################################
##  test_utils.py
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
Tests various utilities
"""

import unittest
from dendropy import get_logger

import dendropy.tests
_LOG = get_logger("Utilities")

### MODULES THAT WE ARE TESTING ###
from dendropy import utils
### MODULES THAT WE ARE TESTING ###

class NormalizedBitmaskDict(unittest.TestCase):

    
    def test_complementing(self):
    
        mask = 0xFF # 1111 1111
        splits = [
            ((0x03, '0000 0011'), (0x03, '0000 0011')),
            ((0x34, '0011 0100'), (0xCB, '1100 1011')),
            ((0x44, '0100 0100'), (0xBB, '1011 1011')),
            ((0x12, '0001 0010'), (0xED, '1110 1101')),
            ((0x75, '0111 0101'), (0x75, '0111 0101')),
            ]
        d = utils.NormalizedBitmaskDict(mask=mask)
        for s in splits:
            d[s[0][0]] = s[0][1]
            
        for s in splits:
            assert s[0][0] in d
            assert s[1][0] in d
            assert d[s[0][0]] == d[s[1][0]]

        for k, v in d.items():
            pass
        
        del d[splits[0][0][0]]
        del d[splits[1][1][0]]
        assert splits[0][0][0] not in d
        assert splits[0][1][0] not in d
        assert splits[1][0][0] not in d
        assert splits[1][1][0] not in d
                               
if __name__ == "__main__":
    unittest.main()
        