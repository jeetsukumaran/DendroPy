#! /usr/bin/env python

############################################################################
##  test_char_struct.py
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
Tests character block manipulations.
"""

import unittest
from copy import deepcopy
from dendropy import get_logger
from dendropy import get_logging_level
import dendropy.tests
from dendropy import datasets
_LOG = get_logger("CharacterBlockStructure")

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
            cb2.append_taxon_sequence(t, state_symbols="CCCCCCCCCC")            
            
        ds1b = deepcopy(ds1)
        cb = ds1b.char_blocks[0]
        ntax_pre = len(cb)
        nchars_pre = len(cb.values()[0])
        cb.extend_characters(ds2.char_blocks[0])
        assert len(cb) == ntax_pre, \
                    "Number of taxa have changed after from %d to %d" % (ntax_pre, len(cb))
        for t in cb:
            _LOG.debug("\n%s: %s" \
                % (str(t), cb[t].values_as_string()))
            assert len(cb[t]) == 20, \
                "Data vector is incorrect length (%d):\n%s: %s" \
                % (len(cb[t]), str(t), cb[t].values_as_string())
            assert cb[t].values_as_string() == "AAAAAAAAAACCCCCCCCCC", \
                "Incorrect sequence:\n%s: %s" % (str(t), cb[t].values_as_string())
                
        ds1b = deepcopy(ds1)
        cb = ds1b.char_blocks[0]
        cb.extend(ds2.char_blocks[0], overwrite_existing=True)
        target_ntax = 20
        assert len(cb) == target_ntax,  \
                    "Number of rows in character block have not changed to %d (%d)" % (target_ntax, len(cb))
        assert len(cb.taxa_block) == target_ntax, \
                    "Number of taxa in taxa block have not changed to %d (%d)" % (target_ntax, len(cb))
                    
        for t in tb2:
            cb_tb_labels = cb.taxa_block.labels()
            assert t.label in cb_tb_labels, \
                "Taxon '%s' not found in taxa block:\n%s" % (str(t), str(cb_tb_labels))
            cb_labels = [t.label for t in cb]
            assert t.label in cb_labels, \
                "Taxon '%s' not found in char block:\n%s" % (str(t), str(cb_labels))
        for t in cb:
            _LOG.debug("\n%s: %s" \
                % (str(t), cb[t].values_as_string()))                    
            assert len(cb[t]) == 10, \
                "Data vector is incorrect length (%d):\n%s: %s" \
                % (len(cb[t]), str(t), cb[t].values_as_string())
            assert cb[t].values_as_string() == "CCCCCCCCCC", \
                "Incorrect sequence:\n%s: %s" % (str(t), cb[t].values_as_string())
                
        ds1b = deepcopy(ds1)
        cb = ds1b.char_blocks[0]
        cb.extend(ds2.char_blocks[0], append_existing=True)
        target_ntax = 20
        assert len(cb) == target_ntax, \
                    "Number of rows in character block have not changed to %d (%d)" % (target_ntax, len(cb))
        assert len(cb.taxa_block) == target_ntax, \
                    "Number of taxa in taxa block have not changed to %d (%d)" % (target_ntax, len(cb))
        for t in tb2:
            cb_tb_labels = cb.taxa_block.labels()
            assert t.label in cb_tb_labels, \
                "Taxon '%s' not found in taxa block:\n%s" % (str(t), str(cb_tb_labels))
            cb_labels = [t.label for t in cb]
            assert t.label in cb_labels, \
                "Taxon '%s' not found in char block:\n%s" % (str(t), str(cb_labels))
        for t in cb:
            _LOG.debug("\n%s: %s" \
                % (str(t), cb[t].values_as_string()))
            tnum = int(t.label[-2:])
            if tnum > 10:
                assert len(cb[t]) == 10, \
                    "Data vector is incorrect length (%d):\n%s: %s" \
                    % (len(cb[t]), str(t), cb[t].values_as_string())
                assert cb[t].values_as_string() == "CCCCCCCCCC", \
                    "Incorrect sequence:\n%s: %s" % (str(t), cb[t].values_as_string())
            else:
                assert len(cb[t]) == 20, \
                    "Data vector is incorrect length (%d):\n%s: %s" \
                    % (len(cb[t]), str(t), cb[t].values_as_string())
                assert cb[t].values_as_string() == "AAAAAAAAAACCCCCCCCCC", \
                    "Incorrect sequence:\n%s: %s" % (str(t), cb[t].values_as_string())
                                                                           
if __name__ == "__main__":
    unittest.main()
