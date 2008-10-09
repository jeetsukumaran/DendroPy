#! /usr/bin/env python

############################################################################
##  test_char_gen.py
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
Test character generation
"""

import unittest
import subprocess
import re
import os
import sys

from dendropy import dataio
from dendropy import get_logger
import dendropy.tests
from dendropy.tests import paup

_LOG = get_logger("Splits")

from dendropy import chargen

model_tree_string = """
#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=5;
    TAXLABELS
        A
        B
        C
        D
        E
  ;
END;
begin trees;
    tree true=(A:0.25,(B:0.25,(C:0.25,(D:0.25,E:0.25):0.25):0.25):0.25):0.25;
end;
"""

class CharGenTest(unittest.TestCase):
    
    def testCharGen(self):
        source_ds = dataio.get_nexus(string=model_tree_string)
        tree_model = source_ds.trees_blocks[0][0]
        output_ds = chargen.generate_hky_dataset(10000, tree_model=tree_model)
        
        mle = dendropy.tests.paup.estimate_char_model(
            model_tree=tree_model,
            char_block=output_ds.char_blocks[0],
            num_states=6,
            unequal_base_freqs=True,
            gamma_rates=True,
            prop_invar=True)
            
        _LOG.info(mle)            
            
                    
if __name__ == "__main__":
    unittest.main()
