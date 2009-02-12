#! /usr/bin/env python

############################################################################
##  test_char_gen.py
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
Test character generation
"""

import unittest
import subprocess
import re
import os
import sys

from dendropy import dataio
from dendropy import get_logger
_LOG = get_logger("CharacterSimulationAndGeneration")

from dendropy.tests import paup, is_test_enabled, TestLevel
from dendropy import chargen

tree_model_string = """
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

if is_test_enabled(TestLevel.SLOW, _LOG, __name__, "skipping all sequence generation frequency checking"):

    class CharGenTest(unittest.TestCase):
    
        def setUp(self):
            source_ds = dataio.from_nexus(string=tree_model_string)
            self.tree_model = source_ds.trees_blocks[0][0]
       
        def estimate_params(self,
            seq_len=10000,
            kappa=1.0,
            base_freqs=[0.25, 0.25, 0.25, 0.25],
            unequal_base_freqs=True,
            gamma_rates=False,
            prop_invar=False):
            
            output_ds = chargen.generate_hky_dataset(seq_len, 
                tree_model=self.tree_model,
                kappa=kappa,
                base_freqs=base_freqs)
            source_ds = dataio.from_nexus(string=tree_model_string)
            tree_model = source_ds.trees_blocks[0][0]
            
            mle = paup.estimate_char_model(
                tree_model=tree_model,
                char_block=output_ds.char_blocks[0],
                num_states=2,
                unequal_base_freqs=unequal_base_freqs,
                gamma_rates=gamma_rates,
                prop_invar=prop_invar) 
                
            return mle
                        
        def testCharGen(self):
            kappas = [float(i)/10 for i in range(10, 200, 20)]
            base_freq_sets = [
                [0.25, 0.25, 0.25, 0.25],
                [0.4, 0.1, 0.4, 0.1]
            ]
            seq_len = 1000
            for kappa in kappas:
                for bf in base_freq_sets:
                    mle = self.estimate_params(kappa=kappa, base_freqs=bf)
                    kappa_comp = "True = %f, Estimated = %f" % (kappa, mle['kappa'])
                    #kappa_threshold = 10000.00 / float(seq_len)
                    kappa_threshold = kappa * 0.2
                    _LOG.info("Kappa: %s [error threshold: %f]" % (kappa_comp, kappa_threshold))
                    assert abs(kappa-mle['kappa']) < kappa_threshold, \
                                "Estimate over threshold (%f): %s" % (kappa_threshold, kappa_comp)
                                
    if __name__ == "__main__":
        #@ Mark, if the slow test skip condition test goes here, 
        # it only gets checked when the module is executed directly; which
        # means the condition does not get checked when invoked by setuptools
        unittest.main()
