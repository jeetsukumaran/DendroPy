#! /usr/bin/env python

############################################################################
##  test_popgensim.py
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
Test character generation under pop gen models.
"""

import unittest

from dendropy import datasets
from dendropy import get_logger
_LOG = get_logger("PopulationGeneticSimulatorTests")

from dendropy import GLOBAL_RNG
from dendropy import utils
from dendropy.tests import paup

### MODULE THAT WE ARE TESTING ###
from dendropy import popgensim
### MODULE THAT WE ARE TESTING ### 

class PopGenSimTest(unittest.TestCase):

    def testFragmentedPopulations(self):
        fp = popgensim.FragmentedPopulations(div_time_gens=400,
                                             num_desc_pops = 2,
                                             mutrate_per_gene_per_generation=10e-5,                  
                                             desc_pop_size=10000,
                                             rng=GLOBAL_RNG)
        ds1 = fp.generate_sequences('A', 10, 4000)
        ds2 = fp.generate_sequences('A', 10, 4000, use_seq_gen=False)
        
if __name__ == "__main__":
    unittest.main()        
    