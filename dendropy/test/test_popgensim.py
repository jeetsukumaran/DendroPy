#! /usr/bin/env python

############################################################################
##  test_tree_io.py
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
Tests population genetic statistic calculation.
"""

import unittest
from dendropy.utility import GLOBAL_RNG
from dendropy import popgensim

class PopGenSimTest(unittest.TestCase):

    def testFragmentedPopulations(self):
        """Generates data only: no verification [TODO: Verification!]"""
        fp = popgensim.FragmentedPopulations(div_time_gens=400,
                                             num_desc_pops = 2,
                                             mutrate_per_site_per_generation=10e-8,
                                             desc_pop_size=10000,
                                             rng=GLOBAL_RNG)
        ds2 = fp.generate_sequences('A', 10, 4000, use_seq_gen=False)

if __name__ == "__main__":
    unittest.main()

