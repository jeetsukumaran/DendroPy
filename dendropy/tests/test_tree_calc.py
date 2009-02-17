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
Tests optional tree modifiers.
"""

import random
import unittest
import itertools
from dendropy import get_logger
from dendropy.datasets import Dataset
from dendropy.tests.debugging_random import DebuggingRandom
from dendropy.tests.util_for_testing import assert_vec_approx_equal, assert_approx_equal
import dendropy.tests
_LOG = get_logger("TreeModification")

### MODULE THAT WE ARE TESTING ###
from dendropy.treecalc import pybus_harvey_gamma
### MODULE THAT WE ARE TESTING ###


class CalcTreeTest(unittest.TestCase):
    def testPHGamma(self):
        newick = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        d = Dataset()
        tree = d.trees_from_string(string=newick, format="NEWICK")[0]
        assert_approx_equal(pybus_harvey_gamma(tree), 0.546276)

if __name__ == "__main__":
    unittest.main()

