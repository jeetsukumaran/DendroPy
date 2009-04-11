############################################################################
##  test_coal.py
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
Tests coalescence calculations.
"""

import unittest
from dendropy import get_logger
from dendropy import datasets
_LOG = get_logger("TreeCoal")

### MODULE THAT WE ARE TESTING ###
from dendropy import coalescent
### MODULE THAT WE ARE TESTING ###

class CalcIntervalsTest(unittest.TestCase):

    def testSimple1(self):
        d = datasets.Dataset()
        t = d.trees_from_string("((((a:1, b:1):1, c:2):1, d:3, e:3):2, (f:4, g:4):1)", "newick")[0]
        intervals = coalescent.coalescence_intervals(t)
        assert intervals == [1.0, 1.0, 1.0, 1.0, 1.0], "intervals found = %s" % ", ".join(intervals)

if __name__ == "__main__":
    unittest.main()

