#! /usr/bin/env python

############################################################################
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
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)
import dendropy
from dendropy import coalescent

class CalcIntervalsTest(unittest.TestCase):

    def testSimple1(self):

        t = dendropy.Tree.get_from_string("((((a:1, b:1):1, c:2):1, d:3, e:3):2, (f:4, g:4):1)", "newick")
        i1 = coalescent.node_waiting_time_pairs(t)
        assert [x[1] for x in i1] == [1.0, 1.0, 1.0, 1.0, 1.0]
        i2 = coalescent.extract_coalescent_frames(t)
        assert i2 == {7: 1.0, 6:1.0, 5:1.0, 3:1.0, 2:1.0}
        check = coalescent.log_probability_of_coalescent_tree(t, 10)


#    if coalescent.de_hoon_statistics:
#        def testKLDiv(self):
#            from dendropy import treegen
#            taxa_block = treegen.random_taxa_block(100)
#            ctrees = [treegen.pure_kingman(taxa_block, 20000) for i in range(10)]
#            _LOG.info("KL divergence from Kingman coalescent of trees generated under pure Kingman process: %s" % coalescent.kl_divergence_coalescent_trees(ctrees, 20000))
#            ytrees = [treegen.uniform_pure_birth(taxa_block) for i in range(10)]
#            _LOG.info("KL divergence from Kingman coalescent of trees generated under Yule process: %s" % coalescent.kl_divergence_coalescent_trees(ytrees, 20000))

if __name__ == "__main__":
    unittest.main()

