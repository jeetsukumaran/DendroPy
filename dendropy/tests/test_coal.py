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
from dendropy import splits
from dendropy import treesum
_LOG = get_logger("TreeCoal")

### MODULE THAT WE ARE TESTING ###
from dendropy import coalescent
### MODULE THAT WE ARE TESTING ###

class CalcIntervalsTest(unittest.TestCase):

    def testSimple1(self):
        d = datasets.Dataset()
        t = d.trees_from_string("((((a:1, b:1):1, c:2):1, d:3, e:3):2, (f:4, g:4):1)", "newick")[0]
        i1 = coalescent.coalescence_intervals(t)
        assert i1 == [1.0, 1.0, 1.0, 1.0, 1.0], "intervals found = %s" % ", ".join(intervals)
        i2 = coalescent.coalescent_frames(t)
        assert i2 == [(7, 1.0), (6, 1.0), (5, 1.0), (3, 1.0), (2, 1.0)]
        check = coalescent.probability_of_coalescent_tree(t, 10)
#         check2 = coalescent.debug_coal_prob(t, 10)
#         print check, check2
        ### TODO: Actually come up a with a decent coalescent tree, calculated the probability,
        ###       and check if it is equal
                  
class DeepCoalTest(unittest.TestCase):
    
    def setUp(self):
        self.dataset = datasets.Dataset()                    
        self.gene_trees = self.dataset.trees_from_string("""
            [&R] (A,(B,(C,D))); [&R] ((A,C),(B,D)); [&R] (C,(A,(B,D)));
            """, "NEWICK")
        self.species_trees = self.dataset.trees_from_string("""
            [&R] (A,(B,(C,D)));
            [&R] (A,(C,(B,D)));
            [&R] (A,(D,(C,B)));
            [&R] (B,(A,(C,D)));
            [&R] (B,(C,(A,D)));
            [&R] (B,(D,(C,A)));
            [&R] (C,(A,(B,D)));
            [&R] (C,(B,(A,D)));
            [&R] (C,(D,(B,A)));
            [&R] (D,(A,(B,C)));
            [&R] (D,(B,(A,C)));
            [&R] (D,(C,(B,A)));
            [&R] ((A,B),(C,D));
            [&R] ((A,C),(B,D));
            [&R] ((A,D),(C,B));
            """, "NEWICK")
                        
        # expected results, for each gene tree / species tree pairing, with
        # cycling through species trees for each gene tree
        self.expected_deep_coalescences = [ 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 1, 2, 2,
                                            2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 2, 0, 2,
                                            2, 1, 2, 3, 3, 3, 0, 1, 1, 3, 3, 3, 2, 1, 2 ]                                                 
        assert len(self.expected_deep_coalescences) == len(self.gene_trees) * len(self.species_trees)       
                
        ## prep trees ##
        assert len(self.dataset.taxa_blocks) == 1      
        tb = self.dataset.taxa_blocks[0]
        for t in self.gene_trees + self.species_trees:
            assert t.taxa_block == tb
            t.is_rooted = True
            splits.encode_splits(t)        

    def testDeepCoalCounting(self):
        idx = 0
        _LOG.info("Species\t\tGene\t\tDC\t\tExp.DC\t\tDiff")
        for gt in self.gene_trees:
            for st in self.species_trees:
                dc = coalescent.num_deep_coalescences(st, gt)        
                _LOG.info("%s\t\t%s\t\t%s\t\t%s\t\t%s" 
                    % (st.compose_newick(),
                       gt.compose_newick(),
                       dc, 
                       self.expected_deep_coalescences[idx], 
                       dc - self.expected_deep_coalescences[idx]))
                assert dc == self.expected_deep_coalescences[idx]                       
                idx += 1          
 
    
if __name__ == "__main__":
    unittest.main()

