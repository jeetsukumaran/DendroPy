#! /usr/bin/env python

############################################################################
##  test_tree_gens.py
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
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Tests tree generation.
"""

import random
import unittest
from dendropy import get_logger
import dendropy.tests
_LOG = get_logger("test_tree_gens")

### MODULE THAT WE ARE TESTING ###
from dendropy import treegens
### MODULE THAT WE ARE TESTING ###

class TreeGenTest(unittest.TestCase):

    def get_species_tree(self, ntax=10):
        """
        Returns a species population tree, testing it in the process.
        """
        ages = [random.randint(1000,10000) for age in range(ntax)]
        ages.sort()
        pop_sizes = [random.randint(1000,10000) for pop in range(2*ntax+1)]
        taxa_block = treegens.random_taxa_block(ntax)
        species_tree = treegens.pop_gen_tree(taxa_block=taxa_block,
                                                 ages=ages,
                                                 num_genes=4,
                                                 pop_sizes=pop_sizes)
        nxfpath = dendropy.tests.test_target_path('species_pop.xml')
        nwfpath = dendropy.tests.test_target_path('species_pop.tre')

        ages2 = []
        for node in species_tree.postorder_node_iter():
            distance_from_tip = node.distance_from_tip()
            if distance_from_tip > 0:
                ages2.append(distance_from_tip)
        ages2.sort()
        _LOG.info("-- AGES --")
        _LOG.info([float(age) for age in ages])
        _LOG.info([float(age) for age in ages2])
        for index in range(len(ages2)):
            self.failIf(ages[index] - ages2[index] > 10e-6)

        pop_sizes2 = []
        for edge in species_tree.postorder_edge_iter():
            pop_sizes2.append(edge.pop_size)
        pop_sizes2.sort()
        _LOG.info("-- POP_SIZES --")
        _LOG.info([float(pop_size) for pop_size in pop_sizes])
        _LOG.info([float(pop_size) for pop_size in pop_sizes2])

        # difficult to test due to uncertainty in order of child visits
#         for index in range(len(pop_sizes2)):
#             self.failIf(pop_sizes[index] - pop_sizes2[index] > 10e-6)

        return species_tree

    def test_kingman_tree_in_species_tree(self, ntax=10):
        species_tree = self.get_species_tree(ntax)
        _LOG.info("Generating 20 gene trees conditional on this species tree ...")
        gene_trees = []
        while len(gene_trees) < 20:
            gene_trees.append(treegens.constrained_kingman(species_tree)[0])


def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(TreeGenTest)


def additional_tests():
    "returns all tests in this file as suite"
    return unittest.TestLoader().loadTestsFromTestCase(TreeGenTest)

def getTestSuite():
    """Alias to the additional_tests().  This is unittest-style.
    `additional_tests` is used by setuptools.
    """
    return additional_tests()

if __name__ == "__main__":
    unittest.main()
