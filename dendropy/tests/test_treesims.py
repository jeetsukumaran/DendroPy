#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

"""
Tests of tree simulation.
"""

import sys
import os
import unittest
import random

from dendropy.utility import messaging
import dendropy
from dendropy import treesim

_LOG = messaging.get_logger(__name__)

class TreeSimTest(unittest.TestCase):

    def get_species_tree(self, ntax=10):
        "Returns a species population tree, testing it in the process."
        ages = [random.randint(1000,10000) for age in range(ntax)]
        ages.sort()
        pop_sizes = [random.randint(1000,10000) for pop in range(2*ntax+1)]
        taxon_set = dendropy.new_taxon_set(ntax)
        species_tree = treesim.pop_gen_tree(taxon_set=taxon_set,
                                                 ages=ages,
                                                 num_genes=4,
                                                 pop_sizes=pop_sizes)
        ages2 = []
        for node in species_tree.postorder_node_iter():
            distance_from_tip = node.distance_from_tip()
            if distance_from_tip > 0:
                ages2.append(distance_from_tip)
        ages2.sort()
        _LOG.debug("-- AGES --")
        _LOG.debug([float(age) for age in ages])
        _LOG.debug([float(age) for age in ages2])
        for index in range(len(ages2)):
            assert (ages[index] - ages2[index]) < 10e-6

        pop_sizes2 = []
        for edge in species_tree.postorder_edge_iter():
            pop_sizes2.append(edge.pop_size)
        pop_sizes2.sort()
        _LOG.debug("-- POP_SIZES --")
        _LOG.debug([float(pop_size) for pop_size in pop_sizes])
        _LOG.debug([float(pop_size) for pop_size in pop_sizes2])

        return species_tree

    def testKingmanTreeInSpeciesTree(self, ntax=10):
        species_tree = self.get_species_tree(ntax)
        _LOG.debug("Generating 20 gene trees conditional on this species tree ...")
        gene_trees = []
        while len(gene_trees) < 20:
            gene_trees.append(treesim.constrained_kingman(species_tree)[0])

    def testPureKingman(self):
        """
        Right now, just making sure we get no errors while running. TODO:
        proper test!
        """
        t = treesim.pure_kingman(dendropy.new_taxon_set(100))
        assert t._debug_tree_is_valid()

if __name__ == "__main__":
    unittest.main()
