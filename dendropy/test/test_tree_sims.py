#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

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
from dendropy.test.support.datagen import RepeatedRandom
_LOG = messaging.get_logger(__name__)
class BirthDeathTreeTest(unittest.TestCase):
    def testGSABD(self):
        """test that the birth-death process produces the correct number of tips with GSA."""
        _RNG = RepeatedRandom()
        for num_leaves in range(2, 15):
            _LOG.debug("Generating tree with %d leaves" % num_leaves)
            t = treesim.birth_death(birth_rate=1.0, death_rate=0.2, ntax=num_leaves, gsa_ntax=3*num_leaves, rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEquals(num_leaves, len(t.leaf_nodes()))

    def testYule(self):
        """test that the pure-birth process produces the correct number of tips."""
        _RNG = RepeatedRandom()
        for num_leaves in range(2, 20):
            t = treesim.birth_death(birth_rate=1.0, death_rate=0.0, ntax=num_leaves, rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEquals(num_leaves, len(t.leaf_nodes()))

    def testGSA(self):
        """test that the pure-birth process produces the correct number of tips with GSA."""
        _RNG = RepeatedRandom()
        for num_leaves in range(2, 20):
            t = treesim.birth_death(birth_rate=1.0, death_rate=0.0, ntax=num_leaves, gsa_ntax=4*num_leaves, rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEquals(num_leaves, len(t.leaf_nodes()))

    def testBDTree(self):
        """PureCoalescentTreeTest -- tree generation without checking [TODO: checks]"""
        _RNG = RepeatedRandom()
        for num_leaves in range(2, 20):
            t = treesim.birth_death(birth_rate=1.0, death_rate=0.2, ntax=num_leaves, rng=_RNG)
            self.assertTrue(t._debug_tree_is_valid())
            self.assertEquals(num_leaves, len(t.leaf_nodes()))

    

class TruncatedCoalescentTreeTest(unittest.TestCase):

    def get_species_tree(self, ntax=10):
        _RNG = RepeatedRandom()
        ages = [random.randint(1000,10000) for age in range(ntax)]
        ages.sort()
        pop_sizes = [random.randint(1000,10000) for pop in range(2*ntax+1)]
        taxon_set = dendropy.new_taxon_set(ntax)
        species_tree = treesim.pop_gen_tree(taxon_set=taxon_set,
                                                 ages=ages,
                                                 num_genes=4,
                                                 pop_sizes=pop_sizes,
                                                 rng=_RNG)
        ages2 = []
        for node in species_tree.postorder_node_iter():
            distance_from_tip = node.distance_from_tip()
            if distance_from_tip > 0:
                ages2.append(distance_from_tip)
        ages2.sort()
        for index in range(len(ages2)):
            assert (ages[index] - ages2[index]) < 10e-6

        pop_sizes2 = []
        for edge in species_tree.postorder_edge_iter():
            pop_sizes2.append(edge.pop_size)
        pop_sizes2.sort()

        return species_tree

    def runTest(self, ntax=10):
        """TruncatedCoalescentTreeTest -- tree generation without checking [TODO: checks]"""
        species_tree = self.get_species_tree(ntax)
        gene_trees = []
        while len(gene_trees) < 20:
            gene_trees.append(treesim.constrained_kingman(species_tree)[0])

class PureCoalescentTreeTest(unittest.TestCase):

    def runTest(self):
        """PureCoalescentTreeTest -- tree generation without checking [TODO: checks]"""
        _RNG = RepeatedRandom()
        t = treesim.pure_kingman(dendropy.new_taxon_set(100), rng=_RNG)
        assert t._debug_tree_is_valid()

if __name__ == "__main__":
    unittest.main()
