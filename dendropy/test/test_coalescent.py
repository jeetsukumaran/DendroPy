#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Tests of birth-death model fitting.
"""

import unittest
import dendropy
from dendropy.test.support.mockrandom import MockRandom
from dendropy.test.support import pathmap
from dendropy.model import coalescent
from dendropy.simulate import popgensim

class TruncatedCoalescentTreeTest(unittest.TestCase):

    def get_species_tree(self, ntax=10):
        _RNG = MockRandom()
        ages = [_RNG.randint(1000,10000) for age in range(ntax)]
        ages.sort()
        pop_sizes = [_RNG.randint(1000,10000) for pop in range(2*ntax+1)]
        taxon_namespace = dendropy.TaxonNamespace(["t{}".format(i+1) for i in range(ntax)])
        species_tree = popgensim.pop_gen_tree(taxon_namespace=taxon_namespace,
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
            gene_trees.append(coalescent.constrained_kingman_tree(species_tree)[0])

class PureCoalescentTreeTest(unittest.TestCase):

    def runTest(self):
        """PureCoalescentTreeTest -- tree generation without checking [TODO: checks]"""
        _RNG = MockRandom()
        tns = dendropy.TaxonNamespace(["t{}".format(i+1) for i in range(100)])
        t = coalescent.pure_kingman_tree(tns, rng=_RNG)
        assert t._debug_tree_is_valid()

if __name__ == "__main__":
    unittest.main()
