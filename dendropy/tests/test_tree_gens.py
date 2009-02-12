#! /usr/bin/env python

############################################################################
##  test_tree_gens.py
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
Tests tree generation.
"""

import random
import unittest
from dendropy import get_logger
from dendropy.tests.debugging_random import DebuggingRandom
import dendropy.dataio as dataio
from dendropy.splits import encode_splits
from dendropy.treedists import symmetric_difference
import dendropy.tests
_LOG = get_logger("TreeGenerationAndSimulation")

### MODULE THAT WE ARE TESTING ###
from dendropy.treegen import *
### MODULE THAT WE ARE TESTING ###

class TreeGenTest(unittest.TestCase):

    def get_species_tree(self, ntax=10):
        "Returns a species population tree, testing it in the process."
        ages = [random.randint(1000,10000) for age in range(ntax)]
        ages.sort()
        pop_sizes = [random.randint(1000,10000) for pop in range(2*ntax+1)]
        taxa_block = random_taxa_block(ntax)
        species_tree = pop_gen_tree(taxa_block=taxa_block,
                                                 ages=ages,
                                                 num_genes=4,
                                                 pop_sizes=pop_sizes)
        nxfpath = dendropy.tests.data_target_path('species_pop.xml')
        nwfpath = dendropy.tests.data_target_path('species_pop.tre')

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
            gene_trees.append(constrained_kingman(species_tree)[0])
            
    def testRandomlyRotate(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        m = [n, n]
        dataset = dataio.trees_from_newick(m)
        trees = [i[0] for i in dataset.trees_blocks]
        ref = trees[0]
        changing = trees[1]
        rng = DebuggingRandom()
        encode_splits(ref)
        encode_splits(changing)
        orig_root = changing.seed_node
        for i in xrange(50):
            randomly_rotate(changing, rng=rng)
            self.assertNotEqual(str(changing), n)
            self.assertEqual(orig_root, changing.seed_node)
            changing.debug_check_tree(logger_obj=_LOG, splits=True)
            if symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

    def testRandomlyReorient(self):
        n = '(Basichlsac,(Lamprothma,Mougeotisp),(((Haplomitr2,Petalaphy),((Angiopteri,(((Azollacaro,((Dennstasam,(Oleandrapi,Polypodapp)),Dicksonant)),Vittarifle),Botrychbit)),(Isoetesmel,((((Agathismac,Agathisova),Pseudotsu),(((Libocedrus,Juniperusc),Callitris),Athrotaxi)),((Liriodchi,Nelumbo),Sagittari))))),Thuidium));'
        m = [n, n]
        dataset = dataio.trees_from_newick(m)
        trees = [i[0] for i in dataset.trees_blocks]
        ref = trees[0]
        changing = trees[1]
        rng = DebuggingRandom()
        encode_splits(ref)
        encode_splits(changing)
        for i in xrange(50):
            randomly_reorient_tree(changing, rng=rng, splits=True)
            self.assertNotEqual(str(changing), n)
            changing.debug_check_tree(logger_obj=_LOG, splits=True)
            if symmetric_difference(ref, changing) != 0:
                self.fail("\n%s\n!=\n%s" % (str(ref), str(changing)))

if __name__ == "__main__":
    unittest.main()
