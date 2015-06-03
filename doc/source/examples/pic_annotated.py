#! /usr/bin/env python

import dendropy
from dendropy.model import continuous

taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(
        path="primates.cc.tre",
        schema="newick",
        taxon_namespace=taxa)
chars = dendropy.ContinuousCharacterMatrix.get(
        path="primates.cc.nex",
        schema="nexus",
        taxon_namespace=taxa)
pic = dendropy.continuous.PhylogeneticIndependentConstrasts(
        tree=tree,
        char_matrix=chars)
pic_trees = dendropy.TreeList(taxon_namespace=taxa)
for cidx in range(chars.vector_size):
    ctree1 = pic.contrasts_tree(character_index=cidx)
    ctree1.label = "PIC %d" % (cidx+1)
    pic_trees.append(ctree1)
print(pic_trees.as_string(schema="nexus"))









