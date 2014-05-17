#! /usr/bin/env python

import dendropy
from dendropy import continuous

taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get_from_path(
        "primates.cc.tre",
        "newick",
        taxon_namespace=taxa)
chars = dendropy.ContinuousCharacterMatrix.get_from_path(
        "primates.cc.nex",
        "nexus",
        taxon_namespace=taxa)
pic = dendropy.continuous.PhylogeneticIndependentConstrasts(
        tree=tree,
        char_matrix=chars)
pic_trees = dendropy.TreeList(taxon_namespace=taxa)
for cidx in range(chars.vector_size):
    ctree1 = pic.contrasts_tree(character_index=cidx)
    ctree1.label = "PIC %d" % (cidx+1)
    pic_trees.append(ctree1)
print(pic_trees.as_string("nexus"))









