#! /usr/bin/env python

import dendropy
from dendropy import continuous

taxa = dendropy.TaxonSet()
tree = dendropy.Tree.get_from_path(
        "primates.cc.tre",
        "newick",
        taxon_set=taxa)
chars = dendropy.ContinuousCharacterMatrix.get_from_path(
        "primates.cc.nex",
        "nexus",
        taxon_set=taxa)
pic = dendropy.continuous.PhylogeneticIndependentConstrasts(
        tree=tree,
        char_matrix=chars)
for cidx in range(chars.vector_size):
    ctree1 = pic.contrasts_tree(character_index=cidx,
            annotate_pic_statistics=True,
            state_values_as_node_labels=True,
            corrected_edge_lengths=False)
    print(ctree1.as_string("newick",
                suppress_leaf_taxon_labels=True,
                suppress_leaf_node_labels=False,
                suppress_internal_taxon_labels=True,
                suppress_internal_node_labels=False))








