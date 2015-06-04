import dendropy
from dendropy.model import continuous

taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(
        path="primates.cc.tre",
        schema="newick",
        taxon_namespace=taxa)
chars = dendropy.ContinuousCharacterMatrix.get_from_path(
        "primates.cc.nex",
        "nexus",
        taxon_namespace=taxa)
pic = continuous.PhylogeneticIndependentConstrasts(
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










