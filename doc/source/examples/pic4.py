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
pic = continuous.PhylogeneticIndependentConstrasts(tree=tree, char_matrix=chars)
ctree1 = pic.contrasts_tree(character_index=0,
    annotate_pic_statistics=True,
    state_values_as_node_labels=False,
    corrected_edge_lengths=False)
for nd in ctree1.postorder_internal_node_iter():
    row = [nd.pic_state_value,
            nd.pic_state_variance,
            nd.pic_contrast_raw,
            nd.pic_edge_length_error]
    row_str = [(("%10.8f") % i) for i in row]
    row_str = "    ".join(row_str)
    label = nd.label.ljust(6)
    print "%s %s" % (label, row_str)
