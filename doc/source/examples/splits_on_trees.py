#! /usr/bin/env python

import dendropy

labels = ["A","B","C","D","E","F","G","H"]
taxa = dendropy.TaxonNamespace(labels)
tree = dendropy.Tree.get(
        data="((A,(B,(C,D))),((E,F),(G,H)));",
        schema="newick",
        taxon_namespace=taxa)
tree.is_rooted = False
tree.encode_bipartitions()
for node in tree:
    node.label = taxa.bitmask_as_bitstring(node.edge.split_bitmask)
print(tree.as_ascii_plot(show_internal_node_labels=True,
    width=40))
