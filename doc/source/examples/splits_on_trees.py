#! /usr/bin/env python

import dendropy

labels = ["A","B","C","D","E","F","G","H"]
taxa = dendropy.TaxonNamespace(labels)
tree = dendropy.Tree.get(
        data="((A,(B,(C,D))),((E,F),(G,H)));",
        schema="newick",
        taxon_namespace=taxa)
tree.is_rooted = False
tree.encode_splits()
for node in tree:
    node.label = taxa.split_as_string(node.edge.split_bitmask)
print(tree.as_ascii_plot(show_internal_node_labels=True,
    width=40))
