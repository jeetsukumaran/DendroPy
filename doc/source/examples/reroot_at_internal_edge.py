#! /usr/bin/env python

import dendropy

tree_str = "[&R] (A, (B, (C, (D, E))));"

tree = dendropy.Tree.get(
        data=tree_str,
        schema="newick")

print("Before:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())
mrca = tree.mrca(taxon_labels=["D", "E"])
tree.reroot_at_edge(mrca.edge, update_bipartitions=False)
print("After:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

