#! /usr/bin/env python

import dendropy

tree_str = "[&R] (A, (B, (C, (D, E))));"

tree = dendropy.Tree.get(
        data=tree_str,
        schema="newick")

print("Before:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())
node_D = tree.find_node_with_taxon_label("D")
tree.reroot_at_edge(node_D.edge, update_bipartitions=False)
print("After:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

