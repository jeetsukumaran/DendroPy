#! /usr/bin/env python

import dendropy

tree_str = "[&R] ((A, (B, (C, (D, E)))),(F, (G, H)));"

tree = dendropy.Tree.get(
        data=tree_str,
        schema="newick")

print("Before:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())
tree.retain_taxa_with_labels(["A", "C", "G"])
print("After:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

