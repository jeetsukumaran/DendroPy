#! /usr/bin/env python

import dendropy

tree_str = "[&R] ((A, (B, (C, (D, E)))),(F, (G, H)));"

tree = dendropy.Tree.get(
        data=tree_str,
        schema="newick")

print("Before:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())
tree.ladderize(ascending=True)
print("Ladderize, ascending=True:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())
tree.ladderize(ascending=False)
print("Ladderize, ascending=False:")
print(tree.as_string(schema='newick'))
print(tree.as_ascii_plot())

