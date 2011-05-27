#! /usr/bin/env python

import dendropy

tree_str = "[&R] ((A, (B, (C, (D, E)))),(F, (G, H)));"

tree = dendropy.Tree.get_from_string(
        tree_str,
        "newick")

print("Before:")
print(tree.as_string('newick'))
print(tree.as_ascii_plot())
tree.ladderize(False)
print("Ladderize Left:")
print(tree.as_string('newick'))
print(tree.as_ascii_plot())
tree.ladderize(True)
print("Ladderize Right:")
print(tree.as_string('newick'))
print(tree.as_ascii_plot())

