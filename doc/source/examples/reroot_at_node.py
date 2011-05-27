#! /usr/bin/env python

import dendropy

tree_str = "[&R] (A, (B, (C, (D, E))));"

tree = dendropy.Tree.get_from_string(
        tree_str,
        "newick")

print("Before:")
print(tree.as_string('newick'))
print(tree.as_ascii_plot())
mrca = tree.mrca(taxon_labels=["D", "E"])
tree.reroot_at_node(mrca, update_splits=False)
print("After:")
print(tree.as_string('newick'))
print(tree.as_ascii_plot())
