import dendropy

tree_str = "[&R] (A, (B, (C, (D, E))));"

tree = dendropy.Tree.get_from_string(
        tree_str,
        "newick")

print("Original:")
print(tree.as_ascii_plot())

tree.is_rooted = False
print("After `is_rooted=False`:")
print(tree.as_ascii_plot())

tree.update_bipartitions()
print("After `update_bipartitions()`:")
print(tree.as_ascii_plot())

tree2 = dendropy.Tree.get_from_string(
        tree_str,
        "newick")
tree2.is_rooted = False
tree2.update_bipartitions(suppress_unifurcations=False)
print("After `update_bipartitions(suppress_unifurcations=False)`:")
print(tree2.as_ascii_plot())
