import dendropy

tree_strs = [
    "     (A, (B, (C, (D, E))));",
    "[&U] (A, (B, (C, (D, E))));",
    "[&R] (A, (B, (C, (D, E))));",
]

rootings = (None, "force-unrooted", "force-rooted", "default-unrooted", "default-rooted",)

for tree_str in tree_strs:
    for rooting in rootings:
        tree = dendropy.Tree.get(
                data=tree_str,
                schema="newick",
                rooting=rooting)
        rooting_token = tree_str[:4]
        rooting_keyword_value = "'{}'".format(rooting) if rooting is not None else "None"
        print("rooting={:20}, token='{}', tree.is_rooted: {}".format(
            rooting_keyword_value, rooting_token, tree.is_rooted))
