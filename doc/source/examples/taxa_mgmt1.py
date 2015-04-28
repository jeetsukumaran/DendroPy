import dendropy

tree_str1 = "((A,B),C);"

tree_list = dendropy.TreeList()
tree_list.read(data=tree_str1, schema="newick")
print(tree_list.taxon_namespace)
tree_list.read(data=tree_str1, schema="newick")
print(tree_list.taxon_namespace)
for nd1, nd2 in zip(tree_list[0], tree_list[1]):
    assert nd1.taxon is nd2.taxon # OK

