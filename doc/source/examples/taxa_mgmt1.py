import dendropy

tree_str1 = "((A,B),C);"
tree_str2 = "((A,B),C);"

tree_list = dendropy.TreeList()
tree_list.read(data=tree_str1, schema="newick")
print(tree_list.taxon_namespace)
tree_list.read(data=tree_str2, schema="newick")
print(tree_list.taxon_namespace)

