import dendropy

tree_str1 = "((A,B),C);"

tree1 = dendropy.Tree.get(data=tree_str1, schema="newick")
tree2 = dendropy.Tree.get(
        data=tree_str1,
        schema="newick",
        taxon_namespace=tree1.taxon_namespace)
print(tree1.taxon_namespace is  tree2.taxon_namespace) # True
for nd1, nd2 in zip(tree1, tree2):
    assert nd1.taxon is nd2.taxon # OK

